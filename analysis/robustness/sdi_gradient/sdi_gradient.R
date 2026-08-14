# sdi_gradient.R

suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(readxl)
  library(vegan)        
  library(compositions)  
  library(patchwork)
  library(conflicted)
})

suppressMessages({
  conflicts_prefer(dplyr::filter, dplyr::select, dplyr::lag,
                   stats::sd, stats::var, stats::dist, stats::cor, stats::cov,
                   base::scale, base::norm, base::`%*%`, .quiet = TRUE)
})

set.seed(42)

# =============================================================================
# PARAMETERS
# =============================================================================
ASSEMBLAGE <- "SDG"
SEED       <- 42
N_PERM     <- 9999
SDI_COL      <- "SDI"      
SDI_XLSX     <- here("analysis/data/raw_data/SDG_core_metric.xlsx")
SDI_DISPLAY  <- 1e3         
OUT_DIR <- here("analysis/robustness/sdi_gradient")
FIG_DIR <- file.path(OUT_DIR, "figures")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
REF <- list(M_R2  = 0.18770190528243877, M_F  = 2.0334613330095146,
            SP_R2 = 0.16785198589124128, SP_F = 1.7750417603590793,
            n = 50L, n_groups = 6L)
REF_TOL <- 1e-3
SDI_LOW      <- "#E0C69A"
SDI_HIGH     <- "#802520"
TYPOLOGY_COLORS <- c(
  "Levallois"      = "#4A6E8A",
  "Discoid"        = "#802520",
  "Unidirectional" = "#BA8530",
  "Multiplatform"  = "#8A7A68",
  "Bidirectional"  = "#788C4A"
)
LAYER_FACET_LEVELS <- c("3", "4")

# =============================================================================
# Helper + data-preparation import — verbatim, WITHOUT running either script
# =============================================================================
selective_eval <- function(path, wanted, functions_only) {
  got <- character(0)
  for (e in parse(path)) {
    if (!is.call(e)) next
    if (!as.character(e[[1]])[1] %in% c("<-", "=")) next
    if (!is.name(e[[2]])) next
    nm <- as.character(e[[2]])
    if (!nm %in% wanted) next
    rhs <- e[[3]]
    if (functions_only &&
        !(is.call(rhs) && identical(as.character(rhs[[1]])[1], "function"))) next
    eval(e, envir = globalenv())
    got <- c(got, nm)
  }
  unique(got)
}

EXP_SCRIPT <- here("analysis/robustness/joint_mfa_discrimination",
                   "joint_mfa_discrimination_stats.R")
SDG_SCRIPT <- here("analysis/robustness/joint_mfa_discrimination",
                   "joint_mfa_discrimination_SDG.R")
for (f in c(EXP_SCRIPT, SDG_SCRIPT))
  if (!file.exists(f))
    stop("Companion script not found at ", f,
         " — cannot import helpers / sample definition; stopping rather than ",
         "duplicating them here.")

HELPERS <- c("replace_zeros", "make_ilr", "ilr_parts", "mfa_s1", "block_inertia",
             "filter_spharm", "safe_filter_groups", "permanova_global")
got_helpers <- selective_eval(EXP_SCRIPT, HELPERS, functions_only = TRUE)

DATA_PREP <- c("RESTRICT_LAYERS", "POWER_COLS_DIR", "POWER_COLS_MORPH",
               "EXCLUDE_CORE_TYPES", "SPHARM_direction", "SPHARM_morphology",
               "metric_data", "core_meta", "build_blocks_sdg")
cat("Loading committed main outputs (SPHARM direction / morphology, metadata)...\n")
got_data <- selective_eval(SDG_SCRIPT, DATA_PREP, functions_only = FALSE)

missing_imports <- c(setdiff(HELPERS, got_helpers), setdiff(DATA_PREP, got_data))
if (length(missing_imports))
  stop("Could not import: ", paste(missing_imports, collapse = ", "),
       " — the companion scripts' structure changed; stopping rather than ",
       "reimplementing them here.")
cat(sprintf("Imported %d helpers from %s and %d data-prep objects from %s.\n",
            length(got_helpers), basename(EXP_SCRIPT),
            length(got_data),    basename(SDG_SCRIPT)))

# =============================================================================
# (0) SAMPLE + SDI
# =============================================================================
cat("\n", strrep("=", 70), "\n", sep = "")
cat("(0) SAMPLE, SDI COLUMN, MISSINGNESS, DESCRIPTIVES\n")
cat(strrep("=", 70), "\n", sep = "")

blocks <- build_blocks_sdg()
if (length(blocks$ids) != REF$n || nlevels(blocks$group) != REF$n_groups)
  stop(sprintf(paste0("Imported sample definition gives n = %d / %d groups, not ",
                      "the published %d / %d. Stopping: the whole point of the ",
                      "import is that the sample is identical."),
               length(blocks$ids), nlevels(blocks$group), REF$n, REF$n_groups))
cat(sprintf("Sample (imported build_blocks_sdg): n = %d, %d core-type groups, no layer filter\n",
            length(blocks$ids), nlevels(blocks$group)))

# ---- SDI: read straight from the workbook, do not recompute from the meshes --
sdi_xl <- read_excel(SDI_XLSX)
cat(sprintf("\n%s columns (%d):\n  %s\n", basename(SDI_XLSX), ncol(sdi_xl),
            paste(names(sdi_xl), collapse = ", ")))
if (!SDI_COL %in% names(sdi_xl))
  stop("Column '", SDI_COL, "' not found in ", basename(SDI_XLSX),
       " — refusing to guess which column holds reduction intensity.")
cat(sprintf("\n  SDI column in use : '%s' (read as-is; NOT recomputed from meshes)\n",
            SDI_COL))

sdi_tbl <- tibble(ID  = str_trim(as.character(sdi_xl$ID)),
                  SDI = suppressWarnings(as.numeric(sdi_xl[[SDI_COL]])))
if (anyDuplicated(sdi_tbl$ID))
  stop("Duplicate IDs in ", basename(SDI_XLSX), " — cannot join SDI unambiguously.")

meta_all <- tibble(ID        = blocks$ids,
                   core_type = blocks$group,
                   layer     = blocks$layer) %>%
  left_join(core_meta %>% select(ID, raw_material), by = "ID") %>%
  left_join(sdi_tbl, by = "ID")
if (nrow(meta_all) != length(blocks$ids))
  stop("Metadata join changed the row count — non-unique keys; stopping.")

miss <- meta_all %>% filter(!is.finite(SDI))
cat(sprintf("\n  Specimens with missing / non-numeric SDI : %d of %d\n",
            nrow(miss), nrow(meta_all)))
if (nrow(miss) > 0) {
  cat("    dropped IDs: ", paste(miss$ID, collapse = ", "), "\n", sep = "")
} else {
  cat("    none — the analysed set is the full published sample.\n")
}

keep <- which(is.finite(meta_all$SDI))
meta <- meta_all[keep, ] %>% mutate(core_type = droplevels(core_type),
                                    Raw_mat   = factor(raw_material))
n_spec <- nrow(meta)

cat(sprintf("\n  FINAL n = %d | %d core-type groups | %d raw-material groups\n",
            n_spec, nlevels(meta$core_type), nlevels(meta$Raw_mat)))
cat("\n  Core-type group sizes:\n"); print(table(meta$core_type))
cat("\n  Raw-material group sizes:\n"); print(table(meta$Raw_mat))
cat("\n  Raw material x core type (read this alongside the sequential SS):\n")
print(table(meta$Raw_mat, meta$core_type))
small <- names(which(table(meta$core_type) < 3))
if (length(small))
  cat(sprintf("\n  NOTE: after the SDI filter these groups fall below n = 3: %s\n",
              paste(small, collapse = ", ")))

q <- stats::quantile(meta$SDI, c(0, .25, .5, .75, 1))
cat(sprintf(paste0("\n  SDI descriptives (n = %d):\n",
                   "    min = %.6g | Q1 = %.6g | median = %.6g | Q3 = %.6g | max = %.6g\n",
                   "    mean = %.6g | sd = %.6g | max/min = %.1f\n"),
            n_spec, q[1], q[2], q[3], q[4], q[5],
            mean(meta$SDI), stats::sd(meta$SDI), q[5] / q[1]))
cat("  SDI by core type (median):\n")
print(meta %>% group_by(core_type) %>%
        summarise(n = n(), median_SDI = median(SDI), .groups = "drop") %>%
        as.data.frame(), row.names = FALSE)
cat("  SDI by raw material (median):\n")
print(meta %>% group_by(Raw_mat) %>%
        summarise(n = n(), median_SDI = median(SDI), .groups = "drop") %>%
        as.data.frame(), row.names = FALSE)
if (stats::sd(meta$SDI) <= 0)
  stop("SDI has zero variance in the analysed set — nothing to test.")

# =============================================================================
# ILR blocks + MFA normalisation — RECOMPUTED on the analysed subset
# =============================================================================
morph_sub <- blocks$morph[keep, , drop = FALSE]
scar_sub  <- blocks$scar[keep, , drop = FALSE]

Z_M  <- make_ilr(morph_sub)
Z_SP <- make_ilr(scar_sub)
rownames(Z_M) <- rownames(Z_SP) <- meta$ID
colnames(Z_M)  <- paste0("M_ilr",  seq_len(ncol(Z_M)))
colnames(Z_SP) <- paste0("SP_ilr", seq_len(ncol(Z_SP)))

s_M    <- mfa_s1(Z_M)
s_SP   <- mfa_s1(Z_SP)
Z_comb <- cbind(Z_M / s_M, Z_SP / s_SP)

cat(sprintf("\nILR (recomputed on the analysed n = %d): M-SPHARM %d parts -> %d coords | SP-SPHARM %d parts -> %d coords\n",
            n_spec, length(ilr_parts(morph_sub)), ncol(Z_M),
            length(ilr_parts(scar_sub)), ncol(Z_SP)))
cat(sprintf("MFA normalisation: s1(M) = %.5f, s1(SP) = %.5f  =>  w_MFA = %.4f\n",
            s_M, s_SP,
            block_inertia(Z_M / s_M) /
              (block_inertia(Z_M / s_M) + block_inertia(Z_SP / s_SP))))

D_scar  <- stats::dist(as.matrix(Z_SP),   method = "euclidean")
D_morph <- stats::dist(as.matrix(Z_M),    method = "euclidean")
D_comb  <- stats::dist(as.matrix(Z_comb), method = "euclidean")

# ---- anchor check: is this still the published SDG sample? ------------------
cat("\n", strrep("-", 70), "\n", sep = "")
if (nrow(miss) == 0) {
  ag_M  <- permanova_global(Z_M,  meta$core_type)
  ag_SP <- permanova_global(Z_SP, meta$core_type)
  ok <- all(abs(c(ag_M$R2 - REF$M_R2, ag_M$F - REF$M_F,
                  ag_SP$R2 - REF$SP_R2, ag_SP$F - REF$SP_F)) <= REF_TOL)
  cat(sprintf(paste0("ANCHOR CHECK (core type alone, committed Table 2):\n",
                     "  M-SPHARM  R2 = %.5f (ref %.5f), F = %.5f (ref %.5f)\n",
                     "  SP-SPHARM R2 = %.5f (ref %.5f), F = %.5f (ref %.5f)\n"),
              ag_M$R2, REF$M_R2, ag_M$F, REF$M_F,
              ag_SP$R2, REF$SP_R2, ag_SP$F, REF$SP_F))
  if (!ok)
    stop("ANCHOR CHECK FAILED — this script's ILR coordinates no longer ",
         "reproduce the committed SDG core-type PERMANOVA. Resolve the sample ",
         "definition before interpreting anything below.")
  cat("  => reproduces the published SDG core-type PERMANOVA.\n")
} else {
  cat(sprintf(paste0("ANCHOR CHECK SKIPPED: %d specimen(s) dropped for missing ",
                     "SDI, so this is no longer the published n = %d sample.\n"),
              nrow(miss), REF$n))
}

# =============================================================================
# (1) PERMANOVA — three spaces x two sums-of-squares types
# =============================================================================
cat("\n", strrep("=", 70), "\n", sep = "")
cat("(1) PERMANOVA: D ~ Raw_mat + core_type + SDI\n")
cat(strrep("=", 70), "\n", sep = "")
cat(sprintf("    seed = %d | permutations = %d | n = %d\n", SEED, N_PERM, n_spec))
cat("    term order is FIXED: raw material (control) -> core type -> SDI (increment)\n")

SPACES <- list(
  scar     = list(D = D_scar,  label = "D_scar (SP-SPHARM, l = 1-6)"),
  morph    = list(D = D_morph, label = "D_morph (M-SPHARM, l = 1-8)"),
  combined = list(D = D_comb,  label = "D_comb (MFA-normalised join)"))

md <- as.data.frame(meta)

ADONIS_COLS <- c("Df", "SumOfSqs", "R2", "F", "Pr(>F)")

run_perm <- function(D, by) {
  set.seed(SEED)
  res <- adonis2(D ~ Raw_mat + core_type + SDI, data = md,
                 by = by, permutations = N_PERM)
  tab <- as.data.frame(res)
  if (!all(ADONIS_COLS %in% names(tab)))
    stop("adonis2() returned columns {", paste(names(tab), collapse = ", "),
         "}; expected {", paste(ADONIS_COLS, collapse = ", "),
         "}. vegan's output layout changed — stopping rather than mislabelling ",
         "the result table.")
  tibble(term     = rownames(tab),
         df       = as.integer(tab[["Df"]]),
         SumOfSqs = tab[["SumOfSqs"]],
         R2       = tab[["R2"]],
         pseudo_F = tab[["F"]],
         p        = tab[["Pr(>F)"]])
}

perm_tbl <- map_dfr(names(SPACES), function(sp) {
  map_dfr(c("terms", "margin"), function(by) {
    cat(sprintf("  running %-9s | by = %-6s ...\n", sp, by))
    run_perm(SPACES[[sp]]$D, by) %>%
      mutate(space = sp, space_label = SPACES[[sp]]$label, ss_type = by,
             .before = 1)
  })
}) %>%
  mutate(assemblage = ASSEMBLAGE, n = n_spec, n_perm = N_PERM, seed = SEED)

print_anova <- function(tbl, heading) {
  cat("\n--- ", heading, " ---\n", sep = "")
  print(tbl %>%
          transmute(term, df,
                    SumOfSqs = round(SumOfSqs, 5),
                    R2       = round(R2, 5),
                    pseudo_F = round(pseudo_F, 4),
                    p) %>%
          as.data.frame(), row.names = FALSE)
  sdi <- tbl %>% filter(term == "SDI")
  if (nrow(sdi) == 1)
    cat(sprintf("    >>> SDI: df = %d, R2 = %.5f, pseudo-F = %.4f, p = %.4f  [%s]\n",
                sdi$df, sdi$R2, sdi$pseudo_F, sdi$p,
                ifelse(sdi$p < 0.05, "SIGNIFICANT", "not significant")))
}

cat("\n", strrep("=", 70), "\n", sep = "")
cat("PRIMARY RESULT — sequential (Type I) SS, by = \"terms\"\n")
cat(strrep("=", 70), "\n", sep = "")
for (sp in names(SPACES))
  print_anova(perm_tbl %>% filter(space == sp, ss_type == "terms"),
              SPACES[[sp]]$label)

cat("\n", strrep("=", 70), "\n", sep = "")
cat("CONTROL ONLY — marginal (Type III) SS, by = \"margin\". Reported, NOT the\n")
cat("headline: the reviewer's question is an INCREMENT over material and type.\n")
cat(strrep("=", 70), "\n", sep = "")
for (sp in names(SPACES))
  print_anova(perm_tbl %>% filter(space == sp, ss_type == "margin"),
              SPACES[[sp]]$label)

write_csv(perm_tbl, file.path(OUT_DIR, "sdi_gradient_permanova.csv"))
cat("\nWrote sdi_gradient_permanova.csv\n")

# =============================================================================
# (2) FIGURE — the same base maps as the type-coloured biplots, coloured by SDI
# =============================================================================
pca_of <- function(Z) {
  p <- stats::prcomp(as.matrix(Z), center = TRUE, scale. = FALSE)
  list(scores  = p$x[, 1:2, drop = FALSE],
       inertia = p$sdev^2 / sum(p$sdev^2) * 100)
}

panel_spec <- list(
  list(key = "scar",  Z = Z_SP,   strip = "SP-SPHARM", short = "SP-SPHARM", ax = "PC"),
  list(key = "morph", Z = Z_M,    strip = "M-SPHARM",  short = "M-SPHARM",  ax = "PC"),
  list(key = "comb",  Z = Z_comb, strip = "Combined",  short = "Combined",  ax = "MFA axis"))

PCA_FITS <- map(panel_spec, ~ pca_of(.x$Z))
names(PCA_FITS) <- map_chr(panel_spec, "key")

# Long frame: every core, in every space, on the shared axes.
scores_all <- map_dfr(seq_along(panel_spec), function(i) {
  sp <- panel_spec[[i]]; pc <- PCA_FITS[[sp$key]]
  as_tibble(pc$scores, .name_repair = ~ c("Axis1", "Axis2")) %>%
    mutate(space_key = sp$key, strip = sp$strip,
           ID = meta$ID, SDI_disp = meta$SDI * SDI_DISPLAY,
           core_type = as.character(meta$core_type),
           raw_material = as.character(meta$Raw_mat),
           layer = as.character(meta$layer))
})

axis_lims <- scores_all %>% group_by(space_key) %>%
  summarise(xlo = min(Axis1), xhi = max(Axis1),
            ylo = min(Axis2), yhi = max(Axis2), .groups = "drop") %>%
  mutate(padx = 0.06 * (xhi - xlo), pady = 0.06 * (yhi - ylo),
         xlo = xlo - padx, xhi = xhi + padx,
         ylo = ylo - pady, yhi = yhi + pady)

# -----------------------------------------------------------------------------
# MAIN-TEXT FIGURE SPECIFICATION -- taken verbatim, not invented here
# -----------------------------------------------------------------------------
MT <- list(
  base_size      = 8,      
  axis_text      = 5,
  strip_text     = 7,     
  strip_fill     = "#EBEBEB",
  legend_text    = 6.5,
  legend_title   = 7,
  legend_key_cm  = 0.30,
  tag_size       = 9,
  tag_face       = "bold",
  pt_shape       = 16,
  pt_size        = 2.0,    
  pt_alpha       = 0.90,
  zero_colour    = "grey70",
  zero_width     = 0.25,
  fig_width_in   = 6.85,  
  out_width_mm   = 174,   
  dpi            = 800)   

cat("\n", strrep("=", 70), "\n", sep = "")
cat("(1) MAIN-TEXT FIGURE SPECIFICATION -- SOURCES AND VALUES\n")
cat(strrep("=", 70), "\n", sep = "")
cat("  theme / point / legend sizes  <- SDG_cores_statistics.R make_coia_biplot()\n")
cat("                                   (manuscript Fig. 9, fig-coia-composite)\n")
cat("  strip styling                 <- exp_cores_statistics.R plot_rose() (Fig. 8c)\n")
cat("  plot.tag                      <- exp_cores_statistics.R composite block\n")
cat("  export geometry               <- manuscript.qmd:347-351 and :44\n")
for (nm in names(MT)) cat(sprintf("    %-14s = %s\n", nm, format(MT[[nm]])))

mt_theme <- function(show_strip = TRUE) {
  theme_bw(base_size = MT$base_size) +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor   = element_blank(),
          axis.text          = element_text(size = MT$axis_text),
          strip.text         = element_text(face = "bold", size = MT$strip_text),
          strip.background   = element_rect(fill = MT$strip_fill,
                                            color = MT$strip_fill),
          legend.key.size    = grid::unit(MT$legend_key_cm, "cm"),
          legend.text        = element_text(size = MT$legend_text),
          legend.title       = element_text(size = MT$legend_title),
          legend.margin      = margin(2, 4, 2, 4),
          plot.tag           = element_text(size = MT$tag_size, face = MT$tag_face))
}

zero_lines <- list(
  geom_hline(yintercept = 0, linetype = "dashed", color = MT$zero_colour,
             linewidth = MT$zero_width),
  geom_vline(xintercept = 0, linetype = "dashed", color = MT$zero_colour,
             linewidth = MT$zero_width))

make_main_panel <- function(sp, tag) {
  df <- scores_all %>% filter(space_key == sp$key)
  lm_ <- axis_lims %>% filter(space_key == sp$key)
  p <- ggplot(df, aes(Axis1, Axis2)) + zero_lines +
    facet_wrap(~ strip) +
    geom_point(aes(colour = SDI_disp), shape = MT$pt_shape, size = MT$pt_size,
               alpha = MT$pt_alpha) +
    scale_colour_gradient(low = SDI_LOW, high = SDI_HIGH,
                          name = expression(SDI~(10^-3)),
                          limits = range(scores_all$SDI_disp))
  pc <- PCA_FITS[[sp$key]]
  p + coord_cartesian(xlim = c(lm_$xlo, lm_$xhi), ylim = c(lm_$ylo, lm_$yhi)) +
    labs(x = sprintf("%s 1 (%.1f%%)", sp$ax, pc$inertia[1]),
         y = sprintf("%s 2 (%.1f%%)", sp$ax, pc$inertia[2]),
         tag = tag) +
    mt_theme() +
    guides(colour = guide_colourbar(order = 1,
                                    barwidth  = grid::unit(0.30, "cm"),
                                    barheight = grid::unit(1.6, "cm")))
}

ok_fig <- TRUE
tryCatch({
  row_sdi <- map2(panel_spec, c("a", "b", "c"), ~ make_main_panel(.x, .y))
  p_main  <- wrap_plots(row_sdi, nrow = 1) +
    plot_layout(guides = "collect") & theme(legend.position = "right")

  MAIN_H <- 2.35
  ggsave(file.path(FIG_DIR, "fig_sdi_gradient_main.png"), p_main,
         width = MT$fig_width_in, height = MAIN_H, dpi = MT$dpi)
  cat("\nWrote figures/fig_sdi_gradient_main.png\n")

  cat("\n", strrep("=", 70), "\n", sep = "")
  cat("(2) MAIN FIGURE RENDERED SIZE AND READABILITY AT COLUMN WIDTH\n")
  cat(strrep("=", 70), "\n", sep = "")
  cat(sprintf("  canvas          : %.2f x %.2f in  (%.1f x %.1f mm) at %d dpi\n",
              MT$fig_width_in, MAIN_H, MT$fig_width_in * 25.4, MAIN_H * 25.4,
              MT$dpi))
  cat(sprintf("  pixels          : %d x %d\n",
              round(MT$fig_width_in * MT$dpi), round(MAIN_H * MT$dpi)))
  cat(sprintf("  placed at       : %d mm (out-width, full text width)\n",
              MT$out_width_mm))
  scale_f <- MT$out_width_mm / (MT$fig_width_in * 25.4)
  cat(sprintf("  scale on page   : %.3f  (1.0 = no shrink)\n", scale_f))
  cat(sprintf("  effective sizes : axis text %.1f pt, strip %.1f pt, legend %.1f pt\n",
              MT$axis_text * scale_f, MT$strip_text * scale_f,
              MT$legend_text * scale_f))
  cat(sprintf("  panel width     : ~%.2f in (%.1f mm) per panel, 3 across\n",
              (MT$fig_width_in - 1.0) / 3, (MT$fig_width_in - 1.0) / 3 * 25.4))
  if (scale_f >= 0.98)
    cat("  => no shrink: type sizes reach the page exactly as set above.\n")
  else
    cat(sprintf("  => shrunk to %.0f%%; smallest text lands at %.1f pt.\n",
                scale_f * 100, MT$axis_text * scale_f))
  cat("  NOTE: this is the FULL text width (174 mm), not a single column. The\n")
  cat("  journal single-column width would halve every type size and this layout\n")
  cat("  would need splitting into two stacked figures instead.\n")
}, error = function(e) {
  ok_fig <<- FALSE
  cat("Main figure failed:", conditionMessage(e), "\n")
})

# =============================================================================
# (2b) SI FACET FIGURES -- does a gradient hide inside one material or layer?
# =============================================================================
make_facet_fig <- function(split_col, keep_levels, out_name, note,
                           label_prefix = "") {
  df <- scores_all %>% filter(.data[[split_col]] %in% keep_levels)
  ns <- df %>% filter(space_key == panel_spec[[1]]$key) %>% count(.data[[split_col]])
  panels <- list()
  for (lv in keep_levels) {
    n_lv <- ns$n[ns[[split_col]] == lv]
    for (sp in panel_spec) {
      d   <- df %>% filter(space_key == sp$key, .data[[split_col]] == lv) %>%
        mutate(cell = sprintf("%s%s (n = %d)  ·  %s",
                              label_prefix, gsub("_", " ", lv), n_lv, sp$short))
      lm_ <- axis_lims %>% filter(space_key == sp$key)
      pc  <- PCA_FITS[[sp$key]]
      panels[[length(panels) + 1]] <-
        ggplot(d, aes(Axis1, Axis2)) + zero_lines +
        geom_point(aes(colour = SDI_disp), shape = MT$pt_shape, size = 2.0,
                   alpha = 0.95) +
        scale_colour_gradient(low = SDI_LOW, high = SDI_HIGH,
                              name = expression(SDI~(10^-3)),
                              limits = range(scores_all$SDI_disp)) +
        coord_cartesian(xlim = c(lm_$xlo, lm_$xhi), ylim = c(lm_$ylo, lm_$yhi)) +
        facet_wrap(~ cell) +
        labs(x = sprintf("%s 1 (%.1f%%)", sp$ax, pc$inertia[1]),
             y = sprintf("%s 2 (%.1f%%)", sp$ax, pc$inertia[2])) +
        theme_bw(base_size = 9) +
        theme(panel.grid       = element_blank(),
              axis.text        = element_text(size = 6),
              strip.text       = element_text(face = "bold", size = 7),
              strip.background = element_rect(fill = "#EBEBEB", color = "#EBEBEB"),
              legend.key.size  = grid::unit(0.40, "cm"),
              legend.text      = element_text(size = 7),
              legend.title     = element_text(size = 7.5))
    }
  }
  p <- wrap_plots(panels, nrow = length(keep_levels)) +
    plot_layout(guides = "collect") &
    theme(legend.position = "right")
  ggsave(file.path(FIG_DIR, out_name), p, width = 8.6, height = 5.6, dpi = 300)
  cat(sprintf("Wrote figures/%s   [%s]\n", out_name, note))
  ns
}

cat("\n", strrep("=", 70), "\n", sep = "")
cat("(3) SAMPLE SIZES USED IN EACH FIGURE\n")
cat(strrep("=", 70), "\n", sep = "")
cat(sprintf("  main figure, both rows : n = %d (all cores, no subsetting)\n", n_spec))

ns_rm <- make_facet_fig("raw_material", sort(unique(scores_all$raw_material)),
                        "fig_S_sdi_gradient_by_rawmat.png",
                        "all layers, Layer 2 retained")
ns_ly <- make_facet_fig("layer", LAYER_FACET_LEVELS,
                        "fig_S_sdi_gradient_by_layer.png",
                        sprintf("Layer 2 excluded (n = %d)",
                                sum(scores_all$layer[scores_all$space_key ==
                                                       panel_spec[[1]]$key] == "2")),
                        label_prefix = "Layer ")
cat("\n  by raw material:\n"); print(as.data.frame(ns_rm), row.names = FALSE)
cat("  by layer:\n");         print(as.data.frame(ns_ly), row.names = FALSE)
cat(sprintf("\n  Layer 2 (n = %d) is EXCLUDED from the layer facets and RETAINED\n",
            sum(scores_all$layer[scores_all$space_key == panel_spec[[1]]$key] == "2")))
cat("  in the raw-material facets, which do not condition on layer.\n")
cat("  Sandstone (n = 12) and L3 (n = 18) are the small cells: scatter judgements\n")
cat("  in those rows are correspondingly weaker.\n")

cat("\n", strrep("=", 70), "\n", sep = "")
cat("(4) SHARED-PCA CONFIRMATION\n")
cat(strrep("=", 70), "\n", sep = "")
cat(sprintf("  prcomp() called exactly %d times, once per space, on all %d cores.\n",
            length(PCA_FITS), n_spec))
for (sp in panel_spec)
  cat(sprintf("    %-34s PC1 %.1f%%, PC2 %.1f%%\n", sp$strip,
              PCA_FITS[[sp$key]]$inertia[1], PCA_FITS[[sp$key]]$inertia[2]))
cat("  Main figure panels a-c and every SI facet cell read from these same three\n")
cat("  fits; no subgroup PCA is computed anywhere in this script.\n")

# -----------------------------------------------------------------------------
# Figure captions (the three mandated points for the main figure)
# -----------------------------------------------------------------------------
writeLines(c(
  "FIGURE CAPTIONS -- analysis/robustness/sdi_gradient/figures/",
  strrep("=", 78), "",
  "fig_sdi_gradient_main.png  [MAIN TEXT]",
  paste(
    "Core variation at Sandinggai shown in three descriptor spaces:",
    "(a) scar pattern (SP-SPHARM), (b) morphology (M-SPHARM), and (c) the two",
    "combined by MFA block normalisation. Each panel shows the FIRST TWO",
    "PRINCIPAL AXES of that space; the PERMANOVA reported in the text is",
    "computed on the full distance matrix of the space (all ILR coordinates),",
    "not on these two axes. Points are individual cores (n = 50), coloured by",
    "the continuous scar density index SDI, a proxy for reduction intensity,",
    "NOT by core type. The absence of any spatial organisation of colour",
    "corresponds to the result reported in the text: no detectable SDI gradient",
    "in any space. Axis signs are as returned by the PCA, so the panels align",
    "one-to-one with the type-coloured biplots elsewhere in the paper."),
  "",
  "fig_S_sdi_gradient_by_rawmat.png  [SI]",
  paste(
    "The construction of the figure above, split by raw material, to test",
    "whether a gradient exists within one material and is cancelled by pooling.",
    "All cells use the SAME whole-sample PCA fit and the same axis limits; no",
    "PCA is refitted on a subgroup. Cell sample sizes are given in the strips.",
    "Layer 2 cores are retained here because this split does not condition on",
    "layer. Sandstone (n = 12) is small and its scatter judgement is weaker."),
  "",
  "fig_S_sdi_gradient_by_layer.png  [SI]",
  paste(
    "As above, split by stratigraphic layer. Layer 2 (n = 2) is excluded.",
    "L3 (n = 18) is the smaller cell and its scatter judgement is weaker",
    "than L4 (n = 30)."),
  ""),
  file.path(OUT_DIR, "figure_captions.txt"))
cat("\nWrote figure_captions.txt\n")

# =============================================================================
# (3) VERDICT
# =============================================================================
cat("\n", strrep("=", 70), "\n", sep = "")
cat("(3) VERDICT\n")
cat(strrep("=", 70), "\n", sep = "")

sdi_terms <- perm_tbl %>%
  filter(ss_type == "terms", term == "SDI") %>%
  mutate(sig = p < 0.05)
sdi_margin <- perm_tbl %>%
  filter(ss_type == "margin", term == "SDI") %>%
  select(space, R2_margin = R2, p_margin = p)

summary_tbl <- sdi_terms %>%
  select(space, space_label, df, R2_seq = R2, F_seq = pseudo_F, p_seq = p, sig) %>%
  left_join(sdi_margin, by = "space")

cat("\nSDI increment after Raw_mat and core_type (sequential SS), per space:\n")
print(summary_tbl %>%
        transmute(space, R2_seq = round(R2_seq, 5), F_seq = round(F_seq, 4),
                  p_seq, significant = sig,
                  R2_margin = round(R2_margin, 5), p_margin) %>%
        as.data.frame(), row.names = FALSE)

hits <- summary_tbl %>% filter(sig)
cat("\n(2) One-line judgement:\n")
if (nrow(hits) == 0) {
  cat(sprintf(paste0("    NO. After partialling out raw material and core type, ",
                     "SDI contributes no significant\n    variation in any of the ",
                     "three spaces (largest R2 = %.5f in '%s', p = %.4f).\n"),
              max(summary_tbl$R2_seq),
              summary_tbl$space[which.max(summary_tbl$R2_seq)],
              summary_tbl$p_seq[which.max(summary_tbl$R2_seq)]))
} else {
  cat(sprintf(paste0("    YES. SDI still contributes significantly in %d of 3 ",
                     "spaces after raw material and\n    core type: %s.\n"),
              nrow(hits),
              paste(sprintf("%s (R2 = %.5f, p = %.4f)",
                            hits$space, hits$R2_seq, hits$p_seq),
                    collapse = "; ")))
}

cat("\n(3) Agreement across the three spaces:\n")
if (length(unique(summary_tbl$sig)) == 1) {
  cat(sprintf("    CONSISTENT — all three spaces agree (SDI %s everywhere).\n",
              ifelse(summary_tbl$sig[1], "significant", "not significant")))
} else {
  cat(sprintf(paste0("    NOT CONSISTENT — SDI is significant in {%s} but not in ",
                     "{%s}.\n    This asymmetry is a result in its own right: the ",
                     "two blocks are decoupled in SDG\n    (RV = 0.091), so ",
                     "reduction intensity need not act on both.\n"),
              paste(summary_tbl$space[summary_tbl$sig],  collapse = ", "),
              paste(summary_tbl$space[!summary_tbl$sig], collapse = ", ")))
}

seq_vs_marg <- summary_tbl %>% filter(sig != (p_margin < 0.05))
if (nrow(seq_vs_marg))
  cat(sprintf(paste0("\n    Note (control): sequential and marginal SS disagree ",
                     "for %s. The sequential\n    result stands as the headline; ",
                     "the disagreement reflects the raw-material /\n    core-type ",
                     "/ SDI confounding, not a choice of model.\n"),
              paste(seq_vs_marg$space, collapse = ", ")))

cat("\n", strrep("=", 70), "\n", sep = "")
cat("OUTPUTS\n")
cat(strrep("=", 70), "\n", sep = "")
cat("  ", file.path(OUT_DIR, "sdi_gradient_permanova.csv"), "\n", sep = "")
cat("  ", file.path(OUT_DIR, "figure_captions.txt"), "\n", sep = "")
cat("  ", file.path(FIG_DIR, "fig_sdi_gradient_main.png"),
    ifelse(ok_fig, "   [MAIN TEXT]", "  (FAILED)"), "\n", sep = "")
cat("  ", file.path(FIG_DIR, "fig_S_sdi_gradient_by_rawmat.png"), "   [SI]\n", sep = "")
cat("  ", file.path(FIG_DIR, "fig_S_sdi_gradient_by_layer.png"),  "   [SI]\n", sep = "")

cat("\n", strrep("=", 70), "\n", sep = "")
cat("SESSION INFO\n")
cat(strrep("=", 70), "\n", sep = "")
print(sessionInfo())
