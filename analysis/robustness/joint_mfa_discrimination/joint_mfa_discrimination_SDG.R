# joint_mfa_discrimination_SDG.R

suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(readxl)
  library(vegan)         
  library(compositions)   
  library(ggrepel)
  library(conflicted)
})

suppressMessages({
  conflicts_prefer(dplyr::filter, dplyr::select, dplyr::lag,
                   stats::sd, stats::var, stats::dist, stats::cor, stats::cov,
                   base::scale, base::norm, base::`%*%`, .quiet = TRUE)
})

set.seed(42)

# =============================================================================
# PARAMETERS (seed / permutations identical to EXP)
# =============================================================================
ASSEMBLAGE <- "SDG"
SEED       <- 42
N_PERM     <- 9999
W_GRID_COARSE   <- seq(0, 1, by = 0.05)
DENSE_STEP      <- 0.01
DENSE_HALFWIDTH <- 0.05   
DOMINANCE_THR  <- 0.65    
DEGENERACY_THR <- 0.10    
RESTRICT_LAYERS <- NULL
OUT_DIR <- here("analysis/robustness/joint_mfa_discrimination/SDG")
FIG_DIR <- file.path(OUT_DIR, "figures")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
POWER_COLS_DIR   <- paste0("power_l", 1:6)   
POWER_COLS_MORPH <- paste0("power_l", 1:8)   
EXCLUDE_CORE_TYPES <- c("Handaxe", "Pick")   
TYPOLOGY_COLORS <- c(
  "Levallois"      = "#4A6E8A",
  "Discoid"        = "#802520",
  "Unidirectional" = "#BA8530",
  "Multiplatform"  = "#8A7A68",
  "Bidirectional"  = "#788C4A"
)
REF <- list(
  M_R2  = 0.18770190528243877, M_F  = 2.0334613330095146, M_p  = 0.0104,
  SP_R2 = 0.16785198589124128, SP_F = 1.7750417603590793, SP_p = 0.0405,
  n = 50, n_groups = 6
)
REF_TOL <- 1e-3   

EXP_REF <- list(
  n = 58L, n_groups = 5L, n_pairs = 10L,
  M_R2 = 0.12334, SP_R2 = 0.30174,
  best_w = 0.25, plateau_end = 0.50, rho_improve = -0.845
)

# =============================================================================
# Helper import — the EXP script's helpers, verbatim, WITHOUT running EXP
# =============================================================================
EXP_SCRIPT <- here("analysis/robustness/joint_mfa_discrimination",
                   "joint_mfa_discrimination_stats.R")
if (!file.exists(EXP_SCRIPT))
  stop("EXP companion script not found at ", EXP_SCRIPT,
       " — cannot import helpers; stopping rather than duplicating them here.")

HELPERS <- c("replace_zeros", "make_ilr", "ilr_parts", "mfa_s1", "block_inertia",
             "permanova_global", "permanova_pairwise", "permdisp", "pair_p",
             "filter_spharm", "split_by_group", "safe_filter_groups")

imported <- character(0)
for (e in parse(EXP_SCRIPT)) {
  if (!is.call(e)) next
  if (!as.character(e[[1]])[1] %in% c("<-", "=")) next
  if (!is.name(e[[2]])) next
  nm <- as.character(e[[2]])
  if (!nm %in% HELPERS) next
  rhs <- e[[3]]
  if (!(is.call(rhs) && identical(as.character(rhs[[1]])[1], "function"))) next
  eval(e, envir = globalenv())
  imported <- c(imported, nm)
}
missing_helpers <- setdiff(HELPERS, imported)
if (length(missing_helpers))
  stop("Could not import these helpers from the EXP script: ",
       paste(missing_helpers, collapse = ", "),
       " — the EXP script's structure changed; stopping rather than reimplementing.")
cat(sprintf("Imported %d helpers verbatim from %s:\n  %s\n",
            length(imported), basename(EXP_SCRIPT),
            paste(sort(imported), collapse = ", ")))

P_FLOOR <- 1 / (N_PERM + 1)   # permutation floor, as in EXP

# =============================================================================
# Inputs
# =============================================================================
cat("\nLoading committed main outputs (SPHARM direction / morphology, metadata)...\n")

SPHARM_direction  <- read_csv(
  here("analysis/data/derived_data/SPHARM_direction.csv"),  show_col_types = FALSE)
SPHARM_morphology <- read_csv(
  here("analysis/data/derived_data/SPHARM_morphology.csv"), show_col_types = FALSE)
SPHARM_morphology <- SPHARM_morphology %>%
  left_join(SPHARM_direction %>% select(ID, Typology), by = "ID")

metric_data <- read_xlsx(here("analysis/data/raw_data/SDG_core_metric.xlsx")) %>%
  select(ID, Layer, Core_type_Li_merged, Raw_mat) %>%
  mutate(Layer = as.factor(Layer), Raw_mat = as.factor(Raw_mat))
core_meta <- read_excel(here("analysis/data/raw_data/SDG_core_metric.xlsx")) %>%
  select(ID = ID, raw_material = Raw_mat, core_type = Core_type_Li_merged) %>%
  mutate(across(everything(), ~ str_trim(as.character(.))))

# =============================================================================
# Data preparation — SDG branch, rewritten to match SDG_cores_statistics.R exactly
# =============================================================================
build_blocks_sdg <- function(restrict_layers = RESTRICT_LAYERS) {
  dir_f   <- filter_spharm(SPHARM_direction,  POWER_COLS_DIR,   metric_data)
  morph_f <- filter_spharm(SPHARM_morphology, POWER_COLS_MORPH, metric_data)

  common   <- intersect(dir_f$ID, morph_f$ID)
  df_dir   <- dir_f   %>% filter(ID %in% common) %>% arrange(ID)
  df_morph <- morph_f %>% filter(ID %in% common) %>% arrange(ID)
  if (!identical(df_dir$ID, df_morph$ID))
    stop("SDG direction / morphology frames are not row-aligned after arrange(ID); ",
         "stopping rather than guessing an alignment.")

  arch_idx  <- !str_starts(df_dir$ID, "IM_") & !str_starts(df_dir$ID, "EXP")
  meta_arch <- tibble(ID = df_dir$ID[arch_idx],
                      Layer = as.character(df_dir$Layer[arch_idx])) %>%
    left_join(core_meta %>% select(ID, core_type), by = "ID") %>%
    filter(!core_type %in% EXCLUDE_CORE_TYPES | is.na(core_type))
  n_arch <- nrow(meta_arch)

  n_before_layer <- n_arch
  if (!is.null(restrict_layers))
    meta_arch <- meta_arch %>% filter(Layer %in% restrict_layers)
  n_after_layer <- nrow(meta_arch)

  counts_before <- sort(table(meta_arch$core_type, useNA = "no"))
  mc <- safe_filter_groups(meta_arch, "core_type")
  if (is.null(mc))
    stop("SDG: fewer than two core-type groups with n >= 3 after filtering.")
  dropped <- setdiff(names(counts_before), unique(mc$core_type))

  keep_idx <- match(mc$ID, df_dir$ID)
  list(ids   = df_dir$ID[keep_idx],
       morph = df_morph[keep_idx, POWER_COLS_MORPH],
       scar  = df_dir[keep_idx,   POWER_COLS_DIR],
       group = factor(mc$core_type),
       layer = mc$Layer,
       n_arch = n_arch, n_before_layer = n_before_layer,
       n_after_layer = n_after_layer,
       counts_before = counts_before, dropped = dropped)
}

blocks <- build_blocks_sdg()
n_spec  <- length(blocks$ids)
n_grp   <- nlevels(blocks$group)
n_pairs <- choose(n_grp, 2)

cat("\n", strrep("=", 70), "\n", sep = "")
cat("(0) SAMPLE / GROUP ALIGNMENT\n")
cat(strrep("=", 70), "\n", sep = "")
cat(sprintf("  Archaeological specimens after dropping %s : %d\n",
            paste(EXCLUDE_CORE_TYPES, collapse = " / "), blocks$n_arch))
cat(sprintf("  Layer restriction                          : %s (%d -> %d)\n",
            ifelse(is.null(RESTRICT_LAYERS), "NONE (published core-type recipe)",
                   paste(RESTRICT_LAYERS, collapse = "+")),
            blocks$n_before_layer, blocks$n_after_layer))
cat("  Core-type counts before the min_n = 3 filter:\n")
print(blocks$counts_before)
cat(sprintf("\n  Groups dropped by safe_filter_groups(min_n = 3): %s\n",
            ifelse(length(blocks$dropped) == 0, "none",
                   paste(blocks$dropped, collapse = ", "))))
cat("  (this exclusion is not declared in the manuscript — recorded here)\n")
cat("\n  Retained group sizes:\n"); print(table(blocks$group))
cat(sprintf("\n  FINAL: n = %d specimens | %d groups | %d pairwise comparisons | df = %d\n",
            n_spec, n_grp, n_pairs, n_grp - 1))
cat(sprintf("  Layer composition of the retained set: %s\n",
            paste(sprintf("L%s = %d", names(table(blocks$layer)),
                          as.integer(table(blocks$layer))), collapse = ", ")))
if (n_spec != REF$n || n_grp != REF$n_groups)
  cat(sprintf("\n  NOTE: n = %d / groups = %d differs from the reference (%d / %d).\n",
              n_spec, n_grp, REF$n, REF$n_groups))

# =============================================================================
# ILR blocks + MFA normalisation
# =============================================================================
parts_M  <- ilr_parts(blocks$morph)
parts_SP <- ilr_parts(blocks$scar)
Z_M  <- make_ilr(blocks$morph)
Z_SP <- make_ilr(blocks$scar)
rownames(Z_M) <- rownames(Z_SP) <- blocks$ids
colnames(Z_M)  <- paste0("M_ilr",  seq_len(ncol(Z_M)))
colnames(Z_SP) <- paste0("SP_ilr", seq_len(ncol(Z_SP)))
cat(sprintf("\nILR blocks: M-SPHARM %d parts -> %d coords | SP-SPHARM %d parts -> %d coords\n",
            length(parts_M), ncol(Z_M), length(parts_SP), ncol(Z_SP)))

s_M  <- mfa_s1(Z_M)
s_SP <- mfa_s1(Z_SP)
Z_comb <- cbind(Z_M / s_M, Z_SP / s_SP)
cat(sprintf("MFA normalisation: s1(M) = %.5f, s1(SP) = %.5f\n", s_M, s_SP))

I_M   <- block_inertia(Z_M  / s_M)
I_SP  <- block_inertia(Z_SP / s_SP)
w_mfa <- I_M / (I_M + I_SP)
cat(sprintf("Inertia after s1 normalisation: I_M = %.4f, I_SP = %.4f  =>  w_MFA = %.4f\n",
            I_M, I_SP, w_mfa))

# =============================================================================
# Three-model PERMANOVA (always on the FULL Z distance matrix, never on axes)
# =============================================================================
cat("\n", strrep("=", 70), "\n", sep = "")
cat("PERMANOVA: M-SPHARM alone / SP-SPHARM alone / MFA-combined (SDG)\n")
cat(strrep("=", 70), "\n", sep = "")

MODELS <- list("M-SPHARM" = Z_M, "SP-SPHARM" = Z_SP, "MFA-combined" = Z_comb)
model_res <- map(names(MODELS), function(nm) {
  cat(sprintf("  running %s (%d coords, %d pairs) ...\n",
              nm, ncol(MODELS[[nm]]), n_pairs))
  list(global   = permanova_global(MODELS[[nm]], blocks$group),
       pairwise = permanova_pairwise(MODELS[[nm]], blocks$group))
})
names(model_res) <- names(MODELS)

cat("  running PERMDISP on the combined distance matrix ...\n")
pd_comb <- permdisp(Z_comb, blocks$group)

comparison_tbl <- map_dfr(names(model_res), function(nm) {
  g <- model_res[[nm]]$global; pw <- model_res[[nm]]$pairwise
  tibble(model = nm, n = n_spec, n_coords = ncol(MODELS[[nm]]),
         R2 = g$R2, pseudo_F = g$F, p = g$p,
         n_sig_pairs = sum(pw$significant_holm), n_pairs = n_pairs,
         sig_pairs = paste(sort(pw$comparison[pw$significant_holm]), collapse = "; "),
         permdisp_F = if (nm == "MFA-combined") pd_comb$F else NA_real_,
         permdisp_p = if (nm == "MFA-combined") pd_comb$p else NA_real_)
})

cat("\n--- Three-model comparison (SDG) ---\n")
print(comparison_tbl %>%
        transmute(model, n, n_coords, R2 = round(R2, 5),
                  pseudo_F = round(pseudo_F, 5), p,
                  resolved = sprintf("%d/%d", n_sig_pairs, n_pairs)) %>%
        as.data.frame())

cat("\n--- PERMDISP (combined distance matrix) ---\n")
print(pd_comb$test)

# ---- anchor check ----------------------------------------------------------
cat("\n", strrep("=", 70), "\n", sep = "")
cat("ANCHOR CHECK vs committed SDG core-type PERMANOVA (L3_permanova.csv)\n")
cat(strrep("=", 70), "\n", sep = "")
anchor_ok <- TRUE
check <- function(name, got, expect, tol, deterministic = TRUE) {
  ok <- is.finite(got) && abs(got - expect) <= tol
  if (!ok && deterministic) anchor_ok <<- FALSE
  cat(sprintf("  %-26s got %10.5f  expected %10.5f  %s\n", name, got, expect,
              ifelse(ok, "OK", ifelse(deterministic, "<-- CHECK",
                                      "<-- CHECK (permutation)"))))
}
mg <- model_res[["M-SPHARM"]]$global; sg <- model_res[["SP-SPHARM"]]$global
check("n specimens",        n_spec,  REF$n,        0)
check("n groups",           n_grp,   REF$n_groups, 0)
check("M-SPHARM R2",        mg$R2,   REF$M_R2,     REF_TOL)
check("M-SPHARM pseudo-F",  mg$F,    REF$M_F,      REF_TOL)
check("M-SPHARM p",         mg$p,    REF$M_p,      0.01, deterministic = FALSE)
check("SP-SPHARM R2",       sg$R2,   REF$SP_R2,    REF_TOL)
check("SP-SPHARM pseudo-F", sg$F,    REF$SP_F,     REF_TOL)
check("SP-SPHARM p",        sg$p,    REF$SP_p,     0.01, deterministic = FALSE)
if (!anchor_ok)
  stop("ANCHOR CHECK FAILED — the recomputed SDG single-block PERMANOVA does not ",
       "match the committed Table 2 values. Stopping before the diagnostics: the ",
       "sample definition must be resolved first.")
cat("\n  => Single-block models reproduce the committed SDG Table 2 core-type PERMANOVA.\n")
cat("  (Permutation p-values jitter by ~+/-0.005 and are informational.)\n")
cat(sprintf("\n  Block ordering: M-SPHARM R2 = %.5f %s SP-SPHARM R2 = %.5f  (EXP: %.5f < %.5f)\n",
            mg$R2, ifelse(mg$R2 > sg$R2, ">", "<"), sg$R2, EXP_REF$M_R2, EXP_REF$SP_R2))

# =============================================================================
# MFA biplot
# =============================================================================
pca_comb <- stats::prcomp(Z_comb, center = TRUE, scale. = FALSE)
inertia  <- pca_comb$sdev^2 / sum(pca_comb$sdev^2) * 100
cat(sprintf("\nMFA global PCA: axis 1 = %.1f%%, axis 2 = %.1f%% of inertia\n",
            inertia[1], inertia[2]))

V_M  <- ilrBase(D = length(parts_M))
V_SP <- ilrBase(D = length(parts_SP))
chk_clr <- function(power_df, parts, Z, V) {
  X   <- replace_zeros(as.matrix(power_df)[, parts, drop = FALSE])
  ref <- unclass(clr(X))
  max(abs(Z %*% t(V) - ref))
}
clr_err <- max(chk_clr(blocks$morph, parts_M,  Z_M,  V_M),
               chk_clr(blocks$scar,  parts_SP, Z_SP, V_SP))
cat(sprintf("  ILR basis check (max |Z %%*%% t(V) - clr|) = %.3e\n", clr_err))
if (clr_err > 1e-8)
  stop("ilrBase() does not reproduce the basis used by ilr(); CLR arrows would be wrong.")

idx_M  <- seq_len(ncol(Z_M))
idx_SP <- ncol(Z_M) + seq_len(ncol(Z_SP))
clr_M  <- V_M  %*% pca_comb$rotation[idx_M,  1:2]
clr_SP <- V_SP %*% pca_comb$rotation[idx_SP, 1:2]

deg_label <- function(cols) paste0("l", str_extract(cols, "[0-9]+$"))
load_df <- bind_rows(
  as_tibble(clr_M,  .name_repair = ~ c("Axis1", "Axis2")) %>%
    mutate(block = "M-SPHARM",  degree = deg_label(parts_M)),
  as_tibble(clr_SP, .name_repair = ~ c("Axis1", "Axis2")) %>%
    mutate(block = "SP-SPHARM", degree = deg_label(parts_SP))
) %>%
  mutate(block = factor(block, levels = c("M-SPHARM", "SP-SPHARM")),
         label = paste(block, degree))

scores_df <- as_tibble(pca_comb$x[, 1:2], .name_repair = ~ c("Axis1", "Axis2")) %>%
  mutate(ID = blocks$ids, group = blocks$group)

arrow_scale <- 0.80 * max(abs(as.matrix(scores_df[, c("Axis1", "Axis2")]))) /
  max(sqrt(load_df$Axis1^2 + load_df$Axis2^2))
load_df <- load_df %>% mutate(x = Axis1 * arrow_scale, y = Axis2 * arrow_scale)

# SDG has 6 core types; extend the 5 base colours (EXP script's non-EXP branch).
group_pal <- setNames(colorRampPalette(unname(TYPOLOGY_COLORS))(n_grp),
                      levels(blocks$group))
pretty_lab <- function(x) str_replace_all(x, "_", " ")

p_biplot <- ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70", linewidth = 0.25) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey70", linewidth = 0.25) +
  stat_ellipse(data = scores_df, aes(Axis1, Axis2, color = group),
               type = "norm", level = 0.95, linewidth = 0.4, alpha = 0.8) +
  geom_point(data = scores_df, aes(Axis1, Axis2, fill = group, color = group),
             shape = 21, size = 2, stroke = 0.4, alpha = 0.90) +
  geom_segment(data = load_df,
               aes(x = 0, y = 0, xend = x, yend = y, linetype = block),
               arrow = grid::arrow(length = grid::unit(0.16, "cm"), type = "closed"),
               color = "grey25", linewidth = 0.4) +
  ggrepel::geom_text_repel(data = load_df, aes(x = x, y = y, label = label),
                           size = 2.2, color = "grey20", min.segment.length = 0.2,
                           segment.size = 0.2, segment.color = "grey60",
                           box.padding = 0.2, max.overlaps = 30) +
  scale_color_manual(values = group_pal, name = NULL, labels = pretty_lab) +
  scale_fill_manual(values  = group_pal, name = NULL, labels = pretty_lab) +
  scale_linetype_manual(values = c("M-SPHARM" = "solid", "SP-SPHARM" = "22"),
                        name = "CLR loading") +
  labs(x = sprintf("MFA axis 1 (%.1f%%)", inertia[1]),
       y = sprintf("MFA axis 2 (%.1f%%)", inertia[2])) +
  coord_fixed() +
  theme_bw(base_size = 9) +
  guides(color = guide_legend(order = 1,
                              override.aes = list(shape = 21, size = 2, linetype = 0)),
         fill  = "none",
         linetype = guide_legend(order = 2, override.aes = list(color = "grey25"))) +
  theme(panel.grid.minor = element_blank(),
        legend.key.size  = grid::unit(0.34, "cm"),
        legend.text      = element_text(size = 7),
        legend.title     = element_text(size = 7.5))

ggsave(file.path(FIG_DIR, "fig_S_joint_mfa_biplot_SDG.png"), p_biplot,
       width = 7, height = 7, dpi = 300)
cat("Wrote figures/fig_S_joint_mfa_biplot_SDG.png\n")

# =============================================================================
# D1. Raw vs Holm-corrected pairwise p
# =============================================================================
cat("\n", strrep("=", 70), "\n", sep = "")
cat("D1. RAW vs HOLM-CORRECTED PAIRWISE p (SDG)\n")
cat(strrep("=", 70), "\n", sep = "")

pairwise_all <- map_dfr(names(model_res),
                        ~ model_res[[.x]]$pairwise %>% mutate(model = .x))

pairwise_raw_holm <- bind_rows(
  pairwise_all %>% transmute(comparison, model, p_type = "raw",  value = p),
  pairwise_all %>% transmute(comparison, model, p_type = "holm", value = p_holm)
) %>%
  pivot_wider(names_from = model, values_from = value) %>%
  arrange(p_type, comparison)

d1_tbl <- pairwise_all %>%
  select(comparison, model, p, p_holm) %>%
  pivot_longer(c(p, p_holm), names_to = "stat", values_to = "value") %>%
  unite("key", model, stat) %>%
  pivot_wider(names_from = key, values_from = value) %>%
  transmute(comparison,
            p_M         = `M-SPHARM_p`,
            p_SP        = `SP-SPHARM_p`,
            p_comb      = `MFA-combined_p`,
            p_holm_comb = `MFA-combined_p_holm`,
            improve_x   = p_SP / p_comb,
            censored    = p_SP <= P_FLOOR | p_comb <= P_FLOOR) %>%
  arrange(p_M)

cat(sprintf("  (raw p floor = %.4f at %d permutations; ratios marked * are censored)\n\n",
            P_FLOOR, N_PERM))
cat(sprintf("  %-38s %8s %8s %8s %12s %10s\n",
            "comparison (sorted by p_M)", "p_M", "p_SP", "p_comb", "p_holm_comb", "SP/comb"))
for (i in seq_len(nrow(d1_tbl))) {
  r <- d1_tbl[i, ]
  cat(sprintf("  %-38s %8.4f %8.4f %8.4f %12.4f %9.2f%s\n",
              r$comparison, r$p_M, r$p_SP, r$p_comb, r$p_holm_comb,
              r$improve_x, ifelse(r$censored, "*", " ")))
}

d1_free <- d1_tbl %>% filter(!censored)
rho_improve <- if (nrow(d1_free) >= 3)
  suppressWarnings(stats::cor(d1_free$p_M, d1_free$improve_x, method = "spearman")) else NA_real_
rho_p <- if (nrow(d1_free) >= 3)
  suppressWarnings(stats::cor.test(d1_free$p_M, d1_free$improve_x,
                                   method = "spearman", exact = FALSE)$p.value) else NA_real_
d1_monotone <- is.finite(rho_improve) && rho_improve <= -0.5
cat(sprintf("\n  Spearman(p_M, improvement factor) on the %d uncensored pair(s) = %s (p = %s)\n",
            nrow(d1_free),
            ifelse(is.finite(rho_improve), sprintf("%.3f", rho_improve), "NA"),
            ifelse(is.finite(rho_p), sprintf("%.4f", rho_p), "NA")))
cat(sprintf("  distinct p_M values: %d / %d (ties inflate the Spearman standard error)\n",
            dplyr::n_distinct(d1_free$p_M), nrow(d1_free)))
cat(sprintf("  => improvement %s monotone in morphological signal (EXP: %.3f)\n",
            ifelse(d1_monotone, "IS", "is NOT"), EXP_REF$rho_improve))

write_csv(pairwise_raw_holm, file.path(OUT_DIR, "pairwise_raw_and_holm_SDG.csv"))
cat("\n  Wrote pairwise_raw_and_holm_SDG.csv\n")

# =============================================================================
# D2. Block structure of the leading axes
# =============================================================================
cat("\n", strrep("=", 70), "\n", sep = "")
cat("D2. BLOCK STRUCTURE OF THE LEADING MFA AXES (SDG)\n")
cat(strrep("=", 70), "\n", sep = "")

n_axes_diag <- min(4, ncol(pca_comb$rotation))
axis_block_tbl <- map_dfr(seq_len(n_axes_diag), function(k) {
  sM <- sum(pca_comb$rotation[idx_M,  k]^2)
  sS <- sum(pca_comb$rotation[idx_SP, k]^2)
  tibble(axis = k, inertia_pct = inertia[k], share_M = sM, share_SP = sS,
         dominant_block = dplyr::case_when(sM >= DOMINANCE_THR ~ "M-SPHARM",
                                           sS >= DOMINANCE_THR ~ "SP-SPHARM",
                                           TRUE                ~ "mixed"))
})
print(axis_block_tbl %>%
        mutate(inertia_pct = round(inertia_pct, 2),
               share_M = round(share_M, 3), share_SP = round(share_SP, 3)) %>%
        as.data.frame())

d2_split <- axis_block_tbl$dominant_block[1] != "mixed" &&
  axis_block_tbl$dominant_block[2] != "mixed" &&
  axis_block_tbl$dominant_block[1] != axis_block_tbl$dominant_block[2]
cat(sprintf("\n  => The first two axes %s a one-axis-per-block split (threshold: share >= %.2f).\n",
            ifelse(d2_split, "DO show", "do NOT show"), DOMINANCE_THR))

rel_gap <- abs(inertia[1] - inertia[2]) / inertia[1]
d2_degenerate <- rel_gap < DEGENERACY_THR
if (d2_degenerate) {
  cat(sprintf("\n  CAUTION: axes 1 and 2 are near-degenerate (%.1f%% vs %.1f%%, relative gap %.1f%%).\n",
              inertia[1], inertia[2], rel_gap * 100))
  cat("  The plane is stable but the individual axis directions are not, so the exact\n")
  cat("  bearing of any single loading arrow must not be over-interpreted.\n")
}

write_csv(axis_block_tbl, file.path(OUT_DIR, "axis_block_contributions_SDG.csv"))
cat("\n  Wrote axis_block_contributions_SDG.csv\n")

# =============================================================================
# D3. Weight scan — coarse pass, then an AUTOMATICALLY located dense window
# =============================================================================
cat("\n", strrep("=", 70), "\n", sep = "")
cat("D3. WEIGHT SCAN (SDG): coarse pass, then adaptive densification\n")
cat(strrep("=", 70), "\n", sep = "")

f_M  <- sqrt(block_inertia(Z_M))
f_SP <- sqrt(block_inertia(Z_SP))
ZM_n <- Z_M  / f_M
ZS_n <- Z_SP / f_SP

scan_at <- function(w) {
  Zw <- cbind(sqrt(w) * ZM_n, sqrt(1 - w) * ZS_n)
  g  <- permanova_global(Zw, blocks$group)
  pw <- permanova_pairwise(Zw, blocks$group)
  n05h <- sum(pw$p_holm < 0.05); n01h <- sum(pw$p_holm < 0.01)
  n05r <- sum(pw$p < 0.05)
  cat(sprintf("  w = %.2f  R2 = %.4f  holm05 = %2d  holm01 = %2d  raw05 = %2d  med(-log10 p) = %.3f\n",
              w, g$R2, n05h, n01h, n05r, stats::median(-log10(pw$p))))
  bind_cols(
    tibble(w = w, R2 = g$R2, pseudo_F = g$F, p = g$p, n_pairs = n_pairs,
           n_sig_pairs_holm_05 = n05h, n_sig_pairs_holm_01 = n01h,
           n_sig_pairs_raw_05  = n05r,
           median_neglog10_p_raw = stats::median(-log10(pw$p)),
           n_pairs_at_p_floor  = sum(pw$p <= P_FLOOR)),
    as_tibble_row(setNames(pw$p,
                           paste0("p_raw_", str_replace_all(pw$comparison, "-", "_")))))
}

cat("\n  -- coarse pass --\n")
coarse_res <- map_dfr(W_GRID_COARSE, scan_at)

# Locate the densification window from the coarse result.
step_change <- diff(coarse_res$n_sig_pairs_holm_05)
if (length(step_change) > 0 && min(step_change) < 0) {
  i_drop  <- which.min(step_change)             # steepest fall
  drop_lo <- coarse_res$w[i_drop]; drop_hi <- coarse_res$w[i_drop + 1]
  dense_lo <- max(0, drop_lo - DENSE_HALFWIDTH)
  dense_hi <- min(1, drop_hi + DENSE_HALFWIDTH)
  dense_reason <- sprintf("steepest Holm(0.05) drop %d -> %d between w = %.2f and %.2f",
                          coarse_res$n_sig_pairs_holm_05[i_drop],
                          coarse_res$n_sig_pairs_holm_05[i_drop + 1], drop_lo, drop_hi)
} else {
  i_pk <- which.max(coarse_res$median_neglog10_p_raw)
  dense_lo <- max(0, coarse_res$w[i_pk] - DENSE_HALFWIDTH)
  dense_hi <- min(1, coarse_res$w[i_pk] + DENSE_HALFWIDTH)
  dense_reason <- sprintf(
    "no Holm(0.05) drop anywhere (counts flat at %d); densifying around the continuous-curve peak at w = %.2f",
    coarse_res$n_sig_pairs_holm_05[1], coarse_res$w[i_pk])
}
cat(sprintf("\n  Densification window: [%.2f, %.2f] at step %.2f\n  reason: %s\n",
            dense_lo, dense_hi, DENSE_STEP, dense_reason))

W_extra <- setdiff(round(seq(dense_lo, dense_hi, by = DENSE_STEP), 10),
                   round(coarse_res$w, 10))
cat(sprintf("  -- dense pass (%d extra points) --\n", length(W_extra)))
dense_res <- if (length(W_extra) > 0)
  bind_rows(coarse_res, map_dfr(W_extra, scan_at)) %>% arrange(w) else
    coarse_res %>% arrange(w)

write_csv(dense_res, file.path(OUT_DIR, "weight_scan_dense_SDG.csv"))
cat(sprintf("\n  Wrote weight_scan_dense_SDG.csv (%d rows)\n", nrow(dense_res)))

# ---- read the curve ---------------------------------------------------------
d3_best_i <- which.max(dense_res$median_neglog10_p_raw)
d3_best_w <- dense_res$w[d3_best_i]
d3_best_v <- dense_res$median_neglog10_p_raw[d3_best_i]
d3_flat   <- diff(range(dense_res$median_neglog10_p_raw)) < 1e-9   # guard which.max

lvl0 <- dense_res$n_sig_pairs_holm_05[1]
d3_plateau_end <- if (lvl0 == 0) NA_real_ else
  max(dense_res$w[dplyr::cumall(dense_res$n_sig_pairs_holm_05 == lvl0)])

zone <- dense_res %>% filter(w >= dense_lo, w <= dense_hi)
zone_gap   <- zone$n_sig_pairs_raw_05 - zone$n_sig_pairs_holm_05
d3_max_gap <- if (nrow(zone) > 0) max(zone_gap) else NA_integer_
d3_gap_w   <- if (nrow(zone) > 0) zone$w[which.max(zone_gap)] else NA_real_

cat(sprintf("\n  Continuous curve peaks at w = %.2f (median -log10 raw p = %.3f)%s;\n",
            d3_best_w, d3_best_v,
            ifelse(d3_flat, "  [curve is flat — peak location not meaningful]", "")))
cat(sprintf("    value at w = 0 (pure SP) = %.3f, at w = 1 (pure M) = %.3f, at w_MFA = %.2f it is %.3f.\n",
            dense_res$median_neglog10_p_raw[dense_res$w == 0],
            dense_res$median_neglog10_p_raw[dense_res$w == 1], w_mfa,
            dense_res$median_neglog10_p_raw[which.min(abs(dense_res$w - w_mfa))]))
cat(sprintf("  Holm(0.05) resolved pairs: min %d, max %d (of %d).\n",
            min(dense_res$n_sig_pairs_holm_05), max(dense_res$n_sig_pairs_holm_05), n_pairs))
if (is.na(d3_plateau_end)) {
  cat("  No plateau: 0 pairs are resolved at w = 0, so a 'plateau right endpoint' is undefined.\n")
} else {
  cat(sprintf("  Plateau at %d/%d pairs runs from w = 0 to w = %.2f.\n",
              lvl0, n_pairs, d3_plateau_end))
}
cat(sprintf("  Largest raw-minus-Holm count gap inside the dense window: %d (at w = %.2f).\n",
            d3_max_gap, d3_gap_w))

# ---- figures ----------------------------------------------------------------
CRIT_COLORS <- c("Holm, alpha = 0.05" = "#4A6E8A",
                 "Holm, alpha = 0.01" = "#788C4A",
                 "Raw, alpha = 0.05"  = "#BA8530")
dense_long <- dense_res %>%
  select(w,
         `Holm, alpha = 0.05` = n_sig_pairs_holm_05,
         `Holm, alpha = 0.01` = n_sig_pairs_holm_01,
         `Raw, alpha = 0.05`  = n_sig_pairs_raw_05) %>%
  pivot_longer(-w, names_to = "criterion", values_to = "n_sig") %>%
  mutate(criterion = factor(criterion, levels = names(CRIT_COLORS)))

y_top <- max(n_pairs, max(dense_long$n_sig))
p_dense <- ggplot(dense_long, aes(w, n_sig, color = criterion)) +
  geom_vline(xintercept = w_mfa, linetype = "dashed",
             color = "#802520", linewidth = 0.4) +
  annotate("text", x = w_mfa, y = y_top, hjust = -0.08, vjust = 1.2, size = 2.6,
           color = "#802520", label = sprintf("MFA (s1) setting: w = %.2f", w_mfa)) +
  geom_step(linewidth = 0.7) +
  scale_color_manual(values = CRIT_COLORS, name = NULL) +
  scale_x_continuous(breaks = seq(0, 1, 0.1), limits = c(0, 1)) +
  scale_y_continuous(breaks = scales::breaks_width(1), limits = c(0, y_top)) +
  labs(x = "Morphology share of squared distance (w)",
       y = sprintf("Core-type pairs resolved (of %d)", n_pairs)) +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom",
        legend.key.size = grid::unit(0.34, "cm"))

ggsave(file.path(FIG_DIR, "fig_S_joint_mfa_weight_scan_dense_SDG.png"), p_dense,
       width = 7, height = 4.8, dpi = 300)
cat("  Wrote figures/fig_S_joint_mfa_weight_scan_dense_SDG.png\n")

p_cont <- ggplot(dense_res, aes(w, median_neglog10_p_raw)) +
  geom_vline(xintercept = w_mfa, linetype = "dashed",
             color = "#802520", linewidth = 0.4) +
  annotate("text", x = w_mfa, y = max(dense_res$median_neglog10_p_raw),
           hjust = -0.08, vjust = 1.2, size = 2.6, color = "#802520",
           label = sprintf("MFA (s1) setting: w = %.2f", w_mfa)) +
  geom_line(linewidth = 0.7, color = "#4A6E8A") +
  geom_point(size = 1.2, color = "#4A6E8A") +
  scale_x_continuous(breaks = seq(0, 1, 0.1), limits = c(0, 1)) +
  labs(x = "Morphology share of squared distance (w)",
       y = "Median -log10(raw p) over all core-type pairs") +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank())

ggsave(file.path(FIG_DIR, "fig_S_joint_mfa_weight_resolution_continuous_SDG.png"),
       p_cont, width = 7, height = 4.5, dpi = 300)
cat("  Wrote figures/fig_S_joint_mfa_weight_resolution_continuous_SDG.png\n")

# =============================================================================
# D4. PERMDISP group dispersions
# =============================================================================
cat("\n", strrep("=", 70), "\n", sep = "")
cat("D4. PERMDISP GROUP DISPERSIONS (SDG: M alone / SP alone / combined)\n")
cat(strrep("=", 70), "\n", sep = "")

group_n <- table(blocks$group)
permdisp_by_group <- map_dfr(names(MODELS), function(nm) {
  pdm <- permdisp(MODELS[[nm]], blocks$group)
  md  <- tapply(pdm$disp$distances, pdm$disp$group, mean)
  tibble(model = nm, group = names(md), n = as.integer(group_n[names(md)]),
         mean_dist_to_centroid = as.numeric(md),
         permdisp_F = pdm$F, permdisp_p = pdm$p) %>%
    mutate(dispersion_ratio_vs_min = mean_dist_to_centroid / min(mean_dist_to_centroid))
})

for (nm in names(MODELS)) {
  sub <- permdisp_by_group %>% filter(model == nm) %>%
    arrange(desc(dispersion_ratio_vs_min))
  cat(sprintf("\n  --- %s (PERMDISP F = %.4f, p = %.4f) ---\n",
              nm, sub$permdisp_F[1], sub$permdisp_p[1]))
  print(sub %>% transmute(group, n,
                          mean_dist = round(mean_dist_to_centroid, 4),
                          ratio_vs_min = round(dispersion_ratio_vs_min, 3)) %>%
          as.data.frame())
}

d4_ratio <- permdisp_by_group %>% group_by(model) %>%
  summarise(max_min_ratio = max(mean_dist_to_centroid) / min(mean_dist_to_centroid),
            most_dispersed = group[which.max(mean_dist_to_centroid)],
            permdisp_p = first(permdisp_p), .groups = "drop")
cat("\n  Scale-invariant max/min dispersion ratio per model:\n")
print(d4_ratio %>% mutate(max_min_ratio = round(max_min_ratio, 3)) %>% as.data.frame())

r_M    <- d4_ratio$max_min_ratio[d4_ratio$model == "M-SPHARM"]
r_sp   <- d4_ratio$max_min_ratio[d4_ratio$model == "SP-SPHARM"]
r_comb <- d4_ratio$max_min_ratio[d4_ratio$model == "MFA-combined"]

tol <- 0.02
d4_bracketed <- is.finite(r_M) && is.finite(r_sp) && is.finite(r_comb) &&
  r_comb >= min(r_M, r_sp) - tol && r_comb <= max(r_M, r_sp) + tol
d4_substantive <- is.finite(r_comb) && r_comb < min(r_M, r_sp) - tol

disp_wide <- permdisp_by_group %>%
  select(model, group, mean_dist_to_centroid) %>%
  pivot_wider(names_from = model, values_from = mean_dist_to_centroid)
rho_disp <- suppressWarnings(stats::cor(disp_wide$`M-SPHARM`, disp_wide$`SP-SPHARM`,
                                        method = "spearman"))
top_M  <- disp_wide$group[which.max(disp_wide$`M-SPHARM`)]
bot_M  <- disp_wide$group[which.min(disp_wide$`M-SPHARM`)]
top_SP <- disp_wide$group[which.max(disp_wide$`SP-SPHARM`)]
bot_SP <- disp_wide$group[which.min(disp_wide$`SP-SPHARM`)]
extremes_swap <- (top_M == bot_SP) || (top_SP == bot_M)

cat(sprintf("\n  M alone = %.3f, SP alone = %.3f, combined = %.3f (combined/SP = %.3f).\n",
            r_M, r_sp, r_comb, r_comb / r_sp))
cat(sprintf("  Spearman between the two blocks' group dispersion orderings = %.3f (%s)\n",
            rho_disp,
            dplyr::case_when(!is.finite(rho_disp) ~ "undefined",
                             rho_disp <= -0.5     ~ "reversed orderings",
                             rho_disp <   0.3     ~ "essentially unrelated orderings",
                             TRUE                 ~ "similar orderings")))
cat(sprintf("  Extremes: M most/least = %s / %s; SP most/least = %s / %s%s\n",
            top_M, bot_M, top_SP, bot_SP,
            ifelse(extremes_swap, "  -> an extreme group swaps ends", "  -> extremes do not swap")))
if (d4_substantive) {
  cat("  => SUBSTANTIVE: the combined ratio falls BELOW both single-block ratios, so the\n")
  cat("     blocks partly cancel each other's dispersion imbalance.\n")
} else if (d4_bracketed) {
  cat("  => DILUTION: the combined ratio sits BETWEEN the two single-block ratios, which is\n")
  cat("     what appending a more evenly dispersed block must produce. No group's dispersion\n")
  cat("     has actually fallen; the same imbalance is spread over more coordinates.\n")
} else {
  cat("  => The combined ratio EXCEEDS both block ratios; inspect before interpreting.\n")
}

write_csv(permdisp_by_group %>%
            select(model, group, n, mean_dist_to_centroid, dispersion_ratio_vs_min),
          file.path(OUT_DIR, "permdisp_by_group_SDG.csv"))
cat("\n  Wrote permdisp_by_group_SDG.csv\n")

# =============================================================================
# PERMANOVA comparison CSV (same tidy layout as EXP)
# =============================================================================
permanova_csv <- bind_rows(
  comparison_tbl %>%
    transmute(model, scope = "global", comparison = "core_type (all groups)",
              n, n_coords, R2, pseudo_F, p, p_holm = NA_real_,
              significant_holm = p < 0.05,
              n_sig_pairs, n_pairs, sig_pairs, permdisp_F, permdisp_p),
  map_dfr(names(model_res), function(nm)
    model_res[[nm]]$pairwise %>%
      transmute(model = nm, scope = "pairwise", comparison,
                n, n_coords = ncol(MODELS[[nm]]), R2, pseudo_F, p, p_holm,
                significant_holm,
                n_sig_pairs = NA_integer_, n_pairs = NA_integer_,
                sig_pairs = NA_character_,
                permdisp_F = NA_real_, permdisp_p = NA_real_))
)
write_csv(permanova_csv, file.path(OUT_DIR, "permanova_comparison_SDG.csv"))
cat("Wrote permanova_comparison_SDG.csv\n")

# =============================================================================
# Final report
# =============================================================================
cat("\n", strrep("=", 70), "\n", sep = "")
cat("SDG DIAGNOSTIC CONCLUSIONS\n")
cat(strrep("=", 70), "\n", sep = "")
cat(sprintf("(0) n = %d | %d groups | %d pairs | dropped by min_n=3: %s | anchors: PASS\n",
            n_spec, n_grp, n_pairs,
            ifelse(length(blocks$dropped) == 0, "none", paste(blocks$dropped, collapse = ", "))))
cat(sprintf("(1) Spearman(p_M, improvement factor) = %s on %d uncensored pairs -> %s monotone.\n",
            ifelse(is.finite(rho_improve), sprintf("%.3f", rho_improve), "NA"),
            nrow(d1_free), ifelse(d1_monotone, "IS", "is NOT")))
cat(sprintf("(2) Leading axes %s one-axis-per-block (axis 1 %s %.0f%%, axis 2 %s %.0f%%)%s.\n",
            ifelse(d2_split, "DO show", "do NOT show"),
            axis_block_tbl$dominant_block[1],
            100 * max(axis_block_tbl$share_M[1], axis_block_tbl$share_SP[1]),
            axis_block_tbl$dominant_block[2],
            100 * max(axis_block_tbl$share_M[2], axis_block_tbl$share_SP[2]),
            ifelse(d2_degenerate, "; axes near-degenerate, arrow bearings not interpretable", "")))
cat(sprintf("(3) Resolved pairs range %d-%d of %d; continuous curve peaks at w = %.2f (%.3f).\n",
            min(dense_res$n_sig_pairs_holm_05), max(dense_res$n_sig_pairs_holm_05),
            n_pairs, d3_best_w, d3_best_v))
cat(sprintf("(4) PERMDISP max/min ratio: %.3f (M), %.3f (SP), %.3f (combined) -> %s.\n",
            r_M, r_sp, r_comb,
            ifelse(d4_substantive, "cancellation",
                   ifelse(d4_bracketed, "dilution", "inspect"))))

cat("\n", strrep("-", 70), "\n", sep = "")
cat("(5) PRE-REGISTERED PREDICTIONS (mechanism test on a reversed-ordering assemblage)\n")
cat(strrep("-", 70), "\n", sep = "")
pa_ok <- is.finite(d3_best_w) && !d3_flat && d3_best_w > EXP_REF$best_w
cat(sprintf("  (a) best w shifts toward morphology: SDG %.2f vs EXP %.2f -> %s\n",
            d3_best_w, EXP_REF$best_w,
            ifelse(pa_ok, "MATCHES PREDICTION", "DOES NOT MATCH")))
if (is.na(d3_plateau_end)) {
  cat(sprintf("  (b) plateau right endpoint > EXP %.2f -> NOT EVALUABLE (no plateau: %d/%d resolved at w = 0)\n",
              EXP_REF$plateau_end, lvl0, n_pairs))
  pb_ok <- NA
} else {
  pb_ok <- d3_plateau_end > EXP_REF$plateau_end
  cat(sprintf("  (b) plateau right endpoint: SDG %.2f vs EXP %.2f -> %s\n",
              d3_plateau_end, EXP_REF$plateau_end,
              ifelse(pb_ok, "MATCHES PREDICTION", "DOES NOT MATCH")))
}
pc_ok <- d1_monotone
cat(sprintf("  (c) Spearman(p_M, improvement) still clearly negative: SDG %s vs EXP %.3f -> %s\n",
            ifelse(is.finite(rho_improve), sprintf("%.3f", rho_improve), "NA"),
            EXP_REF$rho_improve, ifelse(pc_ok, "MATCHES PREDICTION", "DOES NOT MATCH")))

cat("\n", strrep("-", 70), "\n", sep = "")
cat("(6) EXP vs SDG COMPARISON  (EXP values hard-coded from its committed outputs)\n")
cat(strrep("-", 70), "\n", sep = "")
cmp <- tibble(
  quantity = c("n", "groups", "pairs", "M-SPHARM R2", "SP-SPHARM R2",
               "stronger block", "best w (continuous peak)",
               "plateau right endpoint", "Spearman(p_M, improvement)"),
  EXP = c(sprintf("%d", EXP_REF$n), sprintf("%d", EXP_REF$n_groups),
          sprintf("%d", EXP_REF$n_pairs), sprintf("%.5f", EXP_REF$M_R2),
          sprintf("%.5f", EXP_REF$SP_R2), "SP-SPHARM",
          sprintf("%.2f", EXP_REF$best_w), sprintf("%.2f", EXP_REF$plateau_end),
          sprintf("%.3f", EXP_REF$rho_improve)),
  SDG = c(sprintf("%d", n_spec), sprintf("%d", n_grp), sprintf("%d", n_pairs),
          sprintf("%.5f", mg$R2), sprintf("%.5f", sg$R2),
          ifelse(mg$R2 > sg$R2, "M-SPHARM", "SP-SPHARM"),
          sprintf("%.2f", d3_best_w),
          ifelse(is.na(d3_plateau_end), "none (0 resolved)",
                 sprintf("%.2f", d3_plateau_end)),
          ifelse(is.finite(rho_improve), sprintf("%.3f", rho_improve), "NA")))
print(as.data.frame(cmp), row.names = FALSE)

cat("\nDone.\n\n")
print(sessionInfo())
