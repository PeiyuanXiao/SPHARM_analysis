# sdi_rawmat_interaction.R

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

SDI_COL  <- "SDI"
SDI_XLSX <- here("analysis/data/raw_data/SDG_core_metric.xlsx")

OUT_DIR <- here("analysis/robustness/sdi_rawmat_interaction")
FIG_DIR <- file.path(OUT_DIR, "figures")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

REF <- list(M_R2  = 0.18770190528243877, M_F  = 2.0334613330095146,
            SP_R2 = 0.16785198589124128, SP_F = 1.7750417603590793,
            n = 50L, n_groups = 6L)
REF_TOL <- 1e-3

RAWMAT_COLORS <- c("chert" = "#4A6E8A", "sandstone" = "#BA8530")

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
# SAMPLE + SDI  (identical construction to sdi_gradient.R:246-300)
# =============================================================================
blocks <- build_blocks_sdg()
if (length(blocks$ids) != REF$n || nlevels(blocks$group) != REF$n_groups)
  stop(sprintf(paste0("Imported sample definition gives n = %d / %d groups, not ",
                      "the published %d / %d. Stopping."),
               length(blocks$ids), nlevels(blocks$group), REF$n, REF$n_groups))

sdi_xl <- read_excel(SDI_XLSX)
if (!SDI_COL %in% names(sdi_xl))
  stop("Column '", SDI_COL, "' not found in ", basename(SDI_XLSX), ".")
sdi_tbl <- tibble(ID  = str_trim(as.character(sdi_xl$ID)),
                  SDI = suppressWarnings(as.numeric(sdi_xl[[SDI_COL]])))
if (anyDuplicated(sdi_tbl$ID))
  stop("Duplicate IDs in ", basename(SDI_XLSX), ".")

meta_all <- tibble(ID = blocks$ids, core_type = blocks$group, layer = blocks$layer) %>%
  left_join(core_meta %>% select(ID, raw_material), by = "ID") %>%
  left_join(sdi_tbl, by = "ID")
if (nrow(meta_all) != length(blocks$ids))
  stop("Metadata join changed the row count — non-unique keys; stopping.")

keep0 <- which(is.finite(meta_all$SDI))
meta0 <- meta_all[keep0, ] %>% mutate(core_type = droplevels(core_type),
                                      Raw_mat   = factor(raw_material))
if (any(meta0$SDI <= 0))
  stop("SDI has non-positive values; log(SDI) versions B/D would be undefined. ",
       "Stopping rather than shifting the variable silently.")

# ---- build the three distance matrices for an arbitrary row subset ----------
build_spaces <- function(row_idx_in_blocks) {
  morph_sub <- blocks$morph[row_idx_in_blocks, , drop = FALSE]
  scar_sub  <- blocks$scar[ row_idx_in_blocks, , drop = FALSE]
  Z_M  <- make_ilr(morph_sub)
  Z_SP <- make_ilr(scar_sub)
  colnames(Z_M)  <- paste0("M_ilr",  seq_len(ncol(Z_M)))
  colnames(Z_SP) <- paste0("SP_ilr", seq_len(ncol(Z_SP)))
  s_M  <- mfa_s1(Z_M); s_SP <- mfa_s1(Z_SP)
  Z_comb <- cbind(Z_M / s_M, Z_SP / s_SP)
  list(Z_M = Z_M, Z_SP = Z_SP, Z_comb = Z_comb,
       D = list(scar     = stats::dist(as.matrix(Z_SP),   method = "euclidean"),
                morph    = stats::dist(as.matrix(Z_M),    method = "euclidean"),
                combined = stats::dist(as.matrix(Z_comb), method = "euclidean")),
       s1 = c(M = s_M, SP = s_SP))
}
sp_full <- build_spaces(keep0)

SPACE_LABELS <- c(scar     = "D_scar (SP-SPHARM, l = 1-6)",
                  morph    = "D_morph (M-SPHARM, l = 1-8)",
                  combined = "D_comb (MFA-normalised join)")

# ---- anchor check -----------------------------------------------------------
cat("\n", strrep("-", 74), "\n", sep = "")
if (length(keep0) == REF$n) {
  ag_M  <- permanova_global(sp_full$Z_M,  meta0$core_type)
  ag_SP <- permanova_global(sp_full$Z_SP, meta0$core_type)
  ok <- all(abs(c(ag_M$R2 - REF$M_R2, ag_M$F - REF$M_F,
                  ag_SP$R2 - REF$SP_R2, ag_SP$F - REF$SP_F)) <= REF_TOL)
  cat(sprintf(paste0("ANCHOR CHECK (core type alone, committed Table 2):\n",
                     "  M-SPHARM  R2 = %.5f (ref %.5f) | SP-SPHARM R2 = %.5f (ref %.5f)\n"),
              ag_M$R2, REF$M_R2, ag_SP$R2, REF$SP_R2))
  if (!ok)
    stop("ANCHOR CHECK FAILED — the ILR coordinates no longer reproduce the ",
         "committed SDG core-type PERMANOVA. Resolve the sample definition first.")
  cat("  => reproduces the published SDG core-type PERMANOVA.\n")
} else {
  cat(sprintf("ANCHOR CHECK SKIPPED: %d of %d specimens dropped for missing SDI.\n",
              REF$n - length(keep0), REF$n))
}

# =============================================================================
# (0) PRE-CHECKS
# =============================================================================
cat("\n", strrep("=", 74), "\n", sep = "")
cat("(0) PRE-CHECKS\n")
cat(strrep("=", 74), "\n", sep = "")
cat(sprintf("  n = %d | core types = %d | raw materials = %d\n",
            nrow(meta0), nlevels(meta0$core_type), nlevels(meta0$Raw_mat)))

# ---- 0a raw material x core type, and how collinear they are ----------------
tab_rc <- table(meta0$Raw_mat, meta0$core_type)
cat("\n(0a) RAW MATERIAL x CORE TYPE\n")
print(tab_rc)
cat("\n  raw-material composition within each core type (%):\n")
print(round(100 * prop.table(tab_rc, margin = 2), 1))

chi <- suppressWarnings(stats::chisq.test(tab_rc))
cramers_v <- sqrt(as.numeric(chi$statistic) /
                  (sum(tab_rc) * (min(dim(tab_rc)) - 1)))
type_conc <- apply(prop.table(tab_rc, margin = 2), 2, max)
pure_types <- names(type_conc)[type_conc == 1]
cat(sprintf(paste0("\n  Cramer's V(Raw_mat, core_type) = %.3f (chi-sq = %.2f, df = %d, p = %.4f)\n",
                   "  median within-type concentration in one material = %.1f%%\n",
                   "  core types confined ENTIRELY to one raw material: %s\n"),
            cramers_v, chi$statistic, chi$parameter, chi$p.value,
            100 * median(type_conc),
            ifelse(length(pure_types) == 0, "none",
                   paste(sprintf("%s (n = %d)", pure_types,
                                 colSums(tab_rc)[pure_types]), collapse = ", "))))
cat("  => the sequential SS removes raw material FIRST, so variance shared\n")
cat("     between material and these types is gone before SDI is reached.\n")

# ---- 0b SDI by raw material, and the common support ------------------------
cat("\n(0b) SDI BY RAW MATERIAL\n")
sdi_desc <- meta0 %>% group_by(Raw_mat) %>%
  summarise(n = n(),
            min = min(SDI), Q1 = quantile(SDI, .25), median = median(SDI),
            Q3 = quantile(SDI, .75), max = max(SDI), .groups = "drop")
print(as.data.frame(sdi_desc %>%
  mutate(across(c(min, Q1, median, Q3, max), ~ signif(.x, 4)))), row.names = FALSE)

lo_cs <- max(tapply(meta0$SDI, meta0$Raw_mat, min))
hi_cs <- min(tapply(meta0$SDI, meta0$Raw_mat, max))
in_cs <- meta0$SDI >= lo_cs & meta0$SDI <= hi_cs
cat(sprintf("\n  median ratio (chert / sandstone) = %.2f\n",
            median(meta0$SDI[meta0$Raw_mat == "chert"]) /
            median(meta0$SDI[meta0$Raw_mat == "sandstone"])))
cat(sprintf("  COMMON SUPPORT = [%.4g, %.4g]\n", lo_cs, hi_cs))
cat("  sample inside the common support, by raw material:\n")
print(as.data.frame(meta0 %>% mutate(in_common = in_cs) %>%
  count(Raw_mat, in_common) %>% pivot_wider(names_from = in_common,
    values_from = n, values_fill = 0, names_prefix = "in_cs_")),
  row.names = FALSE)
cat("  core-type sizes inside the common support:\n")
print(table(meta0$core_type[in_cs]))
cat("  => versions C and D are estimated on this subset; outside it, comparing\n")
cat("     the two slopes relies on linear extrapolation.\n")

# ---- 0c skewness, and whether to log ---------------------------------------
skewness <- function(x) mean((x - mean(x))^3) / stats::sd(x)^3
cat("\n(0c) SDI SKEWNESS\n")
cat(sprintf("  raw SDI      skewness = %+.3f | max/min = %.1f\n",
            skewness(meta0$SDI), max(meta0$SDI) / min(meta0$SDI)))
cat(sprintf("  log(SDI)     skewness = %+.3f\n", skewness(log(meta0$SDI))))
cat(sprintf("  Shapiro-Wilk raw p = %.3g | log p = %.3g\n",
            stats::shapiro.test(meta0$SDI)$p.value,
            stats::shapiro.test(log(meta0$SDI))$p.value))
cat("  => both transforms are carried through as versions A/C (raw) and B/D (log);\n")
cat("     the transform is NOT selected on the outcome.\n")

# =============================================================================
# MAIN TEST
# =============================================================================
cat("\n", strrep("=", 74), "\n", sep = "")
cat("MAIN TEST: D ~ Raw_mat + core_type + SDI + Raw_mat:SDI, by = \"terms\"\n")
cat(strrep("=", 74), "\n", sep = "")
cat(sprintf("  seed = %d | permutations = %d\n", SEED, N_PERM))
cat("  term order FIXED: material -> type -> SDI -> interaction (increment)\n")

ADONIS_COLS <- c("Df", "SumOfSqs", "R2", "F", "Pr(>F)")

run_one <- function(D, md, perm_ctrl) {
  set.seed(SEED)
  res <- adonis2(D ~ Raw_mat + core_type + SDI + Raw_mat:SDI, data = md,
                 by = "terms", permutations = perm_ctrl)
  tab <- as.data.frame(res)
  if (!all(ADONIS_COLS %in% names(tab)))
    stop("adonis2() returned columns {", paste(names(tab), collapse = ", "),
         "}; vegan's output layout changed — stopping rather than mislabelling.")
  tibble(term = rownames(tab), df = as.integer(tab[["Df"]]),
         SumOfSqs = tab[["SumOfSqs"]], R2 = tab[["R2"]],
         pseudo_F = tab[["F"]], p = tab[["Pr(>F)"]])
}

# The four versions. `rows` indexes meta0 / the analysed subset; `sdi_fun` is the
# transform applied to SDI before it enters the model.
VERSIONS <- list(
  A = list(label = "A raw SDI, full sample",       rows = rep(TRUE, nrow(meta0)), tf = identity),
  B = list(label = "B log(SDI), full sample",      rows = rep(TRUE, nrow(meta0)), tf = log),
  C = list(label = "C raw SDI, common support",    rows = in_cs,                  tf = identity),
  D = list(label = "D log(SDI), common support",   rows = in_cs,                  tf = log))

results <- list()
for (v in names(VERSIONS)) {
  V   <- VERSIONS[[v]]
  idx <- which(V$rows)
  md  <- meta0[idx, ] %>%
    mutate(core_type = droplevels(core_type), Raw_mat = droplevels(Raw_mat),
           SDI = V$tf(SDI)) %>%
    as.data.frame()
  sp  <- build_spaces(keep0[idx])

  if (nlevels(md$Raw_mat) < 2) {
    cat(sprintf("\n  [%s] SKIPPED: only one raw material survives (n = %d).\n",
                v, nrow(md)))
    next
  }
  cat(sprintf("\n  [%s] n = %d | %s | core types = %d | %s\n",
              v, nrow(md), paste(sprintf("%s %d", levels(md$Raw_mat),
                                         as.integer(table(md$Raw_mat))),
                                 collapse = " / "),
              nlevels(md$core_type),
              sprintf("ILR: M %d, SP %d coords", ncol(sp$Z_M), ncol(sp$Z_SP))))

  for (spn in names(SPACE_LABELS)) {
    # (i) free permutation
    free <- run_one(sp$D[[spn]], md, N_PERM)
    # (ii) restricted to within raw material. Raw_mat is constant inside a block,
    # so its own p is not interpretable here; kept in the table and flagged.
    ctrl <- permute::how(blocks = md$Raw_mat, nperm = N_PERM)
    rest <- run_one(sp$D[[spn]], md, ctrl)
    results[[length(results) + 1]] <- bind_rows(
      free %>% mutate(permutation = "free"),
      rest %>% mutate(permutation = "restricted_within_rawmat")) %>%
      mutate(version = v, version_label = V$label, space = spn,
             space_label = SPACE_LABELS[[spn]], n = nrow(md),
             n_perm = N_PERM, seed = SEED, .before = 1)
  }
}
perm_tbl <- bind_rows(results)

# ---- (1) full ANOVA tables --------------------------------------------------
cat("\n", strrep("=", 74), "\n", sep = "")
cat("(1) FULL ANOVA TABLES\n")
cat(strrep("=", 74), "\n", sep = "")
for (v in unique(perm_tbl$version)) for (spn in names(SPACE_LABELS)) {
  for (pm in c("free", "restricted_within_rawmat")) {
    t1 <- perm_tbl %>% filter(version == v, space == spn, permutation == pm)
    if (!nrow(t1)) next
    cat(sprintf("\n--- %s | %s | permutation = %s ---\n",
                VERSIONS[[v]]$label, SPACE_LABELS[[spn]], pm))
    print(t1 %>% transmute(term, df, R2 = round(R2, 5),
                           pseudo_F = round(pseudo_F, 4), p) %>%
            as.data.frame(), row.names = FALSE)
    ia <- t1 %>% filter(term == "Raw_mat:SDI")
    if (nrow(ia) == 1)
      cat(sprintf("    >>> INTERACTION Raw_mat:SDI  df = %d, R2 = %.5f, F = %.4f, p = %.4f  [%s]\n",
                  ia$df, ia$R2, ia$pseudo_F, ia$p,
                  ifelse(ia$p < 0.05, "SIGNIFICANT", "not significant")))
    if (pm == "restricted_within_rawmat")
      cat("        (Raw_mat's own p is not interpretable under blocking)\n")
  }
}

# ---- (2) free vs restricted -------------------------------------------------
cat("\n", strrep("=", 74), "\n", sep = "")
cat("(2) FREE vs RESTRICTED PERMUTATION — interaction term only\n")
cat(strrep("=", 74), "\n", sep = "")
ia_tbl <- perm_tbl %>% filter(term == "Raw_mat:SDI") %>%
  select(version, space, permutation, R2, pseudo_F, p) %>%
  pivot_wider(names_from = permutation, values_from = c(R2, pseudo_F, p))
print(as.data.frame(ia_tbl %>%
  transmute(version, space,
            R2 = round(R2_free, 5),
            F_free = round(pseudo_F_free, 3), p_free = p_free,
            F_rest = round(pseudo_F_restricted_within_rawmat, 3),
            p_rest = p_restricted_within_rawmat,
            agree = ifelse((p_free < 0.05) == (p_restricted_within_rawmat < 0.05),
                           "yes", "NO"))), row.names = FALSE)
n_disagree <- sum((ia_tbl$p_free < 0.05) !=
                  (ia_tbl$p_restricted_within_rawmat < 0.05))
cat(sprintf("\n  schemes disagree in %d of %d space x version combinations\n",
            n_disagree, nrow(ia_tbl)))
cat("  R2 is identical between schemes by construction: permutation changes the\n")
cat("  null distribution, not the observed partition of variance.\n")

# ---- (3) agreement across the four versions --------------------------------
cat("\n", strrep("=", 74), "\n", sep = "")
cat("(3) AGREEMENT ACROSS THE FOUR VERSIONS (restricted permutation = headline)\n")
cat(strrep("=", 74), "\n", sep = "")
head_tbl <- perm_tbl %>%
  filter(term == "Raw_mat:SDI", permutation == "restricted_within_rawmat") %>%
  select(version, space, df, R2, pseudo_F, p) %>%
  arrange(space, version)
print(as.data.frame(head_tbl %>% mutate(R2 = round(R2, 5),
                                        pseudo_F = round(pseudo_F, 4),
                                        sig = ifelse(p < 0.05, "*", ""))),
      row.names = FALSE)
n_sig <- sum(head_tbl$p < 0.05)
cat(sprintf("\n  significant in %d of %d version x space combinations\n",
            n_sig, nrow(head_tbl)))
if (n_sig > 0 && n_sig < nrow(head_tbl))
  cat("  => NOT ROBUST: the interaction appears only in some versions. Listed above.\n")

# =============================================================================
# FIGURE — presentation only, NOT the basis of any inference
# =============================================================================
pc1_tbl <- map_dfr(names(SPACE_LABELS), function(spn) {
  Z <- switch(spn, scar = sp_full$Z_SP, morph = sp_full$Z_M, combined = sp_full$Z_comb)
  pc <- stats::prcomp(Z, center = TRUE, scale. = FALSE)
  ev <- pc$sdev^2 / sum(pc$sdev^2)
  tibble(space = spn,
         space_label = sprintf("%s\nPC1 = %.1f%% of inertia", SPACE_LABELS[[spn]],
                               100 * ev[1]),
         ID = meta0$ID, Raw_mat = meta0$Raw_mat, SDI = meta0$SDI,
         PC1 = pc$x[, 1])
})

p_fig <- ggplot(pc1_tbl, aes(SDI, PC1, colour = Raw_mat, fill = Raw_mat)) +
  geom_point(size = 1.6, alpha = .85) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, linewidth = .6,
              alpha = .15) +
  geom_vline(xintercept = c(lo_cs, hi_cs), linetype = "dotted", linewidth = .3,
             colour = "grey40") +
  facet_wrap(~ space_label, nrow = 1, scales = "free") +
  scale_x_log10() +
  scale_colour_manual(values = RAWMAT_COLORS, name = "raw material") +
  scale_fill_manual(values = RAWMAT_COLORS, name = "raw material") +
  labs(title = "SDI against PC1 of each space, by raw material",
       subtitle = paste0("presentation only — PC1 is not known to be the SDI direction ",
                         "(no generalisable linear SDI direction exists in these spaces);\n",
                         "dotted lines bound the common support; fits span only each ",
                         "material's own SDI range"),
       x = "SDI (log scale)", y = "PC1 score") +
  theme_bw(base_size = 9) + theme(legend.position = "bottom")

ggsave(file.path(FIG_DIR, "fig_sdi_rawmat_interaction.png"), p_fig,
       width = 10.5, height = 4.4, dpi = 300)

write_csv(perm_tbl, file.path(OUT_DIR, "sdi_rawmat_interaction.csv"))
cat("\nWrote sdi_rawmat_interaction.csv and figures/fig_sdi_rawmat_interaction.png\n")

# =============================================================================
# (4) VERDICT
# =============================================================================
cat("\n", strrep("=", 74), "\n", sep = "")
cat("(4) VERDICT\n")
cat(strrep("=", 74), "\n", sep = "")
if (n_sig == 0) {
  cat("  No raw-material-specific SDI trend is detectable in any space, under\n")
  cat("  either permutation scheme, in any of the four versions.\n")
  cat("  => The null SDI main effect reported in sdi_gradient.R is NOT the result\n")
  cat("     of opposite trends in chert and sandstone cancelling each other out.\n")
  cat("     That explanation is removed, which strengthens the existing null.\n")
  cat("  Power caveat stands: interactions need far more data than main effects,\n")
  cat(sprintf("     and the interaction R2 here is %.4f-%.4f. A small material-specific\n",
              min(head_tbl$R2), max(head_tbl$R2)))
  cat("     difference cannot be excluded; only a large one can.\n")
} else if (n_sig == nrow(head_tbl)) {
  cat("  A raw-material-specific SDI trend is detected consistently, in every\n")
  cat("  space and version under restricted permutation. See the tables above.\n")
} else {
  cat(sprintf("  MIXED: the interaction reaches p < 0.05 in %d of %d version x space\n",
              n_sig, nrow(head_tbl)))
  cat("  combinations and not the rest. This is reported as an unstable result,\n")
  cat("  not as a positive finding. The combinations involved:\n")
  print(as.data.frame(head_tbl %>% filter(p < 0.05) %>%
                        transmute(version, space, R2 = round(R2, 5), p)),
        row.names = FALSE)
}

cat("\n"); print(sessionInfo())
