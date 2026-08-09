# sdi_layer_interaction.R
# =============================================================================
# DOES THE REDUCTION-INTENSITY (SDI) GRADIENT RUN DIFFERENTLY BETWEEN LAYERS?
# — the layer twin of sdi_rawmat_interaction.R, same folder, same machinery.
#
# NEW, self-contained. Does NOT modify the main pipeline, the cached _targets
# store, derived_data, the manuscript, sdi_gradient.R, or its raw-material
# sibling in this folder. Writes only sdi_layer_interaction.* here.
#
# THE QUESTION
#   sdi_gradient.R found no SDI main effect once raw material and core type were
#   removed (R2 = 0.017-0.027, p = 0.19-0.41). The raw-material twin then removed
#   the "opposite slopes in the two materials cancel" explanation. This script
#   asks the same question of stratigraphy: do Layer 3 and Layer 4 cores follow
#   different techno-morphological trajectories along the reduction gradient?
#   A layer-specific gradient would mean the two layers reached similar end forms
#   by different reduction paths — a diachronic reading of the same reviewer-4
#   question about what continuous quantification buys.
#
# TWO DIFFERENCES FROM THE RAW-MATERIAL TWIN, BOTH FORCED BY THE DATA
#
#   1. LAYER 2 IS DROPPED. It holds 2 cores in this sample, which cannot support
#      a slope. This matches both the repo precedent (sdi_gradient.R's layer
#      facets use Layers 3 and 4 only) and the published Table 2 layer test,
#      where safe_filter_groups(min_n = 3) drops Layer 2. The analysed factor is
#      therefore 2-level, exactly parallel to raw material.
#
#   2. A SECOND MODEL CONTROLS FOR RAW MATERIAL. Layer and raw material are not
#      independent here (Layer 3 is 17% sandstone, Layer 4 is 30%), and raw
#      material is itself strongly confounded with SDI (4x median difference).
#      A Layer:SDI signal could therefore be a Raw_mat:SDI signal in disguise.
#      Model 2 removes raw material first so that the layer interaction is an
#      increment over it:
#        M1  D ~ Layer + core_type + SDI + Layer:SDI          (direct analogue)
#        M2  D ~ Raw_mat + Layer + core_type + SDI + Layer:SDI (material-controlled)
#      M1 is the analogue of the raw-material twin and is reported first; M2 is
#      what decides whether any M1 signal is about stratigraphy at all. Neither
#      is chosen on the outcome; both are always reported.
#
# Term order is FIXED in both models. Sequential (Type I) SS throughout: the
# question is what the interaction ADDS after the grouping factors and the SDI
# main effect, not what it would explain alone.
#
# TWO PERMUTATION SCHEMES, as in the twin and for the same reason
#   (i)  free permutation
#   (ii) restricted WITHIN layer, how(blocks = Layer)
#   Free permutation dissolves the layer/SDI association that exists in the data
#   and can understate p; restricted permutation holds layer fixed and tests
#   "given layer, does the SDI effect differ", which is the question asked. Both
#   are reported; where they disagree the RESTRICTED result is the headline.
#   Under blocking, Layer is constant within blocks and its own p is not
#   interpretable — expected, not a failure.
#
# FOUR VERSIONS, as in the twin
#   A raw SDI full | B log(SDI) full | C raw SDI common support | D log common support
#   Unlike raw material, the two layers' SDI ranges overlap almost completely, so
#   C and D are expected to retain nearly the whole sample and to be far less
#   informative as an extrapolation check than they were for raw material. That
#   is reported rather than treated as a stronger replication than it is.
#
# EXPECTATION MANAGEMENT, declared before the results were seen
#   Same as the twin: interactions need far more data than main effects, the SDI
#   main effect is R2 = 0.017-0.027, and Layer 3 has 18 cores. A null is the
#   high-probability outcome and is informative. No term is reordered, no
#   marginal-SS result is promoted, no specimen is dropped and no transform is
#   picked to manufacture significance. If only some version x space combinations
#   are significant, that instability is the finding.
#
# Writes: sdi_layer_interaction.csv
#         figures/fig_sdi_layer_interaction.png
#
# HOW TO RUN (Docker spharm_analysis):
#   Rscript analysis/robustness/sdi_rawmat_interaction/sdi_layer_interaction.R
# =============================================================================

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

# Layers retained. Layer 2 (n = 2 here) cannot support a slope; same restriction
# as sdi_gradient.R's layer facets and the published Table 2 layer test.
LAYER_LEVELS <- c("3", "4")

OUT_DIR <- here("analysis/robustness/sdi_rawmat_interaction")
FIG_DIR <- file.path(OUT_DIR, "figures")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

REF <- list(M_R2  = 0.18770190528243877, M_F  = 2.0334613330095146,
            SP_R2 = 0.16785198589124128, SP_F = 1.7750417603590793,
            n = 50L, n_groups = 6L)
REF_TOL <- 1e-3

# Layer colours from the repo's earthy palette; two levels, hue plus lightness.
LAYER_COLORS <- c("3" = "#788C4A", "4" = "#802520")

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
    stop("Companion script not found at ", f, " — stopping.")

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
  stop("Could not import: ", paste(missing_imports, collapse = ", "), " — stopping.")
cat(sprintf("Imported %d helpers from %s and %d data-prep objects from %s.\n",
            length(got_helpers), basename(EXP_SCRIPT),
            length(got_data),    basename(SDG_SCRIPT)))

# =============================================================================
# SAMPLE + SDI  (identical construction to the raw-material twin)
# =============================================================================
blocks <- build_blocks_sdg()
if (length(blocks$ids) != REF$n || nlevels(blocks$group) != REF$n_groups)
  stop(sprintf("Imported sample gives n = %d / %d groups, not %d / %d. Stopping.",
               length(blocks$ids), nlevels(blocks$group), REF$n, REF$n_groups))

sdi_xl <- read_excel(SDI_XLSX)
if (!SDI_COL %in% names(sdi_xl))
  stop("Column '", SDI_COL, "' not found in ", basename(SDI_XLSX), ".")
sdi_tbl <- tibble(ID  = str_trim(as.character(sdi_xl$ID)),
                  SDI = suppressWarnings(as.numeric(sdi_xl[[SDI_COL]])))
if (anyDuplicated(sdi_tbl$ID)) stop("Duplicate IDs in ", basename(SDI_XLSX), ".")

meta_all <- tibble(ID = blocks$ids, core_type = blocks$group, layer = blocks$layer) %>%
  left_join(core_meta %>% select(ID, raw_material), by = "ID") %>%
  left_join(sdi_tbl, by = "ID")
if (nrow(meta_all) != length(blocks$ids))
  stop("Metadata join changed the row count; stopping.")

build_spaces <- function(row_idx_in_blocks) {
  morph_sub <- blocks$morph[row_idx_in_blocks, , drop = FALSE]
  scar_sub  <- blocks$scar[ row_idx_in_blocks, , drop = FALSE]
  Z_M  <- make_ilr(morph_sub); Z_SP <- make_ilr(scar_sub)
  colnames(Z_M)  <- paste0("M_ilr",  seq_len(ncol(Z_M)))
  colnames(Z_SP) <- paste0("SP_ilr", seq_len(ncol(Z_SP)))
  s_M <- mfa_s1(Z_M); s_SP <- mfa_s1(Z_SP)
  Z_comb <- cbind(Z_M / s_M, Z_SP / s_SP)
  list(Z_M = Z_M, Z_SP = Z_SP, Z_comb = Z_comb,
       D = list(scar     = stats::dist(as.matrix(Z_SP),   method = "euclidean"),
                morph    = stats::dist(as.matrix(Z_M),    method = "euclidean"),
                combined = stats::dist(as.matrix(Z_comb), method = "euclidean")))
}
SPACE_LABELS <- c(scar     = "D_scar (SP-SPHARM, l = 1-6)",
                  morph    = "D_morph (M-SPHARM, l = 1-8)",
                  combined = "D_comb (MFA-normalised join)")

# ---- anchor check on the FULL published sample, before any layer restriction --
keep_full <- which(is.finite(meta_all$SDI))
cat("\n", strrep("-", 74), "\n", sep = "")
if (length(keep_full) == REF$n) {
  sp0 <- build_spaces(keep_full)
  ag_M  <- permanova_global(sp0$Z_M,  droplevels(meta_all$core_type[keep_full]))
  ag_SP <- permanova_global(sp0$Z_SP, droplevels(meta_all$core_type[keep_full]))
  if (!all(abs(c(ag_M$R2 - REF$M_R2, ag_SP$R2 - REF$SP_R2)) <= REF_TOL))
    stop("ANCHOR CHECK FAILED — resolve the sample definition first.")
  cat(sprintf(paste0("ANCHOR CHECK (core type alone, n = %d, committed Table 2):\n",
                     "  M-SPHARM R2 = %.5f (ref %.5f) | SP-SPHARM R2 = %.5f (ref %.5f)\n",
                     "  => reproduces the published SDG core-type PERMANOVA.\n"),
              REF$n, ag_M$R2, REF$M_R2, ag_SP$R2, REF$SP_R2))
} else {
  cat(sprintf("ANCHOR CHECK SKIPPED: %d of %d dropped for missing SDI.\n",
              REF$n - length(keep_full), REF$n))
}

# =============================================================================
# (0) PRE-CHECKS
# =============================================================================
cat("\n", strrep("=", 74), "\n", sep = "")
cat("(0) PRE-CHECKS\n")
cat(strrep("=", 74), "\n", sep = "")

cat("\n  layer composition of the full published sample (before restriction):\n")
print(table(meta_all$layer[keep_full]))
n_l2 <- sum(meta_all$layer[keep_full] == "2")
cat(sprintf("  Layer 2 holds %d core(s) and is dropped: a 2-core group cannot\n", n_l2))
cat("  support a slope. Same restriction as sdi_gradient.R's layer facets and\n")
cat("  the published Table 2 layer test (safe_filter_groups drops it).\n")

keep0 <- which(is.finite(meta_all$SDI) & meta_all$layer %in% LAYER_LEVELS)
meta0 <- meta_all[keep0, ] %>%
  mutate(core_type = droplevels(core_type),
         Layer     = factor(as.character(layer), levels = LAYER_LEVELS),
         Raw_mat   = factor(raw_material))
if (any(meta0$SDI <= 0))
  stop("SDI has non-positive values; log versions would be undefined. Stopping.")
cat(sprintf("\n  ANALYSED n = %d | layers = %s | core types = %d\n",
            nrow(meta0), paste(levels(meta0$Layer), collapse = "/"),
            nlevels(meta0$core_type)))

# ---- 0a layer x core type, and layer x raw material ------------------------
tab_lc <- table(meta0$Layer, meta0$core_type)
cat("\n(0a) LAYER x CORE TYPE\n"); print(tab_lc)
cat("\n  layer composition within each core type (%):\n")
print(round(100 * prop.table(tab_lc, margin = 2), 1))

cramers_v <- function(tb) {
  ch <- suppressWarnings(stats::chisq.test(tb))
  list(v = sqrt(as.numeric(ch$statistic) / (sum(tb) * (min(dim(tb)) - 1))),
       chi = as.numeric(ch$statistic), df = as.integer(ch$parameter),
       p = ch$p.value)
}
cv_lc <- cramers_v(tab_lc)
conc_lc <- apply(prop.table(tab_lc, margin = 2), 2, max)
pure_lc <- names(conc_lc)[conc_lc == 1]
cat(sprintf(paste0("\n  Cramer's V(Layer, core_type) = %.3f (chi-sq = %.2f, df = %d, p = %.4f)\n",
                   "  median within-type concentration in one layer = %.1f%%\n",
                   "  core types confined ENTIRELY to one layer: %s\n"),
            cv_lc$v, cv_lc$chi, cv_lc$df, cv_lc$p, 100 * median(conc_lc),
            ifelse(length(pure_lc) == 0, "none",
                   paste(sprintf("%s (n = %d)", pure_lc,
                                 colSums(tab_lc)[pure_lc]), collapse = ", "))))

tab_lr <- table(meta0$Layer, meta0$Raw_mat)
cat("\n  LAYER x RAW MATERIAL (why model M2 exists):\n"); print(tab_lr)
cat("  raw-material composition within each layer (%):\n")
print(round(100 * prop.table(tab_lr, margin = 1), 1))
cv_lr <- cramers_v(tab_lr)
cat(sprintf(paste0("  Cramer's V(Layer, Raw_mat) = %.3f (chi-sq = %.2f, df = %d, p = %.4f)\n",
                   "  => layer and raw material are%s independent, and raw material is\n",
                   "     itself strongly confounded with SDI, so a Layer:SDI signal could\n",
                   "     be a Raw_mat:SDI signal in disguise. Model M2 removes material first.\n"),
            cv_lr$v, cv_lr$chi, cv_lr$df, cv_lr$p,
            ifelse(cv_lr$p < 0.05, " NOT", " approximately")))

# ---- 0b SDI by layer, and the common support -------------------------------
cat("\n(0b) SDI BY LAYER\n")
print(as.data.frame(meta0 %>% group_by(Layer) %>%
  summarise(n = n(), min = min(SDI), Q1 = quantile(SDI, .25),
            median = median(SDI), Q3 = quantile(SDI, .75), max = max(SDI),
            .groups = "drop") %>%
  mutate(across(c(min, Q1, median, Q3, max), ~ signif(.x, 4)))), row.names = FALSE)

lo_cs <- max(tapply(meta0$SDI, meta0$Layer, min))
hi_cs <- min(tapply(meta0$SDI, meta0$Layer, max))
in_cs <- meta0$SDI >= lo_cs & meta0$SDI <= hi_cs
med_ratio <- median(meta0$SDI[meta0$Layer == "3"]) /
             median(meta0$SDI[meta0$Layer == "4"])
cat(sprintf("\n  median ratio (L3 / L4) = %.2f\n", med_ratio))
cat(sprintf("  COMMON SUPPORT = [%.4g, %.4g]  (retains %d of %d cores, %.0f%%)\n",
            lo_cs, hi_cs, sum(in_cs), nrow(meta0), 100 * mean(in_cs)))
cat("  sample inside the common support, by layer:\n")
print(as.data.frame(meta0 %>% mutate(in_common = in_cs) %>%
  count(Layer, in_common) %>% pivot_wider(names_from = in_common,
    values_from = n, values_fill = 0, names_prefix = "in_cs_")), row.names = FALSE)
cat("  => the two layers' SDI ranges overlap almost completely, so unlike the\n")
cat("     raw-material twin, versions C and D are NOT a strong extrapolation\n")
cat("     check here: they re-test nearly the same sample.\n")

# ---- 0c skewness ------------------------------------------------------------
skewness <- function(x) mean((x - mean(x))^3) / stats::sd(x)^3
cat("\n(0c) SDI SKEWNESS\n")
cat(sprintf("  raw SDI  skewness = %+.3f | max/min = %.1f | Shapiro p = %.3g\n",
            skewness(meta0$SDI), max(meta0$SDI) / min(meta0$SDI),
            stats::shapiro.test(meta0$SDI)$p.value))
cat(sprintf("  log(SDI) skewness = %+.3f | Shapiro p = %.3g\n",
            skewness(log(meta0$SDI)), stats::shapiro.test(log(meta0$SDI))$p.value))
cat("  => both transforms carried through; the transform is NOT selected on the outcome.\n")

# =============================================================================
# MAIN TEST — two models x four versions x three spaces x two permutations
# =============================================================================
cat("\n", strrep("=", 74), "\n", sep = "")
cat("MAIN TEST (sequential SS, by = \"terms\")\n")
cat("  M1  D ~ Layer + core_type + SDI + Layer:SDI\n")
cat("  M2  D ~ Raw_mat + Layer + core_type + SDI + Layer:SDI\n")
cat(strrep("=", 74), "\n", sep = "")
cat(sprintf("  seed = %d | permutations = %d\n", SEED, N_PERM))

ADONIS_COLS <- c("Df", "SumOfSqs", "R2", "F", "Pr(>F)")
MODELS <- list(
  M1 = list(label = "M1 Layer first (direct analogue)",
            fml = D ~ Layer + core_type + SDI + Layer:SDI),
  M2 = list(label = "M2 raw material controlled",
            fml = D ~ Raw_mat + Layer + core_type + SDI + Layer:SDI))

run_one <- function(D, md, fml, perm_ctrl) {
  set.seed(SEED)
  env <- list2env(list(D = D), parent = environment())
  res <- adonis2(fml, data = md, by = "terms", permutations = perm_ctrl,
                 environment = env)
  tab <- as.data.frame(res)
  if (!all(ADONIS_COLS %in% names(tab)))
    stop("adonis2() returned columns {", paste(names(tab), collapse = ", "),
         "}; vegan's output layout changed — stopping.")
  tibble(term = rownames(tab), df = as.integer(tab[["Df"]]),
         SumOfSqs = tab[["SumOfSqs"]], R2 = tab[["R2"]],
         pseudo_F = tab[["F"]], p = tab[["Pr(>F)"]])
}

VERSIONS <- list(
  A = list(label = "A raw SDI, full sample",     rows = rep(TRUE, nrow(meta0)), tf = identity),
  B = list(label = "B log(SDI), full sample",    rows = rep(TRUE, nrow(meta0)), tf = log),
  C = list(label = "C raw SDI, common support",  rows = in_cs,                  tf = identity),
  D = list(label = "D log(SDI), common support", rows = in_cs,                  tf = log))

results <- list()
for (v in names(VERSIONS)) {
  V   <- VERSIONS[[v]]
  idx <- which(V$rows)
  md  <- meta0[idx, ] %>%
    mutate(core_type = droplevels(core_type), Layer = droplevels(Layer),
           Raw_mat = droplevels(Raw_mat), SDI = V$tf(SDI)) %>%
    as.data.frame()
  sp  <- build_spaces(keep0[idx])
  if (nlevels(md$Layer) < 2) {
    cat(sprintf("\n  [%s] SKIPPED: only one layer survives (n = %d).\n", v, nrow(md)))
    next
  }
  cat(sprintf("\n  [%s] n = %d | %s | core types = %d\n", v, nrow(md),
              paste(sprintf("L%s %d", levels(md$Layer),
                            as.integer(table(md$Layer))), collapse = " / "),
              nlevels(md$core_type)))
  for (mn in names(MODELS)) for (spn in names(SPACE_LABELS)) {
    free <- run_one(sp$D[[spn]], md, MODELS[[mn]]$fml, N_PERM)
    ctrl <- permute::how(blocks = md$Layer, nperm = N_PERM)
    rest <- run_one(sp$D[[spn]], md, MODELS[[mn]]$fml, ctrl)
    results[[length(results) + 1]] <- bind_rows(
      free %>% mutate(permutation = "free"),
      rest %>% mutate(permutation = "restricted_within_layer")) %>%
      mutate(model = mn, model_label = MODELS[[mn]]$label,
             version = v, version_label = V$label, space = spn,
             space_label = SPACE_LABELS[[spn]], n = nrow(md),
             n_perm = N_PERM, seed = SEED, .before = 1)
  }
}
perm_tbl <- bind_rows(results)
IA <- "Layer:SDI"

# ---- (1) full ANOVA tables --------------------------------------------------
cat("\n", strrep("=", 74), "\n", sep = "")
cat("(1) FULL ANOVA TABLES\n")
cat(strrep("=", 74), "\n", sep = "")
for (mn in names(MODELS)) for (v in unique(perm_tbl$version))
  for (spn in names(SPACE_LABELS)) for (pm in c("free", "restricted_within_layer")) {
    t1 <- perm_tbl %>% filter(model == mn, version == v, space == spn,
                              permutation == pm)
    if (!nrow(t1)) next
    cat(sprintf("\n--- %s | %s | %s | permutation = %s ---\n",
                MODELS[[mn]]$label, VERSIONS[[v]]$label, SPACE_LABELS[[spn]], pm))
    print(t1 %>% transmute(term, df, R2 = round(R2, 5),
                           pseudo_F = round(pseudo_F, 4), p) %>%
            as.data.frame(), row.names = FALSE)
    ia <- t1 %>% filter(term == IA)
    if (nrow(ia) == 1)
      cat(sprintf("    >>> INTERACTION %s  df = %d, R2 = %.5f, F = %.4f, p = %.4f  [%s]\n",
                  IA, ia$df, ia$R2, ia$pseudo_F, ia$p,
                  ifelse(ia$p < 0.05, "SIGNIFICANT", "not significant")))
    if (pm == "restricted_within_layer")
      cat("        (Layer's own p is not interpretable under blocking)\n")
  }

# ---- (2) free vs restricted -------------------------------------------------
cat("\n", strrep("=", 74), "\n", sep = "")
cat("(2) FREE vs RESTRICTED PERMUTATION — interaction term only\n")
cat(strrep("=", 74), "\n", sep = "")
ia_tbl <- perm_tbl %>% filter(term == IA) %>%
  select(model, version, space, permutation, R2, pseudo_F, p) %>%
  pivot_wider(names_from = permutation, values_from = c(R2, pseudo_F, p))
print(as.data.frame(ia_tbl %>%
  transmute(model, version, space, R2 = round(R2_free, 5),
            F_free = round(pseudo_F_free, 3), p_free = p_free,
            p_rest = p_restricted_within_layer,
            agree = ifelse((p_free < 0.05) == (p_restricted_within_layer < 0.05),
                           "yes", "NO"))), row.names = FALSE)
n_disagree <- sum((ia_tbl$p_free < 0.05) !=
                  (ia_tbl$p_restricted_within_layer < 0.05))
cat(sprintf("\n  schemes disagree in %d of %d model x version x space combinations\n",
            n_disagree, nrow(ia_tbl)))

# ---- (3) agreement across versions and models ------------------------------
cat("\n", strrep("=", 74), "\n", sep = "")
cat("(3) AGREEMENT ACROSS VERSIONS AND MODELS (restricted permutation = headline)\n")
cat(strrep("=", 74), "\n", sep = "")
head_tbl <- perm_tbl %>%
  filter(term == IA, permutation == "restricted_within_layer") %>%
  select(model, version, space, df, R2, pseudo_F, p) %>%
  arrange(model, space, version)
print(as.data.frame(head_tbl %>% mutate(R2 = round(R2, 5),
                                        pseudo_F = round(pseudo_F, 4),
                                        sig = ifelse(p < 0.05, "*", ""))),
      row.names = FALSE)
n_sig <- sum(head_tbl$p < 0.05)
cat(sprintf("\n  significant in %d of %d model x version x space combinations\n",
            n_sig, nrow(head_tbl)))
if (n_sig > 0 && n_sig < nrow(head_tbl))
  cat("  => NOT ROBUST: the interaction appears only in some combinations.\n")

# =============================================================================
# FIGURE — presentation only, NOT the basis of any inference
# =============================================================================
# As in the raw-material twin: PC1 is NOT known to be the SDI-related direction.
# sdi_gradient.R's candidate search found no generalisable linear SDI direction in
# any of the three spaces (29 non-circular candidates, all |rho| < 0.30, LOO rho
# all <= 0). A flat panel is the expected picture and does not contradict the
# PERMANOVA. Fits are per layer, so each line spans only its own SDI range.
sp_an <- build_spaces(keep0)
pc1_tbl <- map_dfr(names(SPACE_LABELS), function(spn) {
  Z  <- switch(spn, scar = sp_an$Z_SP, morph = sp_an$Z_M, combined = sp_an$Z_comb)
  pc <- stats::prcomp(Z, center = TRUE, scale. = FALSE)
  ev <- pc$sdev^2 / sum(pc$sdev^2)
  tibble(space = spn,
         space_label = sprintf("%s\nPC1 = %.1f%% of inertia",
                               SPACE_LABELS[[spn]], 100 * ev[1]),
         ID = meta0$ID, Layer = meta0$Layer, SDI = meta0$SDI, PC1 = pc$x[, 1])
})

p_fig <- ggplot(pc1_tbl, aes(SDI, PC1, colour = Layer, fill = Layer)) +
  geom_point(size = 1.6, alpha = .85) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, linewidth = .6,
              alpha = .15) +
  geom_vline(xintercept = c(lo_cs, hi_cs), linetype = "dotted", linewidth = .3,
             colour = "grey40") +
  facet_wrap(~ space_label, nrow = 1, scales = "free") +
  scale_x_log10() +
  scale_colour_manual(values = LAYER_COLORS, name = "layer") +
  scale_fill_manual(values = LAYER_COLORS, name = "layer") +
  labs(title = "SDI against PC1 of each space, by layer (Layers 3 and 4)",
       subtitle = paste0("presentation only — PC1 is not known to be the SDI direction ",
                         "(no generalisable linear SDI direction exists in these spaces);\n",
                         "dotted lines bound the common support, which here covers ",
                         "almost the whole range; fits span only each layer's own SDI range"),
       x = "SDI (log scale)", y = "PC1 score") +
  theme_bw(base_size = 9) + theme(legend.position = "bottom")

ggsave(file.path(FIG_DIR, "fig_sdi_layer_interaction.png"), p_fig,
       width = 10.5, height = 4.4, dpi = 300)

write_csv(perm_tbl, file.path(OUT_DIR, "sdi_layer_interaction.csv"))
cat("\nWrote sdi_layer_interaction.csv and figures/fig_sdi_layer_interaction.png\n")

# =============================================================================
# (4) VERDICT
# =============================================================================
cat("\n", strrep("=", 74), "\n", sep = "")
cat("(4) VERDICT\n")
cat(strrep("=", 74), "\n", sep = "")
if (n_sig == 0) {
  cat("  No layer-specific SDI trend is detectable in any space, under either\n")
  cat("  permutation scheme, in any version, in either model.\n")
  cat("  => The null SDI main effect is NOT the result of opposite trends in\n")
  cat("     Layer 3 and Layer 4 cancelling out. That explanation is removed.\n")
  cat(sprintf("  Power caveat: interaction R2 spans %.4f-%.4f, L3 has %d cores.\n",
              min(head_tbl$R2), max(head_tbl$R2), sum(meta0$Layer == "3")))
  cat("     A small layer-specific difference cannot be excluded; a large one can.\n")
} else if (n_sig == nrow(head_tbl)) {
  cat("  A layer-specific SDI trend is detected consistently, in every space,\n")
  cat("  version and model under restricted permutation.\n")
} else {
  cat(sprintf("  MIXED: the interaction reaches p < 0.05 in %d of %d combinations.\n",
              n_sig, nrow(head_tbl)))
  cat("  Reported as unstable, not as a positive finding. Combinations involved:\n")
  print(as.data.frame(head_tbl %>% filter(p < 0.05) %>%
          transmute(model, version, space, R2 = round(R2, 5), p)), row.names = FALSE)
  m1 <- head_tbl %>% filter(model == "M1", p < 0.05)
  m2 <- head_tbl %>% filter(model == "M2", p < 0.05)
  cat(sprintf("\n  significant under M1 (layer first): %d | under M2 (material controlled): %d\n",
              nrow(m1), nrow(m2)))
  if (nrow(m1) > 0 && nrow(m2) == 0)
    cat("  => every M1 signal disappears once raw material is removed first, so it\n     is attributable to raw material rather than to stratigraphy.\n")
}

cat("\n"); print(sessionInfo())
