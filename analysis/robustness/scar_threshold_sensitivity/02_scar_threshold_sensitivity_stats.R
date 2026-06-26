# 02_scar_threshold_sensitivity_stats.R
# =============================================================================
# Scar minimum-SIZE-THRESHOLD sensitivity analysis for the SP-SPHARM pipeline —
# SI add-on. NEW, self-contained file. Does NOT modify the main pipeline, the
# cached _targets store, the derived_data cache, or the manuscript. It only READS
# the committed production outputs (for a sanity check) and the per-threshold
# SP-SPHARM power spectra produced by 01_sweep_spharm_threshold.py, plus the
# attrition tables from 00_scar_attrition.py, and WRITES new outputs under
# analysis/robustness/scar_threshold_sensitivity/.
#
# It re-uses the project's existing statistical machinery (the same package
# functions the main pipeline calls — vegan::adonis2 / mantel, ade4::coinertia /
# randtest, compositions::ilr) and replicates, verbatim, the data-prep steps from
# the main scripts (same source-line attributions as 02_bandwidth_sensitivity_stats.R):
#   - r_spharm/power_degree_selection.R     (degree selection)
#   - r_spharm/spharm_analysis.R           (EXP PERMANOVA `perm_dir`)
#   - r_statistics/exp_cores_statistics.R  (EXP Mantel + RV)
#   - r_statistics/SDG_cores_statistics.R  (SDG Mantel + RV + scar~core-type PERMANOVA)
#
# For every threshold T it evaluates whether the paper's conclusions hold:
#   (a) order selection : SP-SPHARM cumulative power (% at l=6) and the degree at
#                         which cross-specimen CV first exceeds 100% (EXP & SDG).
#   (c) experimental    : PERMANOVA core-type R2 / pseudo-F / p and the set of
#                         significant pairwise distinctions (resolution profile).
#   (d) decoupling      : global Mantel r and RV (EXP and SDG); plus the SDG
#                         scar~core-type PERMANOVA.
# Part (b) — ideal-core discriminability — is UNCHANGED BY CONSTRUCTION: the size
# threshold is not applied to the synthetic ideal cores (their scar lengths are
# fixed, non-physical values; see 00_scar_attrition.py / README.md), so the ideal
# cores are held at production values and the ideal-core distance structure is
# identical at every threshold. It is reported once, as the production reference.
#
# Outputs (all NEW):
#   threshold_sensitivity_metrics.csv         tidy: one row per threshold
#   threshold_orderselection_by_degree.csv    per-degree cumulative power & CV per T
#   figures/fig_S_threshold_orderselection.png
#   figures/fig_S_threshold_scarcounts.png
#
# HOW TO RUN (canonical environment, R 4.4 + renv):
#   Rscript analysis/robustness/scar_threshold_sensitivity/02_scar_threshold_sensitivity_stats.R
#   # prerequisite: run 00_scar_attrition.py and 01_sweep_spharm_threshold.py first.
# =============================================================================

suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(readxl)
  library(vegan)          # adonis2, betadisper, mantel
  library(ade4)           # dudi.pca, coinertia, randtest (RV)
  library(compositions)   # ilr
  library(RVAideMemoire)  # pairwise.perm.manova
  library(patchwork)
  library(conflicted)
})

# Resolve namespace clashes exactly as the main pipeline does (_targets.R).
suppressMessages({
  conflicts_prefer(dplyr::filter, dplyr::select, dplyr::lag,
                   stats::sd, stats::var, stats::dist, stats::cor, stats::cov,
                   base::scale, base::norm, .quiet = TRUE)
})

set.seed(42)

# =============================================================================
# PARAMETERS  (edit here; must match the Python sweep 01_sweep_spharm_threshold.py)
# =============================================================================
# Minimum-size cutoffs in mm. T = 0 is the production anchor (all recorded scars).
THRESHOLDS <- c(0.0, 5.0, 10.0)
T_REF      <- 0.0

OUT_DIR     <- here("analysis/robustness/scar_threshold_sensitivity")
SPECTRA_DIR <- file.path(OUT_DIR, "spectra")
FIG_DIR     <- file.path(OUT_DIR, "figures")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

spectra_path <- function(T) file.path(SPECTRA_DIR, sprintf("SPHARM_direction_t%04.1f.csv", T))

# Power-column conventions (identical to the main analysis).
POWER_COLS_ALL   <- paste0("power_l", 1:20)  # order selection
POWER_COLS_DIR   <- paste0("power_l", 1:6)   # scar (SP-SPHARM) descriptor
POWER_COLS_MORPH <- paste0("power_l", 1:8)   # morphology (M-SPHARM) descriptor
POWER_COLS_IM    <- paste0("power_l", 1:4)   # ideal-core discriminability (methods_comparison_IM.R:254)

# Typology handling (spharm_analysis.R:46-62)
EXCLUDE_TYPES   <- c("Biface")
LEVALLOIS_MERGE <- c("Levallois convergent", "Levallois laminar",
                     "Levallois preferential", "Levallois recurrent")
TYPOLOGY_ORDER  <- c("Unidirectional", "Bidirectional", "Levallois",
                     "Discoid", "Multiplatform")

# SDG metadata handling (SDG_cores_statistics.R:38-49)
EXCLUDE_CORE_TYPES <- c("Handaxe", "Pick")

# Committed production reference values used for the sanity check at T = 0 (which is
# identical to the published h = 0.35, all-scar pipeline). Deterministic statistics
# must reproduce these; permutation p-values carry ~+/-0.005 Monte-Carlo jitter and
# are excluded from the check. Values from bandwidth_sensitivity_metrics.csv (h=0.35)
# / the cached derived_data CSVs.
REF <- list(
  exp_cumpower_l6 = 99.92, exp_cv_cross_l = 9,
  exp_perm_R2 = 0.30174, exp_perm_F = 5.72570,  # ILR / Aitchison anchor (was z-score 0.32624 / 6.41589)
  exp_mantel_r = -0.06126, exp_RV = 0.10095,
  sdg_mantel_r = 0.01189,  sdg_RV = 0.09066,
  sdg_perm_scar_coretype_R2 = 0.16785, sdg_perm_scar_coretype_F = 1.77504
)

# =============================================================================
# Helpers — copied VERBATIM from the source scripts (attribution in comments)
# =============================================================================

# exp_cores_statistics.R:65-79 / SDG_cores_statistics.R:76-90
replace_zeros <- function(X, delta = NULL) {
  X <- as.matrix(X)
  for (i in seq_len(nrow(X))) {
    row_i    <- X[i, ]
    zero_idx <- row_i == 0
    if (!any(zero_idx)) next
    nonzero_min      <- min(row_i[!zero_idx])
    d                <- ifelse(is.null(delta), nonzero_min * 0.65, delta)
    n_zero           <- sum(zero_idx)
    row_i[zero_idx]  <- d
    row_i[!zero_idx] <- row_i[!zero_idx] * (1 - n_zero * d)
    X[i, ] <- row_i
  }
  X
}

# ILR / Aitchison helper (exp/SDG_cores_statistics.R convention): drop zero-variance
# columns, multiplicative zero replacement, then ilr into Euclidean space.
make_ilr <- function(power_df) {
  X    <- as.matrix(power_df)
  keep <- apply(X, 2, function(v) sd(v, na.rm = TRUE) > 0)
  as.matrix(ilr(replace_zeros(X[, keep, drop = FALSE])))
}

# exp_cores_statistics.R:81-83
extract_subdist <- function(D_full, ids) as.dist(as.matrix(D_full)[ids, ids])

# SDG_cores_statistics.R:96-110
safe_filter_groups <- function(meta_df, group_col, min_n = 3) {
  counts <- table(meta_df[[group_col]], useNA = "no")
  valid  <- names(counts[counts >= min_n])
  if (length(valid) < 2) return(NULL)
  meta_df %>% filter(!is.na(.data[[group_col]]), .data[[group_col]] %in% valid)
}

# spharm_analysis.R:82-100
filter_spharm <- function(df, power_cols, meta = NULL) {
  result <- df %>% select(ID, Typology, all_of(power_cols))
  if (!is.null(meta)) result <- left_join(result, meta, by = "ID")
  result
}
split_by_group <- function(df) {
  list(
    exp_im = df %>% filter(str_starts(ID, "EXP") | str_starts(ID, "IM_")),
    sdg_im = df %>% filter(str_starts(ID, "SDG") | str_starts(ID, "IM_"))
  )
}

# spharm_analysis.R:109-118 — z-score using the non-IM EXP reference mean/sd
scale_features <- function(df_target, cols) {
  ref_mat  <- df_target %>% filter(!str_starts(ID, "IM_")) %>%
    select(all_of(cols)) %>% as.matrix()
  col_mean <- colMeans(ref_mat)
  col_sd   <- apply(ref_mat, 2, sd)
  mat      <- df_target %>% select(all_of(cols)) %>% as.matrix()
  base::scale(mat, center = col_mean, scale = col_sd)
}

# power_degree_selection.R:74-124 (trimmed to the quantities used here)
compute_order_stats <- function(df, cols) {
  mat       <- df %>% select(all_of(cols)) %>% as.matrix()
  n_orders  <- length(cols)
  col_means <- colMeans(mat, na.rm = TRUE)
  col_sds   <- apply(mat, 2, sd, na.rm = TRUE)
  col_cvs   <- col_sds / col_means * 100
  row_sums  <- rowSums(mat, na.rm = TRUE)
  total_mean<- mean(row_sums)
  cumul_pct <- cumsum(col_means) / total_mean * 100
  tibble(order = seq_len(n_orders), cv_pct = col_cvs, cumul_pct = cumul_pct)
}

# spharm_analysis.R:302-329 — EXP PERMANOVA on z-scored power (global + Holm pairwise)
run_permanova_dir <- function(X, group_vec) {
  df_grp <- data.frame(Typology = group_vec)
  d      <- stats::dist(X, method = "euclidean")
  set.seed(42)
  perm_global <- adonis2(X ~ Typology, data = df_grp,
                         method = "euclidean", permutations = 999)
  set.seed(42)
  pairwise_res <- pairwise.perm.manova(d, group_vec, nperm = 999, p.method = "holm")
  list(
    R2 = perm_global$R2[1], F = perm_global$`F`[1], p = perm_global$`Pr(>F)`[1],
    pairwise = pairwise_res$p.value
  )
}

# =============================================================================
# Static (threshold-INDEPENDENT) inputs — loaded once
# Morphology (M-SPHARM) is mesh-derived and does not involve scars, so it is reused
# unchanged at every threshold.
# =============================================================================
cat("Loading threshold-independent inputs (morphology, metadata)...\n")

SPHARM_morphology <- read_csv(
  here("analysis/data/derived_data/SPHARM_morphology.csv"), show_col_types = FALSE)

metric_data <- read_xlsx(here("analysis/data/raw_data/SDG_core_metric.xlsx")) %>%
  select(ID, Layer, Core_type_Li_merged, Raw_mat) %>%
  mutate(Layer = as.factor(Layer), Raw_mat = as.factor(Raw_mat))

core_meta <- read_excel(here("analysis/data/raw_data/SDG_core_metric.xlsx")) %>%
  select(ID = ID, raw_material = Raw_mat, core_type = Core_type_Li_merged) %>%
  mutate(across(everything(), ~ str_trim(as.character(.))))

# =============================================================================
# Per-threshold analysis blocks (identical estimators to the main pipeline)
# =============================================================================

# ---- (a) order selection: EXP & SDG -----------------------------------------
order_selection_block <- function(dir_df) {
  dir_exp <- dir_df %>% filter(str_starts(ID, "EXP"))
  dir_sdg <- dir_df %>% filter(str_starts(ID, "SDG"), !str_starts(ID, "IM_"))
  cols    <- intersect(POWER_COLS_ALL, colnames(dir_df))
  list(EXP = compute_order_stats(dir_exp, cols),
       SDG = compute_order_stats(dir_sdg, cols))
}
first_cv_cross <- function(os) {
  idx <- which(os$cv_pct > 100)
  if (length(idx) == 0) NA_integer_ else os$order[idx[1]]
}
cumpower_at <- function(os, l) os$cumul_pct[os$order == l]

# ---- (b) ideal-core distance — computed ONCE (threshold-invariant) ----------
im_distance_matrix <- function(dir_df) {
  IM <- dir_df %>% filter(str_starts(ID, "IM_")) %>% arrange(ID)
  X  <- IM %>% select(all_of(POWER_COLS_IM)) %>% as.matrix()
  d  <- as.matrix(stats::dist(base::scale(X), method = "euclidean"))
  rownames(d) <- colnames(d) <- IM$ID
  d
}

# ---- (c) EXP PERMANOVA `perm_dir` (spharm_analysis.R:106-138, 386-388) -------
exp_permanova_block <- function(dir_df) {
  SPHARM_direction_filter <- filter_spharm(dir_df, POWER_COLS_DIR, metric_data)
  df_exp_dir <- split_by_group(SPHARM_direction_filter)$exp_im
  df_exp_only <- df_exp_dir %>%
    filter(!str_starts(ID, "IM_"), !Typology %in% EXCLUDE_TYPES) %>%
    mutate(Typology = case_when(Typology %in% LEVALLOIS_MERGE ~ "Levallois",
                                TRUE ~ Typology),
           Typology = droplevels(as.factor(Typology)))
  non_im_idx <- !str_starts(df_exp_dir$ID, "IM_") &
    !df_exp_dir$Typology %in% EXCLUDE_TYPES
  y_typology <- df_exp_only$Typology
  ilr_dir <- make_ilr(df_exp_dir[non_im_idx, POWER_COLS_DIR])   # ILR / Aitchison (was z-score)
  pm <- run_permanova_dir(ilr_dir, y_typology)

  pv  <- pm$pairwise
  sig <- which(pv < 0.05, arr.ind = TRUE)
  sig_pairs <- if (nrow(sig) == 0) character(0) else
    apply(sig, 1, function(rc) paste(rownames(pv)[rc[1]], colnames(pv)[rc[2]], sep = "-"))
  n_pairs <- sum(!is.na(pv))
  list(R2 = pm$R2, F = pm$F, p = pm$p,
       n_sig = length(sig_pairs), n_pairs = n_pairs,
       sig_pairs = paste(sort(sig_pairs), collapse = "; "))
}

# ---- (d) EXP decoupling: global Mantel + RV (exp_cores_statistics.R:157-403) -
exp_decoupling_block <- function(dir_df) {
  morph_typ <- SPHARM_morphology %>%
    left_join(dir_df %>% select(ID, Typology), by = "ID")
  dir_f   <- filter_spharm(dir_df,    POWER_COLS_DIR,   metric_data)
  morph_f <- filter_spharm(morph_typ, POWER_COLS_MORPH, metric_data)

  df_scar_all  <- split_by_group(dir_f)$exp_im
  df_morph_all <- split_by_group(morph_f)$exp_im
  common <- intersect(df_morph_all$ID, df_scar_all$ID)
  df_morph_all <- df_morph_all %>% filter(ID %in% common) %>% arrange(ID)
  df_scar_all  <- df_scar_all  %>% filter(ID %in% common) %>% arrange(ID)

  morph_power <- df_morph_all %>% select(all_of(POWER_COLS_MORPH)) %>%
    rename_with(~ paste0("M", seq_along(.))) %>% as.data.frame()
  scar_power  <- df_scar_all %>% select(all_of(POWER_COLS_DIR)) %>%
    rename_with(~ paste0("S", seq_along(.))) %>% as.data.frame()
  rownames(morph_power) <- df_morph_all$ID
  rownames(scar_power)  <- df_scar_all$ID
  morph_power <- morph_power[, sapply(morph_power, sd, na.rm = TRUE) > 0]
  scar_power  <- scar_power[,  sapply(scar_power,  sd, na.rm = TRUE) > 0]

  morph_ilr <- as.data.frame(ilr(replace_zeros(as.matrix(morph_power))))
  scar_ilr  <- as.data.frame(ilr(replace_zeros(as.matrix(scar_power))))
  rownames(morph_ilr) <- rownames(morph_power)
  rownames(scar_ilr)  <- rownames(scar_power)
  D_morph_all <- stats::dist(morph_ilr)
  D_scar_all  <- stats::dist(scar_ilr)

  exp_ids <- rownames(morph_power)[!str_starts(rownames(morph_power), "IM_") &
                                     rownames(morph_power) != "EXP43_Biface"]
  D_morph_exp <- extract_subdist(D_morph_all, exp_ids)
  D_scar_exp  <- extract_subdist(D_scar_all,  exp_ids)

  mantel_global <- mantel(D_morph_exp, D_scar_exp, method = "spearman",
                          permutations = 9999)

  morph_ilr_exp <- morph_ilr[exp_ids, ]; scar_ilr_exp <- scar_ilr[exp_ids, ]
  colnames(morph_ilr_exp) <- paste0("M_ilr", seq_len(ncol(morph_ilr_exp)))
  colnames(scar_ilr_exp)  <- paste0("S_ilr", seq_len(ncol(scar_ilr_exp)))
  dudi_morph <- dudi.pca(morph_ilr_exp, center = TRUE, scale = TRUE,
                         scannf = FALSE, nf = ncol(morph_ilr_exp))
  dudi_scar  <- dudi.pca(scar_ilr_exp,  center = TRUE, scale = TRUE,
                         scannf = FALSE, nf = ncol(scar_ilr_exp))
  coin_exp <- coinertia(dudi_morph, dudi_scar, scannf = FALSE, nf = 2)
  set.seed(42)
  rv_test <- randtest(coin_exp, nrepet = 9999)
  list(mantel_r = mantel_global$statistic, mantel_p = mantel_global$signif,
       RV = coin_exp$RV, RV_p = rv_test$pvalue, n = length(exp_ids))
}

# ---- (d) SDG decoupling + scar~core-type PERMANOVA --------------------------
sdg_block <- function(dir_df) {
  morph_typ <- SPHARM_morphology %>%
    left_join(dir_df %>% select(ID, Typology), by = "ID")
  df_scar_raw  <- filter_spharm(dir_df,    POWER_COLS_DIR,   metric_data)
  df_morph_raw <- filter_spharm(morph_typ, POWER_COLS_MORPH, metric_data)
  common <- intersect(df_morph_raw$ID, df_scar_raw$ID)
  df_morph_raw <- df_morph_raw %>% filter(ID %in% common) %>% arrange(ID)
  df_scar_raw  <- df_scar_raw  %>% filter(ID %in% common) %>% arrange(ID)

  morph_power <- df_morph_raw %>% select(all_of(POWER_COLS_MORPH)) %>%
    rename_with(~ paste0("M", seq_along(.))) %>% as.data.frame()
  scar_power  <- df_scar_raw %>% select(all_of(POWER_COLS_DIR)) %>%
    rename_with(~ paste0("S", seq_along(.))) %>% as.data.frame()
  rownames(morph_power) <- df_morph_raw$ID
  rownames(scar_power)  <- df_scar_raw$ID
  morph_power <- morph_power[, sapply(morph_power, sd, na.rm = TRUE) > 0]
  scar_power  <- scar_power[,  sapply(scar_power,  sd, na.rm = TRUE) > 0]

  morph_ilr <- as.data.frame(ilr(replace_zeros(as.matrix(morph_power))))
  scar_ilr  <- as.data.frame(ilr(replace_zeros(as.matrix(scar_power))))
  rownames(morph_ilr) <- rownames(morph_power)
  rownames(scar_ilr)  <- rownames(scar_power)
  D_morph_all <- stats::dist(morph_ilr)
  D_scar_all  <- stats::dist(scar_ilr)

  arch_ids <- rownames(morph_power)[!str_starts(rownames(morph_power), "IM_") &
                                      !str_starts(rownames(morph_power), "EXP")]
  meta_arch <- tibble(ID = arch_ids) %>%
    mutate(layer = case_when(str_detect(ID, "L2") ~ "Layer 2",
                             str_detect(ID, "L3") ~ "Layer 3",
                             str_detect(ID, "L4") ~ "Layer 4",
                             TRUE ~ "Other")) %>%
    left_join(core_meta, by = "ID") %>%
    filter(!core_type %in% EXCLUDE_CORE_TYPES | is.na(core_type))
  arch_ids <- meta_arch$ID
  D_morph_arch <- extract_subdist(D_morph_all, arch_ids)
  D_scar_arch  <- extract_subdist(D_scar_all,  arch_ids)

  mantel_global <- mantel(D_morph_arch, D_scar_arch, method = "spearman",
                          permutations = 9999)

  morph_ilr_arch <- morph_ilr[arch_ids, ]; scar_ilr_arch <- scar_ilr[arch_ids, ]
  colnames(morph_ilr_arch) <- paste0("M_ilr", seq_len(ncol(morph_ilr_arch)))
  colnames(scar_ilr_arch)  <- paste0("S_ilr", seq_len(ncol(scar_ilr_arch)))
  dudi_morph <- dudi.pca(morph_ilr_arch, center = TRUE, scale = TRUE,
                         scannf = FALSE, nf = ncol(morph_ilr_arch))
  dudi_scar  <- dudi.pca(scar_ilr_arch,  center = TRUE, scale = TRUE,
                         scannf = FALSE, nf = ncol(scar_ilr_arch))
  coin_arch <- coinertia(dudi_morph, dudi_scar, scannf = FALSE, nf = 2)
  set.seed(42)
  rv_test <- randtest(coin_arch, nrepet = 9999)

  meta_coretype <- safe_filter_groups(meta_arch, "core_type")
  perm <- list(R2 = NA_real_, F = NA_real_, p = NA_real_)
  if (!is.null(meta_coretype)) {
    D_scar_ct <- extract_subdist(D_scar_arch, meta_coretype$ID)
    set.seed(42)
    res <- adonis2(D_scar_ct ~ group,
                   data = data.frame(group = factor(meta_coretype$core_type)),
                   permutations = 9999, add = "lingoes")
    perm <- list(R2 = res$R2[1], F = res$`F`[1], p = res$`Pr(>F)`[1])
  }
  list(mantel_r = mantel_global$statistic, mantel_p = mantel_global$signif,
       RV = coin_arch$RV, RV_p = rv_test$pvalue, n = length(arch_ids),
       perm_R2 = perm$R2, perm_F = perm$F, perm_p = perm$p)
}

# =============================================================================
# Sweep over thresholds
# =============================================================================
metrics    <- list()
order_long <- list()
im_ref     <- NULL          # ideal-core distance matrix (threshold-invariant)

for (T in THRESHOLDS) {
  csv <- spectra_path(T)
  if (!file.exists(csv)) {
    warning(sprintf("Missing spectra for T=%.1f (%s); run 01_sweep_spharm_threshold.py first. Skipping.",
                    T, basename(csv)))
    next
  }
  cat(sprintf("\n=== T = %.1f mm ===\n", T))
  dir_df <- read_csv(csv, show_col_types = FALSE)
  dir_df$ID <- str_trim(dir_df$ID)

  os  <- order_selection_block(dir_df)
  ep  <- exp_permanova_block(dir_df)
  ed  <- exp_decoupling_block(dir_df)
  sd_ <- sdg_block(dir_df)
  if (abs(T - T_REF) < 1e-9) im_ref <- im_distance_matrix(dir_df)

  # retained scar counts from the spectra n_scars column (cross-check w/ attrition)
  ns <- if ("n_scars" %in% names(dir_df)) dir_df %>% mutate(
          grp = case_when(str_starts(ID, "EXP") ~ "EXP",
                          str_starts(ID, "SDG") ~ "SDG",
                          str_starts(ID, "IM_") ~ "IM", TRUE ~ "OTHER")) else NULL
  scars_exp <- if (!is.null(ns)) sum(ns$n_scars[ns$grp == "EXP"], na.rm = TRUE) else NA
  scars_sdg <- if (!is.null(ns)) sum(ns$n_scars[ns$grp == "SDG"], na.rm = TRUE) else NA

  order_long[[sprintf("%.1f", T)]] <- bind_rows(
    os$EXP %>% mutate(dataset = "EXP"),
    os$SDG %>% mutate(dataset = "SDG")
  ) %>% mutate(threshold_mm = T)

  metrics[[sprintf("%.1f", T)]] <- tibble(
    threshold_mm = T,
    scars_kept_EXP = scars_exp, scars_kept_SDG = scars_sdg,
    exp_cumpower_l6 = cumpower_at(os$EXP, 6),
    exp_cv_cross_l  = first_cv_cross(os$EXP),
    sdg_cumpower_l6 = cumpower_at(os$SDG, 6),
    sdg_cv_cross_l  = first_cv_cross(os$SDG),
    exp_perm_R2 = ep$R2, exp_perm_F = ep$F, exp_perm_p = ep$p,
    exp_perm_nsig = ep$n_sig, exp_perm_npairs = ep$n_pairs,
    exp_perm_sig_pairs = ep$sig_pairs,
    exp_mantel_r = ed$mantel_r, exp_mantel_p = ed$mantel_p,
    exp_RV = ed$RV, exp_RV_p = ed$RV_p, exp_n = ed$n,
    sdg_mantel_r = sd_$mantel_r, sdg_mantel_p = sd_$mantel_p,
    sdg_RV = sd_$RV, sdg_RV_p = sd_$RV_p, sdg_n = sd_$n,
    sdg_perm_scar_coretype_R2 = sd_$perm_R2,
    sdg_perm_scar_coretype_F  = sd_$perm_F,
    sdg_perm_scar_coretype_p  = sd_$perm_p
  )
}

metrics_df <- bind_rows(metrics)
if (nrow(metrics_df) == 0) stop("No spectra found. Run 01_sweep_spharm_threshold.py first.")
order_long_df <- bind_rows(order_long)

# Ideal-core reference separations (threshold-invariant; reported once).
im_summary <- NULL
if (!is.null(im_ref)) {
  pick <- function(a, b) if (all(c(a, b) %in% rownames(im_ref))) im_ref[a, b] else NA_real_
  im_summary <- list(
    mean_dist       = mean(im_ref[upper.tri(im_ref)]),
    disc_lev        = pick("IM_discoid", "IM_Levallois_preferential"),
    biface_discuni  = pick("IM_biface",  "IM_discoid_unifacial"))
}

# =============================================================================
# Sanity check at T = T_REF vs the committed production outputs
# =============================================================================
cat("\n", strrep("=", 70), "\n", sep = "")
cat("SANITY CHECK at T =", T_REF, "mm vs committed production outputs (= h=0.35, all scars)\n")
cat(strrep("=", 70), "\n", sep = "")
ref_row <- metrics_df %>% filter(abs(threshold_mm - T_REF) < 1e-9)
sanity_ok <- TRUE
if (nrow(ref_row) == 1) {
  check <- function(name, got, exp, tol, kind = "deterministic") {
    ok <- is.finite(got) && abs(got - exp) <= tol
    if (!ok && kind == "deterministic") sanity_ok <<- FALSE
    cat(sprintf("  %-30s got %10.5f  expected %10.5f  %s\n",
                name, got, exp, ifelse(ok, "OK", paste0("<-- CHECK (", kind, ")"))))
  }
  check("EXP cumpower @ l6 (%)",  ref_row$exp_cumpower_l6, REF$exp_cumpower_l6, 0.05)
  check("EXP CV-cross degree",    ref_row$exp_cv_cross_l,  REF$exp_cv_cross_l,  0)
  check("EXP PERMANOVA R2",       ref_row$exp_perm_R2,     REF$exp_perm_R2,     1e-3)
  check("EXP PERMANOVA pseudo-F", ref_row$exp_perm_F,      REF$exp_perm_F,      5e-3)
  check("EXP Mantel r",           ref_row$exp_mantel_r,    REF$exp_mantel_r,    1e-3)
  check("EXP RV",                 ref_row$exp_RV,          REF$exp_RV,          1e-3)
  check("SDG Mantel r",           ref_row$sdg_mantel_r,    REF$sdg_mantel_r,    1e-3)
  check("SDG RV",                 ref_row$sdg_RV,          REF$sdg_RV,          1e-3)
  check("SDG scar~coretype R2",   ref_row$sdg_perm_scar_coretype_R2, REF$sdg_perm_scar_coretype_R2, 1e-3)
  check("SDG scar~coretype F",    ref_row$sdg_perm_scar_coretype_F,  REF$sdg_perm_scar_coretype_F,  1e-3)
  cat(sprintf("\n  => %s\n", ifelse(sanity_ok,
      "Deterministic statistics reproduce the cached production (T=0) values.",
      "MISMATCH: recomputed T=0 differs from cache — investigate before trusting the sweep.")))
  cat("  (Permutation p-values jitter by ~+/-0.005 and are not part of the check.)\n")
} else {
  cat("  T=0 not in the grid; cannot run the sanity check.\n")
}

# =============================================================================
# Outputs: tidy CSVs
# =============================================================================
write_csv(metrics_df,    file.path(OUT_DIR, "threshold_sensitivity_metrics.csv"))
write_csv(order_long_df, file.path(OUT_DIR, "threshold_orderselection_by_degree.csv"))
cat("\nWrote threshold_sensitivity_metrics.csv and threshold_orderselection_by_degree.csv\n")

# =============================================================================
# Figures
# =============================================================================
ok_fig <- TRUE
tryCatch({
  tv  <- sort(unique(order_long_df$threshold_mm))
  pal <- setNames(c("#4A6E8A", "#BA8530", "#802520")[seq_along(tv)],
                  sprintf("%.0f", tv))
  tlab <- function(t) ifelse(t == 0, "all (~2 mm)", paste0("> ", t, " mm"))
  order_long_df <- order_long_df %>%
    mutate(tf = factor(sprintf("%.0f", threshold_mm), levels = names(pal)))

  # --- Figure 1: order selection (cumulative power & CV by degree), EXP & SDG ---
  p_cum <- ggplot(order_long_df %>% filter(order <= 12),
                  aes(order, cumul_pct, color = tf, group = tf)) +
    geom_hline(yintercept = c(95, 99), linetype = "dashed", color = "grey75", linewidth = 0.3) +
    geom_vline(xintercept = 6, linetype = "dotted", color = "grey50", linewidth = 0.3) +
    geom_line(linewidth = 0.7) + geom_point(size = 1.3) +
    facet_wrap(~ dataset, labeller = as_labeller(
      c(EXP = "Experimentally knapped cores", SDG = "Sandinggai cores"))) +
    scale_color_manual(values = pal, name = "min. size",
                       labels = tlab(tv)) +
    scale_x_continuous(breaks = 1:12) +
    labs(x = "SPHARM degree (l)", y = "Cumulative power (%)") +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank(),
          strip.text = element_text(face = "bold"))

  p_cv <- ggplot(order_long_df %>% filter(order <= 12),
                 aes(order, cv_pct, color = tf, group = tf)) +
    geom_hline(yintercept = 100, linetype = "dashed", color = "grey50", linewidth = 0.4) +
    geom_vline(xintercept = 9, linetype = "dotted", color = "grey50", linewidth = 0.3) +
    geom_line(linewidth = 0.7) + geom_point(size = 1.3) +
    facet_wrap(~ dataset, labeller = as_labeller(
      c(EXP = "Experimentally knapped cores", SDG = "Sandinggai cores"))) +
    scale_color_manual(values = pal, name = "min. size", labels = tlab(tv)) +
    scale_x_continuous(breaks = 1:12) +
    labs(x = "SPHARM degree (l)", y = "Cross-specimen CV (%)") +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank(),
          strip.text = element_text(face = "bold"))

  ggsave(file.path(FIG_DIR, "fig_S_threshold_orderselection.png"),
         p_cum / p_cv, width = 9, height = 8, dpi = 300)

  # --- Figure 2: per-specimen scar retention (from 00_scar_attrition.py) -------
  att_csv <- file.path(OUT_DIR, "scar_attrition_by_specimen.csv")
  if (file.exists(att_csv)) {
    att <- read_csv(att_csv, show_col_types = FALSE) %>%
      filter(group %in% c("EXP", "SDG"))
    base_rank <- att %>% filter(threshold_mm == 0) %>%
      group_by(group) %>% arrange(desc(n_retained), .by_group = TRUE) %>%
      mutate(rank = row_number()) %>% select(group, ID, rank)
    att <- att %>% left_join(base_rank, by = c("group", "ID")) %>%
      mutate(tf = factor(sprintf("%.0f", threshold_mm), levels = names(pal)))
    p_sc <- ggplot(att, aes(rank, n_retained, color = tf, group = tf)) +
      geom_hline(yintercept = 3, linetype = "dashed", color = "grey60", linewidth = 0.4) +
      geom_line(linewidth = 0.5) + geom_point(size = 1) +
      facet_wrap(~ group, scales = "free_x", labeller = as_labeller(
        c(EXP = "Experimentally knapped cores", SDG = "Sandinggai cores"))) +
      scale_color_manual(values = pal, name = "min. size", labels = tlab(tv)) +
      labs(x = "Specimen (ranked by baseline scar count)", y = "Scars retained") +
      theme_bw(base_size = 11) +
      theme(panel.grid.minor = element_blank(),
            strip.text = element_text(face = "bold"))
    ggsave(file.path(FIG_DIR, "fig_S_threshold_scarcounts.png"),
           p_sc, width = 9, height = 4.5, dpi = 300)
    cat("Wrote 2 figures to", FIG_DIR, "\n")
  } else {
    cat("Wrote 1 figure to", FIG_DIR,
        "(scar-count figure skipped: run 00_scar_attrition.py first)\n")
  }
}, error = function(e) {
  ok_fig <<- FALSE
  cat("Figure generation failed (non-critical):", conditionMessage(e), "\n")
})

# =============================================================================
# Stability flags (printed by the SENSITIVITY FLAGS section below)
# =============================================================================
exp_perm_all_sig <- all(metrics_df$exp_perm_p < 0.05, na.rm = TRUE)
sdg_perm_all_sig <- all(metrics_df$sdg_perm_scar_coretype_p < 0.05, na.rm = TRUE)
mantel_all_ns    <- all(metrics_df$exp_mantel_p >= 0.05, na.rm = TRUE) &&
                    all(metrics_df$sdg_mantel_p >= 0.05, na.rm = TRUE)
rv_all_ns        <- all(metrics_df$exp_RV_p >= 0.05, na.rm = TRUE) &&
                    all(metrics_df$sdg_RV_p >= 0.05, na.rm = TRUE)
cv_cross_stable  <- length(unique(na.omit(metrics_df$exp_cv_cross_l))) == 1

# ---- flag any metric that is sensitive to the threshold (printed for the analyst)
cat("\n", strrep("=", 70), "\n", sep = "")
cat("SENSITIVITY FLAGS (metrics whose conclusion changes across thresholds)\n")
cat(strrep("=", 70), "\n", sep = "")
flag_any <- FALSE
if (!exp_perm_all_sig) { flag_any <- TRUE; cat("  * EXP scar PERMANOVA: p crosses 0.05 across thresholds.\n") }
if (!sdg_perm_all_sig) { flag_any <- TRUE; cat("  * SDG scar~core-type PERMANOVA: p crosses 0.05 across thresholds (report range).\n") }
if (!mantel_all_ns)    { flag_any <- TRUE; cat("  * Global Mantel: becomes significant at some threshold.\n") }
if (!rv_all_ns)        { flag_any <- TRUE; cat("  * RV: becomes significant at some threshold.\n") }
if (!cv_cross_stable)  { flag_any <- TRUE; cat(sprintf("  * EXP CV>100%% crossing degree varies: %s.\n",
                                                       paste(sort(unique(na.omit(metrics_df$exp_cv_cross_l))), collapse = ", "))) }
if (!flag_any) cat("  None: every headline conclusion is stable across the tested thresholds.\n")

cat("\nDone.\n")
