# mesh_sensitivity_stats.R
# =============================================================================
# Mesh-PREPROCESSING sensitivity analysis for the M-SPHARM (morphology) pipeline —
# SI add-on. NEW, self-contained file. Does NOT modify the main pipeline, the cached
# _targets store, the derived_data cache, or the manuscript. It only READS the
# committed production outputs (for the anchor check) + the fixed SP-SPHARM spectra,
# and the per-setting M-SPHARM spectra produced by sweep_mesh_spharm.py, and
# WRITES new outputs under analysis/robustness/mesh_decimation_sensitivity/.
#
# It re-uses the project's statistical machinery (the same package functions the
# main pipeline calls — vegan::adonis2 / mantel, ade4::coinertia / randtest,
# compositions::ilr) and replicates, verbatim, the data-prep from:
#   - r_spharm/power_degree_selection.R     (morphology degree selection)
#   - r_spharm/spharm_analysis.R           (EXP morphology PERMANOVA `perm_morph`)
#   - r_statistics/exp_cores_statistics.R  (EXP Mantel + RV decoupling)
#   - r_statistics/SDG_cores_statistics.R  (SDG morph PERMANOVA: core type / raw
#                                           material / layer; SDG Mantel + RV)
#
# Only the MORPHOLOGY (M-SPHARM) side is perturbed; the scar (SP-SPHARM) side is the
# FIXED committed SPHARM_direction.csv (scar vectors are digitised endpoint
# coordinates, not derived from these meshes). For every (face_target, smooth_iters)
# it evaluates whether the paper's morphology conclusions hold:
#   (a) order selection : M-SPHARM cumulative power (% at l=8) and whether cross-
#                         specimen CV stays < 90% through l=8 (the l=1-8 criterion).
#   (c) experimental    : EXP morphology PERMANOVA core-type R2/F/p + resolution profile.
#   (d) archaeological  : SDG morphology PERMANOVA by core type (should dominate),
#                         raw material and layer (should stay small / non-significant).
#   (e) decoupling      : global Mantel r and RV (EXP and SDG), NEW morphology x FIXED
#                         scar — does the no-covariation conclusion hold?
#
# Outputs (all NEW):
#   mesh_sensitivity_metrics.csv             tidy: one row per setting
#   mesh_orderselection_by_degree.csv        per-degree cumulative power & CV per setting
#   figures/fig_S_mesh_orderselection.png
#
# HOW TO RUN (canonical environment, R 4.4 + renv):
#   Rscript analysis/robustness/mesh_decimation_sensitivity/mesh_sensitivity_stats.R
#   # prerequisite: run sweep_mesh_spharm.py first.
# =============================================================================

suppressPackageStartupMessages({
  library(here); library(tidyverse); library(readxl)
  library(vegan); library(ade4); library(compositions)
  library(RVAideMemoire); library(patchwork); library(conflicted)
})
suppressMessages({
  conflicts_prefer(dplyr::filter, dplyr::select, dplyr::lag,
                   stats::sd, stats::var, stats::dist, stats::cor, stats::cov,
                   base::scale, base::norm, .quiet = TRUE)
})
set.seed(42)

# =============================================================================
# PARAMETERS  (must match the Python sweep sweep_mesh_spharm.py)
# =============================================================================
FACE_TARGETS <- c(10000, 20000, 50000)
SMOOTH_ITERS <- c(0, 3, 6)
PROD_FACES   <- 20000
PROD_SMOOTH  <- 3
FULL_CROSS   <- FALSE

# Two 1-D sweeps through production; de-duplicated by an explicit key (the shared
# production setting otherwise appears in both sweeps).
build_settings <- function() {
  s <- c(map(FACE_TARGETS, ~ list(face = .x, smooth = PROD_SMOOTH)),
         map(SMOOTH_ITERS, ~ list(face = PROD_FACES, smooth = .x)))
  if (FULL_CROSS) {
    g <- expand.grid(face = FACE_TARGETS, smooth = SMOOTH_ITERS)
    s <- c(s, pmap(g, function(face, smooth) list(face = face, smooth = smooth)))
  }
  keys <- map_chr(s, ~ sprintf("f%d_s%d", as.integer(.x$face), as.integer(.x$smooth)))
  s[!duplicated(keys)]
}
SETTINGS <- build_settings()

OUT_DIR     <- here("analysis/robustness/mesh_decimation_sensitivity")
SPECTRA_DIR <- file.path(OUT_DIR, "spectra")
FIG_DIR     <- file.path(OUT_DIR, "figures")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
spectra_path <- function(face, smooth)
  file.path(SPECTRA_DIR, sprintf("SPHARM_morphology_f%d_s%d.csv", face, smooth))

# Power-column conventions (identical to the main analysis).
POWER_COLS_ALL   <- paste0("power_l", 1:20)
POWER_COLS_DIR   <- paste0("power_l", 1:6)   # scar descriptor (FIXED side)
POWER_COLS_MORPH <- paste0("power_l", 1:8)   # morphology descriptor (perturbed side)

EXCLUDE_TYPES      <- c("Biface")
LEVALLOIS_MERGE    <- c("Levallois convergent", "Levallois laminar",
                        "Levallois preferential", "Levallois recurrent")
EXCLUDE_CORE_TYPES <- c("Handaxe", "Pick")

# Committed production reference values for the anchor check at (20000, 3). All are
# deterministic; permutation p-values jitter ~+/-0.005 and are excluded. Sources:
# DegreeSelection_stats_morphology_{EXP,SDG}.csv, L3_permanova.csv, EXP_L1_results.csv,
# L1_results.csv.
REF <- list(
  exp_morph_cumpower_l8 = 98.371, exp_morph_maxcv_l8 = 83.05,
  sdg_morph_cumpower_l8 = 98.291, sdg_morph_maxcv_l8 = 54.64,
  sdg_morph_coretype_R2 = 0.18770, sdg_morph_coretype_F = 2.03346,
  sdg_morph_rawmat_R2 = 0.04066, sdg_morph_layer_R2 = 0.01759,
  exp_mantel_r = -0.06126, exp_RV = 0.10095,
  sdg_mantel_r = 0.01189,  sdg_RV = 0.09066
)

# =============================================================================
# Helpers — copied VERBATIM from the source scripts (same as the bandwidth/scar
# sensitivity add-ons; attribution in comments).
# =============================================================================
replace_zeros <- function(X, delta = NULL) {            # exp/SDG_cores_statistics.R
  X <- as.matrix(X)
  for (i in seq_len(nrow(X))) {
    row_i <- X[i, ]; zero_idx <- row_i == 0
    if (!any(zero_idx)) next
    nonzero_min <- min(row_i[!zero_idx])
    d <- ifelse(is.null(delta), nonzero_min * 0.65, delta)
    n_zero <- sum(zero_idx)
    row_i[zero_idx] <- d; row_i[!zero_idx] <- row_i[!zero_idx] * (1 - n_zero * d)
    X[i, ] <- row_i
  }
  X
}
extract_subdist <- function(D_full, ids) as.dist(as.matrix(D_full)[ids, ids])
# ILR / Aitchison helper (exp/SDG_cores_statistics.R convention): drop zero-variance
# columns, multiplicative zero replacement, then ilr into Euclidean space.
make_ilr <- function(power_df) {
  X    <- as.matrix(power_df)
  keep <- apply(X, 2, function(v) sd(v, na.rm = TRUE) > 0)
  as.matrix(ilr(replace_zeros(X[, keep, drop = FALSE])))
}
safe_filter_groups <- function(meta_df, group_col, min_n = 3) {
  counts <- table(meta_df[[group_col]], useNA = "no")
  valid  <- names(counts[counts >= min_n])
  if (length(valid) < 2) return(NULL)
  meta_df %>% filter(!is.na(.data[[group_col]]), .data[[group_col]] %in% valid)
}
filter_spharm <- function(df, power_cols, meta = NULL) {
  result <- df %>% select(ID, Typology, all_of(power_cols))
  if (!is.null(meta)) result <- left_join(result, meta, by = "ID")
  result
}
split_by_group <- function(df) list(
  exp_im = df %>% filter(str_starts(ID, "EXP") | str_starts(ID, "IM_")),
  sdg_im = df %>% filter(str_starts(ID, "SDG") | str_starts(ID, "IM_")))
scale_features <- function(df_target, cols) {           # spharm_analysis.R:109-118
  ref_mat  <- df_target %>% filter(!str_starts(ID, "IM_")) %>%
    select(all_of(cols)) %>% as.matrix()
  col_mean <- colMeans(ref_mat); col_sd <- apply(ref_mat, 2, sd)
  mat <- df_target %>% select(all_of(cols)) %>% as.matrix()
  base::scale(mat, center = col_mean, scale = col_sd)
}
compute_order_stats <- function(df, cols) {             # power_degree_selection.R:74-124
  mat <- df %>% select(all_of(cols)) %>% as.matrix()
  col_means <- colMeans(mat, na.rm = TRUE)
  col_sds   <- apply(mat, 2, sd, na.rm = TRUE)
  cumul_pct <- cumsum(col_means) / mean(rowSums(mat, na.rm = TRUE)) * 100
  tibble(order = seq_along(cols), cv_pct = col_sds / col_means * 100, cumul_pct = cumul_pct)
}
run_permanova_dir <- function(X, group_vec) {           # spharm_analysis.R:302-329
  d <- stats::dist(X, method = "euclidean")
  set.seed(42)
  perm_global <- adonis2(X ~ Typology, data = data.frame(Typology = group_vec),
                         method = "euclidean", permutations = 999)
  set.seed(42)
  pairwise_res <- pairwise.perm.manova(d, group_vec, nperm = 999, p.method = "holm")
  list(R2 = perm_global$R2[1], F = perm_global$`F`[1], p = perm_global$`Pr(>F)`[1],
       pairwise = pairwise_res$p.value)
}
# SDG_cores_statistics.R:1515-1534 — distance-matrix PERMANOVA (ILR, lingoes).
run_permanova_dist <- function(dist_mat, group_vec) {
  group_vec <- as.character(group_vec)
  if (length(unique(group_vec)) < 2) return(list(R2 = NA_real_, F = NA_real_, p = NA_real_))
  set.seed(42)
  res <- adonis2(dist_mat ~ group, data = data.frame(group = factor(group_vec)),
                 permutations = 9999, add = "lingoes")
  list(R2 = res$R2[1], F = res$`F`[1], p = res$`Pr(>F)`[1])
}

# =============================================================================
# Static (perturbation-INDEPENDENT) inputs — the FIXED scar side + metadata
# =============================================================================
cat("Loading fixed inputs (committed SP-SPHARM scar spectra, metadata)...\n")
SPHARM_direction <- read_csv(here("analysis/data/derived_data/SPHARM_direction.csv"),
                             show_col_types = FALSE)
SPHARM_direction$ID <- str_trim(SPHARM_direction$ID)
TYPOLOGY_MAP <- SPHARM_direction %>% select(ID, Typology)

metric_data <- read_xlsx(here("analysis/data/raw_data/SDG_core_metric.xlsx")) %>%
  select(ID, Layer, Core_type_Li_merged, Raw_mat) %>%
  mutate(Layer = as.factor(Layer), Raw_mat = as.factor(Raw_mat))
core_meta <- read_excel(here("analysis/data/raw_data/SDG_core_metric.xlsx")) %>%
  select(ID = ID, raw_material = Raw_mat, core_type = Core_type_Li_merged) %>%
  mutate(across(everything(), ~ str_trim(as.character(.))))

# =============================================================================
# Per-setting analysis blocks (morphology = swept; scar = fixed committed)
# =============================================================================

# ---- (a) morphology degree selection: EXP & SDG ------------------------------
order_block <- function(morph_df) {
  m_exp <- morph_df %>% filter(str_starts(ID, "EXP"))
  m_sdg <- morph_df %>% filter(str_starts(ID, "SDG"), !str_starts(ID, "IM_"))
  cols  <- intersect(POWER_COLS_ALL, colnames(morph_df))
  list(EXP = compute_order_stats(m_exp, cols), SDG = compute_order_stats(m_sdg, cols))
}
cumpower_at <- function(os, l) os$cumul_pct[os$order == l]
maxcv_through <- function(os, l) max(os$cv_pct[os$order <= l], na.rm = TRUE)

# ---- (c) EXP morphology PERMANOVA `perm_morph` (spharm_analysis.R:121,387) ----
exp_permanova_block <- function(morph_df, power_cols = POWER_COLS_MORPH) {
  filt <- filter_spharm(morph_df, power_cols, metric_data)
  df_exp <- split_by_group(filt)$exp_im
  df_exp_only <- df_exp %>%
    filter(!str_starts(ID, "IM_"), !Typology %in% EXCLUDE_TYPES) %>%
    mutate(Typology = case_when(Typology %in% LEVALLOIS_MERGE ~ "Levallois", TRUE ~ Typology),
           Typology = droplevels(as.factor(Typology)))
  non_im_idx <- !str_starts(df_exp$ID, "IM_") & !df_exp$Typology %in% EXCLUDE_TYPES
  ilr_morph <- make_ilr(df_exp[non_im_idx, power_cols])   # ILR / Aitchison (was z-score)
  pm <- run_permanova_dir(ilr_morph, df_exp_only$Typology)
  pv  <- pm$pairwise; sig <- which(pv < 0.05, arr.ind = TRUE)
  sig_pairs <- if (nrow(sig) == 0) character(0) else
    apply(sig, 1, function(rc) paste(rownames(pv)[rc[1]], colnames(pv)[rc[2]], sep = "-"))
  list(R2 = pm$R2, F = pm$F, p = pm$p, n_sig = length(sig_pairs),
       n_pairs = sum(!is.na(pv)), sig_pairs = paste(sort(sig_pairs), collapse = "; "))
}

# ---- Build ILR morph/scar distances. `restrict` selects the per-endpoint specimen
# ---- set BEFORE the zero-variance column drop + ILR, exactly as the main pipeline:
# ---- EXP decoupling uses the exp_im split (exp_cores_statistics.R:162-206); SDG uses
# ---- the full EXP+SDG+IM set (SDG_cores_statistics.R). Sharing one full-set ILR for
# ---- both endpoints changes the EXP column drop and breaks the EXP Mantel/RV anchor.
build_ilr_dists <- function(morph_df, restrict = NULL) {
  # morph_df already carries Typology (joined once in the sweep loop).
  morph_f <- filter_spharm(morph_df,         POWER_COLS_MORPH, metric_data)
  scar_f  <- filter_spharm(SPHARM_direction,  POWER_COLS_DIR,   metric_data)
  if (!is.null(restrict)) {
    morph_f <- split_by_group(morph_f)[[restrict]]
    scar_f  <- split_by_group(scar_f)[[restrict]]
  }
  common <- intersect(morph_f$ID, scar_f$ID)
  morph_f <- morph_f %>% filter(ID %in% common) %>% arrange(ID)
  scar_f  <- scar_f  %>% filter(ID %in% common) %>% arrange(ID)
  mp <- morph_f %>% select(all_of(POWER_COLS_MORPH)) %>%
    rename_with(~ paste0("M", seq_along(.))) %>% as.data.frame()
  sp <- scar_f %>% select(all_of(POWER_COLS_DIR)) %>%
    rename_with(~ paste0("S", seq_along(.))) %>% as.data.frame()
  rownames(mp) <- morph_f$ID; rownames(sp) <- scar_f$ID
  mp <- mp[, sapply(mp, sd, na.rm = TRUE) > 0]; sp <- sp[, sapply(sp, sd, na.rm = TRUE) > 0]
  morph_ilr <- as.data.frame(ilr(replace_zeros(as.matrix(mp))))
  scar_ilr  <- as.data.frame(ilr(replace_zeros(as.matrix(sp))))
  rownames(morph_ilr) <- rownames(mp); rownames(scar_ilr) <- rownames(sp)
  list(morph_ilr = morph_ilr, scar_ilr = scar_ilr,
       D_morph = stats::dist(morph_ilr), D_scar = stats::dist(scar_ilr))
}

mantel_rv <- function(morph_ilr, scar_ilr, D_morph, D_scar, ids) {
  Dm <- extract_subdist(D_morph, ids); Ds <- extract_subdist(D_scar, ids)
  mg <- mantel(Dm, Ds, method = "spearman", permutations = 9999)
  mi <- morph_ilr[ids, ]; si <- scar_ilr[ids, ]
  colnames(mi) <- paste0("M_ilr", seq_len(ncol(mi))); colnames(si) <- paste0("S_ilr", seq_len(ncol(si)))
  dm <- dudi.pca(mi, center = TRUE, scale = TRUE, scannf = FALSE, nf = ncol(mi))
  ds <- dudi.pca(si, center = TRUE, scale = TRUE, scannf = FALSE, nf = ncol(si))
  co <- coinertia(dm, ds, scannf = FALSE, nf = 2)
  set.seed(42); rv <- randtest(co, nrepet = 9999)
  list(mantel_r = mg$statistic, mantel_p = mg$signif, RV = co$RV, RV_p = rv$pvalue, n = length(ids))
}

# ---- (e) EXP decoupling -------------------------------------------------------
exp_decoupling_block <- function(ilr) {
  exp_ids <- rownames(ilr$morph_ilr)[!str_starts(rownames(ilr$morph_ilr), "IM_") &
                                       rownames(ilr$morph_ilr) != "EXP43_Biface"]
  mantel_rv(ilr$morph_ilr, ilr$scar_ilr, ilr$D_morph, ilr$D_scar, exp_ids)
}

# ---- (d)+(e) SDG decoupling + morphology PERMANOVA terms ---------------------
sdg_block <- function(ilr) {
  arch_ids <- rownames(ilr$morph_ilr)[!str_starts(rownames(ilr$morph_ilr), "IM_") &
                                        !str_starts(rownames(ilr$morph_ilr), "EXP")]
  meta_arch <- tibble(ID = arch_ids) %>%
    mutate(layer = case_when(str_detect(ID, "L2") ~ "Layer 2", str_detect(ID, "L3") ~ "Layer 3",
                             str_detect(ID, "L4") ~ "Layer 4", TRUE ~ "Other")) %>%
    left_join(core_meta, by = "ID") %>%
    filter(!core_type %in% EXCLUDE_CORE_TYPES | is.na(core_type))
  arch_ids <- meta_arch$ID
  D_morph_arch <- extract_subdist(ilr$D_morph, arch_ids)

  dec <- mantel_rv(ilr$morph_ilr, ilr$scar_ilr, ilr$D_morph, ilr$D_scar, arch_ids)

  meta_layer    <- safe_filter_groups(meta_arch %>% filter(layer != "Other"), "layer")
  meta_rawmat   <- safe_filter_groups(meta_arch, "raw_material")
  meta_coretype <- safe_filter_groups(meta_arch, "core_type")
  pj <- function(meta, col) {
    if (is.null(meta)) return(list(R2 = NA_real_, F = NA_real_, p = NA_real_))
    run_permanova_dist(extract_subdist(D_morph_arch, meta$ID), meta[[col]])
  }
  ct <- pj(meta_coretype, "core_type"); rm_ <- pj(meta_rawmat, "raw_material"); ly <- pj(meta_layer, "layer")
  list(mantel_r = dec$mantel_r, mantel_p = dec$mantel_p, RV = dec$RV, RV_p = dec$RV_p, n = dec$n,
       coretype_R2 = ct$R2, coretype_F = ct$F, coretype_p = ct$p,
       rawmat_R2 = rm_$R2, rawmat_p = rm_$p, layer_R2 = ly$R2, layer_p = ly$p)
}

# =============================================================================
# Sweep over settings
# =============================================================================
metrics <- list(); order_long <- list()

for (st in SETTINGS) {
  face <- st$face; smooth <- st$smooth
  key  <- sprintf("f%d_s%d", face, smooth)
  csv  <- spectra_path(face, smooth)
  if (!file.exists(csv)) {
    warning(sprintf("Missing spectra %s; run sweep_mesh_spharm.py first. Skipping.", basename(csv)))
    next
  }
  cat(sprintf("\n=== faces=%d, smooth=%d ===\n", face, smooth))
  morph_df <- read_csv(csv, show_col_types = FALSE); morph_df$ID <- str_trim(morph_df$ID)
  # Attach Typology (a fixed attribute) from the committed scar file once, so the
  # blocks that need it (filter_spharm) work; the spectra CSV has no Typology column.
  morph_df <- morph_df %>% left_join(TYPOLOGY_MAP, by = "ID")

  os  <- order_block(morph_df)
  ep  <- exp_permanova_block(morph_df)
  ilr_exp <- build_ilr_dists(morph_df, restrict = "exp_im")  # EXP set (exp_cores_statistics.R)
  ilr_all <- build_ilr_dists(morph_df)                       # SDG set (full EXP+SDG+IM)
  ed  <- exp_decoupling_block(ilr_exp)
  sb  <- sdg_block(ilr_all)

  order_long[[key]] <- bind_rows(os$EXP %>% mutate(dataset = "EXP"),
                                 os$SDG %>% mutate(dataset = "SDG")) %>%
    mutate(face_target = face, smooth_iters = smooth, setting = key)

  metrics[[key]] <- tibble(
    face_target = face, smooth_iters = smooth, setting = key,
    is_production = (face == PROD_FACES & smooth == PROD_SMOOTH),
    exp_morph_cumpower_l8 = cumpower_at(os$EXP, 8), exp_morph_maxcv_l8 = maxcv_through(os$EXP, 8),
    sdg_morph_cumpower_l8 = cumpower_at(os$SDG, 8), sdg_morph_maxcv_l8 = maxcv_through(os$SDG, 8),
    exp_morph_perm_R2 = ep$R2, exp_morph_perm_F = ep$F, exp_morph_perm_p = ep$p,
    exp_morph_nsig = ep$n_sig, exp_morph_npairs = ep$n_pairs, exp_morph_sig_pairs = ep$sig_pairs,
    sdg_morph_coretype_R2 = sb$coretype_R2, sdg_morph_coretype_F = sb$coretype_F, sdg_morph_coretype_p = sb$coretype_p,
    sdg_morph_rawmat_R2 = sb$rawmat_R2, sdg_morph_rawmat_p = sb$rawmat_p,
    sdg_morph_layer_R2 = sb$layer_R2, sdg_morph_layer_p = sb$layer_p,
    exp_mantel_r = ed$mantel_r, exp_mantel_p = ed$mantel_p, exp_RV = ed$RV, exp_RV_p = ed$RV_p, exp_n = ed$n,
    sdg_mantel_r = sb$mantel_r, sdg_mantel_p = sb$mantel_p, sdg_RV = sb$RV, sdg_RV_p = sb$RV_p, sdg_n = sb$n
  )
}

metrics_df <- bind_rows(metrics)
if (nrow(metrics_df) == 0) stop("No spectra found. Run sweep_mesh_spharm.py first.")
order_long_df <- bind_rows(order_long)

# =============================================================================
# Anchor check at production (20000, 3) vs committed values
# =============================================================================
cat("\n", strrep("=", 70), "\n", sep = "")
cat("ANCHOR CHECK at production (faces=", PROD_FACES, ", smooth=", PROD_SMOOTH,
    ") vs committed outputs\n", sep = "")
cat(strrep("=", 70), "\n", sep = "")
pr <- metrics_df %>% filter(is_production)
anchor_ok <- TRUE
if (nrow(pr) == 1) {
  ck <- function(name, got, exp, tol) {
    ok <- is.finite(got) && abs(got - exp) <= tol
    if (!ok) anchor_ok <<- FALSE
    cat(sprintf("  %-32s got %10.5f  expected %10.5f  %s\n", name, got, exp,
                ifelse(ok, "OK", "<-- CHECK")))
  }
  ck("EXP morph cumpower @ l8 (%)", pr$exp_morph_cumpower_l8, REF$exp_morph_cumpower_l8, 0.05)
  ck("EXP morph max CV (l1-8) (%)", pr$exp_morph_maxcv_l8,    REF$exp_morph_maxcv_l8,    0.05)
  ck("SDG morph cumpower @ l8 (%)", pr$sdg_morph_cumpower_l8, REF$sdg_morph_cumpower_l8, 0.05)
  ck("SDG morph max CV (l1-8) (%)", pr$sdg_morph_maxcv_l8,    REF$sdg_morph_maxcv_l8,    0.05)
  ck("SDG morph~coretype R2",      pr$sdg_morph_coretype_R2,  REF$sdg_morph_coretype_R2, 1e-3)
  ck("SDG morph~coretype F",       pr$sdg_morph_coretype_F,   REF$sdg_morph_coretype_F,  5e-3)
  ck("SDG morph~rawmaterial R2",   pr$sdg_morph_rawmat_R2,    REF$sdg_morph_rawmat_R2,   1e-3)
  ck("SDG morph~layer R2",         pr$sdg_morph_layer_R2,     REF$sdg_morph_layer_R2,    1e-3)
  ck("EXP Mantel r",               pr$exp_mantel_r,           REF$exp_mantel_r,          1e-3)
  ck("EXP RV",                     pr$exp_RV,                 REF$exp_RV,                1e-3)
  ck("SDG Mantel r",               pr$sdg_mantel_r,           REF$sdg_mantel_r,          1e-3)
  ck("SDG RV",                     pr$sdg_RV,                 REF$sdg_RV,                1e-3)
  cat(sprintf("\n  => %s\n", ifelse(anchor_ok,
      "Deterministic statistics reproduce the cached production values.",
      "MISMATCH: production differs from cache — investigate before trusting the sweep.")))
  cat("  (EXP morph PERMANOVA R2 is reported but has no committed scalar; it is anchored\n")
  cat("   via the Python-side spectra reproduction of SPHARM_morphology.csv.)\n")
  cat("  (Permutation p-values jitter ~+/-0.005 and are not part of the check.)\n")
} else cat("  Production setting not in the grid; cannot run the anchor check.\n")

# =============================================================================
# Outputs: tidy CSVs
# =============================================================================
write_csv(metrics_df,    file.path(OUT_DIR, "mesh_sensitivity_metrics.csv"))
write_csv(order_long_df, file.path(OUT_DIR, "mesh_orderselection_by_degree.csv"))
cat("\nWrote mesh_sensitivity_metrics.csv and mesh_orderselection_by_degree.csv\n")

# =============================================================================
# Figures
# =============================================================================
tryCatch({
  setting_lab <- function(face, smooth) sprintf("%dk / s%d", face %/% 1000, smooth)
  order_long_df <- order_long_df %>%
    mutate(lab = setting_lab(face_target, smooth_iters),
           axis = ifelse(smooth_iters == PROD_SMOOTH, "vary faces (smooth = prod)",
                         ifelse(face_target == PROD_FACES, "vary smoothing (faces = prod)", "cross")))
  pal <- setNames(colorRampPalette(c("#4A6E8A", "#5C7F71", "#BA8530", "#802520"))(
                    length(unique(order_long_df$lab))), sort(unique(order_long_df$lab)))

  # Figure 1: order selection (cumulative power & CV by degree), EXP & SDG.
  p_cum <- ggplot(order_long_df %>% filter(order <= 12),
                  aes(order, cumul_pct, color = lab, group = lab)) +
    geom_hline(yintercept = c(95, 98), linetype = "dashed", color = "grey75", linewidth = 0.3) +
    geom_vline(xintercept = 8, linetype = "dotted", color = "grey50", linewidth = 0.3) +
    geom_line(linewidth = 0.6) + geom_point(size = 1.1) +
    facet_wrap(~ dataset, labeller = as_labeller(
      c(EXP = "Experimentally knapped cores", SDG = "Sandinggai cores"))) +
    scale_color_manual(values = pal, name = "faces/smooth") + scale_x_continuous(breaks = 1:12) +
    labs(x = "SPHARM degree (l)", y = "Cumulative power (%)") +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank(),
          strip.text = element_text(face = "bold"))
  p_cv <- ggplot(order_long_df %>% filter(order <= 12),
                 aes(order, cv_pct, color = lab, group = lab)) +
    geom_hline(yintercept = 90, linetype = "dashed", color = "grey50", linewidth = 0.4) +
    geom_vline(xintercept = 8, linetype = "dotted", color = "grey50", linewidth = 0.3) +
    geom_line(linewidth = 0.6) + geom_point(size = 1.1) +
    facet_wrap(~ dataset, labeller = as_labeller(
      c(EXP = "Experimentally knapped cores", SDG = "Sandinggai cores"))) +
    scale_color_manual(values = pal, name = "faces/smooth") + scale_x_continuous(breaks = 1:12) +
    labs(x = "SPHARM degree (l)", y = "Cross-specimen CV (%)") +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank(),
          strip.text = element_text(face = "bold"))
  ggsave(file.path(FIG_DIR, "fig_S_mesh_orderselection.png"), p_cum / p_cv,
         width = 9, height = 8, dpi = 300)

  cat("Wrote 1 figure to", FIG_DIR, "\n")
}, error = function(e) cat("Figure generation failed (non-critical):", conditionMessage(e), "\n"))

# =============================================================================
# Stability flags (printed by the SENSITIVITY FLAGS section below)
# =============================================================================
exp_perm_all_sig <- all(metrics_df$exp_morph_perm_p < 0.05, na.rm = TRUE)
sdg_ct_all_sig   <- all(metrics_df$sdg_morph_coretype_p < 0.05, na.rm = TRUE)
rawmat_all_ns    <- all(metrics_df$sdg_morph_rawmat_p >= 0.05, na.rm = TRUE)
layer_all_ns     <- all(metrics_df$sdg_morph_layer_p  >= 0.05, na.rm = TRUE)
mantel_all_ns    <- all(metrics_df$exp_mantel_p >= 0.05, na.rm = TRUE) &&
                    all(metrics_df$sdg_mantel_p >= 0.05, na.rm = TRUE)
rv_all_ns        <- all(metrics_df$exp_RV_p >= 0.05, na.rm = TRUE) &&
                    all(metrics_df$sdg_RV_p >= 0.05, na.rm = TRUE)
cv_ok            <- all(metrics_df$exp_morph_maxcv_l8 < 90, na.rm = TRUE) &&
                    all(metrics_df$sdg_morph_maxcv_l8 < 90, na.rm = TRUE)

# ---- flag any metric sensitive to decimation/smoothing ----------------------
cat("\n", strrep("=", 70), "\n", sep = "")
cat("SENSITIVITY FLAGS (metrics whose conclusion changes across settings)\n")
cat(strrep("=", 70), "\n", sep = "")
flag_any <- FALSE
if (!exp_perm_all_sig) { flag_any <- TRUE; cat("  * EXP morph PERMANOVA: p crosses 0.05 across settings.\n") }
if (!sdg_ct_all_sig)   { flag_any <- TRUE; cat("  * SDG morph core-type PERMANOVA: p crosses 0.05.\n") }
if (!rawmat_all_ns)    { flag_any <- TRUE; cat("  * SDG raw-material term: becomes significant at some setting.\n") }
if (!layer_all_ns)     { flag_any <- TRUE; cat("  * SDG layer term: becomes significant at some setting.\n") }
if (!mantel_all_ns)    { flag_any <- TRUE; cat("  * Global Mantel: significant at some setting.\n") }
if (!rv_all_ns)        { flag_any <- TRUE; cat("  * RV: significant at some setting.\n") }
if (!cv_ok)            { flag_any <- TRUE; cat("  * Morph CV exceeds 90% within l1-8 at some setting.\n") }
if (!flag_any) cat("  None: every morphology conclusion is stable across the tested settings.\n")

cat("\nDone.\n")
