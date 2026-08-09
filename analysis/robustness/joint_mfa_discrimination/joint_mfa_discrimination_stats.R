# joint_mfa_discrimination_stats.R
# =============================================================================
# JOINT (M-SPHARM + SP-SPHARM) DISCRIMINANT ANALYSIS via MFA block normalisation
# — SI add-on. NEW, self-contained file. Does NOT modify the main pipeline, the
# cached _targets store, the derived_data cache, or the manuscript. It only READS
# the committed main outputs (SPHARM_direction.csv / SPHARM_morphology.csv +
# raw metadata) and WRITES new outputs under
# analysis/robustness/joint_mfa_discrimination/.
#
# The main analysis tests morphology (M-SPHARM, l = 1-8) and scar patterning
# (SP-SPHARM, l = 1-6) as two SEPARATE descriptor blocks. This script asks the
# obvious follow-up: does merging the two blocks into a single MFA-balanced
# feature space buy any discriminant power that neither block has alone —
# specifically for the two core-type pairs the scar block fails to separate
# (discoid-bidirectional, discoid-multiplatform)?
#
# It re-uses the project's statistical machinery (the same package functions the
# main pipeline calls — vegan::adonis2 / betadisper / permutest,
# compositions::ilr) and replicates, verbatim, the data-prep from:
#   - r_spharm/spharm_analysis.R           (EXP PERMANOVA `perm_dir` / `perm_morph`)
#   - r_statistics/SDG_cores_statistics.R  (SDG core-type PERMANOVA; parameterised
#                                           branch, not exercised by this run)
# Source line numbers are cited next to each replicated block.
#
# Steps:
#   (1) ILR each power spectrum separately  -> Z_M (n x 7), Z_SP (n x 5).
#   (2) MFA block normalisation: divide each centred block by its own first
#       singular value -> Z_comb = cbind(Z_M / s_M, Z_SP / s_SP).
#   (3) PERMANOVA (global + all pairwise, Holm) for M alone / SP alone / combined,
#       always on the FULL Z distance matrix, never on truncated MFA axes: with
#       all axes retained the MFA is an orthogonal rotation and distances are
#       unchanged, so the three models are directly comparable; truncating would
#       silently change the test.
#   (4) PERMDISP (betadisper + permutest) on the combined distance matrix.
#   (5) MFA biplot: global PCA of Z_comb, one point per core, 95% ellipses, and
#       variable arrows mapped back from ILR to CLR (clr_loading = V %*% ilr_loading,
#       V = the ILR contrast basis) so each arrow corresponds to one harmonic
#       degree instead of an uninterpretable between-degree contrast.
#   (6) Weight scan: w = 0 (pure SP) .. 1 (pure M) in steps of 0.05 under TOTAL-
#       INERTIA (Frobenius) normalisation, so w is exactly the morphology share of
#       squared distance; coordinates carry sqrt(w) because PERMANOVA partitions
#       squared distances. The curve reports the number of Holm-significant pairs,
#       NOT R^2: under Euclidean geometry squared distances are additive, so R^2(w)
#       is a linear-fractional function of w, strictly monotone on [0, 1] and pinned
#       at both ends by the two single-block values — the whole curve is determined
#       analytically by two numbers already in the table. For the same reason the
#       combined R^2 can never exceed max(R2_M, R2_SP); a drop is arithmetic, not
#       evidence that merging failed.
#
# No LOOCV / classification accuracy is attempted here, so the scan reports
# resolved pairs only.
#
# Outputs (all NEW):
#   joint_mfa_permanova_comparison.csv   global + pairwise results, 3 models
#   joint_mfa_weight_scan.csv            one row per w
#   figures/fig_joint_mfa_main.png       MAIN TEXT composite (a biplot | b scan / c curve)
#   figures/fig_S_joint_mfa_biplot.png                    source panel a
#   figures/fig_S_joint_mfa_weight_scan.png               source panel b
#   figures/fig_S_joint_mfa_weight_resolution_continuous.png  source panel c
#
# HOW TO RUN (canonical environment, Docker spharm_analysis, R 4.4 + renv):
#   Rscript analysis/robustness/joint_mfa_discrimination/joint_mfa_discrimination_stats.R
#   # or in RStudio: source(here::here("analysis/robustness/joint_mfa_discrimination/joint_mfa_discrimination_stats.R"))
# =============================================================================

suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(readxl)
  library(vegan)          # adonis2, betadisper, permutest
  library(compositions)   # ilr, ilrBase, clr
  library(ggrepel)
  library(patchwork)      # main-text composite
  library(conflicted)
})

# Resolve namespace clashes exactly as the main pipeline does (_targets.R).
# compositions masks several base/stats generics (var, cor, dist, scale, ...).
suppressMessages({
  conflicts_prefer(dplyr::filter, dplyr::select, dplyr::lag,
                   stats::sd, stats::var, stats::dist, stats::cor, stats::cov,
                   base::scale, base::norm, base::`%*%`, .quiet = TRUE)
})

set.seed(42)

# =============================================================================
# PARAMETERS
# =============================================================================
ASSEMBLAGE <- "EXP"        # "EXP" (this run) or "SDG"; see build_blocks()
SEED       <- 42
N_PERM     <- 9999         # global + pairwise PERMANOVA, PERMDISP
W_GRID     <- seq(0, 1, by = 0.05)   # weight scan; w = 1 pure M, w = 0 pure SP
# Dense grid for diagnostic D3: 0.01 resolution across 0.61-0.84, where the coarse
# scan drops 4 pairs in a single 0.05 step. A superset of W_GRID, so the coarse
# rows are reproduced exactly (results at a given w do not depend on the grid).
W_GRID_DENSE <- sort(unique(round(c(seq(0, 0.60, by = 0.05),
                                    seq(0.61, 0.84, by = 0.01),
                                    seq(0.85, 1.00, by = 0.05)), 10)))
DOMINANCE_THR <- 0.65   # D2: block share at which an axis counts as block-dominated
DEGENERACY_THR <- 0.10  # D2: relative inertia gap below which axes 1-2 are degenerate

OUT_DIR <- here("analysis/robustness/joint_mfa_discrimination")
FIG_DIR <- file.path(OUT_DIR, "figures")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

# Power-column conventions (identical to the main analysis).
POWER_COLS_DIR   <- paste0("power_l", 1:6)   # scar       (SP-SPHARM) l = 1-6
POWER_COLS_MORPH <- paste0("power_l", 1:8)   # morphology (M-SPHARM)  l = 1-8

# Typology handling (spharm_analysis.R:53-67)
EXCLUDE_TYPES   <- c("Biface")
LEVALLOIS_MERGE <- c("Levallois convergent", "Levallois laminar",
                     "Levallois preferential", "Levallois recurrent")
TYPOLOGY_ORDER  <- c("Unidirectional", "Bidirectional", "Levallois",
                     "Discoid", "Multiplatform")
TYPOLOGY_COLORS <- c(
  "Levallois"      = "#4A6E8A",
  "Discoid"        = "#802520",
  "Unidirectional" = "#BA8530",
  "Multiplatform"  = "#8A7A68",
  "Bidirectional"  = "#788C4A"
)

# SDG metadata handling (SDG_cores_statistics.R:38-49), for the parameterised branch.
EXCLUDE_CORE_TYPES <- c("Handaxe", "Pick")

# The two pairs the scar block alone fails to separate (manuscript.qmd:319) — the
# reason this analysis exists. Matched on membership, order-insensitive.
FOCUS_PAIRS <- list(c("Discoid", "Bidirectional"),
                    c("Discoid", "Multiplatform"))

# Committed production reference values for the EXP anchor check. R2 / pseudo-F are
# deterministic; permutation p-values jitter by ~+/-0.005 and are reported, not
# checked. Sources: mesh_sensitivity_metrics.csv (production row 20000/3) and
# threshold_sensitivity_metrics.csv (threshold 0 = all scars), which are themselves
# anchored to `perm_morph` / `perm_dir` in spharm_analysis.R.
REF <- list(
  M_R2  = 0.12334, M_F  = 1.86418, M_p  = 0.025, M_nsig  = 0,
  SP_R2 = 0.30174, SP_F = 5.72568, SP_p = 0.001, SP_nsig = 8,
  n     = 58
)

# =============================================================================
# Helpers — copied VERBATIM from the source scripts (attribution in comments)
# =============================================================================

# spharm_analysis.R:119-133 / exp_cores_statistics.R:65-79
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

# spharm_analysis.R:134-138 — drop zero-variance parts, multiplicative zero
# replacement, then ilr into Euclidean space.
make_ilr <- function(power_df) {
  X    <- as.matrix(power_df)
  keep <- apply(X, 2, function(v) sd(v, na.rm = TRUE) > 0)
  as.matrix(ilr(replace_zeros(X[, keep, drop = FALSE])))
}

# spharm_analysis.R:87-92
filter_spharm <- function(df, power_cols, meta = NULL) {
  result <- df %>% select(ID, Typology, all_of(power_cols))
  if (!is.null(meta)) result <- left_join(result, meta, by = "ID")
  result
}

# spharm_analysis.R:97-102
split_by_group <- function(df) list(
  exp_im = df %>% filter(str_starts(ID, "EXP") | str_starts(ID, "IM_")),
  sdg_im = df %>% filter(str_starts(ID, "SDG") | str_starts(ID, "IM_")))

# SDG_cores_statistics.R:96-110
safe_filter_groups <- function(meta_df, group_col, min_n = 3) {
  counts <- table(meta_df[[group_col]], useNA = "no")
  valid  <- names(counts[counts >= min_n])
  if (length(valid) < 2) return(NULL)
  meta_df %>% filter(!is.na(.data[[group_col]]), .data[[group_col]] %in% valid)
}

# =============================================================================
# New helpers (no existing equivalent in the repo)
# =============================================================================

# Which parts survive make_ilr()'s zero-variance drop — needed to build the
# matching ILR contrast basis for the CLR back-mapping.
ilr_parts <- function(power_df) {
  X <- as.matrix(power_df)
  colnames(X)[apply(X, 2, function(v) sd(v, na.rm = TRUE) > 0)]
}

# MFA block weight: first singular value of the centred block.
mfa_s1 <- function(Z) svd(base::scale(as.matrix(Z), scale = FALSE))$d[1]

# Total inertia of a centred block = squared Frobenius norm.
block_inertia <- function(Z) sum(base::scale(as.matrix(Z), scale = FALSE)^2)

# =============================================================================
# Data preparation — one block pair + a grouping factor, per assemblage
# =============================================================================
cat("Loading committed main outputs (SPHARM direction / morphology, metadata)...\n")

SPHARM_direction  <- read_csv(
  here("analysis/data/derived_data/SPHARM_direction.csv"),  show_col_types = FALSE)
SPHARM_morphology <- read_csv(
  here("analysis/data/derived_data/SPHARM_morphology.csv"), show_col_types = FALSE)
SPHARM_morphology <- SPHARM_morphology %>%                     # spharm_analysis.R:80-81
  left_join(SPHARM_direction %>% select(ID, Typology), by = "ID")

metric_data <- read_xlsx(here("analysis/data/raw_data/SDG_core_metric.xlsx")) %>%
  select(ID, Layer, Core_type_Li_merged, Raw_mat) %>%          # spharm_analysis.R:76-78
  mutate(Layer = as.factor(Layer), Raw_mat = as.factor(Raw_mat))
core_meta <- read_excel(here("analysis/data/raw_data/SDG_core_metric.xlsx")) %>%
  select(ID = ID, raw_material = Raw_mat, core_type = Core_type_Li_merged) %>%
  mutate(across(everything(), ~ str_trim(as.character(.))))    # SDG_cores_statistics.R:171-177

# Returns: morph / scar POWER matrices (rows = specimens, aligned), the grouping
# factor, and the specimen IDs. Everything downstream is assemblage-agnostic.
build_blocks <- function(assemblage = ASSEMBLAGE) {
  dir_f   <- filter_spharm(SPHARM_direction,  POWER_COLS_DIR,   metric_data)
  morph_f <- filter_spharm(SPHARM_morphology, POWER_COLS_MORPH, metric_data)

  if (assemblage == "EXP") {
    # spharm_analysis.R:113-157 — EXP+IM split, drop IM_ and Biface, merge Levallois.
    df_dir   <- split_by_group(dir_f)$exp_im
    df_morph <- split_by_group(morph_f)$exp_im
    if (!identical(df_dir$ID, df_morph$ID))
      stop("EXP direction / morphology frames are not row-aligned; ",
           "the main pipeline indexes both with one `non_im_idx` (spharm_analysis.R:150). ",
           "Stopping rather than guessing an alignment.")

    df_only <- df_dir %>%
      filter(!str_starts(ID, "IM_"), !Typology %in% EXCLUDE_TYPES) %>%
      mutate(Typology = case_when(Typology %in% LEVALLOIS_MERGE ~ "Levallois",
                                  TRUE ~ Typology),
             Typology = droplevels(as.factor(Typology)))
    keep_idx <- !str_starts(df_dir$ID, "IM_") & !df_dir$Typology %in% EXCLUDE_TYPES

    list(ids   = df_dir$ID[keep_idx],
         morph = df_morph[keep_idx, POWER_COLS_MORPH],
         scar  = df_dir[keep_idx,   POWER_COLS_DIR],
         group = factor(as.character(df_only$Typology),
                        levels = TYPOLOGY_ORDER[TYPOLOGY_ORDER %in%
                                                  unique(as.character(df_only$Typology))]))

  } else if (assemblage == "SDG") {
    # SDG_cores_statistics.R:148-230 — the UN-split filtered frames (full
    # EXP+SDG+IM set), then restrict to archaeological specimens and drop the
    # excluded core types. NOTE: parameterised for a future run; not exercised
    # by this script's EXP-only execution.
    common   <- intersect(dir_f$ID, morph_f$ID)
    df_dir   <- dir_f   %>% filter(ID %in% common) %>% arrange(ID)
    df_morph <- morph_f %>% filter(ID %in% common) %>% arrange(ID)
    stopifnot(identical(df_dir$ID, df_morph$ID))

    arch_idx  <- !str_starts(df_dir$ID, "IM_") & !str_starts(df_dir$ID, "EXP")
    meta_arch <- tibble(ID = df_dir$ID[arch_idx]) %>%
      left_join(core_meta, by = "ID") %>%
      filter(!core_type %in% EXCLUDE_CORE_TYPES | is.na(core_type))
    meta_arch <- safe_filter_groups(meta_arch, "core_type")
    if (is.null(meta_arch))
      stop("SDG: fewer than two core-type groups with n >= 3 (safe_filter_groups).")

    keep_idx <- match(meta_arch$ID, df_dir$ID)
    list(ids   = df_dir$ID[keep_idx],
         morph = df_morph[keep_idx, POWER_COLS_MORPH],
         scar  = df_dir[keep_idx,   POWER_COLS_DIR],
         group = factor(meta_arch$core_type))

  } else {
    stop("ASSEMBLAGE must be 'EXP' or 'SDG'; got ", assemblage)
  }
}

blocks <- build_blocks()
n_spec <- length(blocks$ids)
cat(sprintf("\nAssemblage: %s | specimens: %d | groups: %d\n",
            ASSEMBLAGE, n_spec, nlevels(blocks$group)))
print(table(blocks$group))
if (n_spec == 0) stop("No specimens retained — check the assemblage filters.")

# =============================================================================
# (1) ILR each block; (2) MFA block normalisation
# =============================================================================
parts_M  <- ilr_parts(blocks$morph)
parts_SP <- ilr_parts(blocks$scar)
Z_M  <- make_ilr(blocks$morph)   # n x (|parts_M|  - 1)
Z_SP <- make_ilr(blocks$scar)    # n x (|parts_SP| - 1)
rownames(Z_M) <- rownames(Z_SP) <- blocks$ids
colnames(Z_M)  <- paste0("M_ilr",  seq_len(ncol(Z_M)))
colnames(Z_SP) <- paste0("SP_ilr", seq_len(ncol(Z_SP)))
cat(sprintf("ILR blocks: M-SPHARM %d parts -> %d coords | SP-SPHARM %d parts -> %d coords\n",
            length(parts_M), ncol(Z_M), length(parts_SP), ncol(Z_SP)))

s_M  <- mfa_s1(Z_M)
s_SP <- mfa_s1(Z_SP)
Z_comb <- cbind(Z_M / s_M, Z_SP / s_SP)
cat(sprintf("MFA normalisation: s1(M) = %.5f, s1(SP) = %.5f\n", s_M, s_SP))

# Morphology share of total inertia under the s1 normalisation — where the MFA
# setting sits on the weight axis of step (6).
I_M    <- block_inertia(Z_M  / s_M)
I_SP   <- block_inertia(Z_SP / s_SP)
w_mfa  <- I_M / (I_M + I_SP)
cat(sprintf("Inertia after s1 normalisation: I_M = %.4f, I_SP = %.4f  =>  w_MFA = %.4f\n",
            I_M, I_SP, w_mfa))

# =============================================================================
# PERMANOVA machinery
# =============================================================================
# Global test on the FULL coordinate matrix (Euclidean distances, all axes kept).
permanova_global <- function(Z, group, seed = SEED, nperm = N_PERM) {
  d <- stats::dist(as.matrix(Z), method = "euclidean")
  set.seed(seed)
  res <- adonis2(d ~ type, data = data.frame(type = group), permutations = nperm)
  list(R2 = res$R2[1], F = res$`F`[1], p = res$`Pr(>F)`[1])
}

# All pairwise comparisons: each pair is re-tested on its own subset, then the
# whole family of p-values is Holm-corrected once (p.adjust), matching the
# main analysis's correction strategy (exp_cores_statistics.R:268).
permanova_pairwise <- function(Z, group, seed = SEED, nperm = N_PERM) {
  Z    <- as.matrix(Z)
  lev  <- levels(droplevels(group))
  prs  <- utils::combn(lev, 2, simplify = FALSE)
  out  <- map_dfr(seq_along(prs), function(k) {
    a <- prs[[k]][1]; b <- prs[[k]][2]
    idx <- which(group %in% c(a, b))
    g   <- droplevels(group[idx])
    d   <- stats::dist(Z[idx, , drop = FALSE], method = "euclidean")
    set.seed(seed + k)
    res <- adonis2(d ~ type, data = data.frame(type = g), permutations = nperm)
    tibble(group1 = a, group2 = b, comparison = paste(a, b, sep = "-"),
           n = length(idx), R2 = res$R2[1], pseudo_F = res$`F`[1],
           p = res$`Pr(>F)`[1])
  })
  out %>% mutate(p_holm = p.adjust(p, method = "holm"),
                 significant_holm = p_holm < 0.05)
}

# PERMDISP: homogeneity of multivariate dispersions (spharm_analysis.R:330-333).
permdisp <- function(Z, group, seed = SEED, nperm = N_PERM) {
  d <- stats::dist(as.matrix(Z), method = "euclidean")
  set.seed(seed)
  disp <- betadisper(d, group)
  set.seed(seed)
  pt <- permutest(disp, permutations = nperm)
  list(disp = disp, test = pt,
       F = pt$tab$`F`[1], p = pt$tab$`Pr(>F)`[1])
}

# p for a pair identified by membership (order-insensitive). `col` selects which
# column to read: "p_holm" (default, used throughout the main analysis) or "p" for
# the raw uncorrected value that diagnostic D1 needs.
pair_p <- function(pw, members, col = "p_holm") {
  hit <- pw %>% filter((group1 == members[1] & group2 == members[2]) |
                         (group1 == members[2] & group2 == members[1]))
  if (nrow(hit) != 1) return(NA_real_)
  hit[[col]][1]
}

# =============================================================================
# (3) Three-model PERMANOVA comparison
# =============================================================================
cat("\n", strrep("=", 70), "\n", sep = "")
cat("PERMANOVA: M-SPHARM alone / SP-SPHARM alone / MFA-combined\n")
cat(strrep("=", 70), "\n", sep = "")

MODELS <- list(
  "M-SPHARM"     = Z_M,
  "SP-SPHARM"    = Z_SP,
  "MFA-combined" = Z_comb
)

model_res <- map(names(MODELS), function(nm) {
  cat(sprintf("  running %s (%d coords) ...\n", nm, ncol(MODELS[[nm]])))
  list(global   = permanova_global(MODELS[[nm]], blocks$group),
       pairwise = permanova_pairwise(MODELS[[nm]], blocks$group))
})
names(model_res) <- names(MODELS)

# PERMDISP is reported for the combined space (the model this analysis proposes).
cat("  running PERMDISP on the combined distance matrix ...\n")
pd_comb <- permdisp(Z_comb, blocks$group)

n_pairs <- nrow(model_res[[1]]$pairwise)

comparison_tbl <- map_dfr(names(model_res), function(nm) {
  g  <- model_res[[nm]]$global
  pw <- model_res[[nm]]$pairwise
  tibble(model = nm, n = n_spec, n_coords = ncol(MODELS[[nm]]),
         R2 = g$R2, pseudo_F = g$F, p = g$p,
         n_sig_pairs = sum(pw$significant_holm), n_pairs = n_pairs,
         sig_pairs = paste(sort(pw$comparison[pw$significant_holm]), collapse = "; "),
         permdisp_F = if (nm == "MFA-combined") pd_comb$F else NA_real_,
         permdisp_p = if (nm == "MFA-combined") pd_comb$p else NA_real_)
})

cat("\n--- Three-model comparison ---\n")
print(comparison_tbl %>%
        transmute(model, n, n_coords,
                  R2 = round(R2, 5), pseudo_F = round(pseudo_F, 5), p,
                  resolved = sprintf("%d/%d", n_sig_pairs, n_pairs)) %>%
        as.data.frame())

cat("\n--- PERMDISP (combined distance matrix, betadisper + permutest) ---\n")
print(pd_comb$test)

cat("\n--- Pairwise Holm-corrected p (all models) ---\n")
pairwise_wide <- map_dfr(names(model_res),
                         ~ model_res[[.x]]$pairwise %>% mutate(model = .x)) %>%
  select(comparison, model, p_holm) %>%
  pivot_wider(names_from = model, values_from = p_holm) %>%
  mutate(across(where(is.numeric), ~ round(.x, 4)))
print(as.data.frame(pairwise_wide))

# ---- the two pairs SP-SPHARM alone cannot separate --------------------------
if (ASSEMBLAGE == "EXP") {
  cat("\n--- Focus pairs (the two SP-SPHARM alone fails to resolve) ---\n")
  focus_tbl <- map_dfr(FOCUS_PAIRS, function(pr) {
    tibble(comparison = paste(pr, collapse = "-"),
           `M-SPHARM`     = pair_p(model_res[["M-SPHARM"]]$pairwise, pr),
           `SP-SPHARM`    = pair_p(model_res[["SP-SPHARM"]]$pairwise, pr),
           `MFA-combined` = pair_p(model_res[["MFA-combined"]]$pairwise, pr))
  })
  print(focus_tbl %>% mutate(across(where(is.numeric), ~ round(.x, 4))) %>%
          as.data.frame())
  cat("  (Holm-corrected p; < 0.05 = pair resolved.)\n")
}

# ---- anchor check against the committed single-block results ----------------
cat("\n", strrep("=", 70), "\n", sep = "")
cat("ANCHOR CHECK vs committed single-block PERMANOVA (EXP)\n")
cat(strrep("=", 70), "\n", sep = "")
anchor_ok <- TRUE
if (ASSEMBLAGE == "EXP") {
  check <- function(name, got, expect, tol, deterministic = TRUE) {
    ok <- is.finite(got) && abs(got - expect) <= tol
    if (!ok && deterministic) anchor_ok <<- FALSE
    cat(sprintf("  %-26s got %10.5f  expected %10.5f  %s\n", name, got, expect,
                ifelse(ok, "OK", ifelse(deterministic, "<-- CHECK",
                                        "<-- CHECK (permutation)"))))
  }
  mg <- model_res[["M-SPHARM"]]$global;  mp <- model_res[["M-SPHARM"]]$pairwise
  sg <- model_res[["SP-SPHARM"]]$global; sp <- model_res[["SP-SPHARM"]]$pairwise
  check("n specimens",          n_spec,    REF$n,     0)
  check("M-SPHARM R2",          mg$R2,     REF$M_R2,  1e-3)
  check("M-SPHARM pseudo-F",    mg$F,      REF$M_F,   1e-3)
  check("M-SPHARM p",           mg$p,      REF$M_p,   0.01, deterministic = FALSE)
  check("M-SPHARM sig pairs",   sum(mp$significant_holm),  REF$M_nsig,  0)
  check("SP-SPHARM R2",         sg$R2,     REF$SP_R2, 1e-3)
  check("SP-SPHARM pseudo-F",   sg$F,      REF$SP_F,  1e-3)
  check("SP-SPHARM p",          sg$p,      REF$SP_p,  0.01, deterministic = FALSE)
  check("SP-SPHARM sig pairs",  sum(sp$significant_holm), REF$SP_nsig, 0)
  cat(sprintf("\n  => %s\n", ifelse(anchor_ok,
      "Single-block models reproduce the committed main-text PERMANOVA.",
      "MISMATCH: investigate before trusting the combined model.")))
  cat("  (Permutation p-values jitter by ~+/-0.005 and are informational.)\n")
} else {
  cat("  Anchors are EXP-specific; skipped for", ASSEMBLAGE, "\n")
}

# =============================================================================
# (5) MFA biplot — global PCA of Z_comb
# =============================================================================
# The MFA global analysis IS a PCA of the block-normalised concatenation, so this
# is computed directly on the same Z_comb the PERMANOVA used (guaranteeing the
# figure and the test describe one geometry). Only centring is applied: scaling
# each ILR coordinate to unit variance (FactoMineR's type = "s") would inflate the
# high-order, low-energy coordinates to the same weight as the low-order signal
# and destroy the l = 1-8 / 1-6 truncation design.
pca_comb  <- stats::prcomp(Z_comb, center = TRUE, scale. = FALSE)
inertia   <- pca_comb$sdev^2 / sum(pca_comb$sdev^2) * 100
cat(sprintf("\nMFA global PCA: axis 1 = %.1f%%, axis 2 = %.1f%% of inertia\n",
            inertia[1], inertia[2]))

# Cross-check against FactoMineR::MFA(type = "c"): same normalisation, so the
# eigenvalue percentages must agree (non-critical, skipped if unavailable).
if (requireNamespace("FactoMineR", quietly = TRUE)) {
  mfa_df  <- as.data.frame(cbind(Z_M, Z_SP))
  mfa_fit <- try(FactoMineR::MFA(mfa_df,
                                 group = c(ncol(Z_M), ncol(Z_SP)),
                                 type  = c("c", "c"),
                                 name.group = c("M-SPHARM", "SP-SPHARM"),
                                 graph = FALSE), silent = TRUE)
  if (!inherits(mfa_fit, "try-error")) {
    fm_pct <- mfa_fit$eig[1:2, 2]
    cat(sprintf("  FactoMineR::MFA(type = 'c') cross-check: %.1f%% / %.1f%%  (max diff %.4f pp)\n",
                fm_pct[1], fm_pct[2], max(abs(fm_pct - inertia[1:2]))))
  }
}

# ---- ILR loadings -> CLR loadings -------------------------------------------
# ILR coordinates are contrasts BETWEEN degrees, so an arrow on an ILR axis has no
# reading. The CLR representation has one coordinate per harmonic degree. With V
# the ILR contrast basis (D x (D-1), orthonormal columns; compositions::ilr uses
# clr(x) %*% V), the CLR image of an ILR loading vector a is V %*% a — and because
# V'V = I the arrow lengths are preserved, so M and SP arrows stay comparable.
V_M  <- ilrBase(D = length(parts_M))
V_SP <- ilrBase(D = length(parts_SP))

# Verify V is the basis compositions::ilr actually used: clr = ilr %*% t(V).
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

# Loadings are read in the analysed (MFA-weighted) space, i.e. they are
# contributions to the plotted axes, not raw-ILR direction cosines.
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

# Common arrow scaling: one factor for both blocks, so relative lengths within and
# between blocks survive the rescale to score units.
arrow_scale <- 0.80 * max(abs(as.matrix(scores_df[, c("Axis1", "Axis2")]))) /
  max(sqrt(load_df$Axis1^2 + load_df$Axis2^2))
load_df <- load_df %>% mutate(x = Axis1 * arrow_scale, y = Axis2 * arrow_scale)

group_pal <- if (ASSEMBLAGE == "EXP")
  TYPOLOGY_COLORS[levels(blocks$group)] else
    setNames(colorRampPalette(unname(TYPOLOGY_COLORS))(nlevels(blocks$group)),
             levels(blocks$group))

# -----------------------------------------------------------------------------
# MAIN-TEXT FIGURE SPECIFICATION -- taken verbatim, not invented here
# -----------------------------------------------------------------------------
# This figure goes in the main text, so its typography follows the manuscript's
# existing main-text panels rather than the looser SI settings this script used
# while it was an SI-only analysis.
#   theme / legend / zero-line    <- SDG_cores_statistics.R make_coia_biplot()
#                                    (manuscript Fig. 9, fig-coia-composite)
#   point size                    <- spharm_analysis.R p_dir_plot (Fig. 7c), the
#                                    closest analogue: one point per core inside
#                                    convex hulls, in a panel of similar width.
#                                    NOT Fig. 8a's 2.0 -- that is a 5.14 in panel,
#                                    this one is 3.71 in (widths 1.18 : 1), so a
#                                    2.0 dot reads nearly twice the area of the
#                                    dots in the equally sized Fig. 9 panels.
#   export geometry               <- manuscript.qmd fig-width 6.85 / out-width
#                                    174mm / fig-dpi 800
MT <- list(base_size = 8, axis_text = 5, legend_text = 6.5, legend_title = 7,
           legend_key_cm = 0.30, tag_size = 9, tag_face = "bold",
           pt_shape = 16, pt_size = 1.4, pt_stroke = 0.25, pt_alpha = 0.90,
           zero_col = "grey70", zero_lwd = 0.25,
           fig_width_in = 6.85, dpi = 800)

mt_theme <- function() {
  theme_bw(base_size = MT$base_size) +
    theme(panel.grid   = element_blank(),
          axis.text    = element_text(size = MT$axis_text),
          legend.key.size = grid::unit(MT$legend_key_cm, "cm"),
          legend.text  = element_text(size = MT$legend_text),
          legend.title = element_text(size = MT$legend_title),
          legend.margin = margin(2, 4, 2, 4),
          plot.tag     = element_text(size = MT$tag_size, face = MT$tag_face))
}

# Reference line on the two weight panels. Placed at exactly 0.5 (equal
# weighting) rather than at the data-derived w_mfa. The s1 normalisation lands
# within a hundredth of 0.5 here, so the two are indistinguishable at this
# scale, and a round 0.5 is the honest position for an UNANNOTATED line. The
# realised gap is printed to the console below; what the line marks must be
# stated in the caption, since the panel no longer says so.
W_REF <- 0.5

# Convex hulls, not normal ellipses: the main-text ordination panels
# (spharm_analysis.R run_lda_plot(), Figure 7c) draw group extents as filled
# chull polygons. An ellipse is a parametric summary that extends well beyond
# the observed points for the small, non-elliptical groups here (Discoid n = 6),
# so the hull is both the honest extent and the house style.
hull_df <- scores_df %>%
  dplyr::group_by(group) %>%
  dplyr::slice(chull(Axis1, Axis2)) %>%
  dplyr::ungroup()

p_biplot <- ggplot() +
  geom_hline(yintercept = 0, color = MT$zero_col, linewidth = MT$zero_lwd,
             linetype = "dashed") +
  geom_vline(xintercept = 0, color = MT$zero_col, linewidth = MT$zero_lwd,
             linetype = "dashed") +
  geom_polygon(data = hull_df, aes(Axis1, Axis2, fill = group, group = group),
               alpha = 0.25, color = NA) +
  geom_polygon(data = hull_df, aes(Axis1, Axis2, color = group, group = group),
               fill = NA, linewidth = 0.01) +
  geom_point(data = scores_df, aes(Axis1, Axis2, color = group),
             shape = MT$pt_shape, size = MT$pt_size, stroke = MT$pt_stroke,
             alpha = MT$pt_alpha) +
  # Thin shafts with small open heads: at main-text scale a 0.4 linewidth with a
  # filled 0.16 cm head reads as a blob and competes with the hulls for
  # attention. The arrows are a secondary layer here.
  geom_segment(data = load_df,
               aes(x = 0, y = 0, xend = x, yend = y, linetype = block),
               arrow = grid::arrow(length = grid::unit(0.075, "cm"),
                                   angle = 22, type = "open"),
               color = "grey30", linewidth = 0.22) +
  ggrepel::geom_text_repel(data = load_df, aes(x = x, y = y, label = label),
                           size = 1.9, color = "grey20", min.segment.length = 0.2,
                           segment.size = 0.2, segment.color = "grey60",
                           box.padding = 0.2, max.overlaps = 30) +
  scale_color_manual(values = group_pal, name = NULL) +
  scale_fill_manual(values  = group_pal, name = NULL) +
  scale_linetype_manual(values = c("M-SPHARM" = "solid", "SP-SPHARM" = "22"),
                        # Untitled: what the arrows are (CLR loadings, back-mapped
                        # from the ILR coordinates used in the analysis) is stated
                        # in the caption instead.
                        name = NULL) +
  labs(x = sprintf("MFA axis 1 (%.1f%%)", inertia[1]),
       y = sprintf("MFA axis 2 (%.1f%%)", inertia[2])) +
  # Equal aspect: the analysis is distance-based and the two axes carry nearly the
  # same inertia, so an unequal scale would misrepresent between-group proximity.
  coord_fixed() +
  guides(color = guide_legend(order = 1, override.aes = list(shape = 16, size = 1.6,
                                                             linetype = 0)),
         fill  = "none",
         linetype = guide_legend(order = 2, title = NULL,
                                 override.aes = list(color = "grey30"))) +
  mt_theme() +
  # Inset legend, boxed on a 75% white ground with a grey80 hairline border --
  # the main-text ordination convention (Fig. 8a, Fig. 9). The translucent fill
  # is what lets the hull polygon underneath still read. Anchored bottom LEFT,
  # not bottom right -- the M-SPHARM l1 arrow points into the lower-right
  # quadrant and its repel label landed on top of the legend there.
  theme(legend.position       = c(0.01, 0.01),
        legend.justification  = c(0, 0),
        legend.background     = element_rect(fill = scales::alpha("white", 0.75),
                                             colour = "grey80", linewidth = 0.25),
        legend.box.background = element_rect(fill = "transparent", colour = NA),
        legend.key            = element_rect(fill = "transparent", colour = NA),
        legend.margin         = margin(2, 4, 2, 4),
        # 0.02 cm was fine while the two legends were unboxed; now that each
        # carries a border they need a visible hairline gap between them.
        legend.spacing.y      = grid::unit(0.10, "cm"))

ggsave(file.path(FIG_DIR, "fig_S_joint_mfa_biplot.png"), p_biplot,
       width = 3.7, height = 3.9, dpi = MT$dpi)
cat("Wrote figures/fig_S_joint_mfa_biplot.png\n")

# =============================================================================
# (6) Weight scan — pairwise resolution as a function of the morphology share
# =============================================================================
# TOTAL-INERTIA (Frobenius) normalisation here, not s1: after dividing each block
# by its own Frobenius norm both blocks contribute inertia 1, so multiplying the
# coordinates by sqrt(w) / sqrt(1 - w) makes w exactly the morphology share of the
# squared distances that PERMANOVA partitions. Under the s1 normalisation the
# x-axis would not have that reading.
cat("\n", strrep("=", 70), "\n", sep = "")
cat(sprintf("WEIGHT SCAN: w = %.2f .. %.2f, %d steps (w = 1 pure M, w = 0 pure SP)\n",
            min(W_GRID), max(W_GRID), length(W_GRID)))
cat(strrep("=", 70), "\n", sep = "")

f_M  <- sqrt(block_inertia(Z_M))
f_SP <- sqrt(block_inertia(Z_SP))
ZM_n <- Z_M  / f_M
ZS_n <- Z_SP / f_SP

scan_res <- map_dfr(seq_along(W_GRID), function(i) {
  w  <- W_GRID[i]
  Zw <- cbind(sqrt(w) * ZM_n, sqrt(1 - w) * ZS_n)
  g  <- permanova_global(Zw, blocks$group)
  pw <- permanova_pairwise(Zw, blocks$group)
  cat(sprintf("  w = %.2f  R2 = %.4f  p = %.4f  resolved = %d/%d\n",
              w, g$R2, g$p, sum(pw$significant_holm), n_pairs))
  tibble(w = w, R2 = g$R2, pseudo_F = g$F, p = g$p,
         n_sig_pairs = sum(pw$significant_holm), n_pairs = n_pairs,
         sig_pairs = paste(sort(pw$comparison[pw$significant_holm]), collapse = "; "),
         p_holm_focus1 = pair_p(pw, FOCUS_PAIRS[[1]]),
         p_holm_focus2 = pair_p(pw, FOCUS_PAIRS[[2]]))
}) %>%
  rename(!!paste0("p_holm_", paste(FOCUS_PAIRS[[1]], collapse = "_")) := p_holm_focus1,
         !!paste0("p_holm_", paste(FOCUS_PAIRS[[2]], collapse = "_")) := p_holm_focus2)

# Endpoint identity check: a global rescale of one block cannot change R2, so the
# w = 0 / w = 1 rows must reproduce the single-block models exactly.
end_ok <- isTRUE(all.equal(scan_res$R2[scan_res$w == 0],
                           model_res[["SP-SPHARM"]]$global$R2, tolerance = 1e-6)) &&
  isTRUE(all.equal(scan_res$R2[scan_res$w == 1],
                   model_res[["M-SPHARM"]]$global$R2, tolerance = 1e-6))
cat(sprintf("\n  Endpoint identity (R2 at w = 0 / 1 vs SP / M alone): %s\n",
            ifelse(end_ok, "OK", "<-- CHECK")))

# The reference line carries no in-panel label: the annotation repeated a value
# the caption states, and at 2.6 pt inside a half-width main-text panel it
# collided with the curve. Axis titles are shortened for the same reason -- the
# Holm correction and the alpha level belong in the caption, not on the axis.
p_scan <- ggplot(scan_res, aes(w, n_sig_pairs)) +
  geom_vline(xintercept = W_REF, linetype = "dashed",
             color = "#802520", linewidth = 0.4) +
  geom_line(linewidth = 0.6, color = "#4A6E8A") +
  # Same marker size as panel c: the two sit stacked at identical physical size,
  # so unequal dots read as an inconsistency rather than as denser sampling.
  geom_point(size = 1.1, color = "#4A6E8A") +
  scale_x_continuous(breaks = seq(0, 1, 0.2), limits = c(0, 1)) +
  scale_y_continuous(breaks = seq(0, n_pairs, 2), limits = c(0, n_pairs)) +
  labs(x = "w", y = "Pairs resolved") +
  mt_theme()

ggsave(file.path(FIG_DIR, "fig_S_joint_mfa_weight_scan.png"), p_scan,
       width = 3.1, height = 1.95, dpi = MT$dpi)
cat("Wrote figures/fig_S_joint_mfa_weight_scan.png\n")

# =============================================================================
# Outputs: tidy CSVs
# =============================================================================
# One tidy table: `scope = "global"` rows carry the model-level statistics
# (R2 / pseudo-F / p / resolved pairs / PERMDISP), `scope = "pairwise"` rows carry
# every pair with its raw and Holm-corrected p.
permanova_csv <- bind_rows(
  comparison_tbl %>%
    transmute(model, scope = "global", comparison = "type (all groups)",
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

write_csv(permanova_csv, file.path(OUT_DIR, "joint_mfa_permanova_comparison.csv"))
write_csv(scan_res,      file.path(OUT_DIR, "joint_mfa_weight_scan.csv"))
cat("\nWrote joint_mfa_permanova_comparison.csv and joint_mfa_weight_scan.csv\n")

# =============================================================================
# Console summary
# =============================================================================
cat("\n", strrep("=", 70), "\n", sep = "")
cat("SUMMARY\n")
cat(strrep("=", 70), "\n", sep = "")
print(comparison_tbl %>%
        transmute(model, R2 = round(R2, 5), pseudo_F = round(pseudo_F, 5), p,
                  resolved = sprintf("%d/%d", n_sig_pairs, n_pairs),
                  permdisp_p) %>%
        as.data.frame())
if (ASSEMBLAGE == "EXP") {
  cat("\nFocus pairs, Holm-corrected p:\n")
  print(as.data.frame(focus_tbl %>% mutate(across(where(is.numeric), ~ round(.x, 4)))))
}
cat(sprintf("\nMFA weight setting w_MFA = %.4f; scan best = %d/%d pairs at w = %s\n",
            w_mfa, max(scan_res$n_sig_pairs), n_pairs,
            paste(sprintf("%.2f", scan_res$w[scan_res$n_sig_pairs ==
                                               max(scan_res$n_sig_pairs)]),
                  collapse = ", ")))
cat("Reminder: the combined R2 cannot exceed max(R2_M, R2_SP) under Euclidean\n")
cat("geometry — a drop is arithmetic, not evidence about discriminant power.\n")

# =============================================================================
# =============================================================================
# SUPPLEMENTARY DIAGNOSTICS (D1-D4)
# =============================================================================
# =============================================================================
# Four add-on diagnostics. None changes the analysis above — they add outputs and
# refine the weight grid only. Each reuses the helpers defined earlier
# (make_ilr / permanova_global / permanova_pairwise / block_inertia / pair_p /
# permdisp) rather than re-implementing them, and writes to new filenames.

# The permutation floor: with N_PERM permutations adonis2 cannot return anything
# below 1/(N_PERM + 1). Several pairs sit exactly on that floor, so their p-ratios
# are censored, not measured. D1 and D3 both flag this rather than reading a ratio
# of 1.00 as "no change".
P_FLOOR <- 1 / (N_PERM + 1)

# =============================================================================
# D1. Raw (uncorrected) p-values — three models and the weight scan
# =============================================================================
# Holm is a step-down procedure: its output is forced monotone in the raw ranking,
# so neighbouring pairs are routinely pushed onto the SAME adjusted value — both
# focus pairs land on exactly 0.3662, which tells us nothing about whether merging
# moved them. The "morphology helps where morphology has signal" reading also has
# to be tested on RAW p: Holm-adjusted values are mutually dependent, so one pair
# moving drags its neighbours along.
cat("\n", strrep("=", 70), "\n", sep = "")
cat("D1. RAW vs HOLM-CORRECTED PAIRWISE p\n")
cat(strrep("=", 70), "\n", sep = "")

pairwise_all <- map_dfr(names(model_res),
                        ~ model_res[[.x]]$pairwise %>% mutate(model = .x))

# Long table with a p_type column, models as wide columns.
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
  arrange(p_M)          # ascending in M-alone raw p — the ordering that matters

cat(sprintf("  (raw p floor = %.4f at %d permutations; ratios marked * are censored)\n\n",
            P_FLOOR, N_PERM))
cat(sprintf("  %-30s %8s %8s %8s %12s %10s\n",
            "comparison (sorted by p_M)", "p_M", "p_SP", "p_comb", "p_holm_comb",
            "SP/comb"))
for (i in seq_len(nrow(d1_tbl))) {
  r <- d1_tbl[i, ]
  cat(sprintf("  %-30s %8.4f %8.4f %8.4f %12.4f %9.2f%s\n",
              r$comparison, r$p_M, r$p_SP, r$p_comb, r$p_holm_comb,
              r$improve_x, ifelse(r$censored, "*", " ")))
}

# Does the improvement factor track morphological signal? If morphology helps
# where morphology has signal, improve_x should be LARGE where p_M is SMALL, i.e.
# a NEGATIVE rank correlation. Censored rows carry no information and are dropped.
d1_free <- d1_tbl %>% filter(!censored)
rho_improve <- if (nrow(d1_free) >= 3)
  suppressWarnings(stats::cor(d1_free$p_M, d1_free$improve_x, method = "spearman")) else NA_real_
d1_monotone <- is.finite(rho_improve) && rho_improve <= -0.5
cat(sprintf("\n  Spearman(p_M, improvement factor) on the %d uncensored pair(s) = %s\n",
            nrow(d1_free), ifelse(is.finite(rho_improve), sprintf("%.3f", rho_improve), "NA")))
cat(sprintf("  => improvement %s monotone in morphological signal%s\n",
            ifelse(d1_monotone, "IS", "is NOT"),
            ifelse(nrow(d1_free) < 5, " (too few uncensored pairs to be conclusive)", "")))

if (ASSEMBLAGE == "EXP") {
  cat("\n  --- Focus pairs, RAW p (did merging move them, or did Holm flatten them?) ---\n")
  focus_raw <- map_dfr(FOCUS_PAIRS, function(pr) {
    tibble(comparison = paste(pr, collapse = "-"),
           p_M_raw    = pair_p(model_res[["M-SPHARM"]]$pairwise,     pr, "p"),
           p_SP_raw   = pair_p(model_res[["SP-SPHARM"]]$pairwise,    pr, "p"),
           p_comb_raw = pair_p(model_res[["MFA-combined"]]$pairwise, pr, "p"),
           p_comb_holm = pair_p(model_res[["MFA-combined"]]$pairwise, pr))
  })
  print(as.data.frame(focus_raw %>% mutate(across(where(is.numeric), ~ round(.x, 4)))))
  focus_moved <- any(abs(focus_raw$p_comb_raw - focus_raw$p_SP_raw) > 0.01)
  focus_tied  <- length(unique(round(focus_raw$p_comb_holm, 6))) == 1
  cat(sprintf("  => raw p %s between SP-alone and combined; the two Holm values are %s.\n",
              ifelse(focus_moved, "DID move", "barely moved (< 0.01)"),
              ifelse(focus_tied, "tied by Holm's monotonicity", "distinct")))
}

write_csv(pairwise_raw_holm, file.path(OUT_DIR, "joint_mfa_pairwise_raw_and_holm.csv"))
cat("\n  Wrote joint_mfa_pairwise_raw_and_holm.csv\n")

# =============================================================================
# D2. Loading block structure of the leading axes
# =============================================================================
# Axis 1 and axis 2 carry almost the same inertia. If axis 1 is driven by M and
# axis 2 by SP, that is the main text's decoupling result (RV ~ 0.10, Mantel ~ 0)
# showing up directly in the ordination: a global PCA of two equally weighted,
# mutually uninformative blocks should produce exactly that — one axis per block,
# with near-equal eigenvalues.
cat("\n", strrep("=", 70), "\n", sep = "")
cat("D2. BLOCK STRUCTURE OF THE LEADING MFA AXES\n")
cat(strrep("=", 70), "\n", sep = "")

n_axes_diag <- min(4, ncol(pca_comb$rotation))
axis_block_tbl <- map_dfr(seq_len(n_axes_diag), function(k) {
  # rotation columns are unit vectors, so the two block shares sum to 1 by
  # construction. idx_M / idx_SP index the MFA-normalised coordinates, i.e. the
  # same space the PERMANOVA ran in.
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
if (d2_split)
  cat("     This is the block decoupling (RV ~ 0.10, Mantel ~ 0) seen in the ordination:\n",
      "     two equally weighted, mutually uninformative blocks each claim their own axis.\n", sep = "")

rel_gap <- abs(inertia[1] - inertia[2]) / inertia[1]
d2_degenerate <- rel_gap < DEGENERACY_THR
if (d2_degenerate) {
  cat(sprintf("\n  CAUTION: axes 1 and 2 are near-degenerate (%.1f%% vs %.1f%%, relative gap %.1f%%).\n",
              inertia[1], inertia[2], rel_gap * 100))
  cat("  Within a near-degenerate 2D eigenspace the rotation is not well determined:\n")
  cat("  the plane is stable but the individual axis directions are not, so the exact\n")
  cat("  bearing of any single loading arrow in the biplot must not be over-interpreted.\n")
}

write_csv(axis_block_tbl, file.path(OUT_DIR, "joint_mfa_axis_block_contributions.csv"))
cat("\n  Wrote joint_mfa_axis_block_contributions.csv\n")

# =============================================================================
# D3. Dense weight grid + continuous discrimination curve
# =============================================================================
# Two problems with the coarse scan: (a) it loses 4 pairs between w = 0.70 and
# 0.75, which smells like a Holm cascade — when several raw p bunch near the
# threshold, step-down monotonicity flips them together; (b) "pairs resolved" is a
# 0-10 integer, far too coarse to show the real margin changes visible in the
# pairwise table. So: refine the grid where the cliff is, and add a continuous
# statistic (median -log10 raw p) alongside counts at two alphas and on raw p.
cat("\n", strrep("=", 70), "\n", sep = "")
cat(sprintf("D3. DENSE WEIGHT SCAN: %d steps (0.01 resolution over 0.61-0.84)\n",
            length(W_GRID_DENSE)))
cat(strrep("=", 70), "\n", sep = "")

dense_res <- map_dfr(seq_along(W_GRID_DENSE), function(i) {
  w  <- W_GRID_DENSE[i]
  Zw <- cbind(sqrt(w) * ZM_n, sqrt(1 - w) * ZS_n)      # same construction as (6)
  g  <- permanova_global(Zw, blocks$group)
  pw <- permanova_pairwise(Zw, blocks$group)
  n05h <- sum(pw$p_holm < 0.05); n01h <- sum(pw$p_holm < 0.01)
  n05r <- sum(pw$p < 0.05)
  cat(sprintf("  w = %.2f  R2 = %.4f  holm05 = %2d  holm01 = %2d  raw05 = %2d  med(-log10 p) = %.3f\n",
              w, g$R2, n05h, n01h, n05r, stats::median(-log10(pw$p))))
  bind_cols(
    tibble(w = w, R2 = g$R2, pseudo_F = g$F, p = g$p, n_pairs = n_pairs,
           n_sig_pairs_holm_05 = n05h,
           n_sig_pairs_holm_01 = n01h,
           n_sig_pairs_raw_05  = n05r,
           median_neglog10_p_raw = stats::median(-log10(pw$p)),
           n_pairs_at_p_floor  = sum(pw$p <= P_FLOOR)),
    # all pairwise RAW p as wide columns (pair order is fixed by combn(levels))
    as_tibble_row(setNames(pw$p,
                           paste0("p_raw_", str_replace_all(pw$comparison, "-", "_")))))
})

write_csv(dense_res, file.path(OUT_DIR, "joint_mfa_weight_scan_dense.csv"))
cat("\n  Wrote joint_mfa_weight_scan_dense.csv\n")

# ---- read the curve ---------------------------------------------------------
d3_best_i <- which.max(dense_res$median_neglog10_p_raw)
d3_best_w <- dense_res$w[d3_best_i]
d3_best_v <- dense_res$median_neglog10_p_raw[d3_best_i]
# Right-hand end of the plateau that starts at w = 0 (contiguous run at the w = 0 level).
lvl0 <- dense_res$n_sig_pairs_holm_05[1]
d3_plateau_end <- max(dense_res$w[dplyr::cumall(dense_res$n_sig_pairs_holm_05 == lvl0)])
# Holm vs raw divergence in the cliff zone.
zone <- dense_res %>% filter(w >= 0.65, w <= 0.84)
zone_gap <- zone$n_sig_pairs_raw_05 - zone$n_sig_pairs_holm_05
d3_max_gap <- if (nrow(zone) > 0) max(zone_gap) else NA_integer_
d3_gap_w   <- if (nrow(zone) > 0) zone$w[which.max(zone_gap)] else NA_real_
d3_cliff_real <- is.finite(d3_max_gap) && d3_max_gap <= 1

cat(sprintf("\n  Continuous curve peaks at w = %.2f (median -log10 raw p = %.3f);\n",
            d3_best_w, d3_best_v))
cat(sprintf("    value at w = 0 (pure SP) = %.3f, at w_MFA = %.2f it is %.3f.\n",
            dense_res$median_neglog10_p_raw[dense_res$w == 0], w_mfa,
            dense_res$median_neglog10_p_raw[which.min(abs(dense_res$w - w_mfa))]))
cat(sprintf("  Plateau at %d/%d pairs runs from w = 0 to w = %.2f (exact right endpoint).\n",
            lvl0, n_pairs, d3_plateau_end))
cat(sprintf("  Over w in [0.65, 0.84] the largest raw-minus-Holm count gap is %d (at w = %.2f).\n",
            d3_max_gap, d3_gap_w))
cat(sprintf("  => the drop across that interval is %s.\n",
            ifelse(d3_cliff_real,
                   "REAL loss of signal, not a Holm artefact (raw and Holm counts track each other)",
                   "amplified by the Holm cascade (raw p keep several pairs that Holm rejects)")))

# ---- figures ----------------------------------------------------------------
# The three-criterion dense step plot (fig_S_joint_mfa_weight_scan_dense.png)
# was dropped: the composite main-text figure carries the coarse scan and the
# continuous curve, and the criterion comparison it showed is reported
# numerically in the D3 console block and in joint_mfa_weight_scan_dense.csv.
# `dense_long` existed only to feed that plot and is gone with it.

# Interior peak of the continuous curve. Taken from the data, not hard-coded, so
# it cannot silently drift out of step with the CSV.
# CAVEAT: the dense 0.01 refinement covers only w = 0.61-0.84 (W_GRID_DENSE); around
# the peak the grid is still 0.05, so this locates the maximum to +/- 0.05 and is
# not a resolved optimum.
i_peak <- which.max(dense_res$median_neglog10_p_raw)
W_PEAK <- dense_res$w[i_peak]
cat(sprintf("\n  Continuous-curve peak: w = %.2f (median -log10 p = %.4f) vs %.4f at w = 0\n",
            W_PEAK, dense_res$median_neglog10_p_raw[i_peak],
            dense_res$median_neglog10_p_raw[dense_res$w == 0]))
cat(sprintf("    grid spacing at the peak = %.2f, so the maximum is located to +/- that.\n",
            min(diff(sort(dense_res$w[dense_res$w >= 0.15 & dense_res$w <= 0.40])))))

p_cont <- ggplot(dense_res, aes(w, median_neglog10_p_raw)) +
  # Two markers, deliberately different colours: dark red = equal weighting
  # (as in panel b), olive = the interior peak. Both unlabelled; the caption
  # names them.
  geom_vline(xintercept = W_REF, linetype = "dashed",
             color = "#802520", linewidth = 0.4) +
  geom_vline(xintercept = W_PEAK, linetype = "dashed",
             color = "#788C4A", linewidth = 0.4) +
  geom_line(linewidth = 0.6, color = "#4A6E8A") +
  geom_point(size = 1.1, color = "#4A6E8A") +
  scale_x_continuous(breaks = seq(0, 1, 0.2), limits = c(0, 1)) +
  labs(x = "w", y = expression("Median" ~ -log[10] * "(" * italic(p) * ")")) +
  mt_theme()

ggsave(file.path(FIG_DIR, "fig_S_joint_mfa_weight_resolution_continuous.png"), p_cont,
       width = 3.1, height = 1.95, dpi = MT$dpi)
cat("  Wrote figures/fig_S_joint_mfa_weight_resolution_continuous.png\n")

# =============================================================================
# MAIN-TEXT COMPOSITE -- biplot left, the two weight panels stacked right
# =============================================================================
# Assembled from the ggplot objects, NOT by stitching the exported PNGs: that
# keeps one typographic definition across all three panels and avoids resampling
# raster panels that were each rendered at a different physical size.
#
# Layout: (a) MFA biplot on the left, (b) coarse weight scan top right,
# (c) continuous resolution curve bottom right -- b above c as requested.
p_composite <- (
  (p_biplot + labs(tag = "a")) |
    ((p_scan + labs(tag = "b")) / (p_cont + labs(tag = "c")))
) + plot_layout(widths = c(1.18, 1))

COMP_H <- 3.9
ggsave(file.path(FIG_DIR, "fig_joint_mfa_main.png"), p_composite,
       width = MT$fig_width_in, height = COMP_H, dpi = MT$dpi)
cat("  Wrote figures/fig_joint_mfa_main.png   [MAIN TEXT]\n")

cat("\n", strrep("=", 70), "\n", sep = "")
cat("MAIN-TEXT COMPOSITE: SPECIFICATION AND RENDERED SIZE\n")
cat(strrep("=", 70), "\n", sep = "")
cat("  spec source: SDG_cores_statistics.R make_coia_biplot() (manuscript Fig. 9)\n")
cat("               + manuscript.qmd fig-width / out-width / fig-dpi\n")
for (nm in names(MT)) cat(sprintf("    %-14s = %s\n", nm, format(MT[[nm]])))
cat(sprintf("  canvas   : %.2f x %.2f in (%.1f x %.1f mm) at %d dpi = %d x %d px\n",
            MT$fig_width_in, COMP_H, MT$fig_width_in * 25.4, COMP_H * 25.4,
            MT$dpi, round(MT$fig_width_in * MT$dpi), round(COMP_H * MT$dpi)))
cat(sprintf("  placed at 174 mm (full text width): scale %.3f, so %.1f pt axis text\n",
            174 / (MT$fig_width_in * 25.4), MT$axis_text))
cat(sprintf("  reference line on (b) and (c) at w = %.2f, UNLABELLED (dark red);\n", W_REF))
cat(sprintf("    the s1 normalisation lands at w_mfa = %.4f, i.e. %.4f away.\n",
            w_mfa, abs(w_mfa - W_REF)))
cat(sprintf("  second line on (c) at the curve's peak, w = %.2f, UNLABELLED (olive).\n",
            W_PEAK))
cat("  Neither line is annotated in-panel, so the caption must name both.\n")

# =============================================================================
# D4. PERMDISP group dispersions
# =============================================================================
# Combined-space PERMDISP is non-significant (p = 0.393) where SP alone is
# significant, which would seem to retire the main text's "interpret the PERMANOVA
# cautiously" caveat. But appending 7 dimensions whose between-group dispersions
# are similar mechanically pulls any dispersion RATIO toward 1. So: did discoid's
# dispersion actually come down, or was it just diluted? Absolute distances are
# not comparable across the three spaces (different overall scale), so everything
# is reported as a ratio to the least-dispersed group of the same model.
cat("\n", strrep("=", 70), "\n", sep = "")
cat("D4. PERMDISP GROUP DISPERSIONS (M alone / SP alone / combined)\n")
cat(strrep("=", 70), "\n", sep = "")

group_n <- table(blocks$group)
permdisp_by_group <- map_dfr(names(MODELS), function(nm) {
  pdm <- permdisp(MODELS[[nm]], blocks$group)     # reuses the helper from (4)
  md  <- tapply(pdm$disp$distances, pdm$disp$group, mean)
  tibble(model = nm, group = names(md), n = as.integer(group_n[names(md)]),
         mean_dist_to_centroid = as.numeric(md),
         permdisp_F = pdm$F, permdisp_p = pdm$p) %>%
    mutate(dispersion_ratio_vs_min = mean_dist_to_centroid / min(mean_dist_to_centroid))
})

for (nm in names(MODELS)) {
  sub <- permdisp_by_group %>% filter(model == nm) %>% arrange(desc(dispersion_ratio_vs_min))
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

# Deciding between "discoid really got less variable" and "the imbalance was spread
# thinner" cannot be done with a drop threshold: mechanical dilution produces a
# large drop too. The discriminating question is whether the combined ratio is a
# BLEND of the two single-block ratios or an active CANCELLATION.
#   - Group centroids in the concatenated space are the concatenated block
#     centroids, so combined dispersion is built from the same per-specimen block
#     dispersions — it carries no new information about any group.
#   - Hence a combined ratio lying between the two block ratios is exactly what
#     adding a more evenly dispersed block must produce: dilution.
#   - Only a combined ratio BELOW both block ratios would mean the blocks actively
#     cancel each other's heterogeneity, which no amount of dilution can achieve.
tol <- 0.02
d4_bracketed <- is.finite(r_M) && is.finite(r_sp) && is.finite(r_comb) &&
  r_comb >= min(r_M, r_sp) - tol && r_comb <= max(r_M, r_sp) + tol
d4_substantive <- is.finite(r_comb) && r_comb < min(r_M, r_sp) - tol
d4_same_top <- length(unique(d4_ratio$most_dispersed)) == 1

# Mechanism behind whatever the verdict turns out to be: do the two blocks rank the
# groups' dispersion the same way? Anti-correlated orderings let the blocks cancel;
# similar orderings mean the combined space can only interpolate (dilute).
disp_wide <- permdisp_by_group %>%
  select(model, group, mean_dist_to_centroid) %>%
  pivot_wider(names_from = model, values_from = mean_dist_to_centroid)
rho_disp <- suppressWarnings(stats::cor(disp_wide$`M-SPHARM`, disp_wide$`SP-SPHARM`,
                                        method = "spearman"))

cat(sprintf("\n  M alone = %.3f, SP alone = %.3f, combined = %.3f (combined/SP = %.3f).\n",
            r_M, r_sp, r_comb, r_comb / r_sp))
cat(sprintf("  Spearman between the two blocks' group dispersion orderings = %.3f (%s)\n",
            rho_disp,
            dplyr::case_when(!is.finite(rho_disp)  ~ "undefined",
                             rho_disp <= -0.5      ~ "reversed orderings",
                             rho_disp <   0.3      ~ "essentially unrelated orderings",
                             TRUE                  ~ "similar orderings")))
# A whole-ordering rank correlation is a blunt summary: cancellation is driven by
# the EXTREMES, so also report whether the top/bottom groups swap ends.
top_M   <- disp_wide$group[which.max(disp_wide$`M-SPHARM`)]
bot_M   <- disp_wide$group[which.min(disp_wide$`M-SPHARM`)]
top_SP  <- disp_wide$group[which.max(disp_wide$`SP-SPHARM`)]
bot_SP  <- disp_wide$group[which.min(disp_wide$`SP-SPHARM`)]
extremes_swap <- (top_M == bot_SP) || (top_SP == bot_M)
cat(sprintf("  Extremes: M most/least = %s / %s; SP most/least = %s / %s%s\n",
            top_M, bot_M, top_SP, bot_SP,
            ifelse(extremes_swap,
                   "  -> an extreme group swaps ends, which is what pulls the combined ratio down",
                   "  -> extremes do not swap")))
cat(sprintf("  Most-dispersed group per model: %s%s\n",
            paste(sprintf("%s = %s", d4_ratio$model, d4_ratio$most_dispersed),
                  collapse = "; "),
            ifelse(d4_same_top, "  (unchanged across all three spaces)", "")))
if (d4_substantive) {
  cat("  => SUBSTANTIVE: the combined ratio falls BELOW both single-block ratios, so the\n")
  cat("     two blocks actively cancel each other's dispersion imbalance — more than\n")
  cat("     dilution can explain.\n")
} else if (d4_bracketed) {
  cat("  => DILUTION: the combined ratio sits BETWEEN the two single-block ratios, which\n")
  cat("     is precisely what appending a more evenly dispersed block must produce. Group\n")
  cat("     centroids in the concatenated space are the concatenated block centroids, so\n")
  cat("     no group's dispersion has actually fallen — the same imbalance is spread over\n")
  cat("     more coordinates and the F statistic drops with it. The main text's\n")
  cat("     'interpret the PERMANOVA cautiously' caveat still applies; word any PERMDISP\n")
  cat("     improvement conservatively.\n")
} else {
  cat("  => The combined ratio EXCEEDS both block ratios; inspect before interpreting.\n")
}

write_csv(permdisp_by_group %>%
            select(model, group, n, mean_dist_to_centroid, dispersion_ratio_vs_min),
          file.path(OUT_DIR, "joint_mfa_permdisp_by_group.csv"))
cat("\n  Wrote joint_mfa_permdisp_by_group.csv\n")

# =============================================================================
# Diagnostic conclusions
# =============================================================================
cat("\n", strrep("=", 70), "\n", sep = "")
cat("DIAGNOSTIC CONCLUSIONS\n")
cat(strrep("=", 70), "\n", sep = "")
if (ASSEMBLAGE == "EXP") {
  cat(sprintf("(1) Focus pairs: raw p %s on merging (SP-alone %s vs combined %s).\n",
              ifelse(focus_moved, "DID shift", "barely shifted"),
              paste(sprintf("%.4f", focus_raw$p_SP_raw),   collapse = " / "),
              paste(sprintf("%.4f", focus_raw$p_comb_raw), collapse = " / ")))
  cat(sprintf("    Their equal Holm values (%.4f) are %s.\n",
              focus_raw$p_comb_holm[1],
              ifelse(focus_tied, "an artefact of Holm's forced monotonicity, not equality of evidence",
                     "genuinely distinct")))
  cat(sprintf("    Improvement factor %s monotone in p_M (Spearman %s, %d uncensored pairs).\n",
              ifelse(d1_monotone, "IS", "is NOT"),
              ifelse(is.finite(rho_improve), sprintf("%.3f", rho_improve), "NA"),
              nrow(d1_free)))
}
cat(sprintf("(2) Leading axes %s one-axis-per-block (axis 1 %s %.0f%%, axis 2 %s %.0f%%)%s.\n",
            ifelse(d2_split, "DO show", "do NOT show"),
            axis_block_tbl$dominant_block[1], 100 * max(axis_block_tbl$share_M[1],
                                                        axis_block_tbl$share_SP[1]),
            axis_block_tbl$dominant_block[2], 100 * max(axis_block_tbl$share_M[2],
                                                        axis_block_tbl$share_SP[2]),
            ifelse(d2_degenerate, "; axes near-degenerate, arrow bearings not interpretable", "")))
cat(sprintf("(3) The 0.70-0.75 drop is %s; dense grid puts the %d/%d plateau's right edge at w = %.2f,\n",
            ifelse(d3_cliff_real, "REAL", "a Holm cascade"), lvl0, n_pairs, d3_plateau_end))
cat(sprintf("    and the continuous curve peaks at w = %.2f (%.3f) vs %.3f at pure SP — %s.\n",
            d3_best_w, d3_best_v, dense_res$median_neglog10_p_raw[dense_res$w == 0],
            ifelse(d3_best_w > 0, "an interior peak exists", "no interior peak")))
cat(sprintf("(4) PERMDISP: max/min dispersion ratio is %.3f (M alone), %.3f (SP alone), %.3f (combined).\n",
            r_M, r_sp, r_comb))
if (d4_substantive) {
  cat(sprintf("    Combined falls below BOTH block ratios, so this is not mere dilution: %s\n",
              ifelse(extremes_swap,
                     sprintf("%s is the most\n    dispersed group in one block and the least in the other, so the extremes cancel",
                             ifelse(top_M == bot_SP, top_M, top_SP)),
                     "the blocks' dispersion\n    patterns partly cancel")))
  cat(sprintf("    (whole-ordering Spearman is only %.3f, so the effect is carried by the extremes,\n",
              rho_disp))
  cat(sprintf("    not by the overall ranking). But %s is\n",
              d4_ratio$most_dispersed[d4_ratio$model == "MFA-combined"]))
  cat(sprintf("    still the most dispersed group and the margin below M alone is small (%.3f vs %.3f),\n",
              r_comb, r_M))
  cat("    so the caveat softens rather than disappears.\n")
} else if (d4_bracketed) {
  cat("    Combined sits between the two block ratios — dilution, not a real reduction in\n")
  cat("    heterogeneity. Word the PERMDISP improvement conservatively.\n")
} else {
  cat("    Combined exceeds both block ratios; inspect before interpreting.\n")
}

cat("\nDone.\n\n")
print(sessionInfo())
