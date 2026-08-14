# coia_power_sensitivity_stats.R

suppressPackageStartupMessages({
  library(here); library(tidyverse); library(readxl)
  library(vegan); library(ade4); library(compositions)
  library(MatrixCorrelation); library(conflicted)
})
suppressMessages({
  conflicts_prefer(dplyr::filter, dplyr::select, dplyr::lag,
                   stats::sd, stats::var, stats::dist, stats::cor, stats::cov,
                   base::scale, base::norm, base::`%*%`, .quiet = TRUE)
})
set.seed(42)

# =============================================================================
# PARAMETERS (kept modest; this is supplementary)
# =============================================================================
N_BOOT     <- 2000
N_SIM      <- 500
N_PERM     <- 199
N_POP      <- 5000
POP_RV_GRID <- seq(0.05, 0.40, by = 0.05)
POWER_THR  <- 0.80
ALPHA      <- 0.05
BOOT_SEED  <- 42
SIM_SEED   <- 2000
CAL_SEED   <- 100

OUT_DIR <- here("analysis/robustness/coia_power_sensitivity")
FIG_DIR <- file.path(OUT_DIR, "figures")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

POWER_COLS_DIR     <- paste0("power_l", 1:6)
POWER_COLS_MORPH   <- paste0("power_l", 1:8)
EXCLUDE_CORE_TYPES <- c("Handaxe", "Pick")

REF <- list(exp_RV = 0.10095, exp_mantel_r = -0.06126,
            sdg_RV = 0.09066, sdg_mantel_r = 0.01189)

# =============================================================================
# Helpers — copied VERBATIM from the source scripts (attribution in comments)
# =============================================================================
replace_zeros <- function(X, delta = NULL) {
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
filter_spharm <- function(df, power_cols, meta = NULL) {
  result <- df %>% select(ID, Typology, all_of(power_cols))
  if (!is.null(meta)) result <- left_join(result, meta, by = "ID")
  result
}
split_by_group <- function(df) list(
  exp_im = df %>% filter(str_starts(ID, "EXP") | str_starts(ID, "IM_")),
  sdg_im = df %>% filter(str_starts(ID, "SDG") | str_starts(ID, "IM_")))

# =============================================================================
# Inputs (committed products of the main pipeline + raw metadata)
# =============================================================================
cat("Loading committed main outputs (SPHARM direction / morphology, metadata)...\n")

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
# Build the EXACT CoIA ILR matrices + reported RV / Mantel
# =============================================================================
build_ilr_tables <- function(dir_df, morph_df) {
  list(dir_f   = filter_spharm(dir_df,   POWER_COLS_DIR,   metric_data),
       morph_f = filter_spharm(morph_df, POWER_COLS_MORPH, metric_data))
}

exp_ilr_tables <- function() {
  bt <- build_ilr_tables(SPHARM_direction, SPHARM_morphology)
  df_scar_all  <- split_by_group(bt$dir_f)$exp_im
  df_morph_all <- split_by_group(bt$morph_f)$exp_im
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

  exp_ids <- rownames(morph_power)[
    !str_starts(rownames(morph_power), "IM_") &
      rownames(morph_power) != "EXP43_Biface"]
  list(M = morph_ilr[exp_ids, , drop = FALSE],
       S = scar_ilr[exp_ids,  , drop = FALSE], ids = exp_ids)
}

sdg_ilr_tables <- function() {
  bt <- build_ilr_tables(SPHARM_direction, SPHARM_morphology)
  df_scar_raw  <- bt$dir_f
  df_morph_raw <- bt$morph_f
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

  arch_ids <- rownames(morph_power)[
    !str_starts(rownames(morph_power), "IM_") &
      !str_starts(rownames(morph_power), "EXP")]
  meta_arch <- tibble(ID = arch_ids) %>%
    left_join(core_meta, by = "ID") %>%
    filter(!core_type %in% EXCLUDE_CORE_TYPES | is.na(core_type))
  arch_ids <- meta_arch$ID
  list(M = morph_ilr[arch_ids, , drop = FALSE],
       S = scar_ilr[arch_ids,  , drop = FALSE], ids = arch_ids)
}

raw_rv_coinertia <- function(M, S) {
  Mi <- as.data.frame(M); Si <- as.data.frame(S)
  colnames(Mi) <- paste0("M_ilr", seq_len(ncol(Mi)))
  colnames(Si) <- paste0("S_ilr", seq_len(ncol(Si)))
  dudi_morph <- dudi.pca(Mi, center = TRUE, scale = TRUE, scannf = FALSE, nf = ncol(Mi))
  dudi_scar  <- dudi.pca(Si, center = TRUE, scale = TRUE, scannf = FALSE, nf = ncol(Si))
  coin <- coinertia(dudi_morph, dudi_scar, scannf = FALSE, nf = 2)
  set.seed(42)
  rt <- randtest(coin, nrepet = 9999)
  list(RV = coin$RV, RV_p = rt$pvalue)
}

raw_mantel <- function(M, S) {
  set.seed(42)
  mt <- mantel(dist(M), dist(S), method = "spearman", permutations = 9999)
  list(r = mt$statistic, p = mt$signif)
}

# =============================================================================
# Statistical primitives
# =============================================================================
rv_std <- function(X, Y) {
  X <- base::scale(as.matrix(X)); Y <- base::scale(as.matrix(Y))
  SXY <- crossprod(X, Y); SXX <- crossprod(X); SYY <- crossprod(Y)
  sum(SXY^2) / sqrt(sum(SXX^2) * sum(SYY^2))
}
rv2_std <- function(X, Y) {
  Xs <- base::scale(as.matrix(X)); Ys <- base::scale(as.matrix(Y))
  if (any(!is.finite(Xs)) || any(!is.finite(Ys))) return(NA_real_)
  MatrixCorrelation::RV2(Xs, Ys, center = TRUE)
}
mantel_stat <- function(M, S) cor(as.vector(dist(M)), as.vector(dist(S)), method = "spearman")
mvrnorm_chol <- function(n, mu, R) sweep(matrix(rnorm(n * length(mu)), n) %*% R, 2, mu, "+")
chol_pd <- function(S) {
  eps <- 1e-8 * mean(diag(S))
  for (k in 0:8) {
    R <- tryCatch(chol(S + diag(eps, nrow(S))), error = function(e) NULL)
    if (!is.null(R)) return(R)
    eps <- eps * 10
  }
  chol(S + diag(eps, nrow(S)))
}

# =============================================================================
# (2) Percentile bootstrap CI for RV2 and Mantel r (rows resampled, blocks paired)
# =============================================================================
bootstrap_ci <- function(M, S, B = N_BOOT, seed = BOOT_SEED) {
  M <- as.matrix(M); S <- as.matrix(S); n <- nrow(M)
  set.seed(seed)
  rv2_b <- numeric(B); man_b <- numeric(B)
  for (b in seq_len(B)) {
    idx <- sample.int(n, n, replace = TRUE)
    Mb  <- M[idx, , drop = FALSE]; Sb <- S[idx, , drop = FALSE]
    rv2_b[b] <- suppressWarnings(rv2_std(Mb, Sb))
    man_b[b] <- suppressWarnings(mantel_stat(Mb, Sb))
  }
  qci <- function(v) stats::quantile(v[is.finite(v)], c(0.025, 0.975), names = FALSE)
  list(rv2 = qci(rv2_b), mantel = qci(man_b))
}

# =============================================================================
# (3) Power curve + MDES
# =============================================================================
block_params <- function(X) {
  X <- as.matrix(X); S <- stats::cov(X)
  list(mu = colMeans(X), R = chol_pd(S), u = eigen(S, symmetric = TRUE)$vectors[, 1])
}
simulate_pair <- function(n, pm, ps, a) {
  s <- rnorm(n)
  list(M = mvrnorm_chol(n, pm$mu, pm$R) + a * outer(s, pm$u),
       S = mvrnorm_chol(n, ps$mu, ps$R) + a * outer(s, ps$u))
}
pop_rv <- function(pm, ps, a, n = N_POP, seed) {
  set.seed(seed); sim <- simulate_pair(n, pm, ps, a); rv_std(sim$M, sim$S)
}
calibrate_a <- function(pm, ps, target_grid = POP_RV_GRID, seed_base = CAL_SEED) {
  lead_sd <- sqrt(max(crossprod(pm$R %*% pm$u), crossprod(ps$R %*% ps$u)))
  a_max   <- 2 * as.numeric(lead_sd)
  for (rep in 1:6) {
    a_dense  <- seq(0, a_max, length.out = 41)
    rv_dense <- vapply(seq_along(a_dense),
                       function(i) pop_rv(pm, ps, a_dense[i], seed = seed_base + i),
                       numeric(1))
    if (max(rv_dense, na.rm = TRUE) >= max(target_grid) + 0.03) break
    a_max <- a_max * 1.6
  }
  rv_mono <- cummax(rv_dense); keep <- !duplicated(rv_mono)
  tibble(target_rv = target_grid,
         a = stats::approx(rv_mono[keep], a_dense[keep], xout = target_grid, rule = 2)$y)
}
power_at_a <- function(n, pm, ps, a, nsim = N_SIM, nperm = N_PERM, seed) {
  set.seed(seed)
  rej_rv <- logical(nsim); rej_man <- logical(nsim)
  for (i in seq_len(nsim)) {
    sim <- simulate_pair(n, pm, ps, a)
    rv_p  <- suppressWarnings(
      ade4::RV.rtest(as.data.frame(sim$M), as.data.frame(sim$S), nrepet = nperm)$pvalue)
    man_p <- suppressWarnings(
      vegan::mantel(dist(sim$M), dist(sim$S), method = "spearman", permutations = nperm)$signif)
    rej_rv[i]  <- is.finite(rv_p)  && rv_p  < ALPHA
    rej_man[i] <- is.finite(man_p) && man_p < ALPHA
  }
  c(power_rv = mean(rej_rv), power_mantel = mean(rej_man))
}
mdes <- function(rv, pow, thr = POWER_THR) {
  o <- order(rv); rv <- rv[o]; pow <- pow[o]
  if (all(pow < thr, na.rm = TRUE)) return(NA_real_)
  if (pow[1] >= thr) return(rv[1])
  for (i in 2:length(rv)) {
    if (isTRUE(pow[i] >= thr) && isTRUE(pow[i - 1] < thr))
      return(rv[i - 1] + (thr - pow[i - 1]) * (rv[i] - rv[i - 1]) / (pow[i] - pow[i - 1]))
  }
  rv[which(pow >= thr)[1]]
}
power_curve <- function(M, S, label) {
  cat(sprintf("\n-- power curve: %s (n = %d) --\n", label, nrow(M)))
  pm <- block_params(M); ps <- block_params(S); n <- nrow(M)
  grid <- bind_rows(tibble(target_rv = 0, a = 0), calibrate_a(pm, ps))
  res <- map_dfr(seq_len(nrow(grid)), function(i) {
    a  <- grid$a[i]
    pw <- power_at_a(n, pm, ps, a, seed = SIM_SEED + i)
    prv <- if (a == 0) 0 else pop_rv(pm, ps, a, seed = CAL_SEED + 1000 + i)
    tibble(assemblage = label, a = a, pop_rv = prv,
           power_rv = pw["power_rv"], power_mantel = pw["power_mantel"])
  })
  list(curve = res, mdes_rv = mdes(res$pop_rv, res$power_rv),
       mdes_mantel = mdes(res$pop_rv, res$power_mantel))
}

# =============================================================================
# Run both assemblages
# =============================================================================
analyse_assemblage <- function(tabs, label) {
  cat("\n", strrep("=", 70), "\n", sep = "")
  cat(sprintf("ASSEMBLAGE: %s (n = %d)\n", label, length(tabs$ids)))
  cat(strrep("=", 70), "\n", sep = "")
  M <- as.matrix(tabs$M); S <- as.matrix(tabs$S)

  rv  <- raw_rv_coinertia(M, S)
  man <- raw_mantel(M, S)
  rv2 <- rv2_std(M, S)
  cat(sprintf("  raw RV (coinertia)      = %.5f (p = %.4f)\n", rv$RV, rv$RV_p))
  cat(sprintf("  rv_std (trace, check)   = %.5f  [should match raw RV]\n", rv_std(M, S)))
  cat(sprintf("  Mantel r (spearman)     = %.5f (p = %.4f)\n", man$r, man$p))
  cat(sprintf("  RV2 (bias-corrected)    = %.5f\n", rv2))

  boot <- bootstrap_ci(M, S)
  cat(sprintf("  RV2 95%% CI              = [%.4f, %.4f]\n", boot$rv2[1], boot$rv2[2]))
  cat(sprintf("  Mantel r 95%% CI         = [%.4f, %.4f]\n", boot$mantel[1], boot$mantel[2]))

  pc <- power_curve(M, S, label)
  cat(sprintf("  MDES@80%% (RV)           = %s\n",
              ifelse(is.na(pc$mdes_rv), ">max", sprintf("%.3f", pc$mdes_rv))))
  cat(sprintf("  MDES@80%% (Mantel)       = %s\n",
              ifelse(is.na(pc$mdes_mantel), ">max", sprintf("%.3f", pc$mdes_mantel))))

  metrics <- tibble(
    assemblage = label, n = length(tabs$ids),
    raw_RV = rv$RV, raw_RV_p = rv$RV_p,
    RV2 = rv2, RV2_CI_low = boot$rv2[1], RV2_CI_high = boot$rv2[2],
    mantel_r = man$r, mantel_p = man$p,
    mantel_CI_low = boot$mantel[1], mantel_CI_high = boot$mantel[2],
    MDES_RV = pc$mdes_rv, MDES_mantel = pc$mdes_mantel)
  list(metrics = metrics, curve = pc$curve)
}

exp_tabs <- exp_ilr_tables()
sdg_tabs <- sdg_ilr_tables()

# ---- (0) ANCHOR check vs the committed reported values ----------------------
cat("\n", strrep("=", 70), "\n", sep = "")
cat("ANCHOR CHECK vs committed CoIA results (exp_cia_analysis / sdg_cia_analysis)\n")
cat(strrep("=", 70), "\n", sep = "")
exp_anchor_rv  <- raw_rv_coinertia(as.matrix(exp_tabs$M), as.matrix(exp_tabs$S))$RV
exp_anchor_man <- raw_mantel(as.matrix(exp_tabs$M), as.matrix(exp_tabs$S))$r
sdg_anchor_rv  <- raw_rv_coinertia(as.matrix(sdg_tabs$M), as.matrix(sdg_tabs$S))$RV
sdg_anchor_man <- raw_mantel(as.matrix(sdg_tabs$M), as.matrix(sdg_tabs$S))$r
anchor_ok <- TRUE
check <- function(name, got, exp, tol = 1e-3) {
  ok <- is.finite(got) && abs(got - exp) <= tol
  if (!ok) anchor_ok <<- FALSE
  cat(sprintf("  %-22s got %10.5f  expected %10.5f  %s\n",
              name, got, exp, ifelse(ok, "OK", "<-- CHECK")))
}
check("EXP raw RV",   exp_anchor_rv,  REF$exp_RV)
check("EXP Mantel r", exp_anchor_man, REF$exp_mantel_r)
check("SDG raw RV",   sdg_anchor_rv,  REF$sdg_RV)
check("SDG Mantel r", sdg_anchor_man, REF$sdg_mantel_r)
cat(sprintf("\n  => %s\n", ifelse(anchor_ok,
    "Recomputed ILR matrices reproduce the reported raw RV and Mantel r.",
    "MISMATCH: investigate before trusting the power/sensitivity numbers.")))

exp_res <- analyse_assemblage(exp_tabs, "EXP")
sdg_res <- analyse_assemblage(sdg_tabs, "SDG")

# =============================================================================
# Assemble outputs
# =============================================================================
metrics_df    <- bind_rows(exp_res$metrics, sdg_res$metrics)
power_long_df <- bind_rows(exp_res$curve, sdg_res$curve) %>%
  pivot_longer(c(power_rv, power_mantel), names_to = "test", values_to = "power") %>%
  mutate(test = recode(test, power_rv = "RV", power_mantel = "Mantel"))

write_csv(metrics_df,    file.path(OUT_DIR, "coia_power_metrics.csv"))
write_csv(power_long_df, file.path(OUT_DIR, "coia_power_curves.csv"))
cat("\nWrote coia_power_metrics.csv and coia_power_curves.csv\n")

# =============================================================================
# Figure: power vs population RV, EXP & SDG, both tests
# =============================================================================
DATASET_LABS <- c(EXP = "Experimentally knapped cores", SDG = "Sandinggai cores")
TEST_COLORS  <- c(RV = "#4A6E8A", Mantel = "#BA8530")
obs_rv_df <- metrics_df %>% transmute(assemblage, raw_RV)

p_power <- ggplot(power_long_df, aes(pop_rv, power, color = test, group = test)) +
  geom_hline(yintercept = POWER_THR, linetype = "dashed", color = "grey50", linewidth = 0.4) +
  geom_vline(data = obs_rv_df, aes(xintercept = raw_RV),
             linetype = "dotted", color = "grey40", linewidth = 0.4) +
  geom_line(linewidth = 0.8) + geom_point(size = 1.6) +
  facet_wrap(~ assemblage, labeller = as_labeller(DATASET_LABS)) +
  scale_color_manual(values = TEST_COLORS, name = "Test") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(x = "Population RV coefficient", y = sprintf("Power (alpha = %.2f)", ALPHA)) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold"), legend.position = "bottom")

ggsave(file.path(FIG_DIR, "fig_S_coia_power.png"), p_power, width = 9, height = 4.5, dpi = 300)
cat("Wrote figures/fig_S_coia_power.png\n")

# =============================================================================
# Console summary
# =============================================================================
cat("\n", strrep("=", 70), "\n", sep = "")
cat("SUMMARY (one row per assemblage)\n")
cat(strrep("=", 70), "\n", sep = "")
print(metrics_df %>% mutate(across(where(is.numeric), ~ round(.x, 4))) %>% as.data.frame())
cat("\nDone.\n")
