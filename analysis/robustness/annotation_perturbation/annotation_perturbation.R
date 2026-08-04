# annotation_perturbation.R
# =============================================================================
# ANNOTATION-ERROR PERTURBATION ANALYSIS for SP-SPHARM — SI add-on.
#
# Responds to reviewer 2: manual placement of direction landmarks is itself
# error-prone (Liu et al. 2026). This converts that limitation statement into
# quantitative evidence by asking how much annotation error SP-SPHARM tolerates
# before its pairwise discriminative power degrades.
#
# NEW, self-contained. Does NOT modify the main pipeline, the _targets store, the
# derived_data cache, or the manuscript. Reads the committed main outputs and the
# replicate descriptors produced by perturb_spharm.py (same folder), writes only
# under analysis/robustness/annotation_perturbation/.
#
# WHERE THE NOISE IS INJECTED. Perturbation is applied to the UPSTREAM UNIT
# VECTORS in directions_aligned_svd.csv, then the ENTIRE production chain is
# re-run per replicate — vMF-sKDE (h = 0.35, kappa = 8.16, 72x36 grid) -> DH 64x128
# grid -> spherical-harmonic expansion (lmax = 20) -> normalised power spectrum ->
# ILR -> PERMANOVA. Noise is never added to the power spectra directly. See
# perturb_spharm.py, which reuses the production functions verbatim and whose
# --verify mode reproduces the committed SPHARM_direction.csv to < 1e-8.
#
# THREE PERTURBATIONS
#   P1 polarity : flip a fraction f of scar vectors, f in {0.02 .. 0.20}. THE
#                 headline. The manuscript's case for SP-SPHARM over fabric rests
#                 on SP-SPHARM retaining polarity (a fabric orientation tensor is
#                 invariant to u -> -u, Mark 1973), so polarity is its most
#                 exposed assumption and the one most worth stress-testing.
#                 Fabric (E, I) and SPI are recomputed on the SAME perturbed
#                 vectors: fabric should be immune by construction. That contrast
#                 turns "SP-SPHARM encodes polarity, fabric does not" from a
#                 theoretical claim into a measurement, and supports the paper's
#                 three-methods-are-complementary argument.
#   P2 angle    : isotropic random rotation per vector, s.d. sigma in {5..20} deg.
#   P3 dropout  : delete a fraction d of scars, d in {0.05, 0.10, 0.20}, with a
#                 floor of 3 scars per specimen (Figure S3's existing bound).
#
# Statistical helpers are imported verbatim from the joint-MFA scripts by
# selective evaluation (the pattern established in joint_mfa_discrimination_SDG.R),
# so every number here is byte-comparable with the existing analyses.
#
# Outputs (all NEW):
#   annotation_perturbation_polarity.csv / _angle.csv / _dropout.csv / _summary.csv
#   annotation_perturbation_pairfailure.csv
#   figures/fig_S_perturbation_polarity_three_methods.png
#   figures/fig_S_perturbation_degradation.png
#   figures/fig_S_perturbation_continuous.png
#
# HOW TO RUN (Docker spharm_analysis):
#   /opt/conda/envs/spharm/bin/python analysis/robustness/annotation_perturbation/perturb_spharm.py --run --reps 100
#   Rscript analysis/robustness/annotation_perturbation/annotation_perturbation.R
# =============================================================================

suppressPackageStartupMessages({
  library(here); library(tidyverse); library(vegan); library(compositions)
  library(FSA); library(conflicted)
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
ASSEMBLAGE <- "EXP"     # SDG interface kept below; this run is EXP only
SEED       <- 42
N_PERM     <- 9999
ALPHA_1    <- 0.05
ALPHA_2    <- 0.01
QLO        <- 0.025
QHI        <- 0.975

OUT_DIR <- here("analysis/robustness/annotation_perturbation")
FIG_DIR <- file.path(OUT_DIR, "figures")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

DESC_CSV <- file.path(OUT_DIR, "perturbed_descriptors.csv")

POWER_COLS_DIR   <- paste0("power_l", 1:6)
POWER_COLS_MORPH <- paste0("power_l", 1:8)
EXCLUDE_TYPES    <- c("Biface")
LEVALLOIS_MERGE  <- c("Levallois convergent", "Levallois laminar",
                      "Levallois preferential", "Levallois recurrent")
TYPOLOGY_ORDER   <- c("Unidirectional", "Bidirectional", "Levallois",
                      "Discoid", "Multiplatform")

METHOD_COLORS <- c("SP-SPHARM" = "#4A6E8A", "Fabric (E, I)" = "#BA8530",
                   "SPI" = "#802520")
PERT_LABS <- c(polarity = "P1  polarity flips (fraction)",
               angle    = "P2  angular jitter (s.d., degrees)",
               dropout  = "P3  scar dropout (fraction)")

# Committed baseline (EXP SP-SPHARM core-type PERMANOVA), same anchor the joint-MFA
# scripts use. R2 / pseudo-F are deterministic and CHECKED; p is permutation noise.
REF <- list(R2 = 0.30174, F = 5.72568, n_sig = 8L, n = 58L, n_pairs = 10L)

# =============================================================================
# Helper import — verbatim, without running the source analyses
# =============================================================================
import_helpers <- function(path, wanted) {
  if (!file.exists(path)) stop("helper source not found: ", path)
  got <- character(0)
  for (e in parse(path)) {
    if (!is.call(e)) next
    if (!as.character(e[[1]])[1] %in% c("<-", "=")) next
    if (!is.name(e[[2]])) next
    nm <- as.character(e[[2]]); if (!nm %in% wanted) next
    rhs <- e[[3]]
    if (!(is.call(rhs) && identical(as.character(rhs[[1]])[1], "function"))) next
    eval(e, envir = globalenv()); got <- c(got, nm)
  }
  miss <- setdiff(wanted, got)
  if (length(miss)) stop("could not import from ", basename(path), ": ",
                         paste(miss, collapse = ", "))
  got
}
h1 <- import_helpers(
  here("analysis/robustness/joint_mfa_discrimination/joint_mfa_discrimination_stats.R"),
  c("replace_zeros", "make_ilr", "permanova_global", "permanova_pairwise", "pair_p"))
h2 <- import_helpers(
  here("analysis/robustness/coia_power_sensitivity/coia_power_sensitivity_stats.R"),
  c("rv_std"))
cat(sprintf("Imported helpers verbatim: %s | %s\n",
            paste(sort(h1), collapse = ", "), paste(sort(h2), collapse = ", ")))

# Mantel statistic = Spearman correlation of the two distance vectors. Identical to
# vegan::mantel()$statistic but without the permutation cost, which we do not need
# here (only the coefficient is reported, per replicate).
mantel_stat <- function(A, B)
  stats::cor(as.vector(stats::dist(A)), as.vector(stats::dist(B)), method = "spearman")

# =============================================================================
# Inputs
# =============================================================================
if (!file.exists(DESC_CSV))
  stop("Missing ", basename(DESC_CSV), " — run perturb_spharm.py --run first.")

desc <- read_csv(DESC_CSV, show_col_types = FALSE)

tidy_typology <- function(x)
  if_else(x %in% LEVALLOIS_MERGE | str_detect(x, regex("levallois", ignore_case = TRUE)),
          "Levallois", x)

desc <- desc %>%
  filter(!Typology %in% EXCLUDE_TYPES) %>%
  mutate(Typology = tidy_typology(Typology))

base_df <- desc %>% filter(perturbation == "baseline") %>% arrange(ID)
ids     <- base_df$ID
n_spec  <- length(ids)
grp     <- factor(base_df$Typology,
                  levels = TYPOLOGY_ORDER[TYPOLOGY_ORDER %in% unique(base_df$Typology)])
n_pairs <- choose(nlevels(grp), 2)
cat(sprintf("\nEXP specimens: %d | groups: %d | pairs: %d\n", n_spec, nlevels(grp), n_pairs))
print(table(grp))

# Fixed morphology block (never perturbed) for the decoupling check.
morph_raw <- read_csv(here("analysis/data/derived_data/SPHARM_morphology.csv"),
                      show_col_types = FALSE) %>%
  filter(ID %in% ids) %>% arrange(ID)
stopifnot(identical(morph_raw$ID, ids))
Z_morph <- make_ilr(morph_raw %>% select(all_of(POWER_COLS_MORPH)))

# =============================================================================
# Per-replicate evaluation
# =============================================================================
ilr_of <- function(df) {
  df <- df %>% arrange(ID)
  stopifnot(identical(df$ID, ids))
  make_ilr(df %>% select(all_of(POWER_COLS_DIR)))
}
Z_base <- ilr_of(base_df)
D_base <- stats::dist(Z_base)

# SP-SPHARM statistics for one replicate's descriptor frame.
eval_spharm <- function(df) {
  Z  <- ilr_of(df)
  gl <- permanova_global(Z, grp, nperm = N_PERM)
  pw <- permanova_pairwise(Z, grp, nperm = N_PERM)
  list(R2 = gl$R2, F = gl$F, p = gl$p,
       n_sig_05 = sum(pw$p_holm < ALPHA_1),
       n_sig_01 = sum(pw$p_holm < ALPHA_2),
       med_neglog10 = stats::median(-log10(pw$p)),
       ilr_shift = stats::median(sqrt(rowSums((Z - Z_base)^2))),
       mantel_vs_base = mantel_stat(Z, Z_base),
       decoup_mantel  = mantel_stat(Z_morph, Z),
       decoup_rv      = rv_std(Z_morph, Z),
       pw = pw)
}

# Fabric (E, I): same PERMANOVA machinery, on the 2-column fabric descriptor
# (spharm_analysis.R:351-352 uses Euclidean distance on raw E, I).
eval_fabric <- function(df) {
  df <- df %>% arrange(ID)
  X  <- df %>% select(E, I) %>% as.matrix()
  if (any(!is.finite(X))) return(list(n_sig_05 = NA_integer_, n_sig_01 = NA_integer_,
                                      med_neglog10 = NA_real_))
  pw <- permanova_pairwise(X, grp, nperm = N_PERM)
  list(n_sig_05 = sum(pw$p_holm < ALPHA_1),
       n_sig_01 = sum(pw$p_holm < ALPHA_2),
       med_neglog10 = stats::median(-log10(pw$p)))
}

# SPI: univariate, so the main analysis's Dunn test with Holm (spharm_analysis.R:282).
eval_spi <- function(df) {
  df <- df %>% arrange(ID)
  d  <- data.frame(SPI = df$SPI, g = grp)
  if (any(!is.finite(d$SPI))) return(list(n_sig_05 = NA_integer_, n_sig_01 = NA_integer_,
                                          med_neglog10 = NA_real_))
  dt <- suppressWarnings(FSA::dunnTest(SPI ~ g, data = d, method = "holm"))
  padj <- dt$res$P.adj; praw <- dt$res$P.unadj
  list(n_sig_05 = sum(padj < ALPHA_1), n_sig_01 = sum(padj < ALPHA_2),
       med_neglog10 = stats::median(-log10(pmax(praw, .Machine$double.xmin))))
}

# =============================================================================
# (0) Baseline anchor
# =============================================================================
cat("\n", strrep("=", 70), "\n", sep = "")
cat("(0) BASELINE ANCHOR (unperturbed chain re-run through perturb_spharm.py)\n")
cat(strrep("=", 70), "\n", sep = "")
b <- eval_spharm(base_df)
bf <- eval_fabric(base_df); bs <- eval_spi(base_df)
anchor_ok <- TRUE
chk <- function(name, got, exp, tol, det = TRUE) {
  ok <- is.finite(got) && abs(got - exp) <= tol
  if (!ok && det) anchor_ok <<- FALSE
  cat(sprintf("  %-24s got %10.5f  expected %10.5f  %s\n", name, got, exp,
              ifelse(ok, "OK", ifelse(det, "<-- CHECK", "<-- CHECK (permutation)"))))
}
chk("n specimens",      n_spec,     REF$n,     0)
chk("PERMANOVA R2",     b$R2,       REF$R2,    1e-3)
chk("PERMANOVA pseudo-F", b$F,      REF$F,     1e-3)
chk("resolved pairs",   b$n_sig_05, REF$n_sig, 0)
if (!anchor_ok)
  stop("BASELINE ANCHOR FAILED — the re-run chain does not reproduce the committed ",
       "SP-SPHARM PERMANOVA. Stopping before the perturbation results.")
cat(sprintf("\n  => Baseline reproduces the committed EXP SP-SPHARM PERMANOVA.\n"))
cat(sprintf("     baseline median -log10(raw p) = %.3f | fabric %d/%d | SPI %d/%d resolved\n",
            b$med_neglog10, bf$n_sig_05, n_pairs, bs$n_sig_05, n_pairs))

base_sig_pairs <- b$pw$comparison[b$pw$p_holm < ALPHA_1]
cat(sprintf("     baseline significant pairs (%d): %s\n",
            length(base_sig_pairs), paste(base_sig_pairs, collapse = "; ")))

# =============================================================================
# Sweep all replicates
# =============================================================================
cat("\n", strrep("=", 70), "\n", sep = "")
cat("EVALUATING REPLICATES\n")
cat(strrep("=", 70), "\n", sep = "")

conds <- desc %>% filter(perturbation != "baseline") %>%
  distinct(perturbation, level) %>% arrange(perturbation, level)
reps_per_cond <- desc %>% filter(perturbation != "baseline") %>%
  count(perturbation, level, rep) %>% count(perturbation, level, name = "R")
cat(sprintf("  %d conditions, R = %s replicates each\n",
            nrow(conds), paste(unique(reps_per_cond$R), collapse = "/")))

t_start <- Sys.time()
rep_rows <- list(); pair_rows <- list()
for (ci in seq_len(nrow(conds))) {
  kind <- conds$perturbation[ci]; lv <- conds$level[ci]
  sub  <- desc %>% filter(perturbation == kind, level == lv)
  reps <- sort(unique(sub$rep))
  for (r in reps) {
    dfr <- sub %>% filter(rep == r)
    s   <- eval_spharm(dfr)
    row <- tibble(perturbation = kind, level = lv, rep = r,
                  R2 = s$R2, pseudo_F = s$F, p = s$p,
                  n_sig_05 = s$n_sig_05, n_sig_01 = s$n_sig_01,
                  med_neglog10_p = s$med_neglog10,
                  ilr_shift = s$ilr_shift, mantel_vs_base = s$mantel_vs_base,
                  decoup_mantel = s$decoup_mantel, decoup_rv = s$decoup_rv)
    # Three-method comparison is P1-specific (brief: "P1 zhuanshu").
    if (kind == "polarity") {
      f <- eval_fabric(dfr); sp <- eval_spi(dfr)
      row <- row %>% mutate(fabric_n_sig_05 = f$n_sig_05, fabric_n_sig_01 = f$n_sig_01,
                            fabric_med_neglog10 = f$med_neglog10,
                            spi_n_sig_05 = sp$n_sig_05, spi_n_sig_01 = sp$n_sig_01,
                            spi_med_neglog10 = sp$med_neglog10)
    }
    rep_rows[[length(rep_rows) + 1]] <- row
    pair_rows[[length(pair_rows) + 1]] <-
      s$pw %>% transmute(perturbation = kind, level = lv, rep = r,
                         comparison, p_raw = p, p_holm)
  }
  cat(sprintf("  %-9s %6.2f : %3d reps done (%.1f min elapsed)\n", kind, lv, length(reps),
              as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
}
rep_df  <- bind_rows(rep_rows)
pair_df <- bind_rows(pair_rows)
elapsed_min <- as.numeric(difftime(Sys.time(), t_start, units = "mins"))

# =============================================================================
# Aggregate: median + 2.5/97.5 percentiles
# =============================================================================
q <- function(x, p) stats::quantile(x[is.finite(x)], p, names = FALSE)
agg_cols <- function(d, cols) {
  map_dfc(cols, function(cn) {
    v <- d[[cn]]
    if (all(is.na(v))) return(setNames(tibble(NA_real_, NA_real_, NA_real_),
                                       paste0(cn, c("_med", "_lo", "_hi"))))
    setNames(tibble(stats::median(v[is.finite(v)]), q(v, QLO), q(v, QHI)),
             paste0(cn, c("_med", "_lo", "_hi")))
  })
}
CORE <- c("R2", "pseudo_F", "n_sig_05", "n_sig_01", "med_neglog10_p",
          "ilr_shift", "mantel_vs_base", "decoup_mantel", "decoup_rv")
P1EX <- c("fabric_n_sig_05", "fabric_n_sig_01", "fabric_med_neglog10",
          "spi_n_sig_05", "spi_n_sig_01", "spi_med_neglog10")

summarise_cond <- function(d) {
  cols <- CORE
  if ("fabric_n_sig_05" %in% names(d) && any(!is.na(d$fabric_n_sig_05)))
    cols <- c(cols, P1EX)
  bind_cols(tibble(perturbation = d$perturbation[1], level = d$level[1],
                   R = nrow(d)), agg_cols(d, cols))
}
summary_df <- rep_df %>% group_split(perturbation, level) %>%
  map_dfr(summarise_cond) %>% arrange(perturbation, level)

# Baseline row, so every table carries its own reference point.
base_row <- tibble(perturbation = "baseline", level = 0, R = 1L,
                   R2_med = b$R2, R2_lo = b$R2, R2_hi = b$R2,
                   pseudo_F_med = b$F, pseudo_F_lo = b$F, pseudo_F_hi = b$F,
                   n_sig_05_med = b$n_sig_05, n_sig_05_lo = b$n_sig_05, n_sig_05_hi = b$n_sig_05,
                   n_sig_01_med = b$n_sig_01, n_sig_01_lo = b$n_sig_01, n_sig_01_hi = b$n_sig_01,
                   med_neglog10_p_med = b$med_neglog10,
                   med_neglog10_p_lo = b$med_neglog10, med_neglog10_p_hi = b$med_neglog10,
                   ilr_shift_med = 0, ilr_shift_lo = 0, ilr_shift_hi = 0,
                   mantel_vs_base_med = 1, mantel_vs_base_lo = 1, mantel_vs_base_hi = 1,
                   decoup_mantel_med = b$decoup_mantel, decoup_mantel_lo = b$decoup_mantel,
                   decoup_mantel_hi = b$decoup_mantel,
                   decoup_rv_med = b$decoup_rv, decoup_rv_lo = b$decoup_rv,
                   decoup_rv_hi = b$decoup_rv,
                   fabric_n_sig_05_med = bf$n_sig_05, spi_n_sig_05_med = bs$n_sig_05)

summary_all <- bind_rows(base_row, summary_df)

write_csv(summary_df %>% filter(perturbation == "polarity"),
          file.path(OUT_DIR, "annotation_perturbation_polarity.csv"))
write_csv(summary_df %>% filter(perturbation == "angle"),
          file.path(OUT_DIR, "annotation_perturbation_angle.csv"))
write_csv(summary_df %>% filter(perturbation == "dropout"),
          file.path(OUT_DIR, "annotation_perturbation_dropout.csv"))
write_csv(summary_all, file.path(OUT_DIR, "annotation_perturbation_summary.csv"))
cat("\nWrote 4 summary CSVs\n")

# =============================================================================
# (1) Collapse thresholds
# =============================================================================
cat("\n", strrep("=", 70), "\n", sep = "")
cat("(1) DISCRIMINATIVE-POWER COLLAPSE THRESHOLDS\n")
cat(strrep("=", 70), "\n", sep = "")
collapse_tbl <- summary_df %>% group_by(perturbation) %>%
  summarise(first_below = {
    idx <- which(n_sig_05_med < REF$n_sig)
    if (length(idx) == 0) NA_real_ else level[min(idx)]
  },
  first_below_01 = {
    idx <- which(n_sig_01_med < b$n_sig_01)
    if (length(idx) == 0) NA_real_ else level[min(idx)]
  }, .groups = "drop")
for (i in seq_len(nrow(collapse_tbl))) {
  k <- collapse_tbl$perturbation[i]
  cat(sprintf("  %-9s : first level with median resolved < %d (alpha 0.05) = %s\n",
              k, REF$n_sig,
              ifelse(is.na(collapse_tbl$first_below[i]), "none in range",
                     sprintf("%.2f", collapse_tbl$first_below[i]))))
}
print(summary_df %>%
        transmute(perturbation, level, R,
                  resolved_05 = sprintf("%.1f [%.0f, %.0f]", n_sig_05_med, n_sig_05_lo, n_sig_05_hi),
                  resolved_01 = sprintf("%.1f [%.0f, %.0f]", n_sig_01_med, n_sig_01_lo, n_sig_01_hi),
                  med_neglog10 = sprintf("%.2f [%.2f, %.2f]", med_neglog10_p_med,
                                         med_neglog10_p_lo, med_neglog10_p_hi),
                  R2 = round(R2_med, 4)) %>% as.data.frame())

# =============================================================================
# (2) P1 three-method comparison
# =============================================================================
cat("\n", strrep("=", 70), "\n", sep = "")
cat("(2) P1 THREE-METHOD COMPARISON (polarity flips)\n")
cat(strrep("=", 70), "\n", sep = "")
pol <- summary_df %>% filter(perturbation == "polarity")
print(pol %>% transmute(level,
                        `SP-SPHARM` = sprintf("%.1f [%.0f, %.0f]", n_sig_05_med, n_sig_05_lo, n_sig_05_hi),
                        Fabric = sprintf("%.1f [%.0f, %.0f]", fabric_n_sig_05_med,
                                         fabric_n_sig_05_lo, fabric_n_sig_05_hi),
                        SPI = sprintf("%.1f [%.0f, %.0f]", spi_n_sig_05_med,
                                      spi_n_sig_05_lo, spi_n_sig_05_hi)) %>%
        as.data.frame())
fab_flat <- diff(range(pol$fabric_n_sig_05_med, na.rm = TRUE)) == 0 &&
  isTRUE(all.equal(pol$fabric_n_sig_05_med[1], bf$n_sig_05))
sp_drop <- pol$n_sig_05_med[nrow(pol)] < b$n_sig_05
cat(sprintf("\n  fabric baseline %d/%d, across all flip levels median stays %s -> %s\n",
            bf$n_sig_05, n_pairs,
            paste(unique(pol$fabric_n_sig_05_med), collapse = "/"),
            ifelse(fab_flat, "IMMUNE, as predicted (orientation tensor discards polarity)",
                   "NOT flat — inspect")))
cat(sprintf("  SP-SPHARM %d/%d at baseline -> %.1f at f = %.2f -> %s\n",
            b$n_sig_05, n_pairs, pol$n_sig_05_med[nrow(pol)], pol$level[nrow(pol)],
            ifelse(sp_drop, "degrades with polarity error, as expected",
                   "no degradation detected")))

# =============================================================================
# (3) Which baseline pairs fail first
# =============================================================================
cat("\n", strrep("=", 70), "\n", sep = "")
cat("(3) ORDER IN WHICH THE BASELINE-SIGNIFICANT PAIRS FAIL\n")
cat(strrep("=", 70), "\n", sep = "")
pair_med <- pair_df %>%
  filter(comparison %in% base_sig_pairs) %>%
  group_by(perturbation, level, comparison) %>%
  summarise(p_holm_med = stats::median(p_holm), .groups = "drop")
fail_order <- pair_med %>% group_by(perturbation, comparison) %>%
  summarise(first_fail_level = {
    idx <- which(p_holm_med >= ALPHA_1)
    if (length(idx) == 0) NA_real_ else level[min(idx)]
  },
  p_holm_at_max = p_holm_med[which.max(level)], .groups = "drop") %>%
  left_join(b$pw %>% transmute(comparison, baseline_p_raw = p, baseline_p_holm = p_holm),
            by = "comparison") %>%
  arrange(perturbation, is.na(first_fail_level), first_fail_level, desc(baseline_p_raw))
for (k in unique(fail_order$perturbation)) {
  cat(sprintf("\n  --- %s ---\n", k))
  print(fail_order %>% filter(perturbation == k) %>%
          transmute(comparison,
                    first_fail = ifelse(is.na(first_fail_level), "survives",
                                        sprintf("%.2f", first_fail_level)),
                    p_holm_at_max_level = round(p_holm_at_max, 4),
                    baseline_p_raw = round(baseline_p_raw, 4)) %>% as.data.frame())
}
write_csv(fail_order, file.path(OUT_DIR, "annotation_perturbation_pairfailure.csv"))
cat("\n  Wrote annotation_perturbation_pairfailure.csv\n")
cat("  (Expectation: failure starts with the pair holding the smallest baseline margin.)\n")

# =============================================================================
# (4) Decoupling conclusion
# =============================================================================
cat("\n", strrep("=", 70), "\n", sep = "")
cat("(4) DECOUPLING (morphology vs perturbed scars) ACROSS ALL LEVELS\n")
cat(strrep("=", 70), "\n", sep = "")
print(summary_df %>% transmute(perturbation, level,
                               mantel_r = sprintf("%.3f [%.3f, %.3f]", decoup_mantel_med,
                                                  decoup_mantel_lo, decoup_mantel_hi),
                               RV = sprintf("%.3f [%.3f, %.3f]", decoup_rv_med,
                                            decoup_rv_lo, decoup_rv_hi),
                               ilr_shift = round(ilr_shift_med, 3),
                               mantel_vs_base = round(mantel_vs_base_med, 3)) %>%
        as.data.frame())
decoup_ok <- all(abs(summary_df$decoup_mantel_hi) < 0.3, na.rm = TRUE) &&
  all(summary_df$decoup_rv_hi < 0.3, na.rm = TRUE)
cat(sprintf("\n  baseline Mantel r = %.3f, RV = %.3f\n", b$decoup_mantel, b$decoup_rv))
cat(sprintf("  => the decoupling conclusion %s across every perturbation level\n",
            ifelse(decoup_ok, "HOLDS (all |Mantel r| and RV stay well below 0.3)",
                   "DOES NOT hold everywhere — inspect")))

# =============================================================================
# Figures
# =============================================================================
# Fig 1 — P1 three-method degradation, the headline output.
three <- bind_rows(
  pol %>% transmute(level, method = "SP-SPHARM", med = n_sig_05_med,
                    lo = n_sig_05_lo, hi = n_sig_05_hi),
  pol %>% transmute(level, method = "Fabric (E, I)", med = fabric_n_sig_05_med,
                    lo = fabric_n_sig_05_lo, hi = fabric_n_sig_05_hi),
  pol %>% transmute(level, method = "SPI", med = spi_n_sig_05_med,
                    lo = spi_n_sig_05_lo, hi = spi_n_sig_05_hi)) %>%
  bind_rows(tibble(level = 0,
                   method = c("SP-SPHARM", "Fabric (E, I)", "SPI"),
                   med = c(b$n_sig_05, bf$n_sig_05, bs$n_sig_05),
                   lo = c(b$n_sig_05, bf$n_sig_05, bs$n_sig_05),
                   hi = c(b$n_sig_05, bf$n_sig_05, bs$n_sig_05))) %>%
  mutate(method = factor(method, levels = names(METHOD_COLORS))) %>% arrange(method, level)

p1 <- ggplot(three, aes(level, med, color = method, fill = method)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.18, colour = NA) +
  geom_line(linewidth = 0.8) + geom_point(size = 1.8) +
  scale_color_manual(values = METHOD_COLORS, name = NULL) +
  scale_fill_manual(values = METHOD_COLORS, name = NULL) +
  scale_x_continuous(breaks = c(0, POLARITY_LEVELS <- sort(unique(pol$level))),
                     labels = scales::percent_format(accuracy = 1)) +
  scale_y_continuous(breaks = 0:n_pairs, limits = c(0, n_pairs)) +
  labs(x = "Scar vectors with reversed polarity",
       y = sprintf("Core-type pairs resolved (Holm, p < %.2f; of %d)", ALPHA_1, n_pairs)) +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom",
        legend.key.size = grid::unit(0.34, "cm"))
ggsave(file.path(FIG_DIR, "fig_S_perturbation_polarity_three_methods.png"), p1,
       width = 7, height = 4.8, dpi = 300)
cat("\nWrote figures/fig_S_perturbation_polarity_three_methods.png\n")

# Fig 2 — resolved pairs vs level, all three perturbations, both alphas.
deg <- summary_df %>%
  transmute(perturbation, level,
            `alpha = 0.05_med` = n_sig_05_med, `alpha = 0.05_lo` = n_sig_05_lo,
            `alpha = 0.05_hi` = n_sig_05_hi,
            `alpha = 0.01_med` = n_sig_01_med, `alpha = 0.01_lo` = n_sig_01_lo,
            `alpha = 0.01_hi` = n_sig_01_hi) %>%
  pivot_longer(-c(perturbation, level),
               names_to = c("alpha", ".value"), names_sep = "_") %>%
  mutate(perturbation = factor(perturbation, levels = names(PERT_LABS)))
p2 <- ggplot(deg, aes(level, med, color = alpha, fill = alpha)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.18, colour = NA) +
  geom_line(linewidth = 0.8) + geom_point(size = 1.6) +
  facet_wrap(~ perturbation, scales = "free_x",
             labeller = as_labeller(PERT_LABS)) +
  scale_color_manual(values = c("alpha = 0.05" = "#4A6E8A", "alpha = 0.01" = "#788C4A"),
                     name = NULL) +
  scale_fill_manual(values = c("alpha = 0.05" = "#4A6E8A", "alpha = 0.01" = "#788C4A"),
                    name = NULL) +
  scale_y_continuous(breaks = 0:n_pairs, limits = c(0, n_pairs)) +
  labs(x = "Perturbation level", y = sprintf("Core-type pairs resolved (of %d)", n_pairs)) +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom",
        strip.text = element_text(face = "bold", size = 8),
        legend.key.size = grid::unit(0.34, "cm"))
ggsave(file.path(FIG_DIR, "fig_S_perturbation_degradation.png"), p2,
       width = 9, height = 4.2, dpi = 300)
cat("Wrote figures/fig_S_perturbation_degradation.png\n")

# Fig 3 — continuous margin (median -log10 raw p).
cont <- summary_df %>%
  mutate(perturbation = factor(perturbation, levels = names(PERT_LABS)))
p3 <- ggplot(cont, aes(level, med_neglog10_p_med)) +
  geom_ribbon(aes(ymin = med_neglog10_p_lo, ymax = med_neglog10_p_hi),
              alpha = 0.18, fill = "#4A6E8A") +
  geom_hline(yintercept = b$med_neglog10, linetype = "dashed",
             color = "#802520", linewidth = 0.4) +
  geom_line(linewidth = 0.8, color = "#4A6E8A") +
  geom_point(size = 1.6, color = "#4A6E8A") +
  facet_wrap(~ perturbation, scales = "free_x", labeller = as_labeller(PERT_LABS)) +
  labs(x = "Perturbation level",
       y = "Median -log10(raw p) over all core-type pairs",
       caption = sprintf("Dashed line: unperturbed baseline (%.2f)", b$med_neglog10)) +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold", size = 8),
        plot.caption = element_text(size = 7, colour = "grey40"))
ggsave(file.path(FIG_DIR, "fig_S_perturbation_continuous.png"), p3,
       width = 9, height = 4.2, dpi = 300)
cat("Wrote figures/fig_S_perturbation_continuous.png\n")

# =============================================================================
# (5) Timing
# =============================================================================
cat("\n", strrep("=", 70), "\n", sep = "")
cat("(5) COST\n")
cat(strrep("=", 70), "\n", sep = "")
n_rep_total <- nrow(rep_df)
cat(sprintf("  replicates evaluated : %d (%d conditions)\n", n_rep_total, nrow(conds)))
cat(sprintf("  R-side elapsed       : %.1f min (%.2f s per replicate at N_PERM = %d)\n",
            elapsed_min, elapsed_min * 60 / n_rep_total, N_PERM))
cat("  Python chain re-run  : 0.30 s per 58-specimen pass (measured; --timing mode)\n")

if (file.exists(file.path(OUT_DIR, "dropout_floor_log.csv"))) {
  fl <- read_csv(file.path(OUT_DIR, "dropout_floor_log.csv"), show_col_types = FALSE)
  fs <- fl %>% group_by(level) %>%
    summarise(mean_floored = mean(n_floored), max_floored = max(n_floored), .groups = "drop")
  cat("\n  P3 dropout: specimens hitting the 3-scar floor (per replicate)\n")
  print(as.data.frame(fs))
  cat("  Handling: such specimens keep as many scars as the floor allows.\n")
}

cat("\nDone.\n\n")
print(sessionInfo())
