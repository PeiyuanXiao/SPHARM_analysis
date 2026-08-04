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
# TWO PERTURBATIONS, both with the SP-SPHARM / fabric / SPI three-method contrast
#   P1 polarity : flip a fraction f of scar vectors, f in {0.02 .. 0.20}. The
#                 manuscript's case for SP-SPHARM over fabric rests on SP-SPHARM
#                 retaining polarity (a fabric orientation tensor is invariant to
#                 u -> -u, Mark 1973), so polarity is its most exposed assumption.
#                 NOTE what this contrast is and is not: fabric's immunity here is a
#                 mathematical identity, not an empirical finding — sum(u u') is
#                 literally unchanged by a sign flip, so E and I come out bit-
#                 identical. It is worth reporting as a check, but it discovers
#                 nothing.
#   P2 angle    : isotropic random rotation per vector, s.d. sigma in {5..20} deg.
#                 THIS is the real experiment. No method has analytic immunity to
#                 small rotations, so the ordering of degradation rates is not known
#                 in advance and is determined purely by how each method encodes
#                 direction. It therefore tests a concrete mechanistic claim
#                 (see PREDICTION below).
#
# A scar-dropout perturbation (P3) was tried and REMOVED: the sparsest EXP specimen
# carries 10 scars, so even 20% dropout leaves 8 and the 3-scar floor never once
# triggered. It only ever probed the data-rich regime — "what if a specimen with
# plenty of scars loses a few?" — whose answer is necessarily "little". The question
# that matters, how few scars suffice, is a downsampling analysis (k = 3/5/8/10/15)
# and belongs in its own script.
#
# PREDICTION (declared before the angle contrast was run, and adjudicated at the end)
#   SP-SPHARM truncates the power spectrum at l = 1-6. Spherical-harmonic angular
#   resolution is roughly 180/l degrees, so l = 6 corresponds to about 30 deg:
#   jitter below that scale falls under the truncation and is filtered out by the
#   descriptor itself. The other two methods have no such explicit low-pass stage —
#   SPI is the resultant-length ratio of the raw unit vectors, with no angular
#   smoothing at all, and fabric's second-moment eigenvalues average somewhat but
#   impose no cutoff. Predicted degradation rate (by retention relative to each
#   method's own baseline):  SPI fastest > fabric > SP-SPHARM slowest.
#
# WHY RETENTION RATIOS. Baseline resolved-pair counts differ by method (SP-SPHARM
# 8/10, fabric 6/10, SPI 4/10), so absolute counts cannot rank degradation rates —
# a drop of one pair means something different at each baseline. Every three-method
# panel therefore reports (i) the absolute count, (ii) the continuous median
# -log10(p), and (iii) retention = value / that method's own baseline. Verdicts are
# adjudicated on (iii), using the continuous statistic, which is free of both the
# baseline mismatch and the integer-threshold steps.
#
# EACH METHOD KEEPS ITS OWN TEST, as in the main analysis: SP-SPHARM and fabric are
# tested by PERMANOVA (adonis2; spharm_analysis.R:352, 407), SPI by Kruskal-Wallis
# with Dunn post-hoc (spharm_analysis.R:282). No test type was changed to make the
# three comparable; the test used is stated in every table and figure caption.
#
# Statistical helpers are imported verbatim from the joint-MFA scripts by
# selective evaluation (the pattern established in joint_mfa_discrimination_SDG.R),
# so every number here is byte-comparable with the existing analyses.
#
# Outputs:
#   annotation_perturbation_polarity.csv / _angle.csv / _summary.csv
#   annotation_perturbation_angle_three_methods.csv
#   annotation_perturbation_pairfailure.csv
#   figures/fig_S_perturbation_polarity_three_methods.png
#   figures/fig_S_perturbation_angle_three_methods.png
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
               angle    = "P2  angular jitter (s.d., degrees)")

# Perturbations that get the SP-SPHARM / fabric / SPI contrast. Both, now: the same
# seeds, the same call path, the same per-replicate vector array for all three
# methods (perturb_spharm.py Engine.run perturbs U once per specimen and derives the
# KDE, SPI and E/I from that one array), so the polarity and angle contrasts are
# structurally identical and can be read side by side.
THREE_METHOD_KINDS <- c("polarity", "angle")

# Spherical-harmonic angular resolution at the truncation degree, for the mechanism
# claim: l = 6 resolves features no finer than about 180/6 = 30 degrees.
LMAX_KEEP    <- 6
SH_RES_DEG   <- 180 / LMAX_KEEP
# Committed bandwidth-sensitivity trend for the R2 cross-check
# (bandwidth_sensitivity_metrics.csv: h = 0.35 -> 0.50).
BW_REF <- list(h_lo = 0.35, h_hi = 0.50, R2_lo = 0.30174, R2_hi = 0.32879,
               nsig_lo = 8L, nsig_hi = 8L)

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
  if (any(!is.finite(X)))
    return(list(R2 = NA_real_, F = NA_real_, p = NA_real_,
                n_sig_05 = NA_integer_, n_sig_01 = NA_integer_, med_neglog10 = NA_real_))
  gl <- permanova_global(X, grp, nperm = N_PERM)
  pw <- permanova_pairwise(X, grp, nperm = N_PERM)
  list(R2 = gl$R2, F = gl$F, p = gl$p,
       n_sig_05 = sum(pw$p_holm < ALPHA_1),
       n_sig_01 = sum(pw$p_holm < ALPHA_2),
       med_neglog10 = stats::median(-log10(pw$p)))
}

# SPI: univariate, so the main analysis's Dunn test with Holm (spharm_analysis.R:282).
# SPI is univariate, so it keeps the main analysis's test: Kruskal-Wallis for the
# global effect, Dunn with Holm for the pairwise comparisons. Deliberately NOT
# converted to a PERMANOVA — the point is to compare each method as the paper
# actually uses it. The "global" columns therefore hold the KW chi-squared, not a
# pseudo-F, and R2 is undefined for it.
eval_spi <- function(df) {
  df <- df %>% arrange(ID)
  d  <- data.frame(SPI = df$SPI, g = grp)
  if (any(!is.finite(d$SPI)))
    return(list(R2 = NA_real_, F = NA_real_, p = NA_real_,
                n_sig_05 = NA_integer_, n_sig_01 = NA_integer_, med_neglog10 = NA_real_))
  kw <- stats::kruskal.test(SPI ~ g, data = d)
  # dunnTest prints the KW test and the comparison table to stderr as a side effect;
  # with one call per replicate that would bury the report, so both streams are sunk.
  dt <- NULL
  invisible(utils::capture.output(
    invisible(utils::capture.output(
      dt <- suppressWarnings(FSA::dunnTest(SPI ~ g, data = d, method = "holm")),
      type = "message")),
    type = "output"))
  padj <- dt$res$P.adj; praw <- dt$res$P.unadj
  list(R2 = NA_real_, F = as.numeric(kw$statistic), p = as.numeric(kw$p.value),
       n_sig_05 = sum(padj < ALPHA_1), n_sig_01 = sum(padj < ALPHA_2),
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
    # Three-method contrast, for BOTH perturbations. fabric and SPI are evaluated on
    # `dfr` — the very same replicate rows, hence the very same perturbed vectors —
    # so no additional sampling noise enters the comparison.
    if (kind %in% THREE_METHOD_KINDS) {
      f <- eval_fabric(dfr); sp <- eval_spi(dfr)
      row <- row %>% mutate(fabric_R2 = f$R2, fabric_F = f$F, fabric_p = f$p,
                            fabric_n_sig_05 = f$n_sig_05, fabric_n_sig_01 = f$n_sig_01,
                            fabric_med_neglog10 = f$med_neglog10,
                            spi_KW_chisq = sp$F, spi_p = sp$p,
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
P1EX <- c("fabric_R2", "fabric_F", "fabric_p",
          "fabric_n_sig_05", "fabric_n_sig_01", "fabric_med_neglog10",
          "spi_KW_chisq", "spi_p",
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
write_csv(summary_all, file.path(OUT_DIR, "annotation_perturbation_summary.csv"))
cat("\nWrote 3 summary CSVs\n")

# =============================================================================
# A. Collapse thresholds
# =============================================================================
cat("\n", strrep("=", 70), "\n", sep = "")
cat("A. DISCRIMINATIVE-POWER COLLAPSE THRESHOLDS\n")
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
# B. P1 three-method comparison
# =============================================================================
cat("\n", strrep("=", 70), "\n", sep = "")
cat("B. P1 THREE-METHOD COMPARISON (polarity flips)\n")
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
# C. Which baseline pairs fail first
# =============================================================================
cat("\n", strrep("=", 70), "\n", sep = "")
cat("C. ORDER IN WHICH THE BASELINE-SIGNIFICANT PAIRS FAIL\n")
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
# D. Decoupling conclusion
# =============================================================================
cat("\n", strrep("=", 70), "\n", sep = "")
cat("D. DECOUPLING (morphology vs perturbed scars) ACROSS ALL LEVELS\n")
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
# ---- shared three-method assembly (identical for P1 and P2) ------------------
# Baselines, one per method, used both as the level-0 point and as the denominator
# of the retention ratio.
BASE_MET <- tibble(method = c("SP-SPHARM", "Fabric (E, I)", "SPI"),
                   base_count   = c(b$n_sig_05, bf$n_sig_05, bs$n_sig_05),
                   base_neglog  = c(b$med_neglog10, bf$med_neglog10, bs$med_neglog10))

# quantity = "count" (resolved pairs) or "neglog" (median -log10 raw p).
three_method_tbl <- function(kind, quantity = c("count", "neglog")) {
  quantity <- match.arg(quantity)
  d <- summary_df %>% filter(perturbation == kind)
  pick <- function(sp_pre, fa_pre, spi_pre)
    bind_rows(
      d %>% transmute(level, method = "SP-SPHARM",
                      med = .data[[paste0(sp_pre, "_med")]],
                      lo  = .data[[paste0(sp_pre, "_lo")]],
                      hi  = .data[[paste0(sp_pre, "_hi")]]),
      d %>% transmute(level, method = "Fabric (E, I)",
                      med = .data[[paste0(fa_pre, "_med")]],
                      lo  = .data[[paste0(fa_pre, "_lo")]],
                      hi  = .data[[paste0(fa_pre, "_hi")]]),
      d %>% transmute(level, method = "SPI",
                      med = .data[[paste0(spi_pre, "_med")]],
                      lo  = .data[[paste0(spi_pre, "_lo")]],
                      hi  = .data[[paste0(spi_pre, "_hi")]]))
  out <- if (quantity == "count")
    pick("n_sig_05", "fabric_n_sig_05", "spi_n_sig_05") else
      pick("med_neglog10_p", "fabric_med_neglog10", "spi_med_neglog10")
  base_col <- if (quantity == "count") "base_count" else "base_neglog"
  out %>%
    bind_rows(BASE_MET %>% transmute(level = 0, method,
                                     med = .data[[base_col]],
                                     lo = .data[[base_col]], hi = .data[[base_col]])) %>%
    left_join(BASE_MET %>% select(method, base = all_of(base_col)), by = "method") %>%
    mutate(method = factor(method, levels = names(METHOD_COLORS)),
           ret_med = med / base, ret_lo = lo / base, ret_hi = hi / base) %>%
    arrange(method, level)
}

# Fig 1 — P1 three-method degradation (counts; unchanged from the previous run).
three <- three_method_tbl("polarity", "count")

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

# =============================================================================
# Fig 1b — P2 angle three-method contrast: counts / continuous / retention
# =============================================================================
# Same colours, same method levels, same seeds and call path as the polarity
# figure, so the two can be read side by side. Three panels rather than one
# because absolute counts alone cannot rank degradation across methods whose
# baselines are 8, 6 and 4 of 10 — see the header note on retention.
ang_cnt <- three_method_tbl("angle", "count")
ang_nlg <- three_method_tbl("angle", "neglog")

ang_panels <- bind_rows(
  ang_cnt %>% transmute(level, method, med, lo, hi,
                        panel = sprintf("(i) pairs resolved (of %d)", n_pairs)),
  ang_nlg %>% transmute(level, method, med, lo, hi,
                        panel = "(ii) median -log10(raw p)"),
  ang_nlg %>% transmute(level, method, med = ret_med, lo = ret_lo, hi = ret_hi,
                        panel = "(iii) retention vs own baseline")) %>%
  mutate(panel = factor(panel, levels = c(sprintf("(i) pairs resolved (of %d)", n_pairs),
                                          "(ii) median -log10(raw p)",
                                          "(iii) retention vs own baseline")))

p1b <- ggplot(ang_panels, aes(level, med, color = method, fill = method)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.18, colour = NA) +
  geom_line(linewidth = 0.8) + geom_point(size = 1.7) +
  facet_wrap(~ panel, nrow = 1, scales = "free_y") +
  scale_color_manual(values = METHOD_COLORS, name = NULL) +
  scale_fill_manual(values = METHOD_COLORS, name = NULL) +
  scale_x_continuous(breaks = c(0, sort(unique(ang_cnt$level[ang_cnt$level > 0])))) +
  labs(x = "Angular jitter (s.d., degrees)", y = NULL,
       caption = paste0("Tests as in the main analysis: SP-SPHARM and fabric by PERMANOVA ",
                        "(adonis2), SPI by Kruskal-Wallis with Dunn post-hoc; Holm throughout.\n",
                        "Panel (iii) divides each method by its own unperturbed baseline ",
                        "(SP-SPHARM 8/10, fabric 6/10, SPI 4/10) and is the only panel on ",
                        "which degradation rates are comparable across methods.")) +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom",
        strip.text = element_text(face = "bold", size = 8),
        legend.key.size = grid::unit(0.34, "cm"),
        plot.caption = element_text(size = 6.5, colour = "grey40", hjust = 0))
ggsave(file.path(FIG_DIR, "fig_S_perturbation_angle_three_methods.png"), p1b,
       width = 10, height = 4.4, dpi = 300)
cat("Wrote figures/fig_S_perturbation_angle_three_methods.png\n")

# Three-method angle table, with retention on both quantities.
angle_three_csv <- ang_cnt %>%
  transmute(perturbation = "angle", level, method,
            n_sig_05_med = med, n_sig_05_lo = lo, n_sig_05_hi = hi,
            n_sig_05_baseline = base, n_sig_05_retention = ret_med) %>%
  left_join(ang_nlg %>%
              transmute(level, method,
                        med_neglog10_med = med, med_neglog10_lo = lo,
                        med_neglog10_hi = hi, med_neglog10_baseline = base,
                        med_neglog10_retention = ret_med),
            by = c("level", "method")) %>%
  mutate(test = dplyr::case_when(
    method == "SPI" ~ "Kruskal-Wallis + Dunn (Holm)",
    TRUE            ~ "PERMANOVA adonis2 + Holm")) %>%
  arrange(method, level)
write_csv(angle_three_csv, file.path(OUT_DIR, "annotation_perturbation_angle_three_methods.csv"))
cat("Wrote annotation_perturbation_angle_three_methods.csv\n")

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
# E. Timing
# =============================================================================
cat("\n", strrep("=", 70), "\n", sep = "")
cat("E. COST\n")
cat(strrep("=", 70), "\n", sep = "")
n_rep_total <- nrow(rep_df)
cat(sprintf("  replicates evaluated : %d (%d conditions)\n", n_rep_total, nrow(conds)))
cat(sprintf("  R-side elapsed       : %.1f min (%.2f s per replicate at N_PERM = %d)\n",
            elapsed_min, elapsed_min * 60 / n_rep_total, N_PERM))
cat("  Python chain re-run  : 0.30 s per 58-specimen pass (measured; --timing mode)\n")

# =============================================================================
# (1) Angle three-method retention table and degradation ranking
# =============================================================================
cat("\n", strrep("=", 70), "\n", sep = "")
cat("(1) P2 ANGLE: THREE-METHOD RETENTION AND DEGRADATION RANKING\n")
cat(strrep("=", 70), "\n", sep = "")
cat(sprintf("  Tests: SP-SPHARM / fabric = PERMANOVA (adonis2); SPI = Kruskal-Wallis + Dunn.\n"))
cat(sprintf("  Baselines: SP-SPHARM %d/%d, fabric %d/%d, SPI %d/%d resolved;\n",
            b$n_sig_05, n_pairs, bf$n_sig_05, n_pairs, bs$n_sig_05, n_pairs))
cat(sprintf("             median -log10(p) %.3f / %.3f / %.3f respectively.\n",
            b$med_neglog10, bf$med_neglog10, bs$med_neglog10))
cat("  (res = pairs resolved; medNL = median -log10 raw p; ret_* = value / own baseline)\n")
print(angle_three_csv %>%
        transmute(method, sigma = level,
                  res = sprintf("%.1f [%.0f,%.0f]", n_sig_05_med, n_sig_05_lo, n_sig_05_hi),
                  ret_n = sprintf("%.2f", n_sig_05_retention),
                  medNL = sprintf("%.2f [%.2f,%.2f]", med_neglog10_med,
                                  med_neglog10_lo, med_neglog10_hi),
                  ret_c = sprintf("%.2f", med_neglog10_retention)) %>%
        as.data.frame(), row.names = FALSE)

SIG_MAX <- max(angle_three_csv$level)
ret_at_max <- angle_three_csv %>% filter(level == SIG_MAX) %>%
  select(method, retention_cont = med_neglog10_retention,
         retention_count = n_sig_05_retention)
# Ranking is adjudicated on the CONTINUOUS retention: free of both the differing
# baselines and the integer-threshold steps.
ord_obs <- ret_at_max %>% arrange(retention_cont) %>% pull(method) %>% as.character()
ord_pred <- c("SPI", "Fabric (E, I)", "SP-SPHARM")   # fastest -> slowest decay
pred_ok  <- identical(ord_obs, ord_pred)

cat(sprintf("\n  Retention at sigma = %.0f deg (continuous statistic):\n", SIG_MAX))
for (i in seq_len(nrow(ret_at_max)))
  cat(sprintf("    %-14s %.3f   (count-based %.3f)\n",
              ret_at_max$method[i], ret_at_max$retention_cont[i],
              ret_at_max$retention_count[i]))
cat(sprintf("\n  observed decay order (fastest -> slowest): %s\n",
            paste(ord_obs, collapse = " > ")))
cat(sprintf("  predicted decay order                    : %s\n",
            paste(ord_pred, collapse = " > ")))

# First sigma at which each method falls below its own baseline (continuous stat).
first_below <- angle_three_csv %>% group_by(method) %>%
  summarise(first_sigma = {
    idx <- which(med_neglog10_retention < 1)
    if (length(idx) == 0) NA_real_ else level[min(idx)]
  }, .groups = "drop")

# =============================================================================
# (2) Prediction verdicts
# =============================================================================
cat("\n", strrep("=", 70), "\n", sep = "")
cat("(2) PREDICTION VERDICTS (declared before the angle contrast was run)\n")
cat(strrep("=", 70), "\n", sep = "")
cat(sprintf("  Mechanism: SP-SPHARM truncates at l = %d, i.e. ~%.0f deg angular resolution;\n",
            LMAX_KEEP, SH_RES_DEG))
cat(sprintf("  jitter below that scale should be filtered out by the truncation itself.\n"))
cat(sprintf("\n  (a) retention at sigma = %.0f deg, ordering vs prediction -> %s\n",
            SIG_MAX, ifelse(pred_ok, "MATCHES PREDICTION", "DOES NOT MATCH")))
cat(sprintf("      %s\n", paste(sprintf("%s %.3f", ret_at_max$method,
                                        ret_at_max$retention_cont), collapse = " | ")))
cat("\n  (b) first sigma at which each method drops below its own baseline:\n")
for (i in seq_len(nrow(first_below)))
  cat(sprintf("      %-14s %s\n", first_below$method[i],
              ifelse(is.na(first_below$first_sigma[i]), "never in range",
                     sprintf("%.0f deg", first_below$first_sigma[i]))))
if (!pred_ok) {
  cat("\n  (c) MECHANISM EXPLANATION DOES NOT HOLD — needs re-examination.\n")
  cat(sprintf("      Observed ordering is %s, not the predicted %s.\n",
              paste(ord_obs, collapse = " > "), paste(ord_pred, collapse = " > ")))
  cat("      The low-pass-truncation account of SP-SPHARM's tolerance is not supported\n")
  cat("      by the observed decay rates and should not be asserted in the manuscript.\n")
} else {
  cat("\n  (c) ordering matches; the low-pass-truncation account is consistent with the data.\n")
  cat("      (Consistent with, not proof of — no alternative mechanism was tested.)\n")
}

# =============================================================================
# (3) SP-SPHARM R2 trend under angular jitter vs the bandwidth sweep
# =============================================================================
cat("\n", strrep("=", 70), "\n", sep = "")
cat("(3) SP-SPHARM R2 TREND UNDER ANGULAR JITTER\n")
cat(strrep("=", 70), "\n", sep = "")
ang_sum <- summary_df %>% filter(perturbation == "angle") %>% arrange(level)
r2_rises <- ang_sum$R2_med[nrow(ang_sum)] > b$R2
cat(sprintf("  baseline R2 = %.5f; by sigma:\n", b$R2))
print(ang_sum %>% transmute(sigma = level, R2 = round(R2_med, 5),
                            resolved = n_sig_05_med) %>% as.data.frame(),
      row.names = FALSE)
cat(sprintf("\n  R2 %s with jitter (%.5f -> %.5f) while resolved pairs go %d -> %.0f.\n",
            ifelse(r2_rises, "RISES slightly", "does not rise"),
            b$R2, ang_sum$R2_med[nrow(ang_sum)], b$n_sig_05,
            ang_sum$n_sig_05_med[nrow(ang_sum)]))
if (r2_rises) {
  cat(sprintf("  Same direction as the bandwidth sweep (Table S1): h %.2f -> %.2f gives\n",
              BW_REF$h_lo, BW_REF$h_hi))
  cat(sprintf("  R2 %.5f -> %.5f. Small isotropic jitter acts like a wider KDE kernel.\n",
              BW_REF$R2_lo, BW_REF$R2_hi))
  cat("  This is NOT an improvement in discriminative power: R2 measures between-group\n")
  cat("  spread against total spread, and smoothing shrinks within-group spread too.\n")
  cat(sprintf("  The pairwise resolution falls over the same range (%d -> %.0f pairs), which\n",
              b$n_sig_05, ang_sum$n_sig_05_med[nrow(ang_sum)]))
  cat("  is the quantity that actually matters.\n")
  cat(sprintf("  One asymmetry worth noting: over h %.2f -> %.2f the resolved count HOLDS at\n",
              BW_REF$h_lo, BW_REF$h_hi))
  cat(sprintf("  %d/%d, whereas jitter costs pairs. Wider smoothing and annotation noise raise\n",
              BW_REF$nsig_hi, n_pairs))
  cat("  R2 alike, but only the noise destroys resolution.\n")
}

cat("\nDone.\n\n")
print(sessionInfo())
