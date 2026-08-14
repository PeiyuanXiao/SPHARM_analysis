# truncation_sensitivity.R

suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(readxl)
  library(vegan)
  library(ade4)
  library(compositions)
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
SEED   <- 42
N_PERM <- 9999
ALPHA_1 <- 0.05
ALPHA_2 <- 0.01

SP_GRID_REQ <- 4:8
M_GRID_REQ  <- 6:10

MAX_DEGREE_COL <- 20

OUT_DIR <- here("analysis/robustness/truncation_sensitivity")
FIG_DIR <- file.path(OUT_DIR, "figures")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

SERIES_COLORS <- c("Resolved pairs (Holm p < 0.05)" = "#4A6E8A",
                   "Resolved pairs (Holm p < 0.01)" = "#802520",
                   "Median -log10(raw p)"           = "#BA8530")
ASM_LABEL <- c(EXP = "Experimentally knapped cores", SDG = "Sandinggai cores")

# ---- Anchors -----------------------------------------------------------------
REF <- list(
  EXP = list(n = 58L, n_groups = 5L,
             M_R2  = 0.12334, M_F  = 1.86418, M_nsig05  = 0L,
             SP_R2 = 0.30174, SP_F = 5.72568, SP_nsig05 = 8L),
  SDG = list(n = 50L, n_groups = 6L,
             M_R2  = 0.18770, M_F  = 2.03346,
             SP_R2 = 0.16785, SP_F = 1.77504)
)
REF_TOL <- 1e-3

REF_DECOUP <- list(
  EXP = list(n = 58L, mantel_r = -0.06126, RV = 0.10095),
  SDG = list(n = 51L, mantel_r =  0.01189, RV = 0.09066)
)

# =============================================================================
# Helper import — verbatim from the committed companion scripts, WITHOUT running
# them
# =============================================================================
EXP_SCRIPT <- here("analysis/robustness/joint_mfa_discrimination",
                   "joint_mfa_discrimination_stats.R")
SDG_SCRIPT <- here("analysis/robustness/joint_mfa_discrimination",
                   "joint_mfa_discrimination_SDG.R")

import_from <- function(path, funs = character(0), consts = character(0)) {
  if (!file.exists(path))
    stop("Companion script not found at ", path,
         " — cannot import; stopping rather than duplicating its code here.")
  got <- character(0)
  for (e in parse(path)) {
    if (!is.call(e)) next
    if (!as.character(e[[1]])[1] %in% c("<-", "=")) next
    if (!is.name(e[[2]])) next
    nm  <- as.character(e[[2]])
    rhs <- e[[3]]
    is_fun <- is.call(rhs) && identical(as.character(rhs[[1]])[1], "function")
    if (nm %in% funs && is_fun) {
      eval(e, envir = globalenv()); got <- c(got, nm)
    } else if (nm %in% consts && !is_fun) {
      eval(e, envir = globalenv()); got <- c(got, nm)
    }
  }
  missing <- setdiff(c(funs, consts), got)
  if (length(missing))
    stop("Could not import from ", basename(path), ": ",
         paste(missing, collapse = ", "),
         " — the companion script's structure changed; stopping rather than ",
         "reimplementing.")
  got
}

EXP_FUNS <- c("replace_zeros", "make_ilr", "ilr_parts",
              "permanova_global", "permanova_pairwise", "permdisp", "pair_p",
              "filter_spharm", "split_by_group", "safe_filter_groups",
              "build_blocks")
EXP_CONSTS <- c("POWER_COLS_DIR", "POWER_COLS_MORPH",
                "EXCLUDE_TYPES", "LEVALLOIS_MERGE", "TYPOLOGY_ORDER",
                "EXCLUDE_CORE_TYPES", "FOCUS_PAIRS")
got_exp <- import_from(EXP_SCRIPT, EXP_FUNS, EXP_CONSTS)

SDG_FUNS   <- c("build_blocks_sdg")
SDG_CONSTS <- c("RESTRICT_LAYERS")
got_sdg <- import_from(SDG_SCRIPT, SDG_FUNS, SDG_CONSTS)

cat(sprintf("Imported %d objects from %s:\n  %s\n",
            length(got_exp), basename(EXP_SCRIPT), paste(sort(got_exp), collapse = ", ")))
cat(sprintf("Imported %d objects from %s:\n  %s\n",
            length(got_sdg), basename(SDG_SCRIPT), paste(sort(got_sdg), collapse = ", ")))

ASSEMBLAGE <- "EXP"

ANCHOR_SP <- length(POWER_COLS_DIR)
ANCHOR_M  <- length(POWER_COLS_MORPH)
stopifnot(identical(POWER_COLS_DIR,   paste0("power_l", seq_len(ANCHOR_SP))),
          identical(POWER_COLS_MORPH, paste0("power_l", seq_len(ANCHOR_M))))
cat(sprintf("\nAnchor truncations from the imported constants: SP l = 1-%d | M l = 1-%d\n",
            ANCHOR_SP, ANCHOR_M))

P_FLOOR <- 1 / (N_PERM + 1)

# ---- local helpers with no equivalent in the imported set --------------------

extract_subdist <- function(D_full, ids) as.dist(as.matrix(D_full)[ids, ids])

compute_degree_stats <- function(df, cols) {
  mat        <- df %>% select(all_of(cols)) %>% as.matrix()
  col_means  <- colMeans(mat, na.rm = TRUE)
  col_sds    <- apply(mat, 2, sd, na.rm = TRUE)
  col_cvs    <- col_sds / col_means * 100
  row_sums   <- rowSums(mat, na.rm = TRUE)
  cumul_pct  <- cumsum(col_means) / mean(row_sums) * 100
  tibble(degree = seq_along(cols), mean = col_means, sd = col_sds,
         cv_pct = col_cvs, cumul_pct = cumul_pct)
}

power_cols_upto <- function(K) paste0("power_l", seq_len(K))

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
# (i) AVAILABLE-DEGREE AUDIT — clamp the sweep to what actually exists
# =============================================================================
cat("\n", strrep("=", 70), "\n", sep = "")
cat("AVAILABLE-DEGREE AUDIT\n")
cat(strrep("=", 70), "\n", sep = "")

available_lmax <- function(df) {
  have <- paste0("power_l", 1:MAX_DEGREE_COL) %in% colnames(df)
  if (!have[1]) return(0L)
  as.integer(max(which(cumsum(!have) == 0)))
}
lmax_dir   <- available_lmax(SPHARM_direction)
lmax_morph <- available_lmax(SPHARM_morphology)
cat(sprintf("  SPHARM_direction.csv  contiguous power columns: power_l1..power_l%d\n", lmax_dir))
cat(sprintf("  SPHARM_morphology.csv contiguous power columns: power_l1..power_l%d\n", lmax_morph))

SP_GRID <- SP_GRID_REQ[SP_GRID_REQ <= lmax_dir]
M_GRID  <- M_GRID_REQ[M_GRID_REQ  <= lmax_morph]
if (!setequal(SP_GRID, SP_GRID_REQ))
  cat(sprintf("  ! SP sweep contracted to the available range: requested 1-%d..1-%d, using 1-%d..1-%d\n",
              min(SP_GRID_REQ), max(SP_GRID_REQ), min(SP_GRID), max(SP_GRID)))
if (!setequal(M_GRID, M_GRID_REQ))
  cat(sprintf("  ! M sweep contracted to the available range: requested 1-%d..1-%d, using 1-%d..1-%d\n",
              min(M_GRID_REQ), max(M_GRID_REQ), min(M_GRID), max(M_GRID)))
if (!ANCHOR_SP %in% SP_GRID || !ANCHOR_M %in% M_GRID)
  stop("The main analysis setting is outside the available sweep range — ",
       "the sweep would have no anchor. Stopping.")
cat(sprintf("  Scan A (SP): %s   [anchor 1-%d]\n",
            paste(sprintf("1-%d", SP_GRID), collapse = ", "), ANCHOR_SP))
cat(sprintf("  Scan B (M) : %s   [anchor 1-%d]\n",
            paste(sprintf("1-%d", M_GRID), collapse = ", "), ANCHOR_M))

# ---- ZERO-ENTRY AUDIT --------------------------------------------------------
zero_audit <- function(df, K, label) {
  cols <- power_cols_upto(K)
  X <- as.matrix(df[, cols, drop = FALSE])
  nz <- sum(X == 0, na.rm = TRUE)
  cat(sprintf("  %-34s zeros in power_l1..l%-2d : %d%s\n", label, K, nz,
              ifelse(nz == 0, "  (replace_zeros is inert; ILR closure-invariant)",
                     "  <-- replace_zeros ACTIVE, row closure is NOT a no-op")))
  nz
}
cat("\n  Zero-entry audit over the swept column ranges:\n")
nz_sp <- zero_audit(SPHARM_direction,  max(SP_GRID), "SPHARM_direction (SP sweep)")
nz_m  <- zero_audit(SPHARM_morphology, max(M_GRID),  "SPHARM_morphology (M sweep)")
CLOSURE_INERT <- (nz_sp == 0L && nz_m == 0L)

closure_dev <- function(df, K) {
  X  <- as.matrix(df[, power_cols_upto(K), drop = FALSE])
  Zr <- make_ilr(X)
  Zc <- make_ilr(X / rowSums(X))
  max(abs(Zr - Zc))
}
dev_sp <- closure_dev(SPHARM_direction,  max(SP_GRID))
dev_m  <- closure_dev(SPHARM_morphology, max(M_GRID))
cat(sprintf("\n  Closure equivalence (max |ILR(raw) - ILR(row-closed)|): SP %.3e | M %.3e\n",
            dev_sp, dev_m))
cat(sprintf("  -> the two formulations are %s; this script uses the pipeline's make_ilr()\n",
            ifelse(max(dev_sp, dev_m) < 1e-12, "numerically interchangeable",
                   "NOT interchangeable <-- investigate before reading any table")))

# =============================================================================
# (ii) SAMPLES — fixed once; truncation cannot change WHICH specimens enter
# =============================================================================
cat("\n", strrep("=", 70), "\n", sep = "")
cat("SAMPLE DEFINITIONS (imported verbatim, evaluated once)\n")
cat(strrep("=", 70), "\n", sep = "")

blocks_exp <- build_blocks("EXP")
blocks_sdg <- build_blocks_sdg()

SAMPLES <- list(
  EXP = list(ids = blocks_exp$ids, group = blocks_exp$group),
  SDG = list(ids = blocks_sdg$ids, group = blocks_sdg$group)
)
for (a in names(SAMPLES)) {
  s <- SAMPLES[[a]]
  s$n       <- length(s$ids)
  s$n_grp   <- nlevels(droplevels(s$group))
  s$n_pairs <- choose(s$n_grp, 2)
  SAMPLES[[a]] <- s
  cat(sprintf("  %s: n = %d | groups = %d | pairs = %d\n",
              a, s$n, s$n_grp, s$n_pairs))
  print(table(s$group))
}
if (SAMPLES$EXP$n != REF$EXP$n || SAMPLES$SDG$n != REF$SDG$n)
  stop(sprintf("Sample sizes (EXP %d, SDG %d) do not match the published recipe (%d, %d). Stopping.",
               SAMPLES$EXP$n, SAMPLES$SDG$n, REF$EXP$n, REF$SDG$n))

power_of <- function(asm, block, K) {
  df  <- if (identical(block, "SP")) SPHARM_direction else SPHARM_morphology
  idx <- match(SAMPLES[[asm]]$ids, df$ID)
  if (anyNA(idx))
    stop("Specimen IDs missing from the source CSV for ", asm, "/", block, ".")
  df[idx, power_cols_upto(K), drop = FALSE]
}

chk_align <- function(asm, block, K, ref_df) {
  got <- as.matrix(power_of(asm, block, K))
  ref <- as.matrix(ref_df)
  dimnames(got) <- dimnames(ref) <- NULL
  isTRUE(all.equal(got, ref, tolerance = 0))
}
align_ok <- c(
  chk_align("EXP", "SP", ANCHOR_SP, blocks_exp$scar),
  chk_align("EXP", "M",  ANCHOR_M,  blocks_exp$morph),
  chk_align("SDG", "SP", ANCHOR_SP, blocks_sdg$scar),
  chk_align("SDG", "M",  ANCHOR_M,  blocks_sdg$morph))
if (!all(align_ok))
  stop("Re-extracted anchor power matrices differ from build_blocks() output — ",
       "row alignment is not what this script assumes. Stopping.")
cat("\n  Re-extraction check: anchor power matrices identical to build_blocks() output (all 4).\n")

# =============================================================================
# (iii) PER-DEGREE DIAGNOSTICS (cumulative power + across-specimen CV)
# =============================================================================
cat("\n", strrep("=", 70), "\n", sep = "")
cat("PER-DEGREE DIAGNOSTICS (cumulative power / across-specimen CV)\n")
cat(strrep("=", 70), "\n", sep = "")

COLS_ALL <- paste0("power_l", 1:min(lmax_dir, lmax_morph, MAX_DEGREE_COL))
deg_sets <- list(
  EXP = list(SP = SPHARM_direction  %>% filter(str_starts(ID, "EXP")),
             M  = SPHARM_morphology %>% filter(str_starts(ID, "EXP"))),
  SDG = list(SP = SPHARM_direction  %>% filter(str_starts(ID, "SDG"), !str_starts(ID, "IM_")),
             M  = SPHARM_morphology %>% filter(str_starts(ID, "SDG"), !str_starts(ID, "IM_")))
)
DEG <- map(deg_sets, ~ map(.x, compute_degree_stats, cols = COLS_ALL))
for (a in names(deg_sets))
  cat(sprintf("  %s degree-selection sets: SP n = %d, M n = %d\n",
              a, nrow(deg_sets[[a]]$SP), nrow(deg_sets[[a]]$M)))

deg_ref_file <- function(block, asm)
  here(sprintf("analysis/data/derived_data/DegreeSelection_stats_%s_%s.csv",
               ifelse(block == "SP", "direction", "morphology"), asm))
deg_dev <- map_dfr(names(DEG), function(a) map_dfr(names(DEG[[a]]), function(b) {
  f <- deg_ref_file(b, a)
  if (!file.exists(f)) return(tibble(set = paste(a, b), d_cv = NA_real_, d_cum = NA_real_))
  ref <- read_csv(f, show_col_types = FALSE)
  j   <- DEG[[a]][[b]] %>% inner_join(ref, by = "degree", suffix = c("", "_ref"))
  tibble(set = paste(a, b),
         d_cv  = max(abs(round(j$cv_pct, 2)    - j$cv_pct_ref)),
         d_cum = max(abs(round(j$cumul_pct, 3) - j$cumul_pct_ref)))
}))
cat("\n  Cross-check vs committed DegreeSelection_stats_*.csv (max abs deviation):\n")
print(as.data.frame(deg_dev %>% mutate(across(where(is.numeric), ~ signif(.x, 3)))),
      row.names = FALSE)
if (any(deg_dev$d_cv > 0.011, na.rm = TRUE) || any(deg_dev$d_cum > 0.0011, na.rm = TRUE))
  cat("  ! Diagnostics do not reproduce the committed degree-selection table — investigate.\n")

cum_at   <- function(a, b, K) DEG[[a]][[b]]$cumul_pct[DEG[[a]][[b]]$degree == K]
cv_at    <- function(a, b, K) DEG[[a]][[b]]$cv_pct[DEG[[a]][[b]]$degree == K]
cv_cross <- function(a, b) {
  i <- which(DEG[[a]][[b]]$cv_pct > 100)
  if (!length(i)) NA_integer_ else as.integer(DEG[[a]][[b]]$degree[i[1]])
}
cv_max   <- function(a, b) max(DEG[[a]][[b]]$cv_pct, na.rm = TRUE)
cv_argmax<- function(a, b) as.integer(DEG[[a]][[b]]$degree[which.max(DEG[[a]][[b]]$cv_pct)])

# =============================================================================
# (iv) PER-SETTING EVALUATION (memoised: each block/degree computed once)
# =============================================================================
MEMO_BLOCK  <- new.env(parent = emptyenv())
MEMO_DECOUP <- new.env(parent = emptyenv())
memo_get <- function(env, key)
  if (exists(key, envir = env, inherits = FALSE)) get(key, envir = env) else NULL

eval_block <- function(asm, block, K) {
  key <- sprintf("%s|%s|%d", asm, block, K)
  hit <- memo_get(MEMO_BLOCK, key); if (!is.null(hit)) return(hit)
  s  <- SAMPLES[[asm]]
  P  <- power_of(asm, block, K)
  Z  <- make_ilr(P)
  cat(sprintf("    [%s %s l=1-%d] %d parts -> %d coords; PERMANOVA (%d pairs) + PERMDISP ...\n",
              asm, block, K, length(ilr_parts(P)), ncol(Z), s$n_pairs))
  gl <- permanova_global(Z, s$group)
  pw <- permanova_pairwise(Z, s$group)
  pd <- permdisp(Z, s$group)
  res <- list(
    n_parts   = length(ilr_parts(P)),
    n_coords  = ncol(Z),
    R2 = gl$R2, F = gl$F, p = gl$p,
    nsig05 = sum(pw$p_holm < ALPHA_1),
    nsig01 = sum(pw$p_holm < ALPHA_2),
    med_neglog10 = stats::median(-log10(pw$p)),
    sig_pairs = paste(sort(pw$comparison[pw$p_holm < ALPHA_1]), collapse = "; "),
    permdisp_F = pd$F, permdisp_p = pd$p,
    pw = pw)
  assign(key, res, envir = MEMO_BLOCK)
  res
}

# ---- decoupling: Mantel + RV between the two blocks --------------------------
rv_and_mantel <- function(Zm, Zs, ids) {
  Dm <- extract_subdist(stats::dist(Zm), ids)
  Ds <- extract_subdist(stats::dist(Zs), ids)
  set.seed(SEED)
  mt <- vegan::mantel(Dm, Ds, method = "spearman", permutations = N_PERM)
  Am <- as.data.frame(Zm[ids, , drop = FALSE])
  As <- as.data.frame(Zs[ids, , drop = FALSE])
  colnames(Am) <- paste0("M_ilr", seq_len(ncol(Am)))
  colnames(As) <- paste0("S_ilr", seq_len(ncol(As)))
  dm <- dudi.pca(Am, center = TRUE, scale = TRUE, scannf = FALSE, nf = ncol(Am))
  ds <- dudi.pca(As, center = TRUE, scale = TRUE, scannf = FALSE, nf = ncol(As))
  co <- coinertia(dm, ds, scannf = FALSE, nf = 2)
  set.seed(SEED)
  rt <- randtest(co, nrepet = N_PERM)
  list(mantel_r = mt$statistic, mantel_p = mt$signif,
       RV = co$RV, RV_p = rt$pvalue, n = length(ids))
}

decoup_exp <- function(K_M, K_SP) {
  cM <- power_cols_upto(K_M); cS <- power_cols_upto(K_SP)
  d_all <- split_by_group(filter_spharm(SPHARM_direction,  cS, metric_data))$exp_im
  m_all <- split_by_group(filter_spharm(SPHARM_morphology, cM, metric_data))$exp_im
  common <- intersect(m_all$ID, d_all$ID)
  d_all <- d_all %>% filter(ID %in% common) %>% arrange(ID)
  m_all <- m_all %>% filter(ID %in% common) %>% arrange(ID)
  Zm <- make_ilr(m_all %>% select(all_of(cM))); rownames(Zm) <- m_all$ID
  Zs <- make_ilr(d_all %>% select(all_of(cS))); rownames(Zs) <- d_all$ID
  rv_and_mantel(Zm, Zs, SAMPLES$EXP$ids)
}

decoup_sdg <- function(K_M, K_SP) {
  cM <- power_cols_upto(K_M); cS <- power_cols_upto(K_SP)
  d_all <- filter_spharm(SPHARM_direction,  cS, metric_data)
  m_all <- filter_spharm(SPHARM_morphology, cM, metric_data)
  common <- intersect(m_all$ID, d_all$ID)
  d_all <- d_all %>% filter(ID %in% common) %>% arrange(ID)
  m_all <- m_all %>% filter(ID %in% common) %>% arrange(ID)
  Zm <- make_ilr(m_all %>% select(all_of(cM))); rownames(Zm) <- m_all$ID
  Zs <- make_ilr(d_all %>% select(all_of(cS))); rownames(Zs) <- d_all$ID
  arch <- tibble(ID = rownames(Zm)) %>%
    filter(!str_starts(ID, "IM_"), !str_starts(ID, "EXP")) %>%
    left_join(core_meta %>% select(ID, core_type), by = "ID") %>%
    filter(!core_type %in% EXCLUDE_CORE_TYPES | is.na(core_type))
  rv_and_mantel(Zm, Zs, arch$ID)
}

eval_decoup <- function(asm, K_M, K_SP) {
  key <- sprintf("%s|%d|%d", asm, K_M, K_SP)
  hit <- memo_get(MEMO_DECOUP, key); if (!is.null(hit)) return(hit)
  cat(sprintf("    [%s decoupling M l=1-%d x SP l=1-%d] Mantel + RV ...\n", asm, K_M, K_SP))
  res <- if (identical(asm, "EXP")) decoup_exp(K_M, K_SP) else decoup_sdg(K_M, K_SP)
  assign(key, res, envir = MEMO_DECOUP)
  res
}

# ---- one row of the sweep ----------------------------------------------------
sweep_row <- function(scan, asm, K_SP, K_M) {
  s   <- SAMPLES[[asm]]
  bSP <- eval_block(asm, "SP", K_SP)
  bM  <- eval_block(asm, "M",  K_M)
  dc  <- eval_decoup(asm, K_M, K_SP)
  varied  <- if (identical(scan, "SP")) "SP-SPHARM" else "M-SPHARM"
  vb      <- if (identical(scan, "SP")) bSP else bM
  varied_l<- if (identical(scan, "SP")) K_SP else K_M
  fp <- if (identical(asm, "EXP"))
    map_dbl(FOCUS_PAIRS, ~ pair_p(vb$pw, .x, col = "p_holm")) else rep(NA_real_, 2)

  tibble(
    scan = scan, assemblage = asm,
    varied_block = varied, varied_lmax = varied_l,
    setting = sprintf("%s l = 1-%d", varied, varied_l),
    lmax_SP = K_SP, lmax_M = K_M,
    is_anchor = (K_SP == ANCHOR_SP && K_M == ANCHOR_M),
    n_coords = vb$n_coords,
    n_permanova = s$n, n_groups = s$n_grp, n_pairs = s$n_pairs,

    SP_n_parts = bSP$n_parts, SP_n_coords = bSP$n_coords,
    SP_R2 = bSP$R2, SP_pseudo_F = bSP$F, SP_p = bSP$p,
    SP_nsig05 = bSP$nsig05, SP_nsig01 = bSP$nsig01,
    SP_med_neglog10_p = bSP$med_neglog10,
    SP_permdisp_F = bSP$permdisp_F, SP_permdisp_p = bSP$permdisp_p,
    SP_sig_pairs = bSP$sig_pairs,

    M_n_parts = bM$n_parts, M_n_coords = bM$n_coords,
    M_R2 = bM$R2, M_pseudo_F = bM$F, M_p = bM$p,
    M_nsig05 = bM$nsig05, M_nsig01 = bM$nsig01,
    M_med_neglog10_p = bM$med_neglog10,
    M_permdisp_F = bM$permdisp_F, M_permdisp_p = bM$permdisp_p,
    M_sig_pairs = bM$sig_pairs,

    mantel_r = dc$mantel_r, mantel_p = dc$mantel_p,
    RV = dc$RV, RV_p = dc$RV_p, n_decoupling = dc$n,

    SP_cumpower_at_lmax = cum_at(asm, "SP", K_SP),
    SP_cv_at_lmax       = cv_at(asm,  "SP", K_SP),
    SP_cv_cross_l       = cv_cross(asm, "SP"),
    M_cumpower_at_lmax  = cum_at(asm, "M", K_M),
    M_cv_at_lmax        = cv_at(asm,  "M", K_M),
    M_cv_cross_l        = cv_cross(asm, "M"),
    M_cv_max_all        = cv_max(asm, "M"),

    focus_p_Discoid_Bidirectional = fp[1],
    focus_p_Discoid_Multiplatform = fp[2])
}

# =============================================================================
# (v) RUN BOTH SCANS
# =============================================================================
cat("\n", strrep("=", 70), "\n", sep = "")
cat(sprintf("RUNNING THE SWEEP  (seed = %d, permutations = %d, p floor = %.5f)\n",
            SEED, N_PERM, P_FLOOR))
cat(strrep("=", 70), "\n", sep = "")

scanA <- map_dfr(c("EXP", "SDG"), function(a) {
  cat(sprintf("\n  -- Scan A | %s | vary SP, hold M at l = 1-%d --\n", a, ANCHOR_M))
  map_dfr(SP_GRID, ~ sweep_row("SP", a, K_SP = .x, K_M = ANCHOR_M))
})
scanB <- map_dfr(c("EXP", "SDG"), function(a) {
  cat(sprintf("\n  -- Scan B | %s | vary M, hold SP at l = 1-%d --\n", a, ANCHOR_SP))
  map_dfr(M_GRID, ~ sweep_row("M", a, K_SP = ANCHOR_SP, K_M = .x))
})
summary_tbl <- bind_rows(scanA, scanB)

anchor_rows <- summary_tbl %>% filter(is_anchor)
consistency_ok <- anchor_rows %>%
  group_by(assemblage) %>%
  summarise(ok = n_distinct(round(SP_R2, 12)) == 1 &&
                 n_distinct(round(M_R2, 12))  == 1 &&
                 n_distinct(SP_nsig05) == 1 && n_distinct(M_nsig05) == 1 &&
                 n_distinct(round(RV, 12)) == 1, .groups = "drop")
cat(sprintf("\n  Anchor appears in both scans and matches: %s\n",
            ifelse(all(consistency_ok$ok), "YES",
                   paste("NO —", paste(consistency_ok$assemblage[!consistency_ok$ok],
                                       collapse = ", ")))))

# =============================================================================
# (vi) OUTPUT CSVs
# =============================================================================
write_csv(scanA,       file.path(OUT_DIR, "truncation_sensitivity_SP.csv"))
write_csv(scanB,       file.path(OUT_DIR, "truncation_sensitivity_M.csv"))
write_csv(summary_tbl, file.path(OUT_DIR, "truncation_sensitivity_summary.csv"))
cat(sprintf("\nWrote truncation_sensitivity_SP.csv (%d rows), truncation_sensitivity_M.csv (%d rows), truncation_sensitivity_summary.csv (%d rows)\n",
            nrow(scanA), nrow(scanB), nrow(summary_tbl)))

# =============================================================================
# (vii) FIGURES
# =============================================================================
make_trunc_fig <- function(df, anchor_l, xlab) {
  key <- df %>% transmute(
    assemblage = factor(assemblage, levels = c("EXP", "SDG")),
    l = varied_lmax, n_pairs,
    nsig05 = if_else(varied_block == "SP-SPHARM", SP_nsig05, M_nsig05),
    nsig01 = if_else(varied_block == "SP-SPHARM", SP_nsig01, M_nsig01),
    med    = if_else(varied_block == "SP-SPHARM", SP_med_neglog10_p, M_med_neglog10_p))

  series_of <- c(nsig05 = "Resolved pairs (Holm p < 0.05)",
                 nsig01 = "Resolved pairs (Holm p < 0.01)",
                 med    = "Median -log10(raw p)")
  long <- key %>%
    pivot_longer(c(nsig05, nsig01, med), names_to = "k", values_to = "value") %>%
    mutate(series = factor(unname(series_of[k]), levels = names(SERIES_COLORS)),
           metric = factor(if_else(k == "med", "Median -log10(raw p)",
                                   "Resolved pairs (of all pairs)"),
                           levels = c("Resolved pairs (of all pairs)",
                                      "Median -log10(raw p)")))

  ceil <- key %>% distinct(assemblage, n_pairs) %>%
    mutate(metric = factor("Resolved pairs (of all pairs)",
                           levels = levels(long$metric)))

  ggplot(long, aes(l, value, color = series, group = series)) +
    geom_hline(data = ceil, aes(yintercept = n_pairs), inherit.aes = FALSE,
               linetype = "dashed", color = "grey60", linewidth = 0.3) +
    geom_vline(xintercept = anchor_l, color = "red",
               linetype = "dashed", linewidth = 0.3) +
    geom_line(aes(linetype = series), linewidth = 0.6) +
    geom_point(aes(shape = series), size = 1.9, stroke = 0.7, fill = "white") +
    facet_grid(metric ~ assemblage, scales = "free_y",
               labeller = labeller(assemblage = ASM_LABEL)) +
    scale_color_manual(values = SERIES_COLORS, name = NULL) +
    scale_linetype_manual(values = c("solid", "22", "solid"), name = NULL) +
    scale_shape_manual(values = c(16, 21, 16), name = NULL) +
    scale_x_continuous(breaks = sort(unique(long$l)),
                       labels = sprintf("1-%d", sort(unique(long$l)))) +
    labs(x = xlab, y = NULL) +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank(),
          strip.text = element_text(face = "bold"),
          legend.position = "bottom")
}

ok_fig <- TRUE
tryCatch({
  p_sp <- make_trunc_fig(scanA, ANCHOR_SP,
                         "SP-SPHARM truncation (retained degrees l)")
  ggsave(file.path(FIG_DIR, "fig_S_truncation_SP.png"), p_sp,
         width = 9, height = 6, dpi = 300)
  p_m <- make_trunc_fig(scanB, ANCHOR_M,
                        "M-SPHARM truncation (retained degrees l)")
  ggsave(file.path(FIG_DIR, "fig_S_truncation_M.png"), p_m,
         width = 9, height = 6, dpi = 300)
  cat("Wrote 2 figures to", FIG_DIR, "\n")
}, error = function(e) {
  ok_fig <<- FALSE
  cat("Figure generation failed (non-critical):", conditionMessage(e), "\n")
})

# =============================================================================
# (0) ANCHOR CHECK
# =============================================================================
cat("\n", strrep("=", 70), "\n", sep = "")
cat("(0) ANCHOR CHECK — main analysis setting (SP l = 1-", ANCHOR_SP,
    ", M l = 1-", ANCHOR_M, ")\n", sep = "")
cat(strrep("=", 70), "\n", sep = "")

anchor_ok  <- TRUE
extra_flag <- FALSE
check <- function(name, got, expect, tol, deterministic = TRUE) {
  ok <- is.finite(got) && abs(got - expect) <= tol
  if (!ok) { if (deterministic) anchor_ok <<- FALSE else extra_flag <<- TRUE }
  cat(sprintf("  %-30s got %12.5f  expected %12.5f  %s\n", name, got, expect,
              ifelse(ok, "OK", ifelse(deterministic, "<-- MISMATCH", "<-- CHECK"))))
}

A <- summary_tbl %>% filter(is_anchor, scan == "SP")
for (a in c("EXP", "SDG")) {
  r <- A %>% filter(assemblage == a)
  R <- REF[[a]]
  cat(sprintf("\n  -- %s --\n", a))
  check("n specimens",        SAMPLES[[a]]$n,  R$n,     0)
  check("groups",             SAMPLES[[a]]$n_grp, R$n_groups, 0)
  check("M-SPHARM R2",        r$M_R2,          R$M_R2,  REF_TOL)
  check("M-SPHARM pseudo-F",  r$M_pseudo_F,    R$M_F,   REF_TOL)
  check("SP-SPHARM R2",       r$SP_R2,         R$SP_R2, REF_TOL)
  check("SP-SPHARM pseudo-F", r$SP_pseudo_F,   R$SP_F,  REF_TOL)
  if (!is.null(R$M_nsig05))
    check("M-SPHARM resolved pairs",  r$M_nsig05,  R$M_nsig05,  0)
  if (!is.null(R$SP_nsig05))
    check("SP-SPHARM resolved pairs", r$SP_nsig05, R$SP_nsig05, 0)
  cat(sprintf("  %-30s M p = %.4f | SP p = %.4f   (permutation; reported, not checked)\n",
              "global p", r$M_p, r$SP_p))
  D <- REF_DECOUP[[a]]
  check("decoupling n",       r$n_decoupling,  D$n,        0,      deterministic = FALSE)
  check("Mantel r",           r$mantel_r,      D$mantel_r, REF_TOL, deterministic = FALSE)
  check("RV",                 r$RV,            D$RV,       REF_TOL, deterministic = FALSE)
  cat(sprintf("  %-30s Mantel p = %.4f | RV p = %.4f   (permutation; reported, not checked)\n",
              "decoupling p", r$mantel_p, r$RV_p))
}

if (!anchor_ok) {
  cat("\n  => STOP: the anchor setting does NOT reproduce the committed values.\n")
  cat("     The sweep is not comparable to the main text; nothing below should be trusted.\n")
  cat("     CSVs and figures were still written for inspection.\n\n")
  print(sessionInfo())
  stop("Anchor check failed — see the table above.")
}
cat("\n  => Anchor reproduces the committed main-text PERMANOVA (R2 and pseudo-F within ",
    format(REF_TOL), ").\n", sep = "")
if (extra_flag)
  cat("  !  One or more decoupling anchors deviate — see '<-- CHECK' above.\n")
cat("  (Permutation p-values carry ~+/-0.005 Monte-Carlo jitter and are reported, not checked.)\n")

# =============================================================================
# Reporting helpers
# =============================================================================
fmt_p   <- function(p) ifelse(p <= P_FLOOR, sprintf("<=%.4f", P_FLOOR),
                              sprintf("%.4f", p))
sig_lab <- function(p) ifelse(p < ALPHA_1, "significant", "n.s.")

scan_table <- function(df, block) {
  v <- if (block == "SP") "SP" else "M"
  df %>% transmute(
    assemblage,
    truncation = sprintf("1-%d%s", varied_lmax, ifelse(is_anchor, "*", " ")),
    n_coords,
    R2       = round(.data[[paste0(v, "_R2")]], 5),
    pseudo_F = round(.data[[paste0(v, "_pseudo_F")]], 4),
    p        = .data[[paste0(v, "_p")]],
    resolved_05 = sprintf("%d/%d", .data[[paste0(v, "_nsig05")]], n_pairs),
    resolved_01 = sprintf("%d/%d", .data[[paste0(v, "_nsig01")]], n_pairs),
    med_neglog10 = round(.data[[paste0(v, "_med_neglog10_p")]], 3),
    mantel = sprintf("%+.3f (%.3f)", mantel_r, mantel_p),
    RV     = sprintf("%.3f (%.3f)", RV, RV_p),
    disp_p = round(.data[[paste0(v, "_permdisp_p")]], 3))
}

stable_window <- function(df, block) {
  v <- if (block == "SP") "SP" else "M"
  d <- df %>% arrange(varied_lmax)
  a <- which(d$is_anchor)
  same <- function(i) {
    d[[paste0(v, "_nsig05")]][i] == d[[paste0(v, "_nsig05")]][a] &&
    d[[paste0(v, "_nsig01")]][i] == d[[paste0(v, "_nsig01")]][a] &&
    (d[[paste0(v, "_p")]][i] < ALPHA_1) == (d[[paste0(v, "_p")]][a] < ALPHA_1) &&
    (d$mantel_p[i] < ALPHA_1) == (d$mantel_p[a] < ALPHA_1) &&
    (d$RV_p[i]     < ALPHA_1) == (d$RV_p[a]     < ALPHA_1)
  }
  lo <- a; while (lo > 1 && same(lo - 1)) lo <- lo - 1
  hi <- a; while (hi < nrow(d) && same(hi + 1)) hi <- hi + 1
  c(lo = d$varied_lmax[lo], hi = d$varied_lmax[hi])
}

compare_two <- function(df, block, l_a, l_b, asm) {
  v <- if (block == "SP") "SP" else "M"
  ra <- df %>% filter(assemblage == asm, varied_lmax == l_a)
  rb <- df %>% filter(assemblage == asm, varied_lmax == l_b)
  if (nrow(ra) != 1 || nrow(rb) != 1) return(NULL)
  g <- function(r, col) r[[paste0(v, "_", col)]]
  cat(sprintf("    %s : l = 1-%d  vs  l = 1-%d\n", asm, l_a, l_b))
  cat(sprintf("      ILR coordinates       %6d          %6d          (dimension differs -> R2 not comparable)\n",
              ra$n_coords, rb$n_coords))
  cat(sprintf("      R2 / pseudo-F         %6.4f / %5.3f  %6.4f / %5.3f   [context only, NOT a criterion]\n",
              g(ra, "R2"), g(ra, "pseudo_F"), g(rb, "R2"), g(rb, "pseudo_F")))
  cat(sprintf("      global p              %8s        %8s\n",
              fmt_p(g(ra, "p")), fmt_p(g(rb, "p"))))
  cat(sprintf("      resolved (Holm .05)   %4d/%-4d       %4d/%-4d       delta = %+d\n",
              g(ra, "nsig05"), ra$n_pairs, g(rb, "nsig05"), rb$n_pairs,
              g(rb, "nsig05") - g(ra, "nsig05")))
  cat(sprintf("      resolved (Holm .01)   %4d/%-4d       %4d/%-4d       delta = %+d\n",
              g(ra, "nsig01"), ra$n_pairs, g(rb, "nsig01"), rb$n_pairs,
              g(rb, "nsig01") - g(ra, "nsig01")))
  cat(sprintf("      median -log10(raw p)  %8.3f        %8.3f        delta = %+.3f\n",
              g(ra, "med_neglog10_p"), g(rb, "med_neglog10_p"),
              g(rb, "med_neglog10_p") - g(ra, "med_neglog10_p")))
  cat(sprintf("      Mantel r (p)          %+.3f (%.3f)   %+.3f (%.3f)   %s -> %s\n",
              ra$mantel_r, ra$mantel_p, rb$mantel_r, rb$mantel_p,
              sig_lab(ra$mantel_p), sig_lab(rb$mantel_p)))
  cat(sprintf("      RV (p)                %6.3f (%.3f)   %6.3f (%.3f)   %s -> %s\n",
              ra$RV, ra$RV_p, rb$RV, rb$RV_p,
              sig_lab(ra$RV_p), sig_lab(rb$RV_p)))
  list(d05 = g(rb, "nsig05") - g(ra, "nsig05"),
       d01 = g(rb, "nsig01") - g(ra, "nsig01"),
       dmed = g(rb, "med_neglog10_p") - g(ra, "med_neglog10_p"),
       flip_global = (g(ra, "p") < ALPHA_1) != (g(rb, "p") < ALPHA_1),
       flip_mantel = (ra$mantel_p < ALPHA_1) != (rb$mantel_p < ALPHA_1),
       flip_rv     = (ra$RV_p     < ALPHA_1) != (rb$RV_p     < ALPHA_1))
}

report_scan <- function(df, block, anchor_l, focus_lo, label) {
  cat("\n", strrep("=", 70), "\n", sep = "")
  cat(label, "\n")
  cat(strrep("=", 70), "\n", sep = "")
  cat("  * marks the main analysis setting.\n")
  cat("  READ THIS TABLE WITH TWO CAVEATS:\n")
  cat("   (1) Truncations are NOT nested models. K parts give K-1 ILR coordinates,\n")
  cat("       and PERMANOVA R2 is mechanically higher in lower dimensions. Part of\n")
  cat("       any R2 difference across rows is geometry, not signal. Stability is\n")
  cat("       judged on resolved-pair counts, median -log10(p) and Mantel/RV.\n")
  cat("   (2) The median and the counts see different things. The median tracks the\n")
  cat("       middle-ranked pair (overall margin); the counts track whether marginal\n")
  cat("       pairs survive. Neither overrules the other.\n\n")
  print(as.data.frame(scan_table(df, block)), row.names = FALSE)

  cat("\n  -- Identity of the pairs that move relative to the anchor (Holm p < 0.05) --\n")
  v <- if (block == "SP") "SP" else "M"
  for (a in c("EXP", "SDG")) {
    d <- df %>% filter(assemblage == a) %>% arrange(varied_lmax)
    base_set <- str_split(d[[paste0(v, "_sig_pairs")]][d$is_anchor], "; ")[[1]]
    base_set <- base_set[nzchar(base_set)]
    cat(sprintf("    %s (anchor l = 1-%d resolves %d): \n", a, anchor_l, length(base_set)))
    for (i in seq_len(nrow(d))) {
      if (d$is_anchor[i]) next
      s <- str_split(d[[paste0(v, "_sig_pairs")]][i], "; ")[[1]]
      s <- s[nzchar(s)]
      lost   <- setdiff(base_set, s)
      gained <- setdiff(s, base_set)
      cat(sprintf("      l = 1-%-2d  lost: %-28s gained: %s\n", d$varied_lmax[i],
                  ifelse(length(lost)   == 0, "none", paste(lost,   collapse = ", ")),
                  ifelse(length(gained) == 0, "none", paste(gained, collapse = ", "))))
    }
  }

  verdicts <- list()
  cat(sprintf("\n  -- The comparison the reviewer asked about: l = 1-%d vs l = 1-%d --\n",
              focus_lo, anchor_l))
  for (a in c("EXP", "SDG")) {
    v <- compare_two(df, block, focus_lo, anchor_l, a)
    if (!is.null(v)) verdicts[[a]] <- v
    cat("\n")
  }
  sensitive <- any(map_lgl(verdicts, ~ .x$d05 != 0 || .x$d01 != 0 ||
                                        .x$flip_global || .x$flip_mantel || .x$flip_rv))
  cat(sprintf("  VERDICT: the conclusions ARE%s sensitive to l = 1-%d vs l = 1-%d.\n",
              ifelse(sensitive, "", " NOT"), focus_lo, anchor_l))
  if (sensitive) {
    cat("    At least one count-based or decoupling verdict changes between the two\n")
    cat("    settings; the differences are itemised above and must be reported as found.\n")
  } else {
    cat("    Resolved-pair counts at both alpha levels are identical, the global test\n")
    cat("    stays on the same side of 0.05, and the Mantel/RV decoupling verdicts are\n")
    cat("    unchanged. The median -log10(p) shifts are listed above as magnitudes only.\n")
    cat("    This is a null result: it says the choice does not affect the conclusions.\n")
    cat("    It does NOT say degree ", anchor_l, " is better than degree ", focus_lo, ".\n", sep = "")
  }
  wins <- list()
  cat("\n  Stable window (contiguous run around the anchor with identical resolved-pair\n")
  cat("  counts at both alpha levels, same side of 0.05 globally, and unchanged Mantel/RV verdicts):\n")
  for (a in c("EXP", "SDG")) {
    w <- stable_window(df %>% filter(assemblage == a), block)
    wins[[a]] <- w
    edge <- c(if (w["lo"] > min(df$varied_lmax)) sprintf("l = 1-%d changes something below", w["lo"] - 1),
              if (w["hi"] < max(df$varied_lmax)) sprintf("l = 1-%d changes something above", w["hi"] + 1))
    cat(sprintf("    %s : l = 1-%d .. 1-%d%s\n", a, w["lo"], w["hi"],
                ifelse(length(edge) == 0, "  (stable across the whole swept range)",
                       paste0("   [", paste(edge, collapse = "; "), "]"))))
  }
  invisible(list(sensitive = sensitive, verdicts = verdicts, windows = wins))
}

# =============================================================================
# (1) SCAN A and (2) SCAN B
# =============================================================================
resA <- report_scan(scanA, "SP", ANCHOR_SP, ANCHOR_SP - 1,
                    sprintf("(1) SCAN A — SP-SPHARM truncation (M held at l = 1-%d)", ANCHOR_M))

cat("\n  -- Tail check (TRAP 2): the two EXP pairs SP-SPHARM alone never resolves --\n")
print(as.data.frame(
  scanA %>% filter(assemblage == "EXP") %>%
    transmute(truncation = sprintf("1-%d%s", varied_lmax, ifelse(is_anchor, "*", " ")),
              resolved = sprintf("%d/%d", SP_nsig05, n_pairs),
              `Discoid-Bidirectional (Holm p)` = round(focus_p_Discoid_Bidirectional, 4),
              `Discoid-Multiplatform (Holm p)` = round(focus_p_Discoid_Multiplatform, 4))),
  row.names = FALSE)
cat("  (These stay unresolved by design — they are the reason the joint-MFA analysis exists.\n")
cat("   What matters here is whether the OTHER pairs' resolution moves with truncation.)\n")

resB <- report_scan(scanB, "M", ANCHOR_M, ANCHOR_M - 1,
                    sprintf("(2) SCAN B — M-SPHARM truncation (SP held at l = 1-%d)", ANCHOR_SP))

# =============================================================================
# (3) DIAGNOSTICS vs THE METHODS STATEMENT
# =============================================================================
cat("\n", strrep("=", 70), "\n", sep = "")
cat("(3) CUMULATIVE POWER AND CROSS-SPECIMEN CV vs THE METHODS STATEMENT\n")
cat(strrep("=", 70), "\n", sep = "")
cat("  Diagnostic sets are the degree-selection sets of power_degree_selection.R\n")
cat("  (EXP = all EXP specimens; SDG = all SDG specimens), not the PERMANOVA samples.\n\n")

diag_tbl <- map_dfr(c("EXP", "SDG"), function(a)
  map_dfr(c("SP", "M"), function(b) {
    ks <- if (b == "SP") SP_GRID else M_GRID
    tibble(assemblage = a,
           descriptor = ifelse(b == "SP", "SP-SPHARM", "M-SPHARM"),
           truncation = sprintf("1-%d%s", ks,
                                ifelse(ks == (if (b == "SP") ANCHOR_SP else ANCHOR_M), "*", " ")),
           `cum power %` = round(map_dbl(ks, ~ cum_at(a, b, .x)), 3),
           `CV % at l`   = round(map_dbl(ks, ~ cv_at(a, b, .x)), 2))
  }))
print(as.data.frame(diag_tbl), row.names = FALSE)

cat("\n  Claims in the Methods, evaluated against these numbers:\n")
for (a in c("EXP", "SDG")) {
  sp5 <- cum_at(a, "SP", 5); sp6 <- cum_at(a, "SP", ANCHOR_SP)
  spx <- cv_cross(a, "SP")
  m8  <- cum_at(a, "M", ANCHOR_M); mmx <- cv_max(a, "M"); marg <- cv_argmax(a, "M")
  mx  <- cv_cross(a, "M")
  cat(sprintf("\n   [%s]\n", a))
  cat(sprintf("    SP power 'effectively exhausted by degree 6 (~99.9%%)' : cum(l=6) = %.3f%%  -> %s\n",
              sp6, ifelse(abs(sp6 - 99.9) < 0.15, "ACCURATE", "CHECK WORDING")))
  cat(sprintf("    SP 'noise-dominated, CV > 100%% from l ~ 9'            : first CV>100%% at l = %s  -> %s\n",
              ifelse(is.na(spx), "never", as.character(spx)),
              ifelse(!is.na(spx) && abs(spx - 9) <= 1, "ACCURATE", "CHECK WORDING")))
  cat(sprintf("    M  'CV < 90%% throughout'                             : max CV = %.2f%% at l = %d  -> %s\n",
              mmx, marg, ifelse(mmx < 90, "ACCURATE", "CHECK WORDING")))
  cat(sprintf("    M  '~98%% of cumulative power by degree 8'            : cum(l=8) = %.3f%%  -> %s\n",
              m8, ifelse(abs(m8 - 98) < 1, "ACCURATE", "CHECK WORDING")))
  cat(sprintf("    M  CV > 100%% first at l = %s (within l = 1-%d)\n",
              ifelse(is.na(mx), "never", as.character(mx)), length(COLS_ALL)))
  cat(sprintf("    Reviewer's premise, quantified: a pure cumulative-power rule already\n"))
  cat(sprintf("      passes 99%% at l = 5 (cum = %.3f%%), so power alone cannot separate\n", sp5))
  cat(sprintf("      l = 5 from l = 6 (cum = %.3f%%). The CV criterion is what places the\n", sp6))
  cat(sprintf("      noise floor at l ~ %s — on that criterion l = 6 is conservative, not aggressive.\n",
              ifelse(is.na(spx), "n/a", as.character(spx))))
}

cat(sprintf("\n  Closure/zero handling: replace_zeros was %s over the swept ranges,\n",
            ifelse(CLOSURE_INERT, "inert (no zero entries)", "ACTIVE (zero entries present)")))
cat(sprintf("  so row closure to sum 1 %s the ILR coordinates.\n",
            ifelse(CLOSURE_INERT, "provably does not change",
                   "COULD change — the anchor check is the arbiter")))

# =============================================================================
# (4) SUMMARY
# =============================================================================
cat("\n", strrep("=", 70), "\n", sep = "")
cat("(4) SUMMARY\n")
cat(strrep("=", 70), "\n", sep = "")
wA <- resA$windows
wB <- resB$windows
covers <- function(w, l) unname(l >= w["lo"] && l <= w["hi"])
alt_SP <- ANCHOR_SP - 1
alt_M  <- ANCHOR_M  - 1
sp_alt_in <- all(map_lgl(wA, covers, l = alt_SP))
m_alt_in  <- all(map_lgl(wB, covers, l = alt_M))

cat(sprintf("  Stable windows around the main analysis setting (SP l = 1-%d, M l = 1-%d):\n",
            ANCHOR_SP, ANCHOR_M))
cat(sprintf("    SP truncation : EXP l = 1-%d..1-%d | SDG l = 1-%d..1-%d\n",
            wA$EXP["lo"], wA$EXP["hi"], wA$SDG["lo"], wA$SDG["hi"]))
cat(sprintf("    M  truncation : EXP l = 1-%d..1-%d | SDG l = 1-%d..1-%d\n",
            wB$EXP["lo"], wB$EXP["hi"], wB$SDG["lo"], wB$SDG["hi"]))
cat(sprintf("\n  ONE-SENTENCE SUMMARY: the main analysis setting sits inside a stable interval\n"))
cat(sprintf("  in both blocks and both assemblages, but that interval does %s the l = 1-%d\n",
            ifelse(sp_alt_in, "include", "NOT extend down to"), alt_SP))
cat(sprintf("  alternative the reviewer proposed for SP-SPHARM%s.\n",
            ifelse(sp_alt_in,
                   ", so 5 and 6 are interchangeable on the evidence here",
                   " — dropping to l = 1-5 changes a count-based or decoupling verdict (see (1))")))
cat(sprintf("  The M-SPHARM setting is %s to the analogous one-degree change (l = 1-%d).\n",
            ifelse(m_alt_in, "insensitive", "sensitive"), alt_M))

cat("\n  Scope, stated precisely so the table is not over-read:\n")
cat("   - Where nothing changes, that is a null result: the choice does not affect the\n")
cat("     conclusions. It is NOT evidence that one degree outperforms another.\n")
cat("   - Where something does change, the change is reported as found, with the\n")
cat("     affected pairs named in (1)/(2). It is not smoothed over to protect the\n")
cat("     published setting.\n")
cat("   - No ranking of truncations by R2 is offered anywhere: R2 is not comparable\n")
cat("     across truncations of differing dimension (TRAP 1).\n")
if (!ok_fig) cat("\n  ! Figures were not written — see the error above.\n")

cat("\nDone.\n\n")
print(sessionInfo())
