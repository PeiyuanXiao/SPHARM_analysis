# power_order_selection.R
# Power-spectrum order-selection diagnostic (morphology & direction, all 20 degrees).
# One-off helper, outside the targets pipeline; produces per-degree descriptive
# statistics and writes them to CSV. Plotting is done separately.
#
# Input:
#   analysis/data/derived_data/SPHARM_direction.csv
#   analysis/data/derived_data/SPHARM_morphology.csv
# Output:
#   analysis/data/derived_data/OrderSelection_stats_direction_EXP.csv
#   analysis/data/derived_data/OrderSelection_stats_morphology_EXP.csv
#   analysis/data/derived_data/OrderSelection_stats_direction_SDG.csv
#   analysis/data/derived_data/OrderSelection_stats_morphology_SDG.csv

library(here)
library(tidyverse)


# ---- Parameters ----

POWER_COLS_ALL <- paste0("power_l", 1:20)
EXP_PREFIX     <- "EXP"
SDG_PREFIX     <- "SDG"


# ---- Load data ----

SPHARM_direction  <- read_csv(
  here("analysis/data/derived_data/SPHARM_direction.csv"),
  show_col_types = FALSE
)
SPHARM_morphology <- read_csv(
  here("analysis/data/derived_data/SPHARM_morphology.csv"),
  show_col_types = FALSE
)

# EXP experimental specimens.
dir_exp   <- SPHARM_direction  %>% filter(str_starts(ID, EXP_PREFIX))
morph_exp <- SPHARM_morphology %>% filter(str_starts(ID, EXP_PREFIX))

# SDG archaeological specimens (excluding IM_ references).
dir_sdg   <- SPHARM_direction  %>%
  filter(str_starts(ID, SDG_PREFIX), !str_starts(ID, "IM_"))
morph_sdg <- SPHARM_morphology %>%
  filter(str_starts(ID, SDG_PREFIX), !str_starts(ID, "IM_"))

cat(sprintf("EXP specimens: direction n=%d, morphology n=%d\n",
            nrow(dir_exp), nrow(morph_exp)))
cat(sprintf("SDG specimens: direction n=%d, morphology n=%d\n",
            nrow(dir_sdg), nrow(morph_sdg)))

check_cols <- function(df, label) {
  missing   <- setdiff(POWER_COLS_ALL, colnames(df))
  available <- intersect(POWER_COLS_ALL, colnames(df))
  if (length(missing) > 0) {
    warning(sprintf("%s missing columns: %s", label, paste(missing, collapse = ", ")))
    cat(sprintf("  %s: %d degrees available (%s ~ %s)\n",
                label, length(available),
                available[1], available[length(available)]))
  } else {
    cat(sprintf("  %s: all 20 degrees available\n", label))
  }
  available
}

dir_cols_exp   <- check_cols(dir_exp,   "EXP direction")
morph_cols_exp <- check_cols(morph_exp, "EXP morphology")
dir_cols_sdg   <- check_cols(dir_sdg,   "SDG direction")
morph_cols_sdg <- check_cols(morph_sdg, "SDG morphology")


# ---- Per-degree descriptive statistics ----

compute_order_stats <- function(df, cols, label) {
  
  mat      <- df %>% select(all_of(cols)) %>% as.matrix()
  n_orders <- length(cols)
  
  col_means  <- colMeans(mat, na.rm = TRUE)
  col_sds    <- apply(mat, 2, sd,  na.rm = TRUE)
  col_vars   <- apply(mat, 2, var, na.rm = TRUE)
  col_cvs    <- col_sds / col_means * 100
  col_snr    <- col_means / col_sds
  row_sums   <- rowSums(mat, na.rm = TRUE)
  total_mean <- mean(row_sums)
  cumul_pct  <- cumsum(col_means) / total_mean * 100
  pct_each   <- col_means / total_mean * 100
  decay_rate <- c(NA, col_means[-1] / col_means[-n_orders])
  
  stats_df <- tibble(
    source      = label,
    order       = seq_len(n_orders),
    order_label = paste0("l=", seq_len(n_orders)),
    mean        = round(col_means,  6),
    sd          = round(col_sds,    6),
    variance    = round(col_vars,   6),
    cv_pct      = round(col_cvs,    2),
    snr         = round(col_snr,    4),
    pct_energy  = round(pct_each,   3),
    cumul_pct   = round(cumul_pct,  3),
    decay_rate  = round(decay_rate, 4)
  )
  
  cat(sprintf("\n======== %s per-degree descriptive statistics ========\n", label))
  print(stats_df %>%
          select(order, mean, sd, variance, cv_pct, snr,
                 pct_energy, cumul_pct, decay_rate) %>%
          as.data.frame())
  
  cat(sprintf("\n  Row sums: min=%.4f, max=%.4f, mean=%.4f, sd=%.4f\n",
              min(row_sums), max(row_sums), mean(row_sums), sd(row_sums)))
  
  snr_drop <- which(diff(col_snr) < -0.3)
  if (length(snr_drop) > 0)
    cat(sprintf("  Marked SNR drop (after this degree): l=%s\n",
                paste(snr_drop, collapse = ", ")))
  
  for (thr in c(90, 95, 99, 99.5)) {
    k <- which(cumul_pct >= thr)[1]
    cat(sprintf("  Cumulative energy >= %5.1f%%: first %d degrees\n", thr, k))
  }
  
  stats_df
}

stats_dir_exp   <- compute_order_stats(dir_exp,   dir_cols_exp,   "Direction (EXP)")
stats_morph_exp <- compute_order_stats(morph_exp, morph_cols_exp, "Morphology (EXP)")
stats_dir_sdg   <- compute_order_stats(dir_sdg,   dir_cols_sdg,   "Direction (SDG)")
stats_morph_sdg <- compute_order_stats(morph_sdg, morph_cols_sdg, "Morphology (SDG)")

write_csv(stats_dir_exp,
          here("analysis/data/derived_data/OrderSelection_stats_direction_EXP.csv"))
write_csv(stats_morph_exp,
          here("analysis/data/derived_data/OrderSelection_stats_morphology_EXP.csv"))
write_csv(stats_dir_sdg,
          here("analysis/data/derived_data/OrderSelection_stats_direction_SDG.csv"))
write_csv(stats_morph_sdg,
          here("analysis/data/derived_data/OrderSelection_stats_morphology_SDG.csv"))
cat("\nSaved: 4 OrderSelection_stats_*.csv files\n")


# ---- Per-degree comparison print ----

cat("\n======== EXP per-degree comparison ========\n")
bind_rows(stats_dir_exp, stats_morph_exp) %>%
  select(source, order, mean, sd, snr, cv_pct, cumul_pct, decay_rate) %>%
  mutate(across(where(is.numeric), ~ round(.x, 4))) %>%
  print(n = 50)

cat("\n======== SDG per-degree comparison ========\n")
bind_rows(stats_dir_sdg, stats_morph_sdg) %>%
  select(source, order, mean, sd, snr, cv_pct, cumul_pct, decay_rate) %>%
  mutate(across(where(is.numeric), ~ round(.x, 4))) %>%
  print(n = 50)

cat("\n========== Order-selection diagnostic complete ==========\n")
