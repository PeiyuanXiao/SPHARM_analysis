# power_degree_selection.R
# Power-spectrum degree-selection diagnostic (morphology & direction, all 20 degrees).
# Sourced by the `degree_selection_diagnostic` target; also runnable standalone.
# Computes per-degree descriptive statistics, writes them to CSV, and builds the
# composite diagnostic figure (per-degree CV and cumulative power, EXP and SDG)
# inserted into the manuscript.
#
# Input:
#   analysis/data/derived_data/SPHARM_direction.csv
#   analysis/data/derived_data/SPHARM_morphology.csv
# Output (CSV):
#   analysis/data/derived_data/DegreeSelection_stats_direction_EXP.csv
#   analysis/data/derived_data/DegreeSelection_stats_morphology_EXP.csv
#   analysis/data/derived_data/DegreeSelection_stats_direction_SDG.csv
#   analysis/data/derived_data/DegreeSelection_stats_morphology_SDG.csv
# Output (object):
#   p_degree_selection  (patchwork composite: row 1 = EXP, row 2 = SDG; no titles)

library(here)
library(tidyverse)
library(patchwork)


# ---- Parameters ----

POWER_COLS_ALL <- paste0("power_l", 1:20)
EXP_PREFIX     <- "EXP"
SDG_PREFIX     <- "SDG"

DESCRIPTOR_COLORS <- c("M-SPHARM" = "#4A6E8A", "SP-SPHARM" = "#BA8530")


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

compute_degree_stats <- function(df, cols, label) {

  mat       <- df %>% select(all_of(cols)) %>% as.matrix()
  n_degrees <- length(cols)

  col_means  <- colMeans(mat, na.rm = TRUE)
  col_sds    <- apply(mat, 2, sd,  na.rm = TRUE)
  col_vars   <- apply(mat, 2, var, na.rm = TRUE)
  col_cvs    <- col_sds / col_means * 100
  row_sums   <- rowSums(mat, na.rm = TRUE)
  total_mean <- mean(row_sums)
  cumul_pct  <- cumsum(col_means) / total_mean * 100
  pct_each   <- col_means / total_mean * 100
  decay_rate <- c(NA, col_means[-1] / col_means[-n_degrees])

  stats_df <- tibble(
    source       = label,
    degree       = seq_len(n_degrees),
    degree_label = paste0("l=", seq_len(n_degrees)),
    mean         = round(col_means,  6),
    sd           = round(col_sds,    6),
    variance     = round(col_vars,   6),
    cv_pct       = round(col_cvs,    2),
    pct_energy   = round(pct_each,   3),
    cumul_pct    = round(cumul_pct,  3),
    decay_rate   = round(decay_rate, 4)
  )

  cat(sprintf("\n======== %s per-degree descriptive statistics ========\n", label))
  print(stats_df %>%
          select(degree, mean, sd, variance, cv_pct,
                 pct_energy, cumul_pct, decay_rate) %>%
          as.data.frame())

  cat(sprintf("\n  Row sums: min=%.4f, max=%.4f, mean=%.4f, sd=%.4f\n",
              min(row_sums), max(row_sums), mean(row_sums), sd(row_sums)))

  for (thr in c(90, 95, 99, 99.5)) {
    k <- which(cumul_pct >= thr)[1]
    cat(sprintf("  Cumulative energy >= %5.1f%%: first %d degrees\n", thr, k))
  }

  stats_df
}

stats_dir_exp   <- compute_degree_stats(dir_exp,   dir_cols_exp,   "Direction (EXP)")
stats_morph_exp <- compute_degree_stats(morph_exp, morph_cols_exp, "Morphology (EXP)")
stats_dir_sdg   <- compute_degree_stats(dir_sdg,   dir_cols_sdg,   "Direction (SDG)")
stats_morph_sdg <- compute_degree_stats(morph_sdg, morph_cols_sdg, "Morphology (SDG)")

write_csv(stats_dir_exp,
          here("analysis/data/derived_data/DegreeSelection_stats_direction_EXP.csv"))
write_csv(stats_morph_exp,
          here("analysis/data/derived_data/DegreeSelection_stats_morphology_EXP.csv"))
write_csv(stats_dir_sdg,
          here("analysis/data/derived_data/DegreeSelection_stats_direction_SDG.csv"))
write_csv(stats_morph_sdg,
          here("analysis/data/derived_data/DegreeSelection_stats_morphology_SDG.csv"))
cat("\nSaved: 4 DegreeSelection_stats_*.csv files\n")


# ---- Per-degree comparison print ----

cat("\n======== EXP per-degree comparison ========\n")
bind_rows(stats_dir_exp, stats_morph_exp) %>%
  select(source, degree, mean, sd, cv_pct, cumul_pct, decay_rate) %>%
  mutate(across(where(is.numeric), ~ round(.x, 4))) %>%
  print(n = 50)

cat("\n======== SDG per-degree comparison ========\n")
bind_rows(stats_dir_sdg, stats_morph_sdg) %>%
  select(source, degree, mean, sd, cv_pct, cumul_pct, decay_rate) %>%
  mutate(across(where(is.numeric), ~ round(.x, 4))) %>%
  print(n = 50)


# ---- Composite diagnostic figure (no titles) ----
# Row 1: EXP across-specimen CV | EXP cumulative power
# Row 2: SDG across-specimen CV | SDG cumulative power
# Each panel overlays the morphology (M-SPHARM) and direction (SP-SPHARM) spectra.

degree_plot_df <- bind_rows(
  stats_dir_exp   %>% mutate(dataset = "EXP", descriptor = "SP-SPHARM"),
  stats_morph_exp %>% mutate(dataset = "EXP", descriptor = "M-SPHARM"),
  stats_dir_sdg   %>% mutate(dataset = "SDG", descriptor = "SP-SPHARM"),
  stats_morph_sdg %>% mutate(dataset = "SDG", descriptor = "M-SPHARM")
) %>%
  mutate(descriptor = factor(descriptor, levels = c("M-SPHARM", "SP-SPHARM")))

# dataset label shown in a panel corner (replaces the EXP/SDG y-axis prefix)
ds_label <- function(ds) if (identical(ds, "EXP")) "Experimental cores" else "Sandinggai cores"

# red dashed vertical lines marking the retained truncation degrees (l = 6, 8)
trunc_lines <- geom_vline(xintercept = c(6, 8), color = "red",
                          linetype = "dashed", linewidth = 0.3)

make_cv_panel <- function(ds, y_nbreaks = NULL) {
  p <- ggplot(filter(degree_plot_df, dataset == ds),
         aes(degree, cv_pct, color = descriptor)) +
    geom_hline(yintercept = 100, linetype = "dashed",
               color = "grey55", linewidth = 0.3) +
    trunc_lines +
    geom_line(linewidth = 0.6) +
    geom_point(size = 1.4) +
    annotate("text", x = -Inf, y = 100, label = "CV = 100%",
             hjust = -0.06, vjust = -0.4, size = 2.6, color = "grey40") +
    annotate("text", x = -Inf, y = Inf, label = ds_label(ds),
             hjust = -0.06, vjust = 1.5, size = 3, fontface = "bold") +
    scale_color_manual(values = DESCRIPTOR_COLORS) +
    scale_x_continuous(breaks = seq(2, 20, 2)) +
    labs(x = "Spherical harmonic degree (l)",
         y = "Across-specimen CV (%)",
         color = NULL) +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank())
  if (!is.null(y_nbreaks)) p <- p + scale_y_continuous(n.breaks = y_nbreaks)
  p
}

make_cumulative_panel <- function(ds) {
  ggplot(filter(degree_plot_df, dataset == ds),
         aes(degree, cumul_pct, color = descriptor)) +
    geom_hline(yintercept = c(95, 99), linetype = "dashed",
               color = "grey55", linewidth = 0.3) +
    trunc_lines +
    geom_line(linewidth = 0.6) +
    geom_point(size = 1.4) +
    annotate("text", x = -Inf, y = 99, label = "Power = 99%",
             hjust = -0.06, vjust = -0.4, size = 2.6, color = "grey40") +
    annotate("text", x = -Inf, y = 95, label = "Power = 95%",
             hjust = -0.06, vjust = 1.4, size = 2.6, color = "grey40") +
    annotate("text", x = Inf, y = -Inf, label = ds_label(ds),
             hjust = 1.06, vjust = -0.7, size = 3, fontface = "bold") +
    scale_color_manual(values = DESCRIPTOR_COLORS) +
    scale_x_continuous(breaks = seq(2, 20, 2)) +
    labs(x = "Spherical harmonic degree (l)",
         y = "Cumulative power (%)",
         color = NULL) +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank())
}

p_degree_selection <- (
  (make_cv_panel("EXP")              | make_cumulative_panel("EXP")) /
    (make_cv_panel("SDG", y_nbreaks = 10) | make_cumulative_panel("SDG"))
) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "a")
# patchwork attaches the tags to the subplots with labs(tag = ), so plot.tag has
# to be set on the subplots with `&`; the plot_annotation() theme only styles the
# tag of the composite itself and never reaches the panels.
p_degree_selection <- p_degree_selection &
  theme(legend.position = "bottom",
        plot.tag = element_text(face = "bold", size = 14))

cat("\n========== Degree-selection diagnostic complete ==========\n")
