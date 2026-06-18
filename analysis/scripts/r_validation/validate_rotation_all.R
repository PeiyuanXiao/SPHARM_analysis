# validate_rotation_all.R
# Rotational-invariance validation: R / E / I + SPHARM power spectrum.
# Method: Bland-Altman analysis (pairwise comparison of alignment conditions).
# Sourced by the `p_rotational_invariance_validity` target.
#
# Prerequisites:
#   1. align_svd.R     -> directions_raw.csv + directions_aligned_svd.csv
#   2. align_lin2024.R -> directions_aligned_lin2024.csv
#   3. python kde_to_spharm_main.py --source all
#
# Output:
#   console: bias / LoA summary table for all metrics
#   returns (object): p_rotational_invariance_validity

library(here)
library(tidyverse)
library(patchwork)
library(glue)

# ==============================================================================
# Helpers
# ==============================================================================

# --- direction metrics ---
compute_R <- function(ux, uy, uz) {
  resultant <- sqrt(sum(ux)^2 + sum(uy)^2 + sum(uz)^2)
  total_len <- sum(sqrt(ux^2 + uy^2 + uz^2))
  resultant / total_len
}

compute_EI <- function(ux, uy, uz) {
  U      <- cbind(ux, uy, uz)
  T_mat  <- (t(U) %*% U) / nrow(U)
  lambda <- sort(eigen(T_mat, symmetric = TRUE)$values, decreasing = TRUE)
  lambda <- pmax(lambda, 0)
  list(
    E = ifelse(lambda[1] > 1e-10, 1 - lambda[2] / lambda[1], NA_real_),
    I = ifelse(lambda[1] > 1e-10,     lambda[3] / lambda[1], NA_real_)
  )
}

# --- Bland-Altman core ---
bland_altman_calc <- function(x, y) {
  mean_xy   <- (x + y) / 2
  diff_xy   <- x - y
  bias      <- mean(diff_xy,  na.rm = TRUE)
  sd_diff   <- sd(diff_xy,    na.rm = TRUE)
  loa_upper <- bias + 1.96 * sd_diff
  loa_lower <- bias - 1.96 * sd_diff
  list(
    df        = data.frame(mean = mean_xy, diff = diff_xy),
    bias      = bias,
    sd_diff   = sd_diff,
    loa_upper = loa_upper,
    loa_lower = loa_lower
  )
}

# --- single Bland-Altman plot ---
plot_ba <- function(ba, title_str, x_label = "Mean", y_label = "Difference") {
  ggplot(ba$df, aes(x = mean, y = diff)) +
    geom_hline(yintercept = ba$bias,      color = "#802520",
               linewidth = 0.45, linetype = "dashed") +
    geom_hline(yintercept = ba$loa_upper, color = "#5C7F71",
               linewidth = 0.45, linetype = "dotted") +
    geom_hline(yintercept = ba$loa_lower, color = "#5C7F71",
               linewidth = 0.45, linetype = "dotted") +
    geom_point(size = 1, alpha = 0.75, color = "#B8B8B8") +
    scale_y_continuous(labels = \(x) formatC(x, format = "e", digits = 1)) +
    annotate("text",
             x = -Inf, y = ba$bias,      hjust = -0.1, vjust = -0.4,
             label = sprintf("Bias = %.2e", ba$bias),
             color = "#802520", size = 2) +
    annotate("text",
             x = -Inf, y = ba$loa_upper, hjust = -0.1, vjust = -0.4,
             label = sprintf("+1.96 SD = %.2e", ba$loa_upper),
             color = "#5C7F71", size = 1.8) +
    annotate("text",
             x = -Inf, y = ba$loa_lower, hjust = -0.1, vjust =  1.4,
             label = sprintf("-1.96 SD = %.2e", ba$loa_lower),
             color = "#5C7F71", size = 1.8) +
    theme_bw(base_size = 6) +
    labs(title = title_str, x = x_label, y = y_label) +
    theme(panel.grid = element_blank(),
          plot.title = element_text(face = "bold", size = 6, hjust = 0.5))
}

# --- summary-table helper ---
summary_row <- function(metric, pair, ba) {
  tibble(
    metric    = metric,
    pair      = pair,
    bias      = ba$bias,
    sd_diff   = ba$sd_diff,
    loa_lower = ba$loa_lower,
    loa_upper = ba$loa_upper
  )
}

# ==============================================================================
# Part A: R / E / I validation
# ==============================================================================

cat("====== Part A: R / E / I ======\n\n")

# --- load direction-vector CSVs ---
read_directions <- function(source) {
  path <- here(glue("analysis/data/derived_data/directions_{source}.csv"))
  if (!file.exists(path))
    stop(glue("Not found: {path}"))
  read_csv(path, show_col_types = FALSE) %>% mutate(source = source)
}

dirs <- bind_rows(
  read_directions("raw"),
  read_directions("aligned_svd"),
  read_directions("aligned_lin2024")
) %>%
  filter(!str_starts(ID, "IM_"))

# Unified source labels
dirs <- dirs %>%
  mutate(source = case_when(
    source == "raw"              ~ "raw",
    source == "aligned_svd"     ~ "svd",
    source == "aligned_lin2024" ~ "lin2024"
  ))

# --- compute R / E / I (single eigen-decomposition) ---
rei <- dirs %>%
  group_by(ID, source) %>%
  summarise(
    R  = compute_R(ux, uy, uz),
    ei = list(compute_EI(ux, uy, uz)),
    .groups = "drop"
  ) %>%
  mutate(
    E = map_dbl(ei, "E"),
    I = map_dbl(ei, "I")
  ) %>%
  select(-ei)

# Wide format (one row per specimen; columns per source)
rei_wide <- rei %>%
  pivot_wider(names_from = source,
              values_from = c(R, E, I),
              names_glue = "{.value}_{source}")

common_ids_rei <- rei_wide %>%
  filter(if_all(everything(), ~ !is.na(.))) %>%
  pull(ID)

rei_wide <- rei_wide %>% filter(ID %in% common_ids_rei)
cat(sprintf("R/E/I validation: %d specimens\n\n", nrow(rei_wide)))

# --- comparison pairs ---
# each pair: c(source_a, source_b, display_label)
pairs_label <- list(
  c("raw",    "svd",     "none-aligned vs techno-aligned"),
  c("raw",    "lin2024", "none-aligned vs morph-aligned"),
  c("svd",    "lin2024", "techno-aligned vs morph-aligned")
)

ba_plots_rei <- list()
summary_table <- tibble()

metric_labels <- c(
  R = "SPI",
  E = "Elongation",
  I = "Isotropy"
)

for (metric in c("R", "E", "I")) {
  for (pair in pairs_label) {
    src_a      <- pair[1]
    src_b      <- pair[2]
    pair_label <- pair[3]
    col_a      <- glue("{metric}_{src_a}")
    col_b      <- glue("{metric}_{src_b}")
    
    metric_name <- metric_labels[[metric]]
    
    ba  <- bland_altman_calc(rei_wide[[col_a]], rei_wide[[col_b]])
    plt <- plot_ba(ba,
                   title_str = glue("{metric_name}: {pair_label}"),
                   x_label   = glue("Mean of measures"),
                   y_label   = glue("Difference of measures"))
    
    ba_plots_rei[[glue("{metric}_{pair_label}")]] <- plt
    summary_table <- bind_rows(summary_table,
                               summary_row(metric, pair_label, ba))
  }
}

# Compose: 3 rows (metrics) x 3 cols (comparison pairs)
p_rei <- wrap_plots(ba_plots_rei, ncol = 3)

# ==============================================================================
# Part B: SPHARM power-spectrum validation (per-degree Bland-Altman, faceted)
# ==============================================================================

cat("====== Part B: SPHARM power spectrum ======\n\n")

# --- load the three power-spectrum CSVs ---
read_spharm <- function(source) {
  path <- here(glue("analysis/data/derived_data/validation/{source}/SPHARM_direction.csv"))
  if (!file.exists(path))
    stop(glue("Not found: {path}\nRun first: python kde_to_spharm_main.py --source all"))
  read_csv(path, show_col_types = FALSE) %>% mutate(source = source)
}

spharm_all <- bind_rows(
  read_spharm("raw"),
  read_spharm("svd"),
  read_spharm("lin2024")
) %>%
  filter(!str_starts(ID, "IM_"))

# Keep the svd frame for Part C
df_svd <- spharm_all %>% filter(source == "svd")

power_cols <- spharm_all %>%
  select(starts_with("power_l")) %>%
  colnames()

# Wide format
spharm_wide <- spharm_all %>%
  select(ID, source, all_of(power_cols)) %>%
  pivot_wider(names_from  = source,
              values_from = all_of(power_cols),
              names_glue  = "{.value}__{source}")

common_ids_spharm <- spharm_wide %>%
  filter(if_all(everything(), ~ !is.na(.))) %>%
  pull(ID)
spharm_wide <- spharm_wide %>% filter(ID %in% common_ids_spharm)
cat(sprintf("SPHARM validation: %d specimens, %d degrees\n\n",
            nrow(spharm_wide), length(power_cols)))

# --- per-degree Bland-Altman, collected into a long table ---
ba_power_summary <- tibble()

for (pair in pairs_label) {
  src_a      <- pair[1]
  src_b      <- pair[2]
  pair_label <- pair[3]
  
  for (pcol in power_cols) {
    col_a <- glue("{pcol}__{src_a}")
    col_b <- glue("{pcol}__{src_b}")
    ba    <- bland_altman_calc(spharm_wide[[col_a]], spharm_wide[[col_b]])
    
    degree <- as.integer(str_remove(pcol, "power_l"))
    ba_power_summary <- bind_rows(ba_power_summary, tibble(
      degree    = degree,
      pair      = pair_label,
      bias      = ba$bias,
      loa_upper = ba$loa_upper,
      loa_lower = ba$loa_lower,
      sd_diff   = ba$sd_diff
    ))
    
    summary_table <- bind_rows(summary_table,
                               summary_row(pcol, pair_label, ba))
  }
}

# ==============================================================================
# Summary table
# ==============================================================================

cat("====== Bland-Altman summary, all metrics ======\n\n")

summary_table %>%
  mutate(across(c(bias, sd_diff, loa_lower, loa_upper),
                \(x) formatC(x, format = "e", digits = 3))) %>%
  arrange(metric, pair) %>%
  print(n = Inf)

write_csv(
  summary_table %>% arrange(metric, pair),
  here("analysis/data/derived_data/validation_ba_summary.csv")
)
cat("\nSaved: validation_ba_summary.csv\n")

# ==============================================================================
# Per-source numeric tables
# ==============================================================================

cat("\n====== Numeric summary: R / E / I ======\n\n")

rei_compare <- rei_wide %>%
  select(ID,
         R_raw, R_svd, R_lin2024,
         E_raw, E_svd, E_lin2024,
         I_raw, I_svd, I_lin2024) %>%
  arrange(ID)

rei_compare %>%
  mutate(across(where(is.numeric), \(x) round(x, 6))) %>%
  print(n = Inf)

write_csv(rei_compare,
          here("analysis/data/derived_data/validation_values_REI.csv"))
cat("Saved: validation_values_REI.csv\n")

cat("\n====== Numeric summary: power spectrum (l1-5) ======\n\n")

power_compare <- spharm_wide %>%
  select(ID, matches("^power_l[1-5]__")) %>%
  arrange(ID)

power_compare <- power_compare %>%
  rename_with(
    ~ str_replace(., "__", "_"),
    matches("^power_l[1-5]__")
  )

power_compare %>%
  mutate(across(where(is.numeric), \(x) round(x, 6))) %>%
  print(n = Inf)

write_csv(power_compare,
          here("analysis/data/derived_data/validation_values_power_l1_5.csv"))
cat("Saved: validation_values_power_l1_5.csv\n")

cat("\n====== Full power-spectrum values saved (not printed) ======\n")

power_all_compare <- spharm_wide %>%
  select(ID, matches("^power_l[0-9]+__")) %>%
  rename_with(~ str_replace(., "__", "_"),
              matches("^power_l[0-9]+__")) %>%
  arrange(ID)

write_csv(power_all_compare,
          here("analysis/data/derived_data/validation_values_power_all.csv"))
cat("Saved: validation_values_power_all.csv\n")

# ==============================================================================
# Part C: empirical check — SVD alignment vs SVD alignment + random Z rotation
# ==============================================================================

cat("\n====== Part C: SVD vs SVD + random Z rotation ======\n\n")

# --- load svd_rotated power spectrum ---
path_rotated <- here("analysis/data/derived_data/validation/svd_rotated/SPHARM_direction.csv")

if (!file.exists(path_rotated)) {
  stop(glue(
    "Not found: {path_rotated}\n",
    "Run first:\n",
    "  1. python rotate_svd_directions.py\n",
    "  2. python kde_to_spharm_main.py --source svd_rotated"
  ))
}

df_svd_rotated <- read_csv(path_rotated, show_col_types = FALSE) %>%
  mutate(source = "svd_rotated") %>%
  filter(!str_starts(ID, "IM_"))

# Specimens shared with svd
common_ids_rot <- intersect(
  df_svd %>% pull(ID),
  df_svd_rotated %>% pull(ID)
)
cat(sprintf("svd vs svd_rotated validation: %d specimens\n\n", length(common_ids_rot)))

# Wide-format merge
spharm_rot_wide <- bind_rows(
  df_svd         %>% filter(ID %in% common_ids_rot),
  df_svd_rotated %>% filter(ID %in% common_ids_rot)
) %>%
  select(ID, source, all_of(power_cols)) %>%
  pivot_wider(
    names_from  = source,
    values_from = all_of(power_cols),
    names_glue  = "{.value}__{source}"
  )

# ------------------------------------------------------------------------------
# C-1: per-degree power-spectrum Bland-Altman
# ------------------------------------------------------------------------------

ba_rot_summary <- tibble()

for (pcol in power_cols) {
  col_svd <- glue("{pcol}__svd")
  col_rot <- glue("{pcol}__svd_rotated")
  ba      <- bland_altman_calc(spharm_rot_wide[[col_svd]],
                               spharm_rot_wide[[col_rot]])
  degree  <- as.integer(str_remove(pcol, "power_l"))
  
  ba_rot_summary <- bind_rows(ba_rot_summary, tibble(
    degree    = degree,
    bias      = ba$bias,
    loa_upper = ba$loa_upper,
    loa_lower = ba$loa_lower,
    sd_diff   = ba$sd_diff
  ))
  
  summary_table <- bind_rows(summary_table,
                             summary_row(pcol, "perturbed vs unperturbed", ba))
}

ba_rot_plot_df <- ba_rot_summary %>%
  mutate(pair = "perturbed vs unperturbed")

ba_combined <- bind_rows(
  ba_power_summary,
  ba_rot_plot_df
) %>%
  mutate(pair = factor(pair, levels = c(
    "none-aligned vs techno-aligned",
    "none-aligned vs morph-aligned",
    "techno-aligned vs morph-aligned",
    "perturbed vs unperturbed"
  )))

pair_colors <- c(
  "none-aligned vs techno-aligned" = "#802520",
  "none-aligned vs morph-aligned"  = "#BA8530",
  "techno-aligned vs morph-aligned"= "#5C7F71",
  "perturbed vs unperturbed"       = "#F5EDDC"
)

p_rot_spharm <- ggplot(ba_combined, aes(x = degree)) +
  geom_ribbon(aes(ymin = loa_lower, ymax = loa_upper,
                  fill = pair),
              alpha = 0.15) +
  geom_line(aes(y = bias, color = pair),
            linewidth = 0.75) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "grey50", linewidth = 0.2) +
  scale_color_manual(values = pair_colors) +
  scale_fill_manual(values  = pair_colors) +
  scale_x_continuous(breaks = seq(0, max(ba_combined$degree), by = 1)) +
  theme_bw(base_size = 6) +
  labs(
    x        = "SP-SPHARM power spectra degree (l)",
    y        = "Difference (bias ± 95% LoA)",
    color    = NULL,
    fill     = NULL
  ) +
  theme(
    panel.grid.major.x = element_line(color = "grey50", linewidth = 0.15, linetype = "dashed"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    plot.title      = element_text(face = "bold", size = 8, hjust = 0.5),
    plot.subtitle   = element_text(size = 7, hjust = 0.5, color = "grey40"),
    legend.key.size = unit(0.32, "cm"),
    legend.text     = element_text(size = 6),
    legend.position      = c(0.87, 0.2),
    legend.justification = c(0.5, 0.5),
    legend.background    = element_rect(fill = "transparent", colour = NA),
    legend.box.background = element_rect(fill = "transparent", colour = NA)
  )

# ------------------------------------------------------------------------------
# C-3: numeric summary and conclusion
# ------------------------------------------------------------------------------

cat("==== Part C numeric summary ====\n")

loa_width_intermethod <- ba_power_summary %>%
  group_by(pair) %>%
  summarise(max_loa_width = max(loa_upper - loa_lower), .groups = "drop")

loa_width_intramethod <- ba_rot_summary %>%
  summarise(max_loa_width = max(loa_upper - loa_lower)) %>%
  mutate(pair = "perturbed vs unperturbed")

cat("\nMax LoA width per comparison (smaller = better agreement):\n")
bind_rows(loa_width_intermethod, loa_width_intramethod) %>%
  mutate(max_loa_width = formatC(max_loa_width, format = "e", digits = 3)) %>%
  print()

max_intra <- max(ba_rot_summary$loa_upper - ba_rot_summary$loa_lower)
max_inter <- max(ba_power_summary$loa_upper - ba_power_summary$loa_lower)

cat(sprintf(
  "\nWithin-frame perturbation LoA width (%.2e) %s between-method LoA width (%.2e)\n",
  max_intra,
  ifelse(max_intra < max_inter, "<", "≥"),
  max_inter
))

if (max_intra < max_inter * 0.1) {
  cat("Conclusion: within-frame random perturbation is far smaller than between-method systematic difference (< 10%),\n")
  cat("            given a fixed SVD frame, results are unaffected by within-frame random error.\n")
} else if (max_intra < max_inter) {
  cat("Conclusion: within-frame random perturbation is smaller than between-method systematic difference,\n")
  cat("            given a fixed SVD frame, results are largely unaffected.\n")
} else {
  cat("Conclusion: within-frame random perturbation is non-negligible; check alignment quality or KDE parameters.\n")
}

# Update summary table
write_csv(
  summary_table %>% arrange(metric, pair),
  here("analysis/data/derived_data/validation_ba_summary.csv")
)
cat("\nUpdated: validation_ba_summary.csv\n")

# ==============================================================================
# Final combined plot
# ==============================================================================

p_rotational_invariance_validity <- p_rei / p_rot_spharm +
  plot_layout(heights = c(3, 1)) +
  plot_annotation(
    tag_levels = "a",
    theme = theme(
      plot.tag = element_text(face = "bold", size = 9)
    )
  )
