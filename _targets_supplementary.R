# ==============================================================================
# _targets_supplementary.R
# Standalone targets pipeline for the Supplementary Materials (robustness
# sensitivity analyses). INDEPENDENT of the main pipeline (_targets.R): it READS
# the committed main outputs (SPHARM_direction/morphology.csv) + the committed
# per-setting spectra as `file` inputs, sources each robustness 02_*.R, and
# returns figures + tables consumed by analysis/paper/supplementary.qmd.
#
# A reader who only wants the main analysis never builds this. To reproduce the
# SI (after the main pipeline has produced the committed derived CSVs):
#   Sys.setenv(TAR_PROJECT = "supplementary"); targets::tar_make()
#   quarto render analysis/paper/supplementary.qmd
# The `supplementary` project (script + its own store) is registered in
# _targets.yaml, so it shares nothing with the main pipeline but the input files.
# ==============================================================================

library(targets)

tar_option_set(
  packages = c("tidyverse", "here", "readxl", "patchwork",
               "vegan", "ade4", "compositions", "RVAideMemoire",
               "MatrixCorrelation", "conflicted")
)

suppressMessages(
  conflicted::conflicts_prefer(
    dplyr::filter, dplyr::select, dplyr::lag,
    stats::sd, stats::var, stats::dist, stats::cor, stats::cov,
    base::scale, base::norm, .quiet = TRUE
  )
)

rob <- function(...) here::here("analysis/robustness", ...)
spectra_files <- function(sub)
  list.files(rob(sub, "spectra"), pattern = "\\.csv$", full.names = TRUE)

list(

  # ---- committed inputs (products of the main pipeline + raw metadata) -------
  tar_target(core_metric_xlsx,
             here::here("analysis/data/raw_data/SDG_core_metric.xlsx"),
             format = "file"),
  tar_target(derived_direction_csv,
             here::here("analysis/data/derived_data/SPHARM_direction.csv"),
             format = "file"),
  tar_target(derived_morphology_csv,
             here::here("analysis/data/derived_data/SPHARM_morphology.csv"),
             format = "file"),

  # ---- bandwidth (h) sensitivity ---------------------------------------------
  tar_target(bandwidth_spectra_csvs,
             spectra_files("bandwidth_sensitivity"), format = "file"),
  tar_target(robustness_bandwidth, {
    force(core_metric_xlsx); force(derived_direction_csv); force(derived_morphology_csv)
    force(bandwidth_spectra_csvs)
    local({
      source(rob("bandwidth_sensitivity/bandwidth_sensitivity_stats.R"), local = TRUE)
      list(metrics = metrics_df, order_long = order_long_df,
           fig_orderselection = p_cum / p_cv,
           fig_im_heatmaps = p_heat)
    })
  }),

  # ---- scar minimum-size threshold sensitivity -------------------------------
  tar_target(scar_spectra_csvs,
             spectra_files("scar_threshold_sensitivity"), format = "file"),
  tar_target(scar_attrition_csvs,
             c(rob("scar_threshold_sensitivity/scar_attrition_by_specimen.csv"),
               rob("scar_threshold_sensitivity/scar_attrition_summary.csv")),
             format = "file"),
  tar_target(robustness_scar_threshold, {
    force(core_metric_xlsx); force(derived_direction_csv); force(derived_morphology_csv)
    force(scar_spectra_csvs); force(scar_attrition_csvs)
    local({
      source(rob("scar_threshold_sensitivity/scar_threshold_sensitivity_stats.R"), local = TRUE)
      list(metrics = metrics_df, order_long = order_long_df,
           fig_orderselection = p_cum / p_cv,
           fig_scarcounts = p_sc)
    })
  }),

  # ---- mesh decimation + smoothing sensitivity -------------------------------
  tar_target(mesh_spectra_csvs,
             spectra_files("mesh_decimation_sensitivity"), format = "file"),
  tar_target(robustness_mesh_decimation, {
    force(core_metric_xlsx); force(derived_direction_csv); force(derived_morphology_csv)
    force(mesh_spectra_csvs)
    local({
      source(rob("mesh_decimation_sensitivity/mesh_sensitivity_stats.R"), local = TRUE)
      list(metrics = metrics_df, order_long = order_long_df,
           fig_orderselection = p_cum / p_cv)
    })
  }),

  # ---- CoIA null-result power / sensitivity (RV2 + bootstrap CI + MDES) -------
  # Reuses ONLY the committed main outputs (no per-setting spectra): it rebuilds
  # the exact EXP/SDG CoIA ILR matrices and asks what coupling the null
  # RV / Mantel result can rule out, for both assemblages.
  tar_target(robustness_power, {
    force(core_metric_xlsx); force(derived_direction_csv); force(derived_morphology_csv)
    local({
      source(rob("coia_power_sensitivity/coia_power_sensitivity_stats.R"), local = TRUE)
      list(metrics = metrics_df, power_long = power_long_df,
           fig_power = p_power)
    })
  })
)
