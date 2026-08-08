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

  # ---- direction-landmark annotation error -----------------------------------
  # READ-ONLY. The sweep itself was run separately by
  # annotation_perturbation/perturb_spharm.py (which re-runs the whole production
  # chain per replicate) and annotation_perturbation.R; this target only loads the
  # committed result tables, counts the scars of the analysed EXP sample from the
  # sweep's own input file, and derives the two interpolated crossing levels
  # quoted in the text. No statistic is recomputed here, and the figure is the
  # committed PNG rather than a re-drawn object.
  tar_target(perturbation_csvs,
             c(rob("annotation_perturbation/annotation_perturbation_summary.csv"),
               rob("annotation_perturbation/annotation_perturbation_polarity_three_methods.csv"),
               rob("annotation_perturbation/annotation_perturbation_angle_three_methods.csv")),
             format = "file"),
  tar_target(perturbation_figure_png,
             rob("annotation_perturbation/figures/fig_S_perturbation_combined.png"),
             format = "file"),
  tar_target(derived_directions_svd_csv,
             here::here("analysis/data/derived_data/directions_aligned_svd.csv"),
             format = "file"),
  tar_target(robustness_annotation, {
    rd <- function(pat)
      readr::read_csv(grep(pat, perturbation_csvs, value = TRUE),
                      show_col_types = FALSE)
    sm_df    <- rd("_summary\\.csv$")
    three_df <- dplyr::bind_rows(rd("_polarity_three_methods\\.csv$"),
                                 rd("_angle_three_methods\\.csv$"))

    # Scar counts of the analysed EXP sample, read from the sweep's own input
    # (directions_aligned_svd.csv) under the same filters the sweep applies: EXP
    # specimens with complete unit vectors, Biface excluded. `flips` is the number
    # of vectors actually reversed at each flip fraction, which the sweep draws
    # per specimen as round(f * n_i) — so the smallest fractions touch only the
    # specimens carrying enough scars for that product to reach one.
    dirs <- readr::read_csv(derived_directions_svd_csv, show_col_types = FALSE)
    dirs <- dirs[grepl("^EXP", dirs$ID) & dirs$Typology != "Biface" &
                 !is.na(dirs$ux) & !is.na(dirs$uy) & !is.na(dirs$uz), ]
    per_spec <- as.integer(table(dirs$ID))
    lv <- sort(unique(sm_df$level[sm_df$perturbation == "polarity"]))
    scars <- list(
      n_spec  = length(per_spec),
      n_scar  = sum(per_spec),
      median  = stats::median(per_spec),
      min     = min(per_spec),
      max     = max(per_spec),
      flips   = stats::setNames(vapply(lv, function(f) sum(round(f * per_spec)),
                                       numeric(1)), lv),
      n_touch = stats::setNames(vapply(lv, function(f) sum(round(f * per_spec) > 0),
                                       numeric(1)), lv))

    # Flip fraction at which the SP-SPHARM curve passes through the (constant)
    # fabric curve, by linear interpolation between the swept levels. Reported
    # for both outcome measures because they cross at different fractions.
    cross <- function(col) {
      d   <- three_df[three_df$perturbation == "polarity", ]
      fab <- unique(d[[col]][d$method == "Fabric (E, I)"])
      sp  <- d[d$method == "SP-SPHARM", ]
      sp  <- sp[order(sp$level), ]
      y   <- sp[[col]] - fab
      i   <- which(y[-1] * y[-length(y)] < 0)[1]
      if (is.na(i)) return(NA_real_)
      sp$level[i] + (sp$level[i + 1] - sp$level[i]) * y[i] / (y[i] - y[i + 1])
    }

    list(summary = sm_df, three = three_df, scars = scars,
         cross_count = cross("n_sig_05_med"),
         cross_med   = cross("med_neglog10_med"),
         fig_png     = perturbation_figure_png)
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

  # ---- truncation degree sensitivity -----------------------------------------
  # READ-ONLY. The sweep itself was run separately by
  # truncation_sensitivity/truncation_sensitivity.R; this target only loads the
  # committed result tables and derives the stable truncation intervals quoted
  # in the text. No statistic is recomputed here.
  tar_target(truncation_csvs,
             c(rob("truncation_sensitivity/truncation_sensitivity_SP.csv"),
               rob("truncation_sensitivity/truncation_sensitivity_M.csv"),
               rob("truncation_sensitivity/truncation_sensitivity_summary.csv")),
             format = "file"),
  tar_target(robustness_truncation, {
    rd <- function(pat)
      readr::read_csv(grep(pat, truncation_csvs, value = TRUE),
                      show_col_types = FALSE)
    sp_df <- rd("_SP\\.csv$")
    m_df  <- rd("_M\\.csv$")
    sm_df <- rd("_summary\\.csv$")

    # Stability rule, applied per scan x assemblage relative to that scan's
    # anchor: a truncation counts as stable when it (i) resolves at least as
    # many Holm-corrected pairwise contrasts as the anchor, (ii) keeps the
    # grouping PERMANOVA significant wherever the anchor is significant, and
    # (iii) preserves the decoupling result (RV p >= 0.05). The reported
    # interval is the contiguous run of stable settings containing the anchor.
    stable_run <- function(d) {
      d    <- d[order(d$varied_lmax), ]
      pre  <- if (d$scan[1] == "SP") "SP_" else "M_"
      nsig <- d[[paste0(pre, "nsig05")]]
      pval <- d[[paste0(pre, "p")]]
      a    <- which(d$is_anchor)[1]
      ok   <- nsig >= nsig[a] &
              (pval < 0.05 | pval[a] >= 0.05) &
              d$RV_p >= 0.05
      lo <- a; while (lo > 1       && isTRUE(ok[lo - 1])) lo <- lo - 1
      hi <- a; while (hi < nrow(d) && isTRUE(ok[hi + 1])) hi <- hi + 1
      data.frame(scan = d$scan[1], assemblage = d$assemblage[1],
                 lo = d$varied_lmax[lo], hi = d$varied_lmax[hi],
                 anchor = d$varied_lmax[a])
    }
    stable_df <- do.call(rbind, lapply(
      split(sm_df, list(sm_df$scan, sm_df$assemblage), drop = TRUE), stable_run))

    list(sp = sp_df, m = m_df, summary = sm_df, stable = stable_df)
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
