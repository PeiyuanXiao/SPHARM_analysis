# ==============================================================================
# _targets.R
# Fully reproducible analysis pipeline: R preprocessing -> Python SPHARM -> R statistics
#
# Reproduction steps:
#   1. docker run ...              # start the container
#   2. targets::tar_make()         # run the whole pipeline
#   3. quarto render paper.qmd     # render the paper
#
# Note:
#   Python targets are executed as separate Python processes instead of via
#   reticulate. This avoids shared-library conflicts between rocker/geospatial's
#   system geospatial stack and conda packages such as PIL/VTK/libtiff.
# ==============================================================================

library(targets)
library(tarchetypes)

# ==============================================================================
# Global options
# ==============================================================================

tar_option_set(
  packages = c("tidyverse", "patchwork", "here",
               "readxl", "vegan", "FSA", "RVAideMemoire",
               "ggrepel", "conflicted",
               "linkET", "compositions", "ade4", "circular",
               "rsvg", "png", "grid")
)

conflicted::conflicts_prefer(dplyr::filter)
conflicted::conflicts_prefer(dplyr::select)
conflicted::conflicts_prefer(dplyr::lag)
conflicted::conflicts_prefer(stats::sd)
conflicted::conflicts_prefer(stats::var)
conflicted::conflicts_prefer(stats::dist)
conflicted::conflicts_prefer(stats::cor)
conflicted::conflicts_prefer(stats::cov)
conflicted::conflicts_prefer(stats::anova)
conflicted::conflicts_prefer(graphics::segments)
conflicted::conflicts_prefer(base::norm)
conflicted::conflicts_prefer(base::scale)
conflicted::conflicts_prefer(base::`%*%`)
conflicted::conflicts_prefer(ggplot2::aes)
conflicted::conflicts_prefer(ggplot2::theme_bw)
conflicted::conflicts_prefer(ggplot2::ggsave)
conflicted::conflicts_prefer(ggplot2::annotate)
conflicted::conflicts_prefer(ggplot2::theme_classic)
conflicted::conflicts_prefer(ggplot2::theme_minimal)
conflicted::conflicts_prefer(ggplot2::theme_void)
conflicted::conflicts_prefer(ggplot2::theme_grey)
conflicted::conflicts_prefer(ggplot2::theme_gray)
conflicted::conflicts_prefer(ggplot2::theme_light)
conflicted::conflicts_prefer(ggplot2::theme_dark)
conflicted::conflicts_prefer(ggplot2::theme_linedraw)

# ==============================================================================
# Python environment (Conda env inside Docker)
# ==============================================================================

PYTHON_BIN <- "/opt/conda/envs/spharm/bin/python"

check_python <- function() {
  if (!file.exists(PYTHON_BIN)) {
    stop(
      "Python environment not found: ", PYTHON_BIN, "\n",
      "Python targets must run inside the Docker container.\n",
      "To run only the R analysis targets, use:\n",
      "  tar_make(spharm_analysis)  # provided the CSVs already exist",
      call. = FALSE
    )
  }
}

python_env <- function() {
  old_ld <- Sys.getenv("LD_LIBRARY_PATH", unset = "")
  ld_paths <- c("/opt/conda/envs/spharm/lib", "/opt/conda/lib")
  if (nzchar(old_ld)) {
    ld_paths <- c(ld_paths, old_ld)
  }
  
  c(
    "MPLBACKEND=Agg",
    "PYTHONPATH=/project/analysis/scripts/SPHARM_main:/project/analysis/scripts",
    paste0("LD_LIBRARY_PATH=", paste(ld_paths, collapse = ":"))
  )
}

run_python_code <- function(code) {
  check_python()
  
  script <- tempfile(fileext = ".py")
  on.exit(unlink(script), add = TRUE)
  
  writeLines(code, script)
  
  status <- system2(
    PYTHON_BIN,
    args = script,
    env = python_env(),
    stdout = "",
    stderr = ""
  )
  
  if (!identical(status, 0L)) {
    stop("Python script failed with exit status: ", status, call. = FALSE)
  }
}

check_output_file <- function(path) {
  if (!file.exists(path)) {
    stop("Expected output file was not created: ", path, call. = FALSE)
  }
  path
}

# ==============================================================================
# Python invocation helpers
# ==============================================================================

#' Run batch_process() from SPHARM_main.py
#' @param input_dir  directory of STL files
#' @param output_dir output directory
#' @return output CSV path (tracked via format = "file")
run_spharm_morphology <- function(input_dir, output_dir) {
  code <- sprintf("
from SPHARM_main import batch_process
batch_process(%s, %s)
", shQuote(input_dir), shQuote(output_dir))
  
  run_python_code(code)
  
  check_output_file(file.path(output_dir, "SPHARM_morphology.csv"))
}

#' Run run_pipeline() from kde_to_spharm_main.py
#' @param source     alignment method: "svd"
#' @param validation whether this is validation mode
#' @return output CSV path
run_spharm_direction <- function(source = "svd", validation = FALSE) {
  val_str <- ifelse(validation, "True", "False")
  
  code <- sprintf("
from kde_to_spharm_main import run_pipeline
run_pipeline(%s, validation=%s)
", shQuote(source), val_str)
  
  run_python_code(code)
  
  if (validation) {
    check_output_file(file.path(
      "/project/analysis/data/derived_data/validation",
      source,
      "SPHARM_direction.csv"
    ))
  } else {
    check_output_file("/project/analysis/data/derived_data/SPHARM_direction.csv")
  }
}

#' Run rotate_svd_directions.py
#' @return output CSV path
run_rotate_svd <- function() {
  check_python()
  
  status <- system2(
    PYTHON_BIN,
    args = c("/project/analysis/scripts/SPHARM_modules/rotate_svd_directions.py"),
    env = python_env(),
    stdout = "",
    stderr = ""
  )
  
  if (!identical(status, 0L)) {
    stop("rotate_svd_directions.py failed with exit status: ", status, call. = FALSE)
  }
  
  check_output_file("/project/analysis/data/derived_data/directions_aligned_svd_rotated.csv")
}

# ==============================================================================
# Output directories
# These live under the git-ignored analysis/output/, so they are absent on a
# fresh clone or CI checkout. ggsave()/saveWidget() do not create parent dirs,
# so create them up front (no-op if they already exist).
# ==============================================================================
for (.out_dir in c("analysis/data/derived_data",
                   "analysis/output/figures",
                   "analysis/output/html")) {
  dir.create(here::here(.out_dir), recursive = TRUE, showWarnings = FALSE)
}

# ==============================================================================
# Target list
# ==============================================================================

list(
  
  # ============================================================================
  # Stage 1: R preprocessing — direction-vector alignment
  # ============================================================================
  
  # align_svd.R: raw scar data -> aligned direction vectors
  # produces: directions_raw.csv + directions_aligned_svd.csv
  tar_target(
    align_svd_csvs,
    local({
      cli::cli_alert_info("Stage 1 - Align scar vectors (SVD)")
      source(here::here("analysis/scripts/r_alignment/align_svd.R"),
             local = TRUE)
      c(here::here("analysis/data/derived_data/directions_raw.csv"),
        here::here("analysis/data/derived_data/directions_aligned_svd.csv"))
    }),
    format = "file",
    description = "Stage 1 - Align scar vectors (SVD)"
  ),
  
  # align_lin2024.R: Lin 2024 alignment (validation pipeline only)
  # produces: directions_aligned_lin2024.csv
  tar_target(
    align_lin2024_csv,
    local({
      cli::cli_alert_info("Stage 1 - Align scar vectors (Lin 2024, validation)")
      source(here::here("analysis/scripts/r_alignment/align_lin2024.R"),
             local = TRUE)
      here::here("analysis/data/derived_data/directions_aligned_lin2024.csv")
    }),
    format = "file",
    description = "Stage 1 - Align scar vectors (Lin 2024, validation)"
  ),
  
  # ============================================================================
  # Stage 2: Python SPHARM analysis
  # ============================================================================
  
  # --- morphology SPHARM: STL -> power spectrum ---
  tar_target(
    spharm_morphology_csv,
    {
      cli::cli_alert_info("Stage 2 - M-SPHARM: core meshes -> power spectra")
      run_spharm_morphology(
        "/project/analysis/data/3D_models_cores",
        "/project/analysis/data/derived_data"
      )
    },
    format = "file",
    description = "Stage 2 - M-SPHARM: core meshes -> power spectra"
  ),
  
  # --- direction SPHARM (production): SVD-aligned vectors -> KDE -> power spectrum ---
  tar_target(
    spharm_direction_csv,
    {
      cli::cli_alert_info("Stage 2 - SP-SPHARM: spherical KDE -> power spectra")
      force(align_svd_csvs)
      run_spharm_direction(source = "svd", validation = FALSE)
    },
    format = "file",
    description = "Stage 2 - SP-SPHARM: spherical KDE -> power spectra"
  ),
  
  # --- generate rotation-perturbed data (validation) ---
  tar_target(
    rotate_svd_csv,
    {
      cli::cli_alert_info("Stage 2 - Rotation-perturbed vectors (invariance test)")
      force(align_svd_csvs)
      run_rotate_svd()
    },
    format = "file",
    description = "Stage 2 - Rotation-perturbed vectors (invariance test)"
  ),
  
  # --- direction SPHARM (validation): four alignments x KDE -> power spectrum ---
  tar_target(
    spharm_direction_validation_csvs,
    {
      cli::cli_alert_info("Stage 2 - SP-SPHARM validation (4 alignments)")
      force(align_svd_csvs)
      force(align_lin2024_csv)
      force(rotate_svd_csv)
      
      csvs <- vapply(
        c("raw", "svd", "lin2024", "svd_rotated"),
        function(src) run_spharm_direction(source = src, validation = TRUE),
        character(1)
      )
      unname(csvs)
    },
    format = "file",
    description = "Stage 2 - SP-SPHARM validation (4 alignments)"
  ),
  
  # ============================================================================
  # Stage 3: R statistical analysis
  # ============================================================================
  
  tar_target(
    spharm_analysis,
    {
      cli::cli_alert_info("Stage 3 - Method comparison + PERMANOVA (EXP cores)")
      force(spharm_morphology_csv)
      force(spharm_direction_csv)
      force(align_svd_csvs)
      
      local({
        source(here::here("analysis/scripts/r_spharm/spharm_analysis.R"),
               local = TRUE)
        list(
          spharm_posthoc = list(
            p_dir_disc_bi  = perm_dir$pairwise$p.value["Discoid",       "Bidirectional"],
            p_fab_uni_bi   = perm_EI$pairwise$p.value["Unidirectional", "Bidirectional"],
            p_fab_lev_disc = perm_EI$pairwise$p.value["Levallois",      "Discoid"],
            p_fab_disc_bi  = perm_EI$pairwise$p.value["Discoid",        "Bidirectional"]
          ),
          perm_morph_r2         = round(perm_morph$global$R2[1], 3),
          perm_morph_f          = round(perm_morph$global$`F`[1], 3),
          perm_morph_p          = round(perm_morph$global$`Pr(>F)`[1], 3),
          p_exp_method_combined = exp_method_compare_combined
        )
      })
    },
    description = "Stage 3 - Method comparison + PERMANOVA (EXP cores)"
  ),
  
  tar_target(
    p_rotational_invariance_validity,
    {
      cli::cli_alert_info("Stage 3 - Rotational-invariance test")
      force(spharm_direction_validation_csvs)
      force(align_svd_csvs)
      force(align_lin2024_csv)
      
      local({
        source(here::here("analysis/scripts/r_validation/validate_rotation_all.R"),
               local = TRUE)
        p_rotational_invariance_validity
      })
    },
    description = "Stage 3 - Rotational-invariance test"
  ),
  
  tar_target(
    im_comparison,
    {
      force(spharm_direction_csv)
      
      local({
        source(here::here("analysis/scripts/r_validation/methods_comparison_IM.R"),
               local = TRUE)
        list(
          dist_spi    = avg_dist$mean_dist[avg_dist$method == "SPI"],
          dist_fabric = avg_dist$mean_dist[avg_dist$method == "Fabric"],
          dist_spharm = avg_dist$mean_dist[avg_dist$method == "SPHARM"],
          p_heatmap   = p_heatmap
        )
      })
    },
    description = "Stage 3 - SPI / fabric / SPHARM comparison (IM)"
  ),
  
  tar_target(
    exp_cia_analysis,
    {
      force(spharm_morphology_csv)
      force(spharm_direction_csv)
      force(align_svd_csvs)
      force(spharm_analysis)
      
      local({
        source(here::here("analysis/scripts/r_statistics/exp_cores_statistics.R"),
               local = TRUE)
        list(
          p_final          = p_final,
          mantel_global_r  = round(mantel_global$statistic, 3),
          mantel_global_p  = round(mantel_global$signif, 3),
          rv               = round(coin_exp$RV, 3),
          rv_p             = round(rv_test$pvalue, 3),
          mantel_discoid_r = round(mantel_by_typology$mantel_r[mantel_by_typology$Typology == "Discoid"], 3),
          mantel_discoid_p = round(mantel_by_typology$p_value[mantel_by_typology$Typology == "Discoid"], 3),
          mantel_lev_r     = round(mantel_by_typology$mantel_r[mantel_by_typology$Typology == "Levallois"], 3),
          mantel_uni_r     = round(mantel_by_typology$mantel_r[mantel_by_typology$Typology == "Unidirectional"], 3),
          mantel_multi_r   = round(mantel_by_typology$mantel_r[mantel_by_typology$Typology == "Multiplatform"], 3),
          mantel_bi_r      = round(mantel_by_typology$mantel_r[mantel_by_typology$Typology == "Bidirectional"], 3),
          rayleigh_bi_p    = round(res_circ_typology$rayleigh$rayleigh_p[res_circ_typology$rayleigh$group == "Bidirectional"], 3),
          rayleigh_disc_p  = round(res_circ_typology$rayleigh$rayleigh_p[res_circ_typology$rayleigh$group == "Discoid"], 3),
          rayleigh_lev_p   = round(res_circ_typology$rayleigh$rayleigh_p[res_circ_typology$rayleigh$group == "Levallois"], 3),
          mean_dir_bi      = round(res_circ_typology$desc$mean_dir_deg[res_circ_typology$desc$group == "Bidirectional"], 0),
          mean_dir_disc    = round(res_circ_typology$desc$mean_dir_deg[res_circ_typology$desc$group == "Discoid"], 0),
          mean_dir_lev     = round(res_circ_typology$desc$mean_dir_deg[res_circ_typology$desc$group == "Levallois"], 0),
          kw_arrow_chi2    = round(res_len_typology$kw$statistic, 2),
          kw_arrow_p       = round(res_len_typology$kw$p.value, 3)
        )
      })
    },
    description = "Stage 3 - Experimental cores: CIA / Mantel / circular"
  ),
  
  tar_target(
    sdg_cia_analysis,
    {
      force(spharm_morphology_csv)
      force(spharm_direction_csv)
      force(spharm_analysis)
      
      local({
        source(here::here("analysis/scripts/r_statistics/SDG_cores_statistics.R"),
               local = TRUE)
        list(
          mantel_global_r          = round(mantel_global$statistic, 3),
          mantel_global_p          = round(mantel_global$signif, 3),
          rv                       = round(coin_arch$RV, 3),
          rv_p                     = round(rv_test$pvalue, 3),
          cia_ax1_pct              = round(cia_inertia[1], 1),
          cia_ax2_pct              = round(cia_inertia[2], 1),
          cia_total_pct            = round(sum(cia_inertia[1:2]), 1),
          kw_arrow_layer_chi2      = round(res_len_layer$kw$statistic, 2),
          kw_arrow_layer_p         = round(res_len_layer$kw$p.value, 3),
          kw_arrow_rawmat_chi2     = round(res_len_rawmat$kw$statistic, 2),
          kw_arrow_rawmat_p        = round(res_len_rawmat$kw$p.value, 3),
          kw_arrow_type_chi2       = round(res_len_coretype$kw$statistic, 2),
          kw_arrow_type_p          = round(res_len_coretype$kw$p.value, 3),
          rayleigh_uni_U           = round(
            res_circ_coretype$rayleigh$rayleigh_U[
              res_circ_coretype$rayleigh$group == "Unifacial_unidirection"], 3),
          rayleigh_uni_p           = round(
            res_circ_coretype$rayleigh$rayleigh_p[
              res_circ_coretype$rayleigh$group == "Unifacial_unidirection"], 3),
          mean_dir_uni             = round(
            res_circ_coretype$desc$mean_dir_deg[
              res_circ_coretype$desc$group == "Unifacial_unidirection"], 0),
          perm_morph_type_r2       = round(permanova_results$R2[
            permanova_results$domain == "Morphology" &
              permanova_results$grouping == "Core Type"], 3),
          perm_morph_type_f        = round(permanova_results$F_value[
            permanova_results$domain == "Morphology" &
              permanova_results$grouping == "Core Type"], 2),
          perm_morph_type_p        = round(permanova_results$p_value[
            permanova_results$domain == "Morphology" &
              permanova_results$grouping == "Core Type"], 3),
          perm_scar_type_r2        = round(permanova_results$R2[
            permanova_results$domain == "Scar Direction" &
              permanova_results$grouping == "Core Type"], 3),
          perm_scar_type_f         = round(permanova_results$F_value[
            permanova_results$domain == "Scar Direction" &
              permanova_results$grouping == "Core Type"], 2),
          perm_scar_type_p         = round(permanova_results$p_value[
            permanova_results$domain == "Scar Direction" &
              permanova_results$grouping == "Core Type"], 3),
          perm_morph_rawmat_r2     = round(permanova_results$R2[
            permanova_results$domain == "Morphology" &
              permanova_results$grouping == "Raw Material"], 3),
          perm_morph_rawmat_f      = round(permanova_results$F_value[
            permanova_results$domain == "Morphology" &
              permanova_results$grouping == "Raw Material"], 2),
          perm_morph_rawmat_p      = round(permanova_results$p_value[
            permanova_results$domain == "Morphology" &
              permanova_results$grouping == "Raw Material"], 3),
          perm_rawmat_pairwise_p   = round(pairwise_results$p_holm[
            pairwise_results$domain == "Morphology" &
              pairwise_results$grouping == "Raw Material" &
              pairwise_results$group1 == "chert" &
              pairwise_results$group2 == "sandstone"], 3),
          perm_scar_rawmat_r2      = round(permanova_results$R2[
            permanova_results$domain == "Scar Direction" &
              permanova_results$grouping == "Raw Material"], 3),
          perm_scar_rawmat_f       = round(permanova_results$F_value[
            permanova_results$domain == "Scar Direction" &
              permanova_results$grouping == "Raw Material"], 3),
          perm_scar_rawmat_p       = round(permanova_results$p_value[
            permanova_results$domain == "Scar Direction" &
              permanova_results$grouping == "Raw Material"], 3),
          perm_morph_layer_r2      = round(permanova_results$R2[
            permanova_results$domain == "Morphology" &
              permanova_results$grouping == "Layer"], 3),
          perm_morph_layer_f       = round(permanova_results$F_value[
            permanova_results$domain == "Morphology" &
              permanova_results$grouping == "Layer"], 3),
          perm_morph_layer_p       = round(permanova_results$p_value[
            permanova_results$domain == "Morphology" &
              permanova_results$grouping == "Layer"], 3),
          perm_scar_layer_r2       = round(permanova_results$R2[
            permanova_results$domain == "Scar Direction" &
              permanova_results$grouping == "Layer"], 3),
          perm_scar_layer_f        = round(permanova_results$F_value[
            permanova_results$domain == "Scar Direction" &
              permanova_results$grouping == "Layer"], 3),
          perm_scar_layer_p        = round(permanova_results$p_value[
            permanova_results$domain == "Scar Direction" &
              permanova_results$grouping == "Layer"], 3),
          fig_coia_composite = p_final
        )
      })
    },
    description = "Stage 3 - SDG cores: CIA / PERMANOVA / circular"
  ),

  # degree-selection diagnostic: per-degree across-specimen CV and cumulative
  # power for EXP and SDG (morphology & direction), used to justify the l
  # truncation. Returns the composite figure inserted into the manuscript.
  tar_target(
    degree_selection_diagnostic,
    {
      force(spharm_morphology_csv)
      force(spharm_direction_csv)

      local({
        source(here::here("analysis/scripts/r_spharm/power_degree_selection.R"),
               local = TRUE)
        p_degree_selection
      })
    },
    description = "Stage 3 - Degree-selection diagnostic (power-order truncation)"
  )
)
