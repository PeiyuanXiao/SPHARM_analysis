# ==============================================================================
# _targets.R
# 完整可复现分析流水线：R 预处理 → Python SPHARM → R 统计
#
# 审稿人复现步骤：
#   1. docker run ...              # 启动容器
#   2. targets::tar_make()         # 运行全部流水线
#   3. quarto render paper.qmd     # 渲染论文
#
# DAG 结构：
#   ┌─ align_svd_csvs ─────────────────────────────────────────────────────┐
#   │  ├─ spharm_direction_csv (Python: KDE → SPHARM)                     │
#   │  │    └─ spharm_analysis → exp_cia_analysis                         │
#   │  │                       → sdg_cia_analysis                         │
#   │  │                       → im_comparison                            │
#   │  ├─ align_lin2024_csv ──┐                                           │
#   │  ├─ rotate_svd_csv ────┐│                                           │
#   │  │                     ││                                           │
#   │  │  spharm_direction_validation_csvs (Python: --source all) ◄───┘   │
#   │  │    └─ p_rotational_invariance_validity                           │
#   │  │                                                                  │
#   ├─ spharm_morphology_csv (Python: STL → SPHARM)                      │
#   │    └─ spharm_analysis (同上)                                        │
#   └──────────────────────────────────────────────────────────────────────┘
# ==============================================================================

library(targets)
library(tarchetypes)

# ==============================================================================
# 全局选项
# ==============================================================================

tar_option_set(
  packages = c("tidyverse", "patchwork", "here",
               "readxl", "vegan", "FSA", "RVAideMemoire",
               "ggrepel", "conflicted",
               "linkET", "compositions", "ade4", "circular",
               "rsvg", "png", "grid",
               "reticulate")
)

# 统一声明所有包冲突
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

# ==============================================================================
# Python 环境配置（Docker 内 Conda 环境）
# ==============================================================================

PYTHON_BIN <- "/opt/conda/envs/spharm/bin/python"

#' 延迟初始化 Python 环境：仅在首次调用时执行 use_python()，
#' 且仅当 Docker 内的 Python 路径存在时才设置。
#' 在 Docker 外运行纯 R target 时不会报错。
ensure_python <- function() {
  if (!file.exists(PYTHON_BIN)) {
    stop(
      "Python 环境未找到：", PYTHON_BIN, "\n",
      "Python targets 需要在 Docker 容器内运行。\n",
      "如果你只需要跑 R 分析 target，可以用：\n",
      "  tar_make(spharm_analysis)  # 前提是 CSV 已存在",
      call. = FALSE
    )
  }
  reticulate::use_python(PYTHON_BIN, required = TRUE)
}

# ==============================================================================
# Python 调用辅助函数
# ==============================================================================

#' 运行 SPHARM_main.py 的 batch_process()
#' @param input_dir  STL 文件目录
#' @param output_dir 输出目录
#' @return 输出 CSV 路径（供 format = "file" 追踪）
run_spharm_morphology <- function(input_dir, output_dir) {
  ensure_python()
  reticulate::py_run_string(sprintf("
import sys, matplotlib
matplotlib.use('Agg')
sys.path.insert(0, '/project/analysis/scripts/SPHARM_main')
sys.path.insert(0, '/project/analysis/scripts')
from SPHARM_main import batch_process
batch_process('%s', '%s')
", input_dir, output_dir))
  file.path(output_dir, "SPHARM_morphology.csv")
}

#' 运行 kde_to_spharm_main.py 的 run_pipeline()
#' @param source     对齐方式："svd"
#' @param validation 是否为验证模式
#' @return 输出 CSV 路径
run_spharm_direction <- function(source = "svd", validation = FALSE) {
  ensure_python()
  val_str <- ifelse(validation, "True", "False")
  reticulate::py_run_string(sprintf("
import sys, matplotlib
matplotlib.use('Agg')
sys.path.insert(0, '/project/analysis/scripts/SPHARM_main')
sys.path.insert(0, '/project/analysis/scripts')
from kde_to_spharm_main import run_pipeline
run_pipeline('%s', validation=%s)
", source, val_str))
  
  if (validation) {
    file.path("/project/analysis/data/derived_data/validation",
              source, "SPHARM_direction.csv")
  } else {
    "/project/analysis/data/derived_data/SPHARM_direction.csv"
  }
}

#' 运行 rotate_svd_directions.py
#' @return 输出 CSV 路径
run_rotate_svd <- function() {
  ensure_python()
  reticulate::py_run_file(
    "/project/analysis/scripts/SPHARM_main/rotate_svd_directions.py"
  )
  "/project/analysis/data/derived_data/directions_aligned_svd_rotated.csv"
}

# ==============================================================================
# targets 列表
# ==============================================================================

list(
  
  # ============================================================================
  # 第一层：R 端预处理 — 方向向量对齐
  # ============================================================================
  
  # align_svd.R：原始刮痕数据 → 对齐后方向向量
  # 产出：directions_raw.csv + directions_aligned_svd.csv
  tar_target(
    align_svd_csvs,
    local({
      source(here::here("analysis/scripts/r_spharm/align_svd.R"),
             local = TRUE)
      c(here::here("analysis/data/derived_data/directions_raw.csv"),
        here::here("analysis/data/derived_data/directions_aligned_svd.csv"))
    }),
    format = "file"
  ),
  
  # align_lin2024.R：Lin 2024 法对齐（仅验证流水线需要）
  # 产出：directions_aligned_lin2024.csv
  tar_target(
    align_lin2024_csv,
    local({
      source(here::here("analysis/scripts/r_spharm/align_lin2024.R"),
             local = TRUE)
      here::here("analysis/data/derived_data/directions_aligned_lin2024.csv")
    }),
    format = "file"
  ),
  
  # ============================================================================
  # 第二层：Python SPHARM 分析
  # ============================================================================
  
  # --- 形态 SPHARM：STL → 功率谱 ---
  tar_target(
    spharm_morphology_csv,
    run_spharm_morphology(
      "/project/analysis/data/3D_models_cores",
      "/project/analysis/data/derived_data"
    ),
    format = "file"
  ),
  
  # --- 方向 SPHARM（生产模式）：SVD 对齐方向向量 → KDE → 功率谱 ---
  tar_target(
    spharm_direction_csv,
    {
      # 显式依赖：确保 align_svd.R 已执行
      force(align_svd_csvs)
      run_spharm_direction(source = "svd", validation = FALSE)
    },
    format = "file"
  ),
  
  # --- 旋转扰动数据生成（验证用）---
  tar_target(
    rotate_svd_csv,
    {
      force(align_svd_csvs)
      run_rotate_svd()
    },
    format = "file"
  ),
  
  # --- 方向 SPHARM（验证模式）：四种对齐 × KDE → 功率谱 ---
  # 产出四份 CSV：validation/{raw,svd,lin2024,svd_rotated}/SPHARM_direction.csv
  tar_target(
    spharm_direction_validation_csvs,
    {
      # 显式依赖：确保所有上游数据就绪
      force(align_svd_csvs)
      force(align_lin2024_csv)
      force(rotate_svd_csv)
      
      # 逐一运行四种 source
      csvs <- vapply(
        c("raw", "svd", "lin2024", "svd_rotated"),
        function(src) run_spharm_direction(source = src, validation = TRUE),
        character(1)
      )
      unname(csvs)
    },
    format = "file"
  ),
  
  # ============================================================================
  # 第三层：R 统计分析（现有 targets，加上显式依赖声明）
  # ============================================================================
  
  tar_target(
    spharm_analysis,
    {
      # 显式依赖 Python 产出的 CSV
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
          perm_morph_p          = round(perm_morph$global$`Pr(>F)`[1], 3),
          p_exp_method_combined = exp_method_compare_combined
        )
      })
    }
  ),
  
  tar_target(
    p_rotational_invariance_validity,
    {
      # 显式依赖：验证模式的 Python 输出 + 对齐数据
      force(spharm_direction_validation_csvs)
      force(align_svd_csvs)
      force(align_lin2024_csv)
      
      local({
        source(here::here("analysis/scripts/r_validation/validate_rotation_all.R"),
               local = TRUE)
        p_rotational_invariance_validity
      })
    }
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
    }
  ),
  
  tar_target(
    exp_cia_analysis,
    {
      force(spharm_morphology_csv)
      force(spharm_direction_csv)
      force(align_svd_csvs)
      
      local({
        source(here::here("analysis/scripts/r_statistics/exp_cores_statistics.R"),
               local = TRUE)
        list(
          p_final          = p_final,
          p_se_combined    = p_se_combined,
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
          kw_se_dir_chi2   = round(res_se_dir$kw$statistic, 2),
          kw_se_dir_p      = round(res_se_dir$kw$p.value, 3),
          p_se_lev_multi   = round(res_se_dir$dunn$p.adj[res_se_dir$dunn$Comparison == "Levallois - Multiplatform"], 3),
          kw_se_morph_chi2 = round(res_se_morph$kw$statistic, 2),
          kw_se_morph_p    = round(res_se_morph$kw$p.value, 3),
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
    }
  ),
  
  tar_target(
    sdg_cia_analysis,
    {
      force(spharm_morphology_csv)
      force(spharm_direction_csv)
      
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
          kw_se_dir_type_chi2      = round(res_se_dir_coretype$kw$statistic, 3),
          kw_se_dir_type_p         = round(res_se_dir_coretype$kw$p.value, 3),
          dunn_se_multi_uni_p      = round(
            res_se_dir_coretype$dunn$p.adj[
              res_se_dir_coretype$dunn$Comparison ==
                "Multifacial - Unifacial_centripetal"], 3),
          kw_se_morph_layer_chi2   = round(res_se_morph_layer$kw$statistic, 3),
          kw_se_morph_layer_p      = round(res_se_morph_layer$kw$p.value, 3),
          kw_se_morph_rawmat_chi2  = round(res_se_morph_rawmat$kw$statistic, 3),
          kw_se_morph_rawmat_p     = round(res_se_morph_rawmat$kw$p.value, 3),
          kw_se_morph_type_chi2    = round(res_se_morph_coretype$kw$statistic, 3),
          kw_se_morph_type_p       = round(res_se_morph_coretype$kw$p.value, 3),
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
          fig_coia_composite = p_final,
          fig_se_composite   = p_se_composite
        )
      })
    }
  )
)