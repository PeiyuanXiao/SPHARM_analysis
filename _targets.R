library(targets)
library(tarchetypes)
tar_option_set(
  packages = c("tidyverse", "patchwork", "here",
               "readxl", "vegan", "FSA", "RVAideMemoire",
               "ggrepel", "conflicted",
               "linkET", "compositions", "ade4", "circular", "rsvg", "png", "grid")
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

list(
  
  tar_target(
    spharm_analysis,
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
  ),
  
  tar_target(
    p_rotational_invariance_validity,
    local({
      source(here::here("analysis/scripts/r_validation/validate_rotation_all.R"),
             local = TRUE)
      p_rotational_invariance_validity
    })
  ),
  
  tar_target(
    im_comparison,
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
  ),
  
  tar_target(
    exp_cia_analysis,
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
  )
  
)