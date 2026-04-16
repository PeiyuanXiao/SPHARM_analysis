library(targets)
library(tarchetypes)
tar_option_set(
  packages = c("tidyverse", "patchwork", "here",
               "readxl", "MASS", "vegan", "FSA", "RVAideMemoire")
)

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
        )
      )
    })
  ),
  
  tar_target(
    p_exp_method_combined,            
    local({
      source(here::here("analysis/scripts/r_spharm/spharm_analysis.R"),
             local = TRUE)
      exp_method_compare_combined
    })
  ),
  
  tar_target(
    p_rotational_invariance_validity,
    {
      source(here::here("analysis/scripts/r_validation/validate_rotation_all.R"))
      p_rotational_invariance_validity
    }
  )
)