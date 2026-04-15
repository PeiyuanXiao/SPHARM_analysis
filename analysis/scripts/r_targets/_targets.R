library(targets)
library(tarchetypes)

tar_option_set(packages = c("tidyverse", "patchwork", "glue", "here"))

list(
  tar_target(
    p_rotational_invariance_validity,
    {
      source(here::here("analysis/scripts/r_validation/validate_rotation_all.R"))
      p_rotational_invariance_validity
    }
  )
)