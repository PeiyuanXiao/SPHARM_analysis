# ==============================================================================
# make.R
# One-command reproduction: run the full targets pipeline, then render the
# reproduction report (the complete statistical output that tar_make() prints,
# as a single self-contained HTML).
#
# Usage (inside the Docker container, from /project):
#   Rscript make.R
#
# This is a thin convenience wrapper. It is exactly equivalent to running, in
# order:
#   targets::tar_make()
#   quarto render analysis/paper/reproduction_report.qmd
#
# Rendering is kept here (outside the targets pipeline) on purpose, matching how
# manuscript.qmd and supplementary.qmd are rendered as post-pipeline steps in CI
# (.github/workflows/reproducibility.yml). If tar_make() fails, the script stops
# before rendering, so the report is never built from a half-finished pipeline.
# ==============================================================================

targets::tar_make()

quarto::quarto_render(here::here("analysis/paper/reproduction_report.qmd"))

cat("\nReproduction report written to:\n  ",
    here::here("analysis/paper/reproduction_report.html"), "\n", sep = "")
