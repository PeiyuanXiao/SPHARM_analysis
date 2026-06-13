# =============================================================================
# preview.R  —  render one manuscript figure at its EXACT print size and view it
# =============================================================================
# Usage (in the container's R / RStudio console, working dir = project root):
#
#   source("analysis/scripts/_preview/preview.R")
#   preview_fig("im")          # edit the SZ block in the script, then re-run this
#   preview_fig("im", dpi=300) # crisper; or open the saved PNG for true size
#
# It re-sources the figure's script (so it picks up your latest SZ edits),
# rebuilds the figure, saves a PNG at the manuscript width/height (proportions
# = exactly what the paper gets), and shows it inline. The saved file path is
# printed; open that file in the Files pane / an image viewer to judge absolute
# physical size.
#
# NOTE: "method" loads ggtern, which patches ggplot2 for the session. If you
# tune `method` and then a different figure, restart R first.
# =============================================================================

suppressMessages({ library(here); library(ggplot2) })

preview_fig <- function(fig = c("rotation", "im", "method", "cia", "coia"),
                        dpi = 150, show = TRUE) {
  fig  <- match.arg(fig)
  spec <- list(
    rotation = list(s = "analysis/scripts/r_validation/validate_rotation_all.R",
                    w = 6.85, h = 8.0, get = function(e) e$p_rotational_invariance_validity),
    im       = list(s = "analysis/scripts/r_validation/methods_comparison_IM.R",
                    w = 6.85, h = 2.85, get = function(e) e$p_heatmap),
    method   = list(s = "analysis/scripts/r_spharm/spharm_analysis.R",
                    w = 6.85, h = 8.2, get = function(e) e$exp_method_compare_combined),
    cia      = list(s = "analysis/scripts/r_statistics/exp_cores_statistics.R",
                    w = 6.85, h = 8.6, get = function(e) e$p_final),
    coia     = list(s = "analysis/scripts/r_statistics/SDG_cores_statistics.R",
                    w = 6.85, h = 9.1, get = function(e) e$p_final)
  )[[fig]]

  e <- new.env(parent = globalenv())
  sys.source(here::here(spec$s), envir = e)
  p <- spec$get(e)

  out <- here::here("analysis/output/figures/tuning", paste0("fig-", fig, ".png"))
  dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
  ggsave(out, p, width = spec$w, height = spec$h, units = "in",
         dpi = dpi, bg = "white", limitsize = FALSE)

  if (isTRUE(show)) {
    if (requireNamespace("magick", quietly = TRUE)) {
      print(magick::image_read(out))            # inline, native pixels
    } else {
      grid::grid.newpage()
      grid::grid.raster(png::readPNG(out))      # inline, scaled to Plots pane
    }
  }
  message(sprintf("[preview] %s  (%.2f x %.2f in @ %d dpi)  ->  %s",
                  fig, spec$w, spec$h, dpi, out))
  invisible(out)
}
