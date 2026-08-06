# annotation_perturbation_figure.R
# =============================================================================
# FIGURE ASSEMBLY for the annotation-error perturbation analysis — SI add-on.
#
# PURE PLOTTING. Every number is read from the CSVs written by
# annotation_perturbation.R; nothing is re-perturbed and no test is re-run. The
# script can therefore be executed on its own, without the ~1 h replicate sweep:
#
#   Rscript analysis/robustness/annotation_perturbation/annotation_perturbation_figure.R
#
# It is also source()d by annotation_perturbation.R, which calls
# make_perturbation_figure() after its summary CSVs are on disk, so the full run
# and the standalone run produce byte-identical output.
#
# ONE FIGURE, replacing the four that this analysis used to emit
# (polarity-three-methods, angle-three-methods, degradation, continuous). The
# earlier set repeated the same replicates in four layouts; a reader had to hold
# two files side by side to see the one contrast that matters. The four are
# deleted here (see OBSOLETE_FIGS) so the folder cannot go stale.
#
# LAYOUT  2 x 2:
#   columns = perturbation   polarity flips (left) | angular jitter (right)
#   rows    = metric         pairs resolved, of 10 (top)
#                            median -log10(raw p)  (bottom)
#   Three methods per panel, median with the 2.5-97.5 percentile band. The rows
#   carry different units, so the y scale is free by row; the x units (% flipped,
#   degrees s.d.) sit on the tick labels because the column strips name only the
#   perturbation. Resolved pairs are Holm-corrected at alpha = 0.05 throughout —
#   stated here and in the CSVs' `test` column rather than on the axis.
#
# WHY BOTH ROWS. The two metrics answer different questions and neither implies
# the other. The median tracks the bulk shift in margin; the count tracks whether
# the marginal pairs survive. This project has been bitten by that gap more than
# once — SP-SPHARM's median -log10(p) is still at its baseline at sigma = 10 deg
# ("no degradation") while the resolved count has already fallen 8 -> 7, and in
# the truncation sweep the median moved 0.113 between l = 5 and l = 6 yet a pair
# was lost. Stacking them is the point of the figure; neither row is a summary of
# the other.
#
# ABSOLUTE AXES — AND WHAT THEY DO NOT SHOW. Both rows plot raw values: resolved
# pairs out of 10 on top, median -log10(raw p) below. The three methods do not
# start level (baselines 8/10, 6/10, 4/10 pairs and 2.841 / 2.430 / 1.559 median
# -log10 p), so the vertical offsets between the lines are baseline differences,
# not degradation differences, and SPI's flatness on the top row is partly just
# having least to lose. Slopes are comparable within a method across levels; they
# are NOT comparable between methods without dividing by each method's own
# baseline. Those retention ratios stay in both three-method CSVs and are the
# only basis on which this analysis ranks anything.
#
# NO CAPTION, NO IN-PANEL TEXT (deliberate, 2026-08-05). Three facts therefore
# have to travel with the figure in the SI text instead of on it:
#   (i)   the zero-width fabric band in the polarity column is the Mark 1973
#         identity — sum(u u') is unchanged by u -> -u, so all 500 polarity
#         replicates return bit-identical fabric results, not a plotting fault;
#   (ii)  each method keeps the test the main analysis gives it — SP-SPHARM and
#         fabric by PERMANOVA (adonis2), SPI by Kruskal-Wallis with Dunn post-hoc,
#         Holm throughout; nothing was harmonised to draw them together;
#   (iii) no ranking of the three methods under angular jitter is claimed: the two
#         metrics order them oppositely and every percentile band still spans the
#         method's own baseline at every level.
# (ii) and (iii) are re-checked against the plotted numbers on every run and
# printed to the console, so the SI text can be trusted against the current data.
#
# Reads : annotation_perturbation_polarity.csv              (wide, per level)
#         annotation_perturbation_angle_three_methods.csv   (tidy, incl. baselines)
# Writes: annotation_perturbation_polarity_three_methods.csv (tidy polarity twin
#           of the angle table — same schema, adds the retention columns, all of
#           it arithmetic on numbers already in the wide CSV)
#         figures/fig_S_perturbation_combined.png
# Deletes: the four superseded PNGs listed in OBSOLETE_FIGS.
# =============================================================================

suppressPackageStartupMessages({
  library(here); library(tidyverse); library(conflicted)
})
suppressMessages({
  conflicts_prefer(dplyr::filter, dplyr::select, dplyr::lag, .quiet = TRUE)
})

# Shared with annotation_perturbation.R: names() also fixes the method ordering.
METHOD_COLORS <- c("SP-SPHARM" = "#4A6E8A", "Fabric (E, I)" = "#BA8530",
                   "SPI" = "#802520")

# Column strips name the perturbation only; the units ride on the axis ticks
# (% for the flipped fraction, deg for the jitter s.d.). "(180 deg)" marks a
# polarity flip as the extreme rotation, which is what makes the two columns
# readable as one scale of annotation error.
PERT_LABS <- c(polarity = "polarity flips (180°)",
               angle    = "angular jitter")

# Superseded by fig_S_perturbation_combined.png.
OBSOLETE_FIGS <- c("fig_S_perturbation_polarity_three_methods.png",
                   "fig_S_perturbation_angle_three_methods.png",
                   "fig_S_perturbation_degradation.png",
                   "fig_S_perturbation_continuous.png")

# =============================================================================
make_perturbation_figure <- function(out_dir, n_pairs = 10L) {

  fig_dir <- file.path(out_dir, "figures")
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

  pol_csv <- file.path(out_dir, "annotation_perturbation_polarity.csv")
  ang_csv <- file.path(out_dir, "annotation_perturbation_angle_three_methods.csv")
  for (f in c(pol_csv, ang_csv))
    if (!file.exists(f))
      stop("missing input CSV: ", basename(f),
           " — run annotation_perturbation.R first (this script never recomputes it).")

  # ---- angle: already tidy, and the only file that records the baselines -----
  ang <- read_csv(ang_csv, show_col_types = FALSE)

  # Baselines are properties of the method, not of the perturbation (both sweeps
  # share the one unperturbed run), so the angle table supplies them for the
  # polarity column too.
  base_met <- ang %>%
    distinct(method, base_count = n_sig_05_baseline, base_neglog = med_neglog10_baseline)
  stopifnot(nrow(base_met) == 3, !anyNA(base_met))

  # ---- polarity: wide -> tidy, in the angle table's schema -------------------
  pol <- read_csv(pol_csv, show_col_types = FALSE)
  grab <- function(method, cnt, nlg)
    pol %>% transmute(perturbation = "polarity", level, method = !!method,
                      n_sig_05_med = .data[[paste0(cnt, "_med")]],
                      n_sig_05_lo  = .data[[paste0(cnt, "_lo")]],
                      n_sig_05_hi  = .data[[paste0(cnt, "_hi")]],
                      med_neglog10_med = .data[[paste0(nlg, "_med")]],
                      med_neglog10_lo  = .data[[paste0(nlg, "_lo")]],
                      med_neglog10_hi  = .data[[paste0(nlg, "_hi")]])
  pol_tidy <- bind_rows(
    grab("SP-SPHARM",     "n_sig_05",        "med_neglog10_p"),
    grab("Fabric (E, I)", "fabric_n_sig_05", "fabric_med_neglog10"),
    grab("SPI",           "spi_n_sig_05",    "spi_med_neglog10")) %>%
    # level 0 = the unperturbed run; zero-width band by construction.
    bind_rows(base_met %>% transmute(perturbation = "polarity", level = 0, method,
                                     n_sig_05_med = base_count, n_sig_05_lo = base_count,
                                     n_sig_05_hi = base_count,
                                     med_neglog10_med = base_neglog,
                                     med_neglog10_lo = base_neglog,
                                     med_neglog10_hi = base_neglog)) %>%
    left_join(base_met, by = "method") %>%
    mutate(n_sig_05_baseline = base_count,
           n_sig_05_retention = n_sig_05_med / base_count,
           med_neglog10_baseline = base_neglog,
           med_neglog10_retention = med_neglog10_med / base_neglog,
           test = if_else(method == "SPI", "Kruskal-Wallis + Dunn (Holm)",
                          "PERMANOVA adonis2 + Holm")) %>%
    select(all_of(names(ang))) %>%
    arrange(match(method, names(METHOD_COLORS)), level)

  pol_out <- file.path(out_dir, "annotation_perturbation_polarity_three_methods.csv")
  write_csv(pol_tidy, pol_out)

  # ---- long form for plotting -----------------------------------------------
  # Raw values, not ratios: the ribbon is the 2.5-97.5 percentile band of the
  # statistic itself. The rows therefore carry different units, which is why the
  # y scale is freed by row below.
  as_panel <- function(d, cnt_or_nlg) {
    pre <- if (cnt_or_nlg == "count") "n_sig_05" else "med_neglog10"
    tibble(perturbation = d$perturbation, level = d$level, method = d$method,
           metric = cnt_or_nlg,
           val_med = d[[paste0(pre, "_med")]],
           val_lo  = d[[paste0(pre, "_lo")]],
           val_hi  = d[[paste0(pre, "_hi")]])
  }
  plot_df <- bind_rows(
    as_panel(pol_tidy, "count"), as_panel(pol_tidy, "neglog"),
    as_panel(ang,      "count"), as_panel(ang,      "neglog")) %>%
    mutate(method = factor(method, levels = names(METHOD_COLORS)),
           perturbation = factor(perturbation, levels = names(PERT_LABS)),
           metric = factor(metric, levels = c("count", "neglog"),
                           labels = c(sprintf("pairs resolved (of %d)", n_pairs),
                                      "median -log10(raw p)"))) %>%
    arrange(perturbation, metric, method, level)
  stopifnot(!anyNA(plot_df$val_med), !anyNA(plot_df$val_lo), !anyNA(plot_df$val_hi))

  # Pins the count row to the full 0-n_pairs range so the reader sees how much
  # headroom each method leaves; the median row keeps a data-driven range.
  span <- tibble(metric = factor(levels(plot_df$metric)[1],
                                 levels = levels(plot_df$metric)),
                 val_med = c(0, n_pairs)) %>%
    crossing(perturbation = factor(names(PERT_LABS), levels = names(PERT_LABS))) %>%
    mutate(level = 0, method = factor("SP-SPHARM", levels = names(METHOD_COLORS)))

  # x is a fraction in the polarity column and degrees in the angle column; the
  # two ranges are an order of magnitude apart, which is what lets one breaks
  # function serve both free_x scales. With the column strips naming only the
  # perturbation, the tick labels are where the units live.
  brk_x <- function(lims) if (max(lims, na.rm = TRUE) <= 1)
    c(0, 0.02, 0.05, 0.10, 0.15, 0.20) else c(0, 5, 10, 15, 20)
  lab_x <- function(b) {
    b <- b[!is.na(b)]
    if (length(b) && max(b) <= 1) paste0(round(b * 100), "%") else paste0(round(b), "°")
  }
  # Same trick on y, keyed on the row's range: counts run to 10, the median to ~4.
  brk_y <- function(lims)
    if (max(lims, na.rm = TRUE) > 5) seq(0, n_pairs, 2) else seq(0, 4, 1)

  p <- ggplot(plot_df, aes(level, val_med, color = method, fill = method)) +
    geom_blank(data = span) +
    geom_ribbon(aes(ymin = val_lo, ymax = val_hi), alpha = 0.18, colour = NA) +
    geom_line(linewidth = 0.8) + geom_point(size = 1.7) +
    facet_grid(metric ~ perturbation, scales = "free",
               labeller = labeller(perturbation = PERT_LABS), switch = "y") +
    scale_color_manual(values = METHOD_COLORS, name = NULL) +
    scale_fill_manual(values = METHOD_COLORS, name = NULL) +
    scale_x_continuous(breaks = brk_x, labels = lab_x, expand = expansion(mult = 0.045)) +
    scale_y_continuous(breaks = brk_y) +
    labs(x = NULL, y = NULL) +
    theme_bw(base_size = 10) +
    theme(panel.grid.minor = element_blank(), legend.position = "bottom",
          strip.text = element_text(face = "bold", size = 9),
          strip.placement = "outside", strip.background.y = element_blank(),
          legend.key.size = grid::unit(0.34, "cm"))
  fig_path <- file.path(fig_dir, "fig_S_perturbation_combined.png")
  ggsave(fig_path, p, width = 8.5, height = 5.6, dpi = 300)

  # ---- retire the superseded figures -----------------------------------------
  # Idempotent: a file already gone is reported as such rather than silently
  # dropped from the list, so the console record names all four either way.
  fate <- vapply(OBSOLETE_FIGS, function(f) {
    fp <- file.path(fig_dir, f)
    if (!file.exists(fp)) "already absent"
    else if (file.remove(fp)) "deleted" else "COULD NOT DELETE"
  }, character(1))
  removed <- names(fate)[fate == "deleted"]

  # ---- report ----------------------------------------------------------------
  cat("\n", strrep("=", 70), "\n", sep = "")
  cat("COMBINED PERTURBATION FIGURE\n")
  cat(strrep("=", 70), "\n", sep = "")
  cat(sprintf("  read   %s  (%d rows)\n", basename(pol_csv), nrow(pol)))
  cat(sprintf("  read   %s  (%d rows)\n", basename(ang_csv), nrow(ang)))
  cat(sprintf("  wrote  %s  (%d rows)\n", basename(pol_out), nrow(pol_tidy)))
  # Panels are not equal in size: the polarity sweep has one more level than the
  # angle sweep, so the per-panel counts are listed rather than divided out.
  panel_n <- plot_df %>% count(perturbation, metric, name = "rows")
  cat(sprintf("  plotted rows: %d in 4 panels\n", nrow(plot_df)))
  for (i in seq_len(nrow(panel_n)))
    cat(sprintf("    %-8s / %-32s %2d rows = %d levels (incl. 0) x 3 methods\n",
                as.character(panel_n$perturbation[i]), as.character(panel_n$metric[i]),
                panel_n$rows[i], panel_n$rows[i] / 3))
  # These are the level-0 points on the figure, and the denominators of the
  # retention columns the CSVs keep. Printed because the figure carries no key.
  cat("\n  Baselines (unperturbed EXP run; the level-0 point in every panel):\n")
  for (i in seq_len(nrow(base_met)))
    cat(sprintf("    %-14s %d/%d pairs resolved   median -log10(p) = %.3f\n",
                base_met$method[i], as.integer(base_met$base_count[i]), n_pairs,
                base_met$base_neglog[i]))

  # The figure no longer carries a caption, so the two claims the SI text has to
  # make on its behalf are re-checked against the plotted numbers on every run.
  spans <- ang %>%
    transmute(level, method,
              cnt = n_sig_05_lo / n_sig_05_baseline <= 1 & n_sig_05_hi / n_sig_05_baseline >= 1,
              nlg = med_neglog10_lo / med_neglog10_baseline <= 1 &
                med_neglog10_hi / med_neglog10_baseline >= 1)
  cat(sprintf("\n  SI-text check - every angle band spans its own baseline: %s\n",
              ifelse(all(spans$cnt) && all(spans$nlg), "TRUE (both metrics, all levels)",
                     "FALSE <-- the no-ranking claim does not hold, inspect")))
  # And that the two metrics disagree on the ordering at the largest jitter.
  worst <- ang %>% filter(level == max(level)) %>%
    transmute(method, ret_cnt = n_sig_05_retention, ret_nlg = med_neglog10_retention)
  ord_cnt <- worst$method[order(worst$ret_cnt)]
  ord_nlg <- worst$method[order(worst$ret_nlg)]
  cat(sprintf("  SI-text check - orderings at sigma = %.0f deg differ: %s\n",
              max(ang$level), ifelse(!identical(ord_cnt, ord_nlg), "TRUE", "FALSE <-- inspect")))
  cat(sprintf("    by count  (worst -> best): %s\n", paste(ord_cnt, collapse = " < ")))
  cat(sprintf("    by median (worst -> best): %s\n", paste(ord_nlg, collapse = " < ")))

  cat(sprintf("\n  Wrote figures/%s\n", basename(fig_path)))
  cat(sprintf("  Superseded figures retired (%d of %d deleted this run):\n",
              length(removed), length(OBSOLETE_FIGS)))
  for (f in OBSOLETE_FIGS)
    cat(sprintf("    - figures/%-48s %s\n", f, fate[[f]]))
  cat(sprintf("  figures/ now holds: %s\n",
              paste(sort(basename(list.files(fig_dir, pattern = "\\.png$"))), collapse = ", ")))

  invisible(list(figure = fig_path, csv = pol_out, removed = removed,
                 baselines = base_met, plot_df = plot_df))
}

# Standalone entry point: sys.nframe() is 0 only when this file is the script
# being run, not when annotation_perturbation.R source()s it.
if (sys.nframe() == 0L) {
  make_perturbation_figure(here("analysis/robustness/annotation_perturbation"))
  cat("\nDone.\n\n")
  print(sessionInfo())
}
