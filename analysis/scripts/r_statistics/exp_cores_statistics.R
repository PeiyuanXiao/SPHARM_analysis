# exp_cores_statistics.R
# Techno-morphological joint analysis of experimental cores. Sourced by `exp_cia_analysis`.
#
# Analysis framework:
#   L1     : global Mantel + CoIA (baseline), with Sankey of the contribution flow.
#   L2-A   : per-type Mantel.
#   L2-B   : CoIA arrow-length group differences.
#   L2-C   : CoIA arrow-direction circular statistics.
#
# Input:
#   - analysis/data/derived_data/SPHARM_direction.csv
#   - analysis/data/derived_data/SPHARM_morphology.csv
#   - analysis/data/raw_data/SDG_core_metric.xlsx
#   - analysis/figures/Axis_trajectory.png (panel D of the composite figure)
#
# Returns (objects): p_final, plus statistics consumed by the paper.

library(here)
library(tidyverse)
library(readxl)
library(vegan)
library(linkET)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(compositions)
library(ade4)
library(circular)
library(FSA)
library(png)
library(grid)

# ==============================================================================
# ---- Global constants: palettes and type order ----
# ==============================================================================

TYPOLOGY_COLORS <- c(
  "Levallois"      = "#4A6E8A",
  "Discoid"        = "#802520",
  "Unidirectional" = "#BA8530",
  "Multiplatform"  = "#8A7A68",
  "Bidirectional"  = "#788C4A"
)

TYPOLOGY_ORDER <- c(
  "Unidirectional",
  "Bidirectional",
  "Levallois",
  "Discoid",
  "Multiplatform"
)

# ==============================================================================
# ---- Global helpers ----
# ==============================================================================

cosine_dist <- function(X) {
  X   <- as.matrix(X)
  sim <- X %*% t(X) /
    (sqrt(rowSums(X^2)) %o% sqrt(rowSums(X^2)))
  sim <- pmin(pmax(sim, -1), 1)
  as.dist(1 - sim)
}

replace_zeros <- function(X, delta = NULL) {
  X <- as.matrix(X)
  for (i in seq_len(nrow(X))) {
    row_i     <- X[i, ]
    zero_idx  <- row_i == 0
    if (!any(zero_idx)) next
    nonzero_min          <- min(row_i[!zero_idx])
    d                    <- ifelse(is.null(delta), nonzero_min * 0.65, delta)
    n_zero               <- sum(zero_idx)
    row_i[zero_idx]      <- d
    row_i[!zero_idx]     <- row_i[!zero_idx] * (1 - n_zero * d)
    X[i, ] <- row_i
  }
  X
}

extract_subdist <- function(D_full, ids) {
  as.dist(as.matrix(D_full)[ids, ids])
}

safe_filter_groups <- function(meta_df, group_col, min_n = 3) {
  counts <- table(meta_df[[group_col]], useNA = "no")
  valid  <- names(counts[counts >= min_n])
  if (length(valid) < 2) {
    cat(sprintf(
      "  [skip] %s: too few valid groups (need >= 2 groups, each >= %d). Current: %s\n",
      group_col, min_n,
      paste(names(counts), counts, sep = "=", collapse = ", ")
    ))
    return(NULL)
  }
  meta_df %>%
    filter(!is.na(.data[[group_col]]),
           .data[[group_col]] %in% valid)
}

circ_stats_one <- function(angles_rad) {
  circ_obj <- circular(angles_rad, type = "angles",
                       units = "radians", modulo = "2pi")
  mean_rad <- as.numeric(mean.circular(circ_obj)) %% (2 * pi)
  list(
    mean_rad = mean_rad,
    mean_deg = mean_rad * 180 / pi,
    rho      = as.numeric(rho.circular(circ_obj))
  )
}

# ==============================================================================
# ---- Data preparation ----
# ==============================================================================

POWER_COLS_DIR   <- paste0("power_l", 1:6)   # direction spectrum l=1-6
POWER_COLS_MORPH <- paste0("power_l", 1:8)  # morphology spectrum l=1-8

SPHARM_direction  <- read_csv(here("analysis/data/derived_data/SPHARM_direction.csv"),
                              show_col_types = FALSE)
SPHARM_morphology <- read_csv(here("analysis/data/derived_data/SPHARM_morphology.csv"),
                              show_col_types = FALSE)
metric_data <- read_xlsx(here("analysis/data/raw_data/SDG_core_metric.xlsx")) %>%
  select(ID, Layer, Core_type_Li_merged, Raw_mat) %>%
  mutate(Layer = as.factor(Layer), Raw_mat = as.factor(Raw_mat))

SPHARM_morphology <- SPHARM_morphology %>%
  left_join(SPHARM_direction %>% select(ID, Typology), by = "ID")

filter_spharm <- function(df, power_cols, meta = NULL) {
  result <- df %>%
    select(ID, Typology, all_of(power_cols))
  if (!is.null(meta)) result <- left_join(result, meta, by = "ID")
  result
}

SPHARM_direction_filter  <- filter_spharm(SPHARM_direction,  POWER_COLS_DIR,   metric_data)
SPHARM_morphology_filter <- filter_spharm(SPHARM_morphology, POWER_COLS_MORPH, metric_data)

split_by_group <- function(df) {
  list(
    exp_im = df %>% filter(str_starts(ID, "EXP") | str_starts(ID, "IM_")),
    sdg_im = df %>% filter(str_starts(ID, "SDG") | str_starts(ID, "IM_"))
  )
}
dir_splits   <- split_by_group(SPHARM_direction_filter)
morph_splits <- split_by_group(SPHARM_morphology_filter)

df_scar_all  <- dir_splits$exp_im
df_morph_all <- morph_splits$exp_im

common_ids   <- intersect(df_morph_all$ID, df_scar_all$ID)
df_morph_all <- df_morph_all %>% filter(ID %in% common_ids) %>% arrange(ID)
df_scar_all  <- df_scar_all  %>% filter(ID %in% common_ids) %>% arrange(ID)
stopifnot(all(df_morph_all$ID == df_scar_all$ID))

cat("==== Data alignment (EXP + IM) ====\n")
cat("Shared specimens:", length(common_ids),
    "; IDs fully matched:", all(df_morph_all$ID == df_scar_all$ID), "\n\n")

morph_power_all <- df_morph_all %>%
  select(all_of(POWER_COLS_MORPH)) %>%
  rename_with(~ paste0("M", seq_along(.))) %>%
  as.data.frame()
scar_power_all  <- df_scar_all %>%
  select(all_of(POWER_COLS_DIR)) %>%
  rename_with(~ paste0("S", seq_along(.))) %>%
  as.data.frame()
rownames(morph_power_all) <- df_morph_all$ID
rownames(scar_power_all)  <- df_scar_all$ID

morph_power_clean <- morph_power_all[, sapply(morph_power_all, sd, na.rm = TRUE) > 0]
scar_power_clean  <- scar_power_all[,  sapply(scar_power_all,  sd, na.rm = TRUE) > 0]

# ILR transform (removes the compositional sum constraint) -> Euclidean distance
morph_ilr_all <- as.data.frame(ilr(replace_zeros(as.matrix(morph_power_clean))))
scar_ilr_all  <- as.data.frame(ilr(replace_zeros(as.matrix(scar_power_clean))))
rownames(morph_ilr_all) <- rownames(morph_power_clean)
rownames(scar_ilr_all)  <- rownames(scar_power_clean)

D_morph_all <- dist(morph_ilr_all)
D_scar_all  <- dist(scar_ilr_all)

exp_ids <- rownames(morph_power_clean)[
  !str_starts(rownames(morph_power_clean), "IM_") &
    rownames(morph_power_clean) != "EXP43_Biface"]
cat("Experimental specimens (excluding IM_):", length(exp_ids), "\n")

morph_exp     <- morph_power_clean[exp_ids, ]
scar_exp      <- scar_power_clean[exp_ids, ]
morph_ilr_exp <- morph_ilr_all[exp_ids, ]
scar_ilr_exp  <- scar_ilr_all[exp_ids, ]

D_morph_exp <- extract_subdist(D_morph_all, exp_ids)
D_scar_exp  <- extract_subdist(D_scar_all,  exp_ids)

meta_exp <- df_morph_all %>%
  filter(ID %in% exp_ids) %>%
  select(ID, Typology) %>%
  left_join(metric_data, by = "ID")

# Unify Levallois spelling
meta_exp <- meta_exp %>%
  mutate(Typology = if_else(
    str_detect(Typology, regex("levallois", ignore_case = TRUE)),
    "Levallois", Typology
  ))

df_morph_all <- df_morph_all %>%
  mutate(Typology = if_else(
    str_detect(Typology, regex("levallois", ignore_case = TRUE)),
    "Levallois", Typology
  ))

cat("\n==== Experimental metadata (Levallois merged) ====\n")
cat("Typology：\n"); print(table(meta_exp$Typology, useNA = "ifany"))

meta_typology <- safe_filter_groups(meta_exp, "Typology")

plot_ids_no_biface <- meta_exp %>% filter(Typology != "Biface") %>% pull(ID)

typology_levels <- TYPOLOGY_ORDER[TYPOLOGY_ORDER %in%
                                    unique(meta_exp$Typology[
                                      meta_exp$Typology != "Biface" &
                                        !is.na(meta_exp$Typology)])]
typology_pal    <- TYPOLOGY_COLORS[typology_levels]

# ==============================================================================
# ========== Level 1: global Mantel + CoIA ==========
# ==============================================================================

cat("\n\n## Level 1: global Mantel + CoIA — baseline (EXP) ##\n")

# ------------------------------------------------------------------------------
# L1-1: global Mantel + linkET network
# ------------------------------------------------------------------------------

cat("\n==== L1-1: global Mantel ====\n")
mantel_global <- mantel(D_morph_exp, D_scar_exp,
                        method = "spearman", permutations = 9999)
print(mantel_global)

run_cross_mantel <- function(X_single, D_target, from_label, var_label,
                             n_perm = 999) {
  if (sd(X_single, na.rm = TRUE) == 0) return(NULL)
  res <- mantel(dist(scale(X_single)), D_target,
                method = "spearman", permutations = n_perm)
  tibble(from = from_label, var = var_label,
         r = res$statistic, p = res$signif)
}

morph_exp_df <- as.data.frame(morph_exp)
scar_exp_df  <- as.data.frame(scar_exp)

mantel_cross_l1 <- bind_rows(
  map_dfr(colnames(morph_exp_df),
          ~ run_cross_mantel(morph_exp_df[[.x]], D_scar_exp,
                             "Scar Direction", .x)),
  map_dfr(colnames(scar_exp_df),
          ~ run_cross_mantel(scar_exp_df[[.x]], D_morph_exp,
                             "Morphology", .x))
) %>%
  mutate(
    p_holm       = p.adjust(p, method = "holm"),
    significance = ifelse(p_holm < 0.05, "P\u22640.05", "P>0.05")
  )

spec_exp_full <- bind_cols(morph_exp_df, scar_exp_df)

p_mantel_net <- qcorrplot(
  correlate(spec_exp_full, method = "spearman"),
  type = "upper", diag = FALSE
) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_couple(
    aes(colour = significance, size = abs(r)),
    data         = mantel_cross_l1,
    curvature    = 0.15,
    label.params = list(color = "transparent")
  ) +
  scale_fill_gradient2(
    low      = "#802520",
    mid      = "#F5EDDC",
    high     = "#4A6E8A",
    midpoint = 0,
    limits   = c(-1, 1),
    name     = "Spearman's rho"
  ) +
  scale_color_manual(
    values = c("P\u22640.05" = "#E6A5A5", "P>0.05" = "#BABABA"),
    name   = "Mantel test\n(Holm corrected)"
  ) +
  scale_size_continuous(range = c(0.5, 2.5), name = "Mantel's |r|") +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid    = element_blank(),
    axis.title    = element_blank(),
    legend.position = "right",
    plot.margin   = margin(20, 20, 20, 20)
  )

cat("Figure built: EXP_L1_Mantel_Network.png\n")

# ------------------------------------------------------------------------------
# L1-2: CoIA + RV permutation test
# ------------------------------------------------------------------------------

cat("\n==== L1-2：CoIA ====\n")

morph_exp_ilr <- morph_ilr_exp
scar_exp_ilr  <- scar_ilr_exp
colnames(morph_exp_ilr) <- paste0("M_ilr", seq_len(ncol(morph_exp_ilr)))
colnames(scar_exp_ilr)  <- paste0("S_ilr", seq_len(ncol(scar_exp_ilr)))

dudi_morph <- dudi.pca(morph_exp_ilr, center = TRUE, scale = TRUE,
                       scannf = FALSE, nf = ncol(morph_exp_ilr))
dudi_scar  <- dudi.pca(scar_exp_ilr,  center = TRUE, scale = TRUE,
                       scannf = FALSE, nf = ncol(scar_exp_ilr))

# Detailed PCA report
report_pca <- function(dudi_obj, label) {
  n_ax <- length(dudi_obj$eig)
  eig  <- dudi_obj$eig
  pct  <- eig / sum(eig) * 100
  cum  <- cumsum(pct)
  
  cat(sprintf("\n====== %s PCA report ======\n", label))
  
  cat("\n-- Eigenvalues and explained variance --\n")
  eig_df <- tibble(
    Axis       = paste0("PC", seq_len(n_ax)),
    Eigenvalue = round(eig, 4),
    Pct_var    = round(pct, 2),
    Cumul_pct  = round(cum, 2)
  )
  print(as.data.frame(eig_df))
  
  cat("\n-- Variable loadings (c1: ILR variables on PCA axes) --\n")
  load_df <- as.data.frame(dudi_obj$c1)
  colnames(load_df) <- paste0("PC", seq_len(ncol(load_df)))
  n_ilr <- nrow(load_df)
  load_df$ILR_meaning <- sapply(seq_len(n_ilr), function(k) {
    num_ids <- paste0("l", seq_len(k))
    den_id  <- paste0("l", k + 1)
    sprintf("log(geomean(%s) / %s)", paste(num_ids, collapse = "+"), den_id)
  })
  print(load_df)
  
  cat("\n-- Score summary (first 2 axes) --\n")
  score_df <- as.data.frame(dudi_obj$li)[, 1:min(2, n_ax), drop = FALSE]
  colnames(score_df) <- paste0("PC", seq_len(ncol(score_df)))
  score_stats <- score_df %>%
    pivot_longer(everything(), names_to = "Axis", values_to = "Score") %>%
    group_by(Axis) %>%
    summarise(
      mean = round(mean(Score), 4), sd  = round(sd(Score),  4),
      min  = round(min(Score),  4), max = round(max(Score), 4),
      .groups = "drop"
    )
  print(as.data.frame(score_stats))
  
  cat("\n-- Dominant variable per axis (max |loading|) --\n")
  load_num <- as.data.frame(dudi_obj$c1)
  for (ax in seq_len(min(2, n_ax))) {
    col <- load_num[[ax]]
    idx <- which.max(abs(col))
    cat(sprintf("  PC%d (%.1f%% var): dominant ILR%d (loading %+.4f) -> %s\n",
                ax, pct[ax], idx, col[idx], load_df$ILR_meaning[idx]))
  }
  invisible(list(eig_df = eig_df, load_df = load_df))
}

pca_report_morph <- report_pca(dudi_morph, "morphology spectrum")
pca_report_scar  <- report_pca(dudi_scar,  "direction spectrum")

coin_exp    <- coinertia(dudi_morph, dudi_scar, scannf = FALSE, nf = 2)
cia_inertia <- coin_exp$eig / sum(coin_exp$eig) * 100

cat("RV coefficient:", round(coin_exp$RV, 4), "\n")

set.seed(42)
rv_test <- randtest(coin_exp, nrepet = 9999)
cat("\nRV permutation test:\n"); print(rv_test)

scores_morph <- as.data.frame(coin_exp$lX) %>% rownames_to_column("ID")
scores_scar  <- as.data.frame(coin_exp$lY) %>% rownames_to_column("ID")

scores_combined <- left_join(
  scores_morph %>% select(ID, Axis1_M = AxcX1, Axis2_M = AxcX2),
  scores_scar  %>% select(ID, Axis1_S = AxcY1, Axis2_S = AxcY2),
  by = "ID"
) %>%
  mutate(
    arrow_length = sqrt((Axis1_M - Axis1_S)^2 + (Axis2_M - Axis2_S)^2),
    arrow_angle  = atan2(Axis2_S - Axis2_M, Axis1_S - Axis1_M)
  ) %>%
  left_join(meta_exp %>% select(ID, Typology), by = "ID")

# CoIA coordinate output
cat("\n==== CoIA specimen coordinates ====\n")
cia_coords <- scores_combined %>%
  select(ID, Typology,
         Morph_Axis1 = Axis1_M, Morph_Axis2 = Axis2_M,
         Scar_Axis1  = Axis1_S, Scar_Axis2  = Axis2_S,
         arrow_length, arrow_angle)
print(cia_coords %>%
        mutate(across(where(is.numeric), ~ round(.x, 4))) %>%
        as.data.frame())
write_csv(cia_coords,
          here("analysis/data/derived_data/EXP_CIA_coords_full.csv"))
cat("Saved: EXP_CIA_coords_full.csv\n")

# PCA-axis contribution to CoIA axes
cat("\n==== Per-side PCA-axis contribution to CoIA axes (weight^2) ====\n")

compute_pca_cia_contribution <- function(a_mat, pct_vec, endpoint_label) {
  df <- as.data.frame(a_mat)
  colnames(df) <- paste0("CoIA_Ax", seq_len(ncol(df)))
  df$PC       <- paste0("PC", seq_len(nrow(df)))
  df$var_pct  <- pct_vec[seq_len(nrow(df))]
  df$endpoint <- endpoint_label
  for (ax in colnames(df)[startsWith(colnames(df), "CoIA_Ax")]) {
    df[[paste0(ax, "_w2")]] <- round(df[[ax]]^2, 4)
  }
  for (ax in paste0("CoIA_Ax", seq_len(ncol(a_mat)))) {
    w2col   <- paste0(ax, "_w2")
    rel_col <- paste0(ax, "_contrib_pct")
    df[[rel_col]] <- round(df[[w2col]] / sum(df[[w2col]]) * 100, 1)
  }
  df %>% select(endpoint, PC, var_pct,
                starts_with("CoIA_Ax1"), starts_with("CoIA_Ax2"))
}

morph_contrib <- compute_pca_cia_contribution(
  coin_exp$aX,
  round(dudi_morph$eig / sum(dudi_morph$eig) * 100, 1),
  "Morphology"
)
scar_contrib <- compute_pca_cia_contribution(
  coin_exp$aY,
  round(dudi_scar$eig / sum(dudi_scar$eig) * 100, 1),
  "Scar direction"
)

pca_cia_contrib <- bind_rows(morph_contrib, scar_contrib)
cat("\nMorphology side:\n"); print(morph_contrib %>% select(-endpoint) %>% as.data.frame())
cat("\nDirection side:\n"); print(scar_contrib  %>% select(-endpoint) %>% as.data.frame())
write_csv(pca_cia_contrib,
          here("analysis/data/derived_data/EXP_PCA_CoIA_contribution.csv"))
cat("Saved: EXP_PCA_CoIA_contribution.csv\n")

# ------------------------------------------------------------------------------
# CoIA diagnostic plots
# ------------------------------------------------------------------------------

eig_df <- tibble(
  axis       = paste0("Axis ", seq_along(coin_exp$eig)),
  eigenvalue = coin_exp$eig,
  pct        = coin_exp$eig / sum(coin_exp$eig) * 100,
  cum_pct    = cumsum(pct)
)

p_scree <- ggplot(eig_df, aes(x = axis, y = pct)) +
  geom_col(fill = "#4A6E8A", alpha = 0.85, width = 0.55) +
  geom_line(aes(y = cum_pct, group = 1), color = "#802520", linewidth = 0.8) +
  geom_point(aes(y = cum_pct), color = "#802520", size = 2.8) +
  geom_text(aes(y = pct + 1.5, label = sprintf("%.1f%%", pct)),
            size = 3, color = "grey30") +
  scale_y_continuous(
    name     = "Explained co-inertia (%)",
    sec.axis = sec_axis(~ ., name = "Cumulative (%)")
  ) +
  theme_bw(base_size = 10) +
  labs(
    title    = "CoIA Scree Plot (EXP)",
    subtitle = sprintf("RV = %.3f, p = %.3f", coin_exp$RV, rv_test$pvalue),
    x = "CoIA Axis"
  ) +
  theme(
    plot.title    = element_text(face = "bold", hjust = 0.5, size = 11),
    plot.subtitle = element_text(hjust = 0.5, size = 8.5, color = "grey50")
  )

morph_pct <- round(dudi_morph$eig / sum(dudi_morph$eig) * 100, 1)
scar_pct  <- round(dudi_scar$eig  / sum(dudi_scar$eig)  * 100, 1)

morph_load <- as.data.frame(coin_exp$aX) %>%
  rownames_to_column("variable") %>%
  rename(Axis1 = AxcX1, Axis2 = AxcX2) %>%
  mutate(
    pct            = morph_pct[as.integer(str_extract(variable, "[0-9]+"))],
    variable_label = sprintf("Morph-PCA%s\n(%.1f%% var)",
                             str_extract(variable, "[0-9]+"), pct),
    endpoint       = "Morphology"
  )

scar_load <- as.data.frame(coin_exp$aY) %>%
  rownames_to_column("variable") %>%
  rename(Axis1 = AxcY1, Axis2 = AxcY2) %>%
  mutate(
    pct            = scar_pct[as.integer(str_extract(variable, "[0-9]+"))],
    variable_label = sprintf("Dir-PCA%s\n(%.1f%% var)",
                             str_extract(variable, "[0-9]+"), pct),
    endpoint       = "Scar Direction"
  )

circle_df <- tibble(angle = seq(0, 2 * pi, length.out = 300),
                    x = cos(angle), y = sin(angle))

make_loading_plot <- function(load_df, title_str, col_fill) {
  ggplot(load_df) +
    geom_path(data = circle_df, aes(x = x, y = y),
              color = "grey80", linewidth = 0.4, linetype = "dashed") +
    geom_hline(yintercept = 0, color = "grey70", linewidth = 0.3) +
    geom_vline(xintercept = 0, color = "grey70", linewidth = 0.3) +
    geom_segment(aes(x = 0, y = 0, xend = Axis1, yend = Axis2),
                 arrow = arrow(length = unit(0.20, "cm"), type = "closed"),
                 color = col_fill, linewidth = 0.9, alpha = 0.85) +
    geom_label(aes(x = Axis1 * 1.12, y = Axis2 * 1.12,
                   label = variable_label),
               size = 2.6, color = "grey20",
               label.size = 0.2, fill = "white", alpha = 0.85,
               lineheight = 0.85) +
    coord_fixed(xlim = c(-1.35, 1.35), ylim = c(-1.35, 1.35)) +
    theme_bw(base_size = 10) +
    labs(
      title = title_str,
      x = sprintf("CoIA Axis 1 (%.1f%%)", cia_inertia[1]),
      y = sprintf("CoIA Axis 2 (%.1f%%)", cia_inertia[2])
    ) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 10),
      panel.grid = element_blank()
    )
}

p_load_morph <- make_loading_plot(morph_load,
                                  "Morphology PCA axes on CoIA space", "#4A6E8A")
p_load_scar  <- make_loading_plot(scar_load,
                                  "Scar direction PCA axes on CoIA space", "#BA8530")

p_cia_diagnostics <- (p_scree | p_load_morph | p_load_scar) +
  plot_annotation(
    title   = "CoIA Axis Diagnostics (EXP)",
    caption = paste(
      "Left: scree plot. Middle: morphology PCA axis loadings on CoIA axes.",
      "Right: scar-direction PCA axis loadings on CoIA axes.",
      "\nBoth panels use PCA-axis projections (aX / aY) for symmetric interpretation.",
      "Arrow length = contribution to CoIA structure.",
      "\nAxis 1 = global regularity (low-freq vs mid-freq energy contrast);",
      "Axis 2 = isotropy vs bipolarity (l1 vs l2 contrast)."
    ),
    theme = theme(
      plot.title   = element_text(face = "bold", hjust = 0.5, size = 12),
      plot.caption = element_text(size = 7.5, color = "grey50", hjust = 0)
    )
  )

cat("Figure built: EXP_L1_CIA_Diagnostics.png\n")

# CoIA main biplot
scores_long_plot <- bind_rows(
  scores_combined %>%
    filter(!is.na(Typology), ID %in% plot_ids_no_biface) %>%
    select(ID, Typology, x = Axis1_M, y = Axis2_M, arrow_length) %>%
    mutate(endpoint = "Morphology"),
  scores_combined %>%
    filter(!is.na(Typology), ID %in% plot_ids_no_biface) %>%
    select(ID, Typology, x = Axis1_S, y = Axis2_S, arrow_length) %>%
    mutate(endpoint = "Scar direction")
) %>%
  mutate(
    endpoint = factor(endpoint, levels = c("Morphology", "Scar direction")),
    Typology = factor(Typology, levels = typology_levels)
  )

endpoint_shapes <- c("Morphology" = 21, "Scar direction" = 24)
endpoint_sizes  <- c("Morphology" = 2, "Scar direction" = 1.7)

p_cia_biplot <-
  ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70", linewidth = 0.25) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey70", linewidth = 0.25) +
  geom_segment(
    data = scores_combined %>%
      filter(!is.na(Typology), ID %in% plot_ids_no_biface) %>%
      mutate(Typology = factor(Typology, levels = typology_levels)),
    aes(x = Axis1_M, y = Axis2_M, xend = Axis1_S, yend = Axis2_S,
        color = Typology),
    linewidth = 0.32, alpha = 0.45, lineend = "round"
  ) +
  geom_point(
    data = scores_long_plot,
    aes(x = x, y = y, fill = Typology, color = Typology,
        shape = endpoint, size = endpoint),
    stroke = 0.4, alpha = 0.90
  ) +
  scale_color_manual(values = typology_pal, name = "Typology",
                     breaks = typology_levels) +
  scale_fill_manual(values  = typology_pal, name = "Typology",
                    breaks = typology_levels) +
  scale_shape_manual(values = endpoint_shapes, name = "Endpoint") +
  scale_size_manual(values  = endpoint_sizes,  name = "Endpoint") +
  theme_bw(base_size = 8) +
  labs(
    x = sprintf("CoIA Axis1 (%.1f%%)", cia_inertia[1]),
    y = sprintf("CoIA Axis2 (%.1f%%)", cia_inertia[2])
  ) +
  guides(
    color = guide_legend(order = 1, override.aes = list(shape = 21, size = 2),
                         title = NULL),
    fill  = "none",
    shape = guide_legend(order = 2,
                         override.aes = list(fill = "grey60", color = "grey30",
                                             size = c(2, 1.7)),
                         title = "Endpoint"),
    size  = "none"
  ) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    legend.position   = c(0.01, 0.99),
    legend.justification = c(0, 1),
    legend.background = element_rect(fill = alpha("white", 0.75),
                                     color = "grey80", linewidth = 0.3),
    legend.key.size   = unit(0.32, "cm"),
    legend.text       = element_text(size = 6.5),
    legend.margin     = margin(2, 4, 2, 4)
  )

cat("Figure built: EXP_L1_CIA_Biplot.png\n")

l1_results <- tibble(
  method  = c("Mantel (ILR, Euclidean, Spearman)", "RV (ILR, Euclidean)"),
  stat    = c(mantel_global$statistic, coin_exp$RV),
  p_value = c(mantel_global$signif,    rv_test$pvalue),
  n       = length(exp_ids)
)
write_csv(l1_results, here("analysis/data/derived_data/EXP_L1_results.csv"))
cat("\nLevel 1 summary:\n"); print(l1_results)

# ==============================================================================
# ---- L1-3: CoIA Sankey (ILR -> PCA -> CoIA contribution flow) ----
# ==============================================================================

cat("\n==== L1-3: CoIA Sankey ====\n")

c1_morph <- as.matrix(dudi_morph$c1)
c1_scar  <- as.matrix(dudi_scar$c1)

colnames(c1_morph) <- paste0("MPC", seq_len(ncol(c1_morph)))
colnames(c1_scar)  <- paste0("DPC", seq_len(ncol(c1_scar)))
rownames(c1_morph) <- paste0("Milr", seq_len(nrow(c1_morph)))
rownames(c1_scar)  <- paste0("Dilr", seq_len(nrow(c1_scar)))

w2_ilr_mpc <- c1_morph^2
w2_ilr_dpc <- c1_scar^2

a_morph <- as.matrix(coin_exp$aX)
a_scar  <- as.matrix(coin_exp$aY)

w2_mpc_cia <- a_morph^2
w2_dpc_cia <- a_scar^2

morph_var_pct_full <- round(dudi_morph$eig / sum(dudi_morph$eig) * 100, 1)
scar_var_pct_full  <- round(dudi_scar$eig  / sum(dudi_scar$eig)  * 100, 1)

# Keep axes with cumulative variance <= 95% (at least axis 1)
keep_mpc <- which(cumsum(morph_var_pct_full) <= 95 | seq_along(morph_var_pct_full) == 1)
keep_dpc <- which(cumsum(scar_var_pct_full)  <= 95 | seq_along(scar_var_pct_full)  == 1)
keep_mpc <- keep_mpc[keep_mpc <= ncol(c1_morph)]
keep_dpc <- keep_dpc[keep_dpc <= ncol(c1_scar)]
n_cia_ax <- min(2, ncol(a_morph))

w2_ilr_mpc <- w2_ilr_mpc[, keep_mpc, drop = FALSE]
w2_ilr_dpc <- w2_ilr_dpc[, keep_dpc, drop = FALSE]
w2_mpc_cia <- w2_mpc_cia[keep_mpc, seq_len(n_cia_ax), drop = FALSE]
w2_dpc_cia <- w2_dpc_cia[keep_dpc, seq_len(n_cia_ax), drop = FALSE]

n_ilr_m <- nrow(w2_ilr_mpc)
n_ilr_d <- nrow(w2_ilr_dpc)
n_mpc   <- ncol(w2_ilr_mpc)
n_dpc   <- ncol(w2_ilr_dpc)

cat(sprintf(
  "  morphology: %d ILR -> %d MorphPC; direction: %d ILR -> %d DirPC; shared %d CoIA axes\n",
  n_ilr_m, n_mpc, n_ilr_d, n_dpc, n_cia_ax
))

# Node height (proportional to total weight^2)
h_ilr_m_raw <- rowSums(w2_ilr_mpc)
h_ilr_d_raw <- rowSums(w2_ilr_dpc)
h_mpc_raw   <- rowSums(w2_mpc_cia)
h_dpc_raw   <- rowSums(w2_dpc_cia)
h_cia_raw   <- colSums(w2_mpc_cia) + colSums(w2_dpc_cia)

MORPH_HEIGHT <- 220
SCAR_HEIGHT  <- 220
NODE_GAP     <- 5

scale_nodes <- function(h_raw, total_h, n_nodes, gap = NODE_GAP) {
  avail <- total_h - gap * (n_nodes - 1)
  h_raw * (avail / sum(h_raw))
}

h_ilr_m <- scale_nodes(h_ilr_m_raw, MORPH_HEIGHT, n_ilr_m)
h_mpc   <- scale_nodes(h_mpc_raw,   MORPH_HEIGHT, n_mpc)
h_ilr_d <- scale_nodes(h_ilr_d_raw, SCAR_HEIGHT,  n_ilr_d)
h_dpc   <- scale_nodes(h_dpc_raw,   SCAR_HEIGHT,  n_dpc)
h_cia   <- scale_nodes(h_cia_raw,
                       MORPH_HEIGHT + SCAR_HEIGHT + NODE_GAP * (n_cia_ax - 1),
                       n_cia_ax)

y_start <- 52

y_ilr_m <- y_start + cumsum(c(0, head(h_ilr_m + NODE_GAP, -1)))
y_mpc   <- y_start + cumsum(c(0, head(h_mpc   + NODE_GAP, -1)))

y_ilr_d <- y_start + MORPH_HEIGHT + NODE_GAP +
  cumsum(c(0, head(h_ilr_d + NODE_GAP, -1)))
y_dpc   <- y_start + MORPH_HEIGHT + NODE_GAP +
  cumsum(c(0, head(h_dpc   + NODE_GAP, -1)))

cia_total_h <- sum(h_cia) + NODE_GAP * (n_cia_ax - 1)
two_ends_h  <- MORPH_HEIGHT + NODE_GAP + SCAR_HEIGHT
cia_offset  <- (two_ends_h - cia_total_h) / 2
y_cia <- y_start + cia_offset +
  cumsum(c(0, head(h_cia + NODE_GAP, -1)))

SVG_W  <- 680
SVG_H  <- ceiling(y_start + MORPH_HEIGHT + NODE_GAP + SCAR_HEIGHT + 30)
NODE_W <- 100
x_col1 <- 20
x_col2 <- 240
x_col3 <- 470

col_morph_fill   <- "#B5D4F4"; col_morph_stroke <- "#185FA5"
col_morph_text   <- "#0C447C"; col_morph_band   <- "#85B7EB"
col_scar_fill    <- "#F5C4B3"; col_scar_stroke  <- "#993C1D"
col_scar_text    <- "#712B13"; col_scar_band    <- "#F0997B"
col_cia_fill     <- "#CECBF6"; col_cia_stroke   <- "#534AB7"
col_cia_text     <- "#3C3489"
col_band_cia_m   <- "#AFA9EC"
col_band_cia_d   <- "#C9A8E0"

svg_rect <- function(x, y, w, h, fill, stroke, rx = 5, sw = 0.8)
  sprintf(
    '<rect x="%.1f" y="%.1f" width="%.1f" height="%.1f" rx="%d" fill="%s" stroke="%s" stroke-width="%.1f"/>',
    x, y, w, h, rx, fill, stroke, sw
  )

svg_band <- function(x1, y1t, y1b, x2, y2t, y2b, fill, opacity = 0.42) {
  xm <- (x1 + x2) / 2
  sprintf(
    paste0('<path d="M%.1f,%.1f C%.1f,%.1f %.1f,%.1f %.1f,%.1f ',
           'L%.1f,%.1f C%.1f,%.1f %.1f,%.1f %.1f,%.1f Z" ',
           'fill="%s" fill-opacity="%.2f"/>'),
    x1, y1t, xm, y1t, xm, y2t, x2, y2t,
    x2, y2b, xm, y2b, xm, y1b, x1, y1b,
    fill, opacity
  )
}

svg_text <- function(x, y, txt, size = 11, weight = 400,
                     anchor = "middle", fill = "#2C2C2A")
  sprintf(
    '<text x="%.1f" y="%.1f" font-size="%d" font-weight="%d" text-anchor="%s" fill="%s" dominant-baseline="central">%s</text>',
    x, y, size, weight, anchor, fill, txt
  )

lines <- character(0)
push  <- function(...) lines <<- c(lines, ...)

push(sprintf(
  '<svg xmlns="http://www.w3.org/2000/svg" width="%d" viewBox="0 0 %d %d">',
  SVG_W, SVG_W, SVG_H
))
push('<title>ILR to PCA to CoIA contribution flow (EXP)</title>')
push('<desc>Sankey-style flow diagram showing ILR spectral variable contributions through PCA axes to CoIA axes.</desc>')
push(sprintf('<rect width="%d" height="%d" fill="white"/>', SVG_W, SVG_H))

push(svg_text(x_col1 + NODE_W / 2, 20, "ILR variables",  12, 500))
push(svg_text(x_col1 + NODE_W / 2, 36, "spectral log-contrasts", 10, 400, fill = "#5F5E5A"))
push(svg_text(x_col2 + NODE_W / 2, 20, "PCA axes",        12, 500))
push(svg_text(x_col2 + NODE_W / 2, 36, "per endpoint",    10, 400, fill = "#5F5E5A"))
push(svg_text(x_col3 + NODE_W / 2, 20, "CoIA axes",       12, 500))
push(svg_text(x_col3 + NODE_W / 2, 36, "shared structure",10, 400, fill = "#5F5E5A"))

outlet_ilr_m <- rep(0, n_ilr_m); inlet_mpc  <- rep(0, n_mpc)
outlet_ilr_d <- rep(0, n_ilr_d); inlet_dpc  <- rep(0, n_dpc)
inlet_cia    <- rep(0, n_cia_ax)

for (i in seq_len(n_ilr_m)) for (j in seq_len(n_mpc)) {
  ww <- w2_ilr_mpc[i, j]; if (ww < 0.004) next
  bw_out <- ww * h_ilr_m[i]; bw_in <- ww * h_mpc[j]
  y1t <- y_ilr_m[i] + outlet_ilr_m[i]; y1b <- y1t + bw_out
  y2t <- y_mpc[j]   + inlet_mpc[j];   y2b <- y2t + bw_in
  push(svg_band(x_col1 + NODE_W, y1t, y1b, x_col2, y2t, y2b,
                col_morph_band, 0.28 + 0.32 * ww))
  outlet_ilr_m[i] <- outlet_ilr_m[i] + bw_out
  inlet_mpc[j]    <- inlet_mpc[j]    + bw_in
}

for (i in seq_len(n_ilr_d)) for (j in seq_len(n_dpc)) {
  ww <- w2_ilr_dpc[i, j]; if (ww < 0.004) next
  bw_out <- ww * h_ilr_d[i]; bw_in <- ww * h_dpc[j]
  y1t <- y_ilr_d[i] + outlet_ilr_d[i]; y1b <- y1t + bw_out
  y2t <- y_dpc[j]   + inlet_dpc[j];   y2b <- y2t + bw_in
  push(svg_band(x_col1 + NODE_W, y1t, y1b, x_col2, y2t, y2b,
                col_scar_band, 0.28 + 0.32 * ww))
  outlet_ilr_d[i] <- outlet_ilr_d[i] + bw_out
  inlet_dpc[j]    <- inlet_dpc[j]    + bw_in
}

outlet_mpc <- rep(0, n_mpc)
for (j in seq_len(n_mpc)) for (k in seq_len(n_cia_ax)) {
  ww <- w2_mpc_cia[j, k]; if (ww < 0.004) next
  tot_k  <- colSums(w2_mpc_cia)[k] + colSums(w2_dpc_cia)[k]
  bw_out <- ww * h_mpc[j]
  bw_in  <- ww * h_cia[k] * (colSums(w2_mpc_cia)[k] / tot_k)
  y1t <- y_mpc[j] + outlet_mpc[j]; y1b <- y1t + bw_out
  y2t <- y_cia[k] + inlet_cia[k];  y2b <- y2t + bw_in
  push(svg_band(x_col2 + NODE_W, y1t, y1b, x_col3, y2t, y2b,
                col_band_cia_m, 0.22 + 0.36 * ww))
  outlet_mpc[j] <- outlet_mpc[j] + bw_out
  inlet_cia[k]  <- inlet_cia[k]  + bw_in
}

outlet_dpc <- rep(0, n_dpc)
for (j in seq_len(n_dpc)) for (k in seq_len(n_cia_ax)) {
  ww <- w2_dpc_cia[j, k]; if (ww < 0.004) next
  tot_k  <- colSums(w2_mpc_cia)[k] + colSums(w2_dpc_cia)[k]
  bw_out <- ww * h_dpc[j]
  bw_in  <- ww * h_cia[k] * (colSums(w2_dpc_cia)[k] / tot_k)
  y1t <- y_dpc[j] + outlet_dpc[j]; y1b <- y1t + bw_out
  y2t <- y_cia[k] + inlet_cia[k];  y2b <- y2t + bw_in
  push(svg_band(x_col2 + NODE_W, y1t, y1b, x_col3, y2t, y2b,
                col_band_cia_d, 0.22 + 0.36 * ww))
  outlet_dpc[j] <- outlet_dpc[j] + bw_out
  inlet_cia[k]  <- inlet_cia[k]  + bw_in
}

for (i in seq_len(n_ilr_m)) {
  push(svg_rect(x_col1, y_ilr_m[i], NODE_W, h_ilr_m[i], col_morph_fill, col_morph_stroke))
  push(svg_text(x_col1 + NODE_W / 2, y_ilr_m[i] + h_ilr_m[i] / 2,
                sprintf("ilr%d (morph)", i), 11, 500, fill = col_morph_text))
}
for (i in seq_len(n_ilr_d)) {
  push(svg_rect(x_col1, y_ilr_d[i], NODE_W, h_ilr_d[i], col_scar_fill, col_scar_stroke))
  push(svg_text(x_col1 + NODE_W / 2, y_ilr_d[i] + h_ilr_d[i] / 2,
                sprintf("ilr%d (dir)", i), 11, 500, fill = col_scar_text))
}
for (j in seq_len(n_mpc)) {
  push(svg_rect(x_col2, y_mpc[j], NODE_W, h_mpc[j], col_morph_fill, col_morph_stroke))
  push(svg_text(x_col2 + NODE_W / 2, y_mpc[j] + h_mpc[j] / 2,
                sprintf("Morph-PC%d", j), 11, 500, fill = col_morph_text))
}
for (j in seq_len(n_dpc)) {
  push(svg_rect(x_col2, y_dpc[j], NODE_W, h_dpc[j], col_scar_fill, col_scar_stroke))
  push(svg_text(x_col2 + NODE_W / 2, y_dpc[j] + h_dpc[j] / 2,
                sprintf("Dir-PC%d", j), 11, 500, fill = col_scar_text))
}
for (k in seq_len(n_cia_ax)) {
  push(svg_rect(x_col3, y_cia[k], NODE_W, h_cia[k], col_cia_fill, col_cia_stroke))
  push(svg_text(x_col3 + NODE_W / 2, y_cia[k] + h_cia[k] / 2,
                c("CoIA Axis 1", "CoIA Axis 2")[k], 11, 500, fill = col_cia_text))
}

push("</svg>")

svg_path <- tempfile(fileext = ".svg")
png_path <- here("analysis/output/figures/EXP_L1_CIA_Sankey.png")
writeLines(lines, svg_path, useBytes = FALSE)

if (requireNamespace("rsvg", quietly = TRUE)) {
  rsvg::rsvg_png(svg_path, png_path, width = SVG_W * 2)
  cat("PNG saved (via rsvg): EXP_L1_CIA_Sankey.png\n")
} else if (nzchar(Sys.which("rsvg-convert"))) {
  system2("rsvg-convert", args = c("-d", "300", "-p", "300", "-o", png_path, svg_path))
  cat("PNG saved (via rsvg-convert): EXP_L1_CIA_Sankey.png\n")
} else if (nzchar(Sys.which("inkscape"))) {
  system2("inkscape", args = c("--export-filename", png_path, "--export-dpi", "300", svg_path))
  cat("PNG saved (via Inkscape): EXP_L1_CIA_Sankey.png\n")
} else {
  cat("[warning] rsvg / rsvg-convert / Inkscape not found; cannot render PNG.\n")
  cat("          Install with install.packages('rsvg') and re-run this section.\n")
  cat("          Intermediate SVG at:", svg_path, "\n")
}
unlink(svg_path)
cat("\n==== L1-3 Sankey done ====\n")

# ==============================================================================
# ========== Level 2: joint evidence ==========
# ==============================================================================

cat("\n\n## Level 2: joint evidence (EXP) ##\n")

# ------------------------------------------------------------------------------
# L2-A: grouped Mantel
# ------------------------------------------------------------------------------

cat("\n---------- L2-A: per-type Mantel ----------\n")

mantel_one_type <- function(type_val, meta_df, D_morph_full, D_scar_full,
                            n_perm = 9999) {
  ids <- meta_df %>% filter(Typology == type_val) %>% pull(ID)
  cat(sprintf("  -> %s (n = %d) ...", type_val, length(ids)))
  if (length(ids) < 5) { cat(" skipped (n < 5)\n"); return(NULL) }
  res <- mantel(extract_subdist(D_morph_full, ids),
                extract_subdist(D_scar_full,  ids),
                method = "spearman", permutations = n_perm)
  cat(sprintf(" r = %.4f, p = %.4f\n", res$statistic, res$signif))
  tibble(
    Typology     = type_val,
    n            = length(ids),
    mantel_r     = res$statistic,
    p_value      = res$signif,
    significance = case_when(
      res$signif < 0.001 ~ "***", res$signif < 0.01 ~ "**",
      res$signif < 0.05  ~ "*",   res$signif < 0.10 ~ ".",
      TRUE               ~ "ns"
    )
  )
}

all_types <- meta_exp %>%
  filter(!is.na(Typology)) %>%
  count(Typology) %>%
  arrange(desc(n)) %>%
  pull(Typology)

mantel_by_typology <- map_dfr(all_types,
                              ~ mantel_one_type(.x, meta_exp,
                                                D_morph_exp, D_scar_exp))
mantel_by_typology <- mantel_by_typology %>%
  mutate(group_var = "Typology", group = Typology,
         p_raw = p_value, p_fdr = p_value)

cat("\n==== Per-type Mantel summary ====\n")
print(mantel_by_typology %>%
        select(Typology, n, mantel_r, p_value, significance) %>%
        mutate(across(c(mantel_r, p_value), ~ round(.x, 4))))

l2_mantel <- mantel_by_typology %>%
  mutate(group_var_label = "Typology",
         group = factor(group, levels = typology_levels))

write_csv(l2_mantel, here("analysis/data/derived_data/EXP_L2_grouped_mantel.csv"))
cat("Saved: EXP_L2_grouped_mantel.csv\n")

# ------------------------------------------------------------------------------
# L2-B: CoIA arrow length
# ------------------------------------------------------------------------------

cat("\n---------- L2-B: arrow length (Typology) ----------\n")

run_arrow_length_analysis <- function(group_col, group_label, palette) {
  valid_groups_stat <- scores_combined %>%
    filter(!is.na(.data[[group_col]])) %>%
    group_by(.data[[group_col]]) %>%
    filter(n() >= 3) %>%
    pull(.data[[group_col]]) %>%
    unique()
  if (length(valid_groups_stat) < 2) {
    cat(sprintf("  [skip] %s arrow-length test: too few valid groups\n", group_label))
    return(invisible(NULL))
  }
  sub_df_stat <- scores_combined %>% filter(.data[[group_col]] %in% valid_groups_stat)
  
  cat(sprintf("\n----- %s x arrow length -----\n", group_label))
  kw <- kruskal.test(reformulate(group_col, "arrow_length"), data = sub_df_stat)
  print(kw)
  pw <- pairwise.wilcox.test(sub_df_stat$arrow_length, sub_df_stat[[group_col]],
                             p.adjust.method = "holm", exact = FALSE)
  print(pw)
  
  valid_groups_plot <- intersect(TYPOLOGY_ORDER,
                                 valid_groups_stat[valid_groups_stat != "Biface"])
  sub_df <- scores_combined %>%
    filter(.data[[group_col]] %in% valid_groups_plot,
           ID %in% plot_ids_no_biface) %>%
    mutate(!!group_col := factor(.data[[group_col]], levels = valid_groups_plot))
  
  p <- ggplot(sub_df,
              aes(x = .data[[group_col]], y = arrow_length,
                  fill = .data[[group_col]], color = .data[[group_col]])) +
    geom_boxplot(outlier.shape = 21, outlier.size = 1.6,
                 alpha = 0.25, linewidth = 0.35) +
    geom_jitter(width = 0.15, size = 1.6, alpha = 0.7, shape = 16) +
    stat_summary(fun = mean, geom = "point",
                 shape = 16, size = 2.4, color = "white") +
    annotate("text", x = Inf, y = Inf,
             label = sprintf("Kruskal-Wallis\nchi² = %.2f, P = %.3f",
                             kw$statistic, kw$p.value),
             hjust = 1.05, vjust = 1.2, size = 2.6, color = "grey40") +
    scale_fill_manual(values  = palette) +
    scale_color_manual(values = palette) +
    scale_x_discrete(
      expand = expansion(add = 0.6),
      labels = c(
        "Unidirectional" = "Uni.",
        "Bidirectional"  = "Bi.",
        "Levallois"      = "Lev.",
        "Discoid"        = "Dis.",
        "Multiplatform"  = "Multi."
      )) +
    theme_bw(base_size = 8) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor   = element_blank(),
      axis.text.x        = element_text(size = 7),
      axis.text.y        = element_text(size = 7),
      legend.position    = "none"
    ) +
    labs(
      x = NULL, y = "CoIA line length"
    )
  
  fname <- sprintf("analysis/output/figures/EXP_L2_Arrow_Length_%s.png",
                   tolower(str_replace_all(group_label, " ", "_")))
  cat(sprintf("Figure built: %s\n", basename(fname)))
  list(sub_df = sub_df, kw = kw, pw = pw, p = p)
}

res_len_typology <- run_arrow_length_analysis("Typology", "Typology", typology_pal)

cat("\n==== Arrow-length summary (Typology) ====\n")
scores_combined %>%
  filter(!is.na(Typology)) %>%
  group_by(Typology) %>%
  summarise(
    n      = n(),
    mean   = round(mean(arrow_length),   4),
    median = round(median(arrow_length), 4),
    sd     = round(sd(arrow_length),     4),
    min    = round(min(arrow_length),    4),
    max    = round(max(arrow_length),    4),
    .groups = "drop"
  ) %>%
  print()

# ------------------------------------------------------------------------------
# L2-C: CoIA arrow-direction circular statistics
# NOTE: plot_rose() must be defined before run_circular_analysis(),
#       otherwise the list stores a function reference instead of a ggplot object.
# ------------------------------------------------------------------------------

cat("\n---------- L2-C: arrow direction (Typology) ----------\n")

# ---- plotting helper ----
plot_rose <- function(res, palette, bw = 40) {
  group_col <- res$group_col
  kde_df    <- res$kde_df
  mean_dirs <- res$mean_dirs
  
  # Map angle from [0, 360) to (-180, 180] for linear display
  kde_linear <- kde_df %>%
    mutate(angle_centered = if_else(angle_deg > 180,
                                    angle_deg - 360,
                                    angle_deg),
           !!group_col := factor(.data[[group_col]],
                                 levels = levels(kde_df[[group_col]])))
  
  mean_linear <- mean_dirs %>%
    mutate(angle_centered = if_else(mean_deg > 180,
                                    mean_deg - 360,
                                    mean_deg),
           !!group_col := factor(.data[[group_col]],
                                 levels = levels(kde_df[[group_col]])))
  
  # Rayleigh label
  rayleigh_labels <- res$rayleigh %>%
    mutate(
      label = case_when(
        rayleigh_p < 0.001 ~ "Rayleigh\nP < 0.001",
        rayleigh_p < 0.05  ~ sprintf("Rayleigh\nP = %.3f", rayleigh_p),
        TRUE               ~ sprintf("Rayleigh\nP = %.3f", rayleigh_p)
      ),
      !!group_col := factor(group, levels = levels(kde_df[[group_col]]))
    )
  
  ggplot(kde_linear,
         aes(x     = angle_centered,
             y     = density,
             fill  = .data[[group_col]],
             color = .data[[group_col]])) +
    geom_area(alpha = 0.4, color = NULL) +
    geom_vline(
      data     = mean_linear,
      aes(xintercept = angle_centered,
          color      = .data[[group_col]]),
      linewidth = 0.4, linetype = "dashed", alpha = 0.75
    ) +
    # Rayleigh significance label
    geom_text(
      data = rayleigh_labels,
      aes(label = label),
      x = 170, y = Inf,
      hjust = 0.8, vjust = 1.4,
      size = 2.2, color = "grey35",
      inherit.aes = FALSE
    ) +
    scale_x_continuous(
      limits = c(-180, 180),
      breaks = seq(-180, 180, by = 45),
      labels = seq(-180, 180, by = 45)
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
    scale_fill_manual(values  = palette) +
    scale_color_manual(values = palette, guide = "none") +
    facet_wrap(reformulate(group_col),
               nrow = 1,
               scales = "free_y") +
    labs(x = "CoIA line direction (°)",
         y = "von Mises KDE",
         fill = group_col) +
    theme_bw(base_size = 8) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_blank(),
      strip.text       = element_text(face = "bold", size = 7),
      strip.background = element_rect(fill = "#EBEBEB", color = "#EBEBEB"),
      axis.text.x      = element_text(size = 6),
      axis.text.y      = element_blank(),
      axis.ticks.y     = element_blank(),
      legend.position  = "none"
    )
}

# ---- analysis helper (defined after plot_rose, which it calls) ----
run_circular_analysis <- function(group_col, group_label, palette) {
  valid_groups <- scores_combined %>%
    filter(!is.na(.data[[group_col]])) %>%
    group_by(.data[[group_col]]) %>%
    filter(n() >= 5) %>%
    pull(.data[[group_col]]) %>%
    unique()
  if (length(valid_groups) < 2) {
    cat(sprintf("  [skip] %s circular statistics: too few valid groups\n", group_label))
    return(invisible(NULL))
  }
  valid_groups_ordered <- intersect(TYPOLOGY_ORDER,
                                    valid_groups[valid_groups != "Biface"])
  sub_df <- scores_combined %>% filter(.data[[group_col]] %in% valid_groups)
  
  cat(sprintf("\n----- %s circular summary -----\n", group_label))
  circ_desc <- map_dfr(valid_groups_ordered, function(g) {
    angles <- sub_df %>% filter(.data[[group_col]] == g) %>% pull(arrow_angle)
    cs     <- circ_stats_one(angles)
    tibble(group_var = group_col, group = g, n = length(angles),
           mean_dir_deg = round(cs$mean_deg, 2),
           concentration_r = round(cs$rho, 4))
  })
  print(circ_desc)
  
  cat(sprintf("\n----- %s Rayleigh test -----\n", group_label))
  rayleigh_res <- map_dfr(valid_groups_ordered, function(g) {
    angles   <- sub_df %>% filter(.data[[group_col]] == g) %>% pull(arrow_angle)
    circ_obj <- circular(angles, type = "angles", units = "radians", modulo = "2pi")
    rt       <- rayleigh.test(circ_obj)
    cat(sprintf("  %s: U = %.4f, p = %.4f -> %s\n",
                g, rt$statistic, rt$p.value,
                ifelse(rt$p.value < 0.05, "concentrated", "dispersed")))
    tibble(group_var = group_col, group = g,
           rayleigh_U = round(rt$statistic, 4),
           rayleigh_p = round(rt$p.value,   4),
           conclusion = ifelse(rt$p.value < 0.05, "concentrated", "uniform"))
  })
  
  # Prepare plotting data
  sub_df_plot <- scores_combined %>%
    filter(.data[[group_col]] %in% valid_groups_ordered,
           ID %in% plot_ids_no_biface)
  
  mean_dirs <- map_dfr(valid_groups_ordered, function(g) {
    angles <- sub_df_plot %>% filter(.data[[group_col]] == g) %>% pull(arrow_angle)
    cs     <- circ_stats_one(angles)
    tibble(!!group_col := g, mean_deg = cs$mean_deg)
  }) %>%
    mutate(!!group_col := factor(.data[[group_col]], levels = valid_groups_ordered))
  
  compute_circular_kde <- function(angles_deg, bw = 25, n = 360) {
    circ <- circular(angles_deg * pi / 180,
                     type = "angles", units = "radians", modulo = "2pi")
    dens <- density(circ, bw = bw, n = n)
    angle_deg <- as.numeric(dens$x) * 180 / pi %% 360
    density   <- as.numeric(dens$y)
    tibble(
      angle_deg = c(angle_deg, angle_deg[1]),
      density   = c(density,   density[1])
    )
  }
  
  kde_df <- sub_df_plot %>%
    mutate(
      angle_deg = arrow_angle * 180 / pi,
      angle_deg = ifelse(angle_deg < 0, angle_deg + 360, angle_deg),
      !!group_col := factor(.data[[group_col]], levels = valid_groups_ordered)
    ) %>%
    group_by(.data[[group_col]]) %>%
    group_modify(~ compute_circular_kde(.x$angle_deg, bw = 25)) %>%
    ungroup() %>%
    mutate(!!group_col := factor(.data[[group_col]], levels = valid_groups_ordered))
  
  # Call plot_rose to build the ggplot object
  p_rose <- plot_rose(
    list(group_col = group_col, kde_df = kde_df, mean_dirs = mean_dirs,
         rayleigh  = rayleigh_res),
    palette
  )
  
  fname <- sprintf("analysis/output/figures/EXP_L2_Arrow_Direction_rose_%s.png",
                   tolower(str_replace_all(group_label, " ", "_")))
  cat(sprintf("Figure built: %s\n", basename(fname)))
  
  list(
    desc      = circ_desc,
    rayleigh  = rayleigh_res,
    kde_df    = kde_df,
    mean_dirs = mean_dirs,
    group_col = group_col,
    p_rose    = p_rose   
  )
}

res_circ_typology <- run_circular_analysis("Typology", "Typology", typology_pal)

p_rose <- plot_rose(res_circ_typology, typology_pal)

# ==============================================================================
# ---- Save derived data ----
# ==============================================================================

if (!is.null(res_circ_typology)) {
  circ_out <- left_join(res_circ_typology$desc, res_circ_typology$rayleigh,
                        by = c("group_var", "group"))
  write_csv(circ_out, here("analysis/data/derived_data/EXP_L2_circular_stats.csv"))
  cat("Saved: EXP_L2_circular_stats.csv\n")
}

scores_combined %>%
  select(ID, Typology, arrow_length, arrow_angle,
         Axis1_M, Axis2_M, Axis1_S, Axis2_S) %>%
  write_csv(here("analysis/data/derived_data/EXP_L2_arrow_stats.csv"))
cat("Saved: EXP_L2_arrow_stats.csv\n")

scores_combined %>%
  write_csv(here("analysis/data/derived_data/EXP_CIA_scores_full.csv"))
cat("Saved: EXP_CIA_scores_full.csv\n")

# Write ILR scores
morph_ilr_exp %>%
  rownames_to_column("ID") %>%
  write_csv(here("analysis/data/derived_data/EXP_morph_ILR_scores.csv"))
cat("Saved: EXP_morph_ILR_scores.csv\n")

scar_ilr_exp %>%
  rownames_to_column("ID") %>%
  write_csv(here("analysis/data/derived_data/EXP_scar_ILR_scores.csv"))
cat("Saved: EXP_scar_ILR_scores.csv\n")

# ==============================================================================
# ---- Summary print ----
# ==============================================================================

cat("\n\n## Analysis summary (EXP) ##\n")

cat("\n[Level 1: baseline]\n")
cat(sprintf("  Mantel r = %.4f, p = %.3f  ->  %s\n",
            mantel_global$statistic, mantel_global$signif,
            ifelse(mantel_global$signif < 0.05, "correlated", "independent (n.s.)")))
cat(sprintf("  RV       = %.4f, p = %.3f  ->  %s\n",
            coin_exp$RV, rv_test$pvalue,
            ifelse(rv_test$pvalue < 0.05, "covarying", "independent (n.s.)")))

cat("\n[Level 2A: per-type Mantel (Typology)]\n")
if (!is.null(l2_mantel) && nrow(l2_mantel) > 0) {
  n_sig <- sum(l2_mantel$p_value < 0.05, na.rm = TRUE)
  cat(sprintf("  %d types tested, p < 0.05: %d\n", nrow(l2_mantel), n_sig))
  if (n_sig > 0) {
    l2_mantel %>% filter(p_value < 0.05) %>%
      select(Typology, n, mantel_r, p_value, significance) %>%
      mutate(across(c(mantel_r, p_value), ~ round(.x, 4))) %>%
      print()
  }
}

cat("\n[Level 2B] see EXP_L2_Arrow_Length_typology.png\n")
cat("[Level 2C] see EXP_L2_Arrow_Direction_rose_typology.png and EXP_L2_circular_stats.csv\n")

cat("\n[Sankey] EXP_L1_CIA_Sankey.png\n")

# ==============================================================================
# ---- Composite figure ----
# ==============================================================================

# Step 1: tag each subplot
p_cia_biplot_tagged <- p_cia_biplot +
  labs(tag = "a") +
  theme(plot.tag = element_text(size = 9, face = "bold"))

p_len_tagged <- res_len_typology$p +
  labs(tag = "b") +
  theme(plot.tag = element_text(size = 9, face = "bold"))

p_rose_tagged <- res_circ_typology$p_rose +
  labs(tag = "c") +
  theme(plot.tag = element_text(size = 9, face = "bold"))

# Step 2: recompose p_composite (no plot_annotation)
p_composite <- (
  ((p_cia_biplot_tagged | p_len_tagged) + plot_layout(widths = c(3, 1))) /
    p_rose_tagged
) +
  plot_layout(heights = c(2.8, 1))

# Step 3: tag the external image as D
external_img <- png::readPNG(here("analysis/figures/Axis_trajectory.png"))
# Place the image so its frame spans the same left/right (43..996 px of a
# 1027 px-wide render) as the panels above: width = 0.949 npc, centred at
# x = 0.5055 npc of the full-width cell. Width-constrained, so the right edge
# stays aligned regardless of the PNG's exact aspect ratio (target W:H ~= 2.95).
grob_img     <- grid::rasterGrob(external_img, interpolate = TRUE,
                                 x     = grid::unit(0.5055, "npc"),
                                 width = grid::unit(0.949,  "npc"))
p_external   <- wrap_elements(full = grob_img) +
  labs(tag = "d") +
  theme(
    plot.tag    = element_text(size = 9, face = "bold"),
    plot.margin = margin(0, 0, 0, 0)
  )

# Step 4: assemble
p_final <- wrap_elements(full = p_composite) / p_external +
  plot_layout(heights = c(3, 1))
