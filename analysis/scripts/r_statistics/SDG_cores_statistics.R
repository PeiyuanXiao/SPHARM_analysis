# SDG_cores_statistics.R
# Techno-morphological joint analysis of the SDG archaeological cores.
# Sourced by `sdg_cia_analysis`.
#
# Three levels:
#   Level 1  global Mantel + CoIA (baseline), with a Sankey of the contribution flow.
#   Level 2  grouped Mantel + arrow-length + arrow-direction differences,
#            by layer / raw material / core type.
#   Level 3  PERMANOVA of group structure (morphology and scar direction separately).
#
# Input:
#   - analysis/data/derived_data/SPHARM_direction.csv
#   - analysis/data/derived_data/SPHARM_morphology.csv
#   - analysis/data/raw_data/SDG_core_metric.xlsx
#   - analysis/figures/source_panels/Axis_trajectory_SDG.png (external panel of the composite figure)
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

# ==============================================================================
# ---- Global constants: palettes and order ----
# ==============================================================================

LAYER_ORDER <- c("Layer 2", "Layer 3", "Layer 4")

CORETYPE_ORDER <- c(
  "Unifacial_unidirection",
  "Unifacial_centripetal",
  "Bifacial_adjacent",
  "Bifacial_independent",
  "Bifacial_centripetal",
  "Multifacial",
  "Core_on_flake"
)

# Excluded core types
EXCLUDE_CORE_TYPES <- c("Handaxe", "Pick")

layer_pal <- c(
  "Layer 2" = "#802520",
  "Layer 3" = "#5C7F71",
  "Layer 4" = "#BA8530"
)

rawmat_pal <- c(
  "chert"     = "#4A6E8A",
  "sandstone" = "#802520"
)

coretype_pal <- c(
  "Unifacial_unidirection" = "#5C7F71",
  "Unifacial_centripetal"  = "#BA8530",
  "Bifacial_adjacent"      = "#802520",
  "Bifacial_independent"   = "#B26538",
  "Bifacial_centripetal"   = "#788C4A",
  "Multifacial"            = "#8A7A68",
  "Core_on_flake"          = "#4A6E8A"
)

# ==============================================================================
# ---- Global helpers ----
# ==============================================================================

replace_zeros <- function(X, delta = NULL) {
  X <- as.matrix(X)
  for (i in seq_len(nrow(X))) {
    row_i      <- X[i, ]
    zero_idx   <- row_i == 0
    if (!any(zero_idx)) next
    nonzero_min <- min(row_i[!zero_idx])
    d           <- ifelse(is.null(delta), nonzero_min * 0.65, delta)
    n_zero      <- sum(zero_idx)
    row_i[zero_idx]  <- d
    row_i[!zero_idx] <- row_i[!zero_idx] * (1 - n_zero * d)
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

# Axial circular statistics: CoIA line directions are undirected (an axis,
# theta equivalent to theta+180). Angles are doubled before averaging and the
# mean resultant is halved back; rho is the axial mean resultant length
# (Mardia & Jupp 2000).
circ_stats_one <- function(angles_rad) {
  circ_obj <- circular(2 * angles_rad, type = "angles",
                       units = "radians", modulo = "2pi")
  mean_rad <- (as.numeric(mean.circular(circ_obj)) %% (2 * pi)) / 2
  list(
    mean_rad = mean_rad,
    mean_deg = mean_rad * 180 / pi,
    rho      = as.numeric(rho.circular(circ_obj))
  )
}

# ==============================================================================
# ---- Data preparation (shared by all levels) ----
# ==============================================================================

POWER_COLS_DIR   <- paste0("power_l", 1:6)   # direction spectrum l=1-6
POWER_COLS_MORPH <- paste0("power_l", 1:8)   # morphology spectrum l=1-8

SPHARM_direction_filter  <- readRDS(here("analysis/data/derived_data/SPHARM_direction_filter.rds"))
SPHARM_morphology_filter <- readRDS(here("analysis/data/derived_data/SPHARM_morphology_filter.rds"))

df_morph_raw <- SPHARM_morphology_filter
df_scar_raw  <- SPHARM_direction_filter

common_ids <- intersect(df_morph_raw$ID, df_scar_raw$ID)
df_morph_raw <- df_morph_raw %>% filter(ID %in% common_ids) %>% arrange(ID)
df_scar_raw  <- df_scar_raw  %>% filter(ID %in% common_ids) %>% arrange(ID)

cat("==== Data alignment ====\n")
cat("Shared specimens:", length(common_ids),
    "; IDs fully matched:", all(df_morph_raw$ID == df_scar_raw$ID), "\n\n")

# External metadata
core_meta_raw <- read_excel(
  here("analysis/data/raw_data/SDG_core_metric.xlsx")
)

cat("\n==== External table diagnostics ====\n")
cat("Columns:"); print(colnames(core_meta_raw))
cat("\nFirst 3 rows:\n"); print(head(core_meta_raw, 3))

core_meta <- core_meta_raw %>%
  select(
    ID           = ID,
    raw_material = Raw_mat,
    core_type    = Core_type_Li_merged
  ) %>%
  mutate(across(everything(), ~ str_trim(as.character(.))))

# Extract power-feature matrices (direction and morphology use their own column counts)
morph_power <- df_morph_raw %>%
  select(all_of(POWER_COLS_MORPH)) %>%
  rename_with(~ paste0("M", seq_along(.))) %>%
  as.data.frame()
scar_power <- df_scar_raw %>%
  select(all_of(POWER_COLS_DIR)) %>%
  rename_with(~ paste0("S", seq_along(.))) %>%
  as.data.frame()
rownames(morph_power) <- df_morph_raw$ID
rownames(scar_power)  <- df_scar_raw$ID

morph_power_clean <- morph_power[, sapply(morph_power, sd, na.rm = TRUE) > 0]
scar_power_clean  <- scar_power[,  sapply(scar_power,  sd, na.rm = TRUE) > 0]

# ILR transform -> Euclidean distance
morph_ilr_all <- as.data.frame(ilr(replace_zeros(as.matrix(morph_power_clean))))
scar_ilr_all  <- as.data.frame(ilr(replace_zeros(as.matrix(scar_power_clean))))
rownames(morph_ilr_all) <- rownames(morph_power_clean)
rownames(scar_ilr_all)  <- rownames(scar_power_clean)

D_morph_all <- dist(morph_ilr_all)
D_scar_all  <- dist(scar_ilr_all)

# Separate archaeological specimens (drop IM_ references)
arch_ids <- rownames(morph_power_clean)[!str_starts(rownames(morph_power_clean), "IM_") &
                                          !str_starts(rownames(morph_power_clean), "EXP")]
cat("Archaeological specimens (excluding IM_):", length(arch_ids), "\n")

morph_arch     <- morph_power_clean[arch_ids, ]
scar_arch      <- scar_power_clean[arch_ids, ]
morph_ilr_arch <- morph_ilr_all[arch_ids, ]
scar_ilr_arch  <- scar_ilr_all[arch_ids, ]

D_morph_arch <- extract_subdist(D_morph_all, arch_ids)
D_scar_arch  <- extract_subdist(D_scar_all,  arch_ids)

# Build metadata
meta_arch <- tibble(ID = arch_ids) %>%
  mutate(
    layer = case_when(
      str_detect(ID, "L2") ~ "Layer 2",
      str_detect(ID, "L3") ~ "Layer 3",
      str_detect(ID, "L4") ~ "Layer 4",
      TRUE                  ~ "Other"
    )
  ) %>%
  left_join(core_meta, by = "ID")

# Drop Handaxe and Pick
meta_arch <- meta_arch %>%
  filter(!core_type %in% EXCLUDE_CORE_TYPES | is.na(core_type))

# Update arch_ids to reflect the filtered set
arch_ids <- meta_arch$ID
morph_arch     <- morph_power_clean[arch_ids, ]
scar_arch      <- scar_power_clean[arch_ids, ]
morph_ilr_arch <- morph_ilr_all[arch_ids, ]
scar_ilr_arch  <- scar_ilr_all[arch_ids, ]
D_morph_arch   <- extract_subdist(D_morph_all, arch_ids)
D_scar_arch    <- extract_subdist(D_scar_all,  arch_ids)

cat("Archaeological specimens after dropping Handaxe/Pick:", length(arch_ids), "\n")

cat("\n==== Metadata merge diagnostics ====\n")
cat("Layer:\n");        print(table(meta_arch$layer,        useNA = "ifany"))
cat("Raw material:\n"); print(table(meta_arch$raw_material,  useNA = "ifany"))
cat("Core type:\n");    print(table(meta_arch$core_type,     useNA = "ifany"))
unmatched <- meta_arch %>% filter(is.na(raw_material)) %>% pull(ID)
if (length(unmatched) > 0) {
  cat("Unmatched IDs (check formatting):\n"); print(unmatched)
}

meta_layer    <- safe_filter_groups(meta_arch %>% filter(layer != "Other"), "layer")
meta_rawmat   <- safe_filter_groups(meta_arch, "raw_material")
meta_coretype <- safe_filter_groups(meta_arch, "core_type")

# ==============================================================================
# ========== Level 1: global Mantel + CoIA ==========
# ==============================================================================

cat("\n\n")
cat("##  Level 1: global Mantel + CoIA — baseline                  ##\n")

# ------------------------------------------------------------------------------
# L1-1: global Mantel + linkET network
# ------------------------------------------------------------------------------

cat("\n==== L1-1: global Mantel ====\n")
mantel_global <- mantel(D_morph_arch, D_scar_arch,
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

morph_arch_df <- as.data.frame(morph_arch)
scar_arch_df  <- as.data.frame(scar_arch)

mantel_cross_l1 <- bind_rows(
  map_dfr(colnames(morph_arch_df),
          ~ run_cross_mantel(morph_arch_df[[.x]], D_scar_arch,
                             "Scar Direction", .x)),
  map_dfr(colnames(scar_arch_df),
          ~ run_cross_mantel(scar_arch_df[[.x]], D_morph_arch,
                             "Morphology", .x))
) %>%
  mutate(
    p_holm       = p.adjust(p, method = "holm"),
    significance = ifelse(p_holm < 0.05, "p\u22640.05", "p>0.05")
  )

spec_arch_full <- bind_cols(morph_arch_df, scar_arch_df)

p_mantel_net <- qcorrplot(
  correlate(spec_arch_full, method = "spearman"),
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
    values = c("p\u22640.05" = "#E6A5A5", "p>0.05" = "#BABABA"),
    name   = "Mantel test\n(Holm corrected)"
  ) +
  scale_size_continuous(range = c(0.5, 2.5), name = "Mantel's |r|") +
  theme_minimal(base_size = 10) +
  theme(panel.grid = element_blank(), axis.title = element_blank(),
        legend.position = "right",
        plot.margin = margin(20, 20, 20, 20))

cat("Figure built: L1_Mantel_Network.png\n")

# ------------------------------------------------------------------------------
# L1-2: CoIA + RV permutation test
# ------------------------------------------------------------------------------

cat("\n==== L1-2：CoIA ====\n")

morph_arch_ilr <- morph_ilr_arch
scar_arch_ilr  <- scar_ilr_arch
colnames(morph_arch_ilr) <- paste0("M_ilr", seq_len(ncol(morph_arch_ilr)))
colnames(scar_arch_ilr)  <- paste0("S_ilr", seq_len(ncol(scar_arch_ilr)))

dudi_morph <- dudi.pca(morph_arch_ilr, center = TRUE, scale = TRUE,
                       scannf = FALSE, nf = ncol(morph_arch_ilr))
dudi_scar  <- dudi.pca(scar_arch_ilr,  center = TRUE, scale = TRUE,
                       scannf = FALSE, nf = ncol(scar_arch_ilr))

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
    sprintf("log(geomean(%s) / %s)",
            paste(num_ids, collapse = "+"), den_id)
  })
  print(load_df)
  
  cat("\n-- Score summary (first 2 axes) --\n")
  score_df <- as.data.frame(dudi_obj$li)[, 1:min(2, n_ax), drop = FALSE]
  colnames(score_df) <- paste0("PC", seq_len(ncol(score_df)))
  score_stats <- score_df %>%
    pivot_longer(everything(), names_to = "Axis", values_to = "Score") %>%
    group_by(Axis) %>%
    summarise(
      mean = round(mean(Score), 4), sd = round(sd(Score), 4),
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

coin_arch   <- coinertia(dudi_morph, dudi_scar, scannf = FALSE, nf = 2)
cia_inertia <- coin_arch$eig / sum(coin_arch$eig) * 100

cat("RV coefficient:", round(coin_arch$RV, 4), "\n")

set.seed(42)
rv_test <- randtest(coin_arch, nrepet = 9999)
cat("\nRV permutation test:\n"); print(rv_test)

scores_morph <- as.data.frame(coin_arch$lX) %>% rownames_to_column("ID")
scores_scar  <- as.data.frame(coin_arch$lY) %>% rownames_to_column("ID")

scores_combined <- left_join(
  scores_morph %>% select(ID, Axis1_M = AxcX1, Axis2_M = AxcX2),
  scores_scar  %>% select(ID, Axis1_S = AxcY1, Axis2_S = AxcY2),
  by = "ID"
) %>%
  mutate(
    arrow_length = sqrt((Axis1_M - Axis1_S)^2 + (Axis2_M - Axis2_S)^2),
    arrow_angle  = atan2(Axis2_S - Axis2_M, Axis1_S - Axis1_M)
  ) %>%
  left_join(meta_arch %>% select(ID, layer, raw_material, core_type),
            by = "ID") %>%
  filter(!str_starts(ID, "EXP"))

# CoIA coordinate output
cat("\n==== CoIA specimen coordinates ====\n")
cia_coords <- scores_combined %>%
  select(ID, layer, raw_material, core_type,
         Morph_Axis1 = Axis1_M, Morph_Axis2 = Axis2_M,
         Scar_Axis1  = Axis1_S, Scar_Axis2  = Axis2_S,
         arrow_length, arrow_angle)
print(cia_coords %>%
        mutate(across(where(is.numeric), ~ round(.x, 4))) %>%
        as.data.frame())
write_csv(cia_coords,
          here("analysis/data/derived_data/CoIA_coords_full.csv"))
cat("Saved: CoIA_coords_full.csv\n")

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
  coin_arch$aX,
  round(dudi_morph$eig / sum(dudi_morph$eig) * 100, 1),
  "Morphology"
)
scar_contrib <- compute_pca_cia_contribution(
  coin_arch$aY,
  round(dudi_scar$eig / sum(dudi_scar$eig) * 100, 1),
  "Scar direction"
)

pca_cia_contrib <- bind_rows(morph_contrib, scar_contrib)
cat("\nMorphology side:\n"); print(morph_contrib %>% select(-endpoint) %>% as.data.frame())
cat("\nDirection side:\n"); print(scar_contrib  %>% select(-endpoint) %>% as.data.frame())
write_csv(pca_cia_contrib,
          here("analysis/data/derived_data/PCA_CoIA_contribution.csv"))
cat("Saved: PCA_CoIA_contribution.csv\n")

# ------------------------------------------------------------------------------
# CoIA diagnostic plots
# ------------------------------------------------------------------------------

eig_df <- tibble(
  axis       = paste0("Axis ", seq_along(coin_arch$eig)),
  eigenvalue = coin_arch$eig,
  pct        = coin_arch$eig / sum(coin_arch$eig) * 100,
  cum_pct    = cumsum(pct)
)

p_scree <- ggplot(eig_df, aes(x = axis, y = pct)) +
  geom_col(fill = "#5C7F71", alpha = 0.85, width = 0.55) +
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
    title    = "CoIA Scree Plot (SDG)",
    subtitle = sprintf("RV = %.3f, p = %.3f", coin_arch$RV, rv_test$pvalue),
    x = "CoIA Axis"
  ) +
  theme(plot.title    = element_text(face = "bold", hjust = 0.5, size = 11),
        plot.subtitle = element_text(hjust = 0.5, size = 8.5, color = "grey50"))

morph_pct <- round(dudi_morph$eig / sum(dudi_morph$eig) * 100, 1)
scar_pct  <- round(dudi_scar$eig  / sum(dudi_scar$eig)  * 100, 1)

morph_load <- as.data.frame(coin_arch$aX) %>%
  rownames_to_column("variable") %>%
  rename(Axis1 = AxcX1, Axis2 = AxcX2) %>%
  mutate(
    pct   = morph_pct[as.integer(str_extract(variable, "[0-9]+"))],
    variable_label = sprintf("Morph-PCA%s\n(%.1f%% var)",
                             str_extract(variable, "[0-9]+"), pct),
    endpoint = "Morphology"
  )

scar_load <- as.data.frame(coin_arch$aY) %>%
  rownames_to_column("variable") %>%
  rename(Axis1 = AxcY1, Axis2 = AxcY2) %>%
  mutate(
    pct   = scar_pct[as.integer(str_extract(variable, "[0-9]+"))],
    variable_label = sprintf("Dir-PCA%s\n(%.1f%% var)",
                             str_extract(variable, "[0-9]+"), pct),
    endpoint = "Scar Direction"
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
    labs(title = title_str,
         x = sprintf("CoIA Axis 1 (%.1f%%)", cia_inertia[1]),
         y = sprintf("CoIA Axis 2 (%.1f%%)", cia_inertia[2])) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 10),
          panel.grid = element_blank())
}

p_load_morph <- make_loading_plot(morph_load,
                                  "Morphology PCA axes on CoIA space", "#5C7F71")
p_load_scar  <- make_loading_plot(scar_load,
                                  "Scar direction PCA axes on CoIA space", "#BA8530")

p_cia_diagnostics <- (p_scree | p_load_morph | p_load_scar) +
  plot_annotation(
    title   = "CoIA Axis Diagnostics (SDG)",
    caption = paste(
      "Left: scree plot. Middle: morphology PCA axis loadings on CoIA axes.",
      "Right: scar-direction PCA axis loadings on CoIA axes.",
      "\nBoth panels use PCA-axis projections (aX / aY) for symmetric interpretation.",
      "Arrow length = contribution to CoIA structure."
    ),
    theme = theme(
      plot.title   = element_text(face = "bold", hjust = 0.5, size = 12),
      plot.caption = element_text(size = 7.5, color = "grey50", hjust = 0)
    )
  )

cat("Figure built: L1_CoIA_Diagnostics.png\n")

# ------------------------------------------------------------------------------
# CoIA biplot helper
# show_color_legend: TRUE = standalone shows colour legend; FALSE = composite use, legend moved to direction plot
# ------------------------------------------------------------------------------
make_coia_biplot <- function(group_col, group_label, palette,
                             group_order = NULL, fname_tag = NULL,
                             show_color_legend = TRUE) {
  valid_groups <- scores_combined %>%
    filter(!is.na(.data[[group_col]]), .data[[group_col]] != "Other") %>%
    pull(.data[[group_col]]) %>% unique()
  if (length(valid_groups) < 1) {
    cat(sprintf("  [skip] %s CoIA biplot: no valid groups\n", group_label))
    return(invisible(NULL))
  }
  if (!is.null(group_order)) {
    lvls <- intersect(group_order, valid_groups)
  } else {
    lvls <- sort(valid_groups)
  }
  
  sub_seg <- scores_combined %>%
    filter(!is.na(.data[[group_col]]), .data[[group_col]] != "Other") %>%
    mutate(!!group_col := factor(.data[[group_col]], levels = lvls))
  
  scores_long <- bind_rows(
    sub_seg %>%
      select(ID, !!group_col, x = Axis1_M, y = Axis2_M, arrow_length) %>%
      mutate(endpoint = "Morphology"),
    sub_seg %>%
      select(ID, !!group_col, x = Axis1_S, y = Axis2_S, arrow_length) %>%
      mutate(endpoint = "Scar direction")
  ) %>%
    mutate(
      endpoint     = factor(endpoint, levels = c("Morphology", "Scar direction")),
      !!group_col := factor(.data[[group_col]], levels = lvls)
    )
  
  endpoint_shapes <- c("Morphology" = 21, "Scar direction" = 24)
  endpoint_sizes  <- c("Morphology" = 1.2, "Scar direction" = 1.3)
  
  p <- ggplot() +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey70", linewidth = 0.25) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey70", linewidth = 0.25) +
    geom_segment(
      data = sub_seg,
      aes(x = Axis1_M, y = Axis2_M, xend = Axis1_S, yend = Axis2_S,
          color = .data[[group_col]]),
      linewidth = 0.32, alpha = 0.45, lineend = "round"
    ) +
    geom_point(
      data = scores_long,
      aes(x = x, y = y,
          fill  = .data[[group_col]],
          color = .data[[group_col]],
          shape = endpoint,
          size  = endpoint),
      stroke = 0.4, alpha = 0.90
    ) +
    scale_color_manual(values = palette, name = group_label, breaks = lvls,
                       labels = function(x) gsub("_", " ", x)) +
    scale_fill_manual(values  = palette, name = group_label, breaks = lvls,
                      labels = function(x) gsub("_", " ", x)) +
    scale_shape_manual(values = endpoint_shapes, name = "Endpoint") +
    scale_size_manual(values  = endpoint_sizes,  name = "Endpoint") +
    theme_bw(base_size = 8) +
    labs(
      x = sprintf("CoIA Axis 1 (%.1f%%)", cia_inertia[1]),
      y = sprintf("CoIA Axis 2 (%.1f%%)", cia_inertia[2])
    ) +
    guides(
      # Colour legend: shown when standalone, hidden in composite (moved to direction plot)
      color = if (show_color_legend) {
        guide_legend(order = 1,
                     override.aes = list(shape = 21, size = 2),
                     title = group_label)
      } else {
        "none"
      },
      fill  = "none",
      # Endpoint-shape legend is always shown
      shape = guide_legend(order = 2,
                           override.aes = list(fill  = "grey60",
                                               color = "grey30",
                                               size  = c(1, 1.3)),
                           title = "Endpoint"),
      size  = "none"
    ) +
    theme(
      panel.grid.major.x   = element_blank(),
      panel.grid.major.y   = element_blank(),
      panel.grid.minor     = element_blank(),
      axis.text = element_text(size = 5),
      legend.position      = c(0.01, 0.99),
      legend.justification = c(0, 1),
      legend.box           = "vertical",
      legend.box.just      = "left",
      legend.background    = element_rect(fill  = alpha("white", 0.75),
                                          color = "grey80", linewidth = 0.25),
      legend.key.size      = unit(0.3, "cm"),
      legend.text          = element_text(size = 6.5),
      legend.title         = element_text(size = 7),
      legend.margin        = margin(2, 4, 2, 4)
    )
  
  tag   <- if (!is.null(fname_tag)) fname_tag else tolower(str_replace_all(group_label, " ", "_"))
  fname <- sprintf("analysis/output/figures/L1_CoIA_Biplot_%s.png", tag)
  cat(sprintf("Figure built: %s\n", basename(fname)))
  p
}

# Standalone version (keeps colour legend)
make_coia_biplot("layer",        "Layer",        layer_pal,    LAYER_ORDER,    "layer")
make_coia_biplot("raw_material", "Raw Material", rawmat_pal,   NULL,           "raw_material")
make_coia_biplot("core_type",    "Core Type",    coretype_pal, CORETYPE_ORDER, "core_type")

# Composite version (colour legend hidden, moved to the right of the direction plot)
p_coia_layer    <- make_coia_biplot("layer",        "Layer",        layer_pal,
                                    LAYER_ORDER,    "layer",        show_color_legend = TRUE)
p_coia_rawmat   <- make_coia_biplot("raw_material", "Raw Material", rawmat_pal,
                                    NULL,           "raw_material", show_color_legend = TRUE)
p_coia_coretype <- make_coia_biplot("core_type",    "Core Type",    coretype_pal,
                                    CORETYPE_ORDER, "core_type",    show_color_legend = TRUE)

l1_results <- tibble(
  method  = c("Mantel (ILR, Euclidean, Spearman)", "RV (ILR, Euclidean)"),
  stat    = c(mantel_global$statistic, coin_arch$RV),
  p_value = c(mantel_global$signif,    rv_test$pvalue),
  n       = length(arch_ids)
)
write_csv(l1_results, here("analysis/data/derived_data/L1_results.csv"))
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

a_morph <- as.matrix(coin_arch$aX)
a_scar  <- as.matrix(coin_arch$aY)

w2_mpc_cia <- a_morph^2
w2_dpc_cia <- a_scar^2

morph_var_pct_full <- round(dudi_morph$eig / sum(dudi_morph$eig) * 100, 1)
scar_var_pct_full  <- round(dudi_scar$eig  / sum(dudi_scar$eig)  * 100, 1)

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
push('<title>ILR to PCA to CoIA contribution flow (SDG)</title>')
push('<desc>Sankey-style flow diagram showing ILR spectral variable contributions through PCA axes to CoIA axes.</desc>')
push(sprintf('<rect width="%d" height="%d" fill="white"/>', SVG_W, SVG_H))

push(svg_text(x_col1 + NODE_W / 2, 20, "ILR variables",    12, 500, fill = "#2C2C2A"))
push(svg_text(x_col1 + NODE_W / 2, 36, "spectral log-contrasts", 10, 400, fill = "#5F5E5A"))
push(svg_text(x_col2 + NODE_W / 2, 20, "PCA axes",          12, 500, fill = "#2C2C2A"))
push(svg_text(x_col2 + NODE_W / 2, 36, "per endpoint",      10, 400, fill = "#5F5E5A"))
push(svg_text(x_col3 + NODE_W / 2, 20, "CoIA axes",         12, 500, fill = "#2C2C2A"))
push(svg_text(x_col3 + NODE_W / 2, 36, "shared structure",  10, 400, fill = "#5F5E5A"))

outlet_ilr_m <- rep(0, n_ilr_m); inlet_mpc  <- rep(0, n_mpc)
outlet_ilr_d <- rep(0, n_ilr_d); inlet_dpc  <- rep(0, n_dpc)
inlet_cia    <- rep(0, n_cia_ax)

for (i in seq_len(n_ilr_m)) for (j in seq_len(n_mpc)) {
  ww <- w2_ilr_mpc[i, j]; if (ww < 0.004) next
  bw_out <- ww * h_ilr_m[i]; bw_in <- ww * h_mpc[j]
  y1t <- y_ilr_m[i] + outlet_ilr_m[i]; y1b <- y1t + bw_out
  y2t <- y_mpc[j]   + inlet_mpc[j];    y2b <- y2t + bw_in
  push(svg_band(x_col1 + NODE_W, y1t, y1b, x_col2, y2t, y2b,
                col_morph_band, 0.28 + 0.32 * ww))
  outlet_ilr_m[i] <- outlet_ilr_m[i] + bw_out
  inlet_mpc[j]    <- inlet_mpc[j]    + bw_in
}

for (i in seq_len(n_ilr_d)) for (j in seq_len(n_dpc)) {
  ww <- w2_ilr_dpc[i, j]; if (ww < 0.004) next
  bw_out <- ww * h_ilr_d[i]; bw_in <- ww * h_dpc[j]
  y1t <- y_ilr_d[i] + outlet_ilr_d[i]; y1b <- y1t + bw_out
  y2t <- y_dpc[j]   + inlet_dpc[j];    y2b <- y2t + bw_in
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
png_path <- here("analysis/output/figures/L1_CoIA_Sankey.png")
writeLines(lines, svg_path, useBytes = FALSE)

if (requireNamespace("rsvg", quietly = TRUE)) {
  rsvg::rsvg_png(svg_path, png_path, width = SVG_W * 2)
  cat("PNG saved (via rsvg): L1_CoIA_Sankey.png\n")
} else if (nzchar(Sys.which("rsvg-convert"))) {
  system2("rsvg-convert",
          args = c("-d", "300", "-p", "300", "-o", png_path, svg_path))
  cat("PNG saved (via rsvg-convert): L1_CoIA_Sankey.png\n")
} else if (nzchar(Sys.which("inkscape"))) {
  system2("inkscape",
          args = c("--export-filename", png_path, "--export-dpi", "300", svg_path))
  cat("PNG saved (via Inkscape): L1_CoIA_Sankey.png\n")
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

cat("\n\n")
cat("##  Level 2: joint evidence                                   ##\n")

# ------------------------------------------------------------------------------
# L2-A: grouped Mantel
# ------------------------------------------------------------------------------

cat("\n---------- L2-A: grouped Mantel ----------\n")

mantel_within_group <- function(group_val, group_col, D_morph_full,
                                D_scar_full, meta_df, n_perm = 9999) {
  ids <- meta_df %>%
    filter(.data[[group_col]] == group_val) %>%
    pull(ID)
  cat(sprintf("  -> %s (n = %d) ...", group_val, length(ids)))
  if (length(ids) < 5) {
    cat(" skipped (n < 5)\n"); return(NULL)
  }
  res <- mantel(extract_subdist(D_morph_full, ids),
                extract_subdist(D_scar_full,  ids),
                method = "spearman", permutations = n_perm)
  cat(sprintf(" r = %.4f, p = %.4f\n", res$statistic, res$signif))
  tibble(group_var = group_col, group = group_val,
         n = length(ids), mantel_r = res$statistic, p_raw = res$signif)
}

run_grouped_mantel <- function(meta_obj, group_col, label) {
  if (is.null(meta_obj)) {
    cat(sprintf("  [skip] %s grouped Mantel: too few groups\n", label))
    return(NULL)
  }
  map_dfr(unique(meta_obj[[group_col]]),
          ~ mantel_within_group(.x, group_col,
                                D_morph_arch, D_scar_arch, meta_arch))
}

mantel_by_layer    <- run_grouped_mantel(meta_layer,    "layer",        "layer")
mantel_by_rawmat   <- run_grouped_mantel(meta_rawmat,   "raw_material", "raw material")
mantel_by_coretype <- run_grouped_mantel(meta_coretype, "core_type",    "core type")

cat("\nMantel by layer:\n")
if (!is.null(mantel_by_layer))
  print(mantel_by_layer %>% mutate(across(c(mantel_r, p_raw), ~ round(.x, 4))))
cat("\nMantel by raw material:\n")
if (!is.null(mantel_by_rawmat))
  print(mantel_by_rawmat %>% mutate(across(c(mantel_r, p_raw), ~ round(.x, 4))))
cat("\nMantel by core type:\n")
if (!is.null(mantel_by_coretype))
  print(mantel_by_coretype %>% mutate(across(c(mantel_r, p_raw), ~ round(.x, 4))))

l2_mantel <- bind_rows(
  mantel_by_layer, mantel_by_rawmat, mantel_by_coretype
) %>%
  mutate(
    p_holm = p.adjust(p_raw, method = "holm"),
    significance = case_when(
      p_holm < 0.001 ~ "***", p_holm < 0.01  ~ "**",
      p_holm < 0.05  ~ "*",   p_holm < 0.10  ~ ".",
      TRUE           ~ "ns"
    ),
    group_var_label = recode(group_var,
                             "layer"        = "Layer",
                             "raw_material" = "Raw Material",
                             "core_type"    = "Core Type"),
    group = as.factor(group)
  ) %>%
  group_by(group_var_label) %>%
  mutate(group = fct_reorder(group, mantel_r)) %>%
  ungroup()

if (nrow(l2_mantel) > 0) {
  write_csv(l2_mantel, here("analysis/data/derived_data/L2_grouped_mantel.csv"))
  cat("Saved: L2_grouped_mantel.csv\n")
}

# ------------------------------------------------------------------------------
# L2-B: CoIA arrow length
# ------------------------------------------------------------------------------

cat("\n---------- L2-B: arrow-length group differences ----------\n")

run_arrow_length_analysis <- function(group_col, group_label, palette,
                                      group_order = NULL) {
  valid_groups <- scores_combined %>%
    filter(!is.na(.data[[group_col]]), .data[[group_col]] != "Other") %>%
    group_by(.data[[group_col]]) %>% filter(n() >= 3) %>%
    pull(.data[[group_col]]) %>% unique()
  if (length(valid_groups) < 2) {
    cat(sprintf("  [skip] %s arrow-length test: too few valid groups\n", group_label))
    return(invisible(NULL))
  }
  
  if (!is.null(group_order)) {
    valid_groups_plot <- intersect(group_order, valid_groups)
  } else {
    valid_groups_plot <- valid_groups
  }
  
  sub_df_stat <- scores_combined %>% filter(.data[[group_col]] %in% valid_groups)
  cat(sprintf("\n----- %s x arrow length -----\n", group_label))
  kw <- kruskal.test(reformulate(group_col, "arrow_length"), data = sub_df_stat)
  print(kw)
  pw <- pairwise.wilcox.test(sub_df_stat$arrow_length, sub_df_stat[[group_col]],
                             p.adjust.method = "holm", exact = FALSE)
  print(pw)
  
  sub_df <- scores_combined %>%
    filter(.data[[group_col]] %in% valid_groups_plot) %>%
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
             label = sprintf("Kruskal-Wallis\nchi\u00b2 = %.2f, p = %.3f",
                             kw$statistic, kw$p.value),
             hjust = 1.05, vjust = 1.2, size = 2.6, color = "grey40") +
    scale_fill_manual(values  = palette) +
    scale_color_manual(values = palette) +
    theme_bw(base_size = 8) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor   = element_blank(),
      axis.text.x        = element_blank(),
      axis.ticks.x       = element_blank(),
      axis.text.y        = element_text(size = 5),
      legend.position    = "none"
    ) +
    labs(x = NULL, y = "CoIA line length")
  
  fname <- sprintf("analysis/output/figures/L2_Arrow_Length_%s.png",
                   tolower(str_replace_all(group_label, " ", "_")))
  cat(sprintf("Figure built: %s\n", basename(fname)))
  list(sub_df = sub_df, kw = kw, pw = pw, p = p)
}

res_len_layer    <- run_arrow_length_analysis("layer",        "Layer",
                                              layer_pal,    LAYER_ORDER)
res_len_rawmat   <- run_arrow_length_analysis("raw_material", "Raw_Material",
                                              rawmat_pal)
res_len_coretype <- run_arrow_length_analysis("core_type",    "Core_Type",
                                              coretype_pal, CORETYPE_ORDER)

cat("\n==== Arrow-length summary ====\n")
for (gc in c("layer", "raw_material", "core_type")) {
  cat(sprintf("\n--- %s ---\n", gc))
  scores_combined %>%
    filter(!is.na(.data[[gc]]), .data[[gc]] != "Other") %>%
    group_by(.data[[gc]]) %>%
    summarise(n      = n(),
              mean   = round(mean(arrow_length),   4),
              median = round(median(arrow_length), 4),
              sd     = round(sd(arrow_length),     4),
              min    = round(min(arrow_length),    4),
              max    = round(max(arrow_length),    4),
              .groups = "drop") %>%
    print()
}

# ------------------------------------------------------------------------------
# L2-C: CoIA arrow-direction circular statistics (linear display)
# plot_rose() accepts show_color_legend:
#   TRUE  = standalone, colour legend on the right (default)
#   FALSE = composite internal call, no colour legend
#   "right_panel" = rightmost composite column, colour legend on the right
# ------------------------------------------------------------------------------

cat("\n---------- L2-C: arrow-direction circular statistics ----------\n")

plot_rose <- function(res, palette, show_color_legend = TRUE) {
  group_col <- res$group_col
  kde_df    <- res$kde_df
  mean_dirs <- res$mean_dirs
  
  kde_linear <- kde_df %>%
    mutate(angle_centered = if_else(angle_deg > 180, angle_deg - 360, angle_deg),
           !!group_col := factor(.data[[group_col]],
                                 levels = levels(kde_df[[group_col]])))
  
  mean_linear <- mean_dirs %>%
    mutate(angle_centered = if_else(mean_deg > 180, mean_deg - 360, mean_deg),
           !!group_col := factor(.data[[group_col]],
                                 levels = levels(kde_df[[group_col]])))
  
  rayleigh_labels <- res$rayleigh %>%
    mutate(
      label = case_when(
        rayleigh_p < 0.001 ~ "Rayleigh p < 0.001",
        rayleigh_p < 0.05  ~ sprintf("Rayleigh p = %.3f", rayleigh_p),
        TRUE               ~ sprintf("Rayleigh p = %.3f", rayleigh_p)
      ),
      !!group_col := factor(group, levels = levels(kde_df[[group_col]]))
    )
  
  # Legend placement strategy
  if (isTRUE(show_color_legend)) {
    # Standalone: colour legend on the right
    legend_pos    <- "right"
    color_guide   <- guide_legend(title = group_col)
  } else if (identical(show_color_legend, "right_panel")) {
    # Rightmost composite column: colour legend on the right
    legend_pos    <- "right"
    color_guide   <- guide_legend(title = group_col)
  } else {
    # Non-final composite column: hide colour legend
    legend_pos    <- "none"
    color_guide   <- "none"
  }
  
  ggplot(kde_linear,
         aes(x     = angle_centered,
             y     = density,
             fill  = .data[[group_col]],
             color = .data[[group_col]])) +
    geom_area(alpha = 0.4, color = NULL) +
    geom_vline(
      data      = mean_linear,
      aes(xintercept = angle_centered, color = .data[[group_col]]),
      linewidth = 0.4, linetype = "dashed", alpha = 0.75
    ) +
    geom_text(
      data = rayleigh_labels,
      aes(label = label),
      x = 160, y = Inf,
      hjust = 1, vjust = 1.4,
      size = 2.4, color = "grey35",
      inherit.aes = FALSE
    ) +
    scale_x_continuous(
      limits = c(-180, 180),
      breaks = seq(-180, 180, by = 45),
      labels = seq(-180, 180, by = 45)
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
    scale_fill_manual(values  = palette, guide = color_guide) +
    scale_color_manual(values = palette, guide = "none") +
    # In composite: ncol=1 vertical facets, facet labels hidden
    facet_wrap(reformulate(group_col), ncol = 1, scales = "free_y") +
    labs(x    = "CoIA line direction (\u00b0)",
         y    = "von Mises KDE") +
    theme_bw(base_size = 8) +
    theme(
      panel.grid.minor    = element_blank(),
      panel.grid.major.x  = element_blank(),
      panel.grid.major.y  = element_blank(),
      # Facet labels: shown standalone, hidden in composite
      strip.text          = if (isTRUE(show_color_legend)) {
        element_text(face = "bold", size = 7)
      } else {
        element_blank()
      },
      strip.background    = if (isTRUE(show_color_legend)) {
        element_rect(fill = "#EBEBEB", color = "#EBEBEB")
      } else {
        element_blank()
      },
      axis.text.x         = element_text(size = 5),
      axis.text.y         = element_blank(),
      axis.ticks.y        = element_blank(),
      legend.position     = legend_pos,
      legend.key.size     = unit(0.32, "cm"),
      legend.text         = element_text(size = 6.5),
      legend.title        = element_text(size = 7, face = "bold")
    )
}

# ---- analysis helper ----
run_circular_analysis <- function(group_col, group_label, palette,
                                  group_order = NULL) {
  valid_groups <- scores_combined %>%
    filter(!is.na(.data[[group_col]]), .data[[group_col]] != "Other") %>%
    group_by(.data[[group_col]]) %>% filter(n() >= 5) %>%
    pull(.data[[group_col]]) %>% unique()
  if (length(valid_groups) < 2) {
    cat(sprintf("  [skip] %s circular statistics: too few valid groups\n", group_label))
    return(invisible(NULL))
  }
  
  if (!is.null(group_order)) {
    valid_groups_ordered <- intersect(group_order, valid_groups)
  } else {
    valid_groups_ordered <- sort(valid_groups)
  }
  
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
    # axial: test for a preferred axis by doubling the undirected line angles
    circ_obj <- circular(2 * angles, type = "angles", units = "radians", modulo = "2pi")
    rt       <- rayleigh.test(circ_obj)
    cat(sprintf("  %s: U = %.4f, p = %.4f -> %s\n",
                g, rt$statistic, rt$p.value,
                ifelse(rt$p.value < 0.05, "concentrated", "dispersed")))
    tibble(group_var = group_col, group = g,
           rayleigh_U = round(rt$statistic, 4),
           rayleigh_p = round(rt$p.value,   4),
           conclusion = ifelse(rt$p.value < 0.05, "concentrated", "uniform"))
  })
  
  mean_dirs <- map_dfr(valid_groups_ordered, function(g) {
    angles <- sub_df %>% filter(.data[[group_col]] == g) %>% pull(arrow_angle)
    cs     <- circ_stats_one(angles)
    # axial mean is an axis: draw it at both ends (mean and mean+180)
    tibble(!!group_col := g,
           mean_deg = c(cs$mean_deg, (cs$mean_deg + 180) %% 360))
  }) %>%
    mutate(!!group_col := factor(.data[[group_col]], levels = valid_groups_ordered))
  
  compute_circular_kde <- function(angles_deg, bw = 25, n = 360) {
    # axial: reflect each undirected line by 180 deg so the KDE is symmetric
    # about the axis (theta equivalent to theta+180)
    angles_deg <- c(angles_deg, (angles_deg + 180) %% 360)
    circ <- circular(angles_deg * pi / 180,
                     type = "angles", units = "radians", modulo = "2pi")
    dens <- density(circ, bw = bw, n = n)
    angle_deg <- as.numeric(dens$x) * 180 / pi %% 360
    density   <- as.numeric(dens$y)
    tibble(angle_deg = c(angle_deg, angle_deg[1]),
           density   = c(density,   density[1]))
  }
  
  kde_df <- sub_df %>%
    mutate(
      angle_deg = arrow_angle * 180 / pi,
      angle_deg = ifelse(angle_deg < 0, angle_deg + 360, angle_deg),
      !!group_col := factor(.data[[group_col]], levels = valid_groups_ordered)
    ) %>%
    group_by(.data[[group_col]]) %>%
    group_modify(~ compute_circular_kde(.x$angle_deg, bw = 25)) %>%
    ungroup() %>%
    mutate(!!group_col := factor(.data[[group_col]], levels = valid_groups_ordered))
  
  # Standalone: default show_color_legend = TRUE (horizontal facets, labelled)
  p_rose_standalone <- plot_rose(
    list(group_col = group_col, kde_df = kde_df,
         mean_dirs = mean_dirs, rayleigh = rayleigh_res),
    palette,
    show_color_legend = TRUE
  )
  
  n_g   <- length(valid_groups_ordered)
  fname <- sprintf("analysis/output/figures/L2_Arrow_Direction_rose_%s.png",
                   tolower(str_replace_all(group_label, " ", "_")))
  cat(sprintf("Figure built: %s\n", basename(fname)))
  
  # Composite use: pass raw data and config for on-demand rebuild
  list(desc      = circ_desc,
       rayleigh  = rayleigh_res,
       kde_df    = kde_df,
       mean_dirs = mean_dirs,
       group_col = group_col,
       palette   = palette,
       p_rose    = p_rose_standalone)
}

res_circ_layer    <- run_circular_analysis("layer",        "Layer",
                                           layer_pal,    LAYER_ORDER)
res_circ_rawmat   <- run_circular_analysis("raw_material", "Raw_Material",
                                           rawmat_pal)
res_circ_coretype <- run_circular_analysis("core_type",    "Core_Type",
                                           coretype_pal, CORETYPE_ORDER)

# ==============================================================================
# ---- Save L2 derived data ----
# ==============================================================================

circ_desc_all <- bind_rows(
  if (!is.null(res_circ_layer))    res_circ_layer$desc,
  if (!is.null(res_circ_rawmat))   res_circ_rawmat$desc,
  if (!is.null(res_circ_coretype)) res_circ_coretype$desc
)
rayleigh_all <- bind_rows(
  if (!is.null(res_circ_layer))    res_circ_layer$rayleigh,
  if (!is.null(res_circ_rawmat))   res_circ_rawmat$rayleigh,
  if (!is.null(res_circ_coretype)) res_circ_coretype$rayleigh
)
if (nrow(circ_desc_all) > 0 && nrow(rayleigh_all) > 0) {
  left_join(circ_desc_all, rayleigh_all, by = c("group_var", "group")) %>%
    write_csv(here("analysis/data/derived_data/L2_circular_stats.csv"))
  cat("Saved: L2_circular_stats.csv\n")
}

scores_combined %>%
  select(ID, layer, raw_material, core_type,
         arrow_length, arrow_angle,
         Axis1_M, Axis2_M, Axis1_S, Axis2_S) %>%
  write_csv(here("analysis/data/derived_data/L2_arrow_stats.csv"))
cat("Saved: L2_arrow_stats.csv\n")

scores_combined %>%
  write_csv(here("analysis/data/derived_data/CoIA_scores_full.csv"))
cat("Saved: CoIA_scores_full.csv\n")

morph_ilr_arch %>%
  rownames_to_column("ID") %>%
  write_csv(here("analysis/data/derived_data/SDG_morph_ILR_scores.csv"))
cat("Saved: SDG_morph_ILR_scores.csv\n")

scar_ilr_arch %>%
  rownames_to_column("ID") %>%
  write_csv(here("analysis/data/derived_data/SDG_scar_ILR_scores.csv"))
cat("Saved: SDG_scar_ILR_scores.csv\n")

# ==============================================================================
# ---- Composite: CoIA x arrow length x direction (3 rows x 3 cols) + external image row ----
# ==============================================================================
cat("\n==== Composite: CoIA / arrow length / direction ====\n")

make_rose_for_composite <- function(res_circ) {
  if (is.null(res_circ)) return(NULL)
  plot_rose(
    list(group_col = res_circ$group_col,
         kde_df    = res_circ$kde_df,
         mean_dirs = res_circ$mean_dirs,
         rayleigh  = res_circ$rayleigh),
    res_circ$palette,
    show_color_legend = FALSE
  ) + theme(plot.margin = margin(4, 4, 4, 4))
}

strip_margin <- function(p) p + theme(plot.margin = margin(4, 4, 4, 4))

get_len_plot <- function(res_len) {
  if (is.null(res_len)) return(NULL)
  res_len$p + theme(plot.margin = margin(4, 4, 4, 4))
}

make_composite_row <- function(p_coia, p_len, res_circ, row_tag) {
  p_rose <- make_rose_for_composite(res_circ)
  if (is.null(p_coia) || is.null(p_len) || is.null(p_rose)) {
    cat(sprintf("  [skip] %s row: incomplete subplots\n", row_tag))
    return(NULL)
  }
  (strip_margin(p_coia) | get_len_plot(p_len) | p_rose) +
    plot_layout(widths = c(5, 2, 2))
}

row_layer    <- make_composite_row(p_coia_layer,    res_len_layer,    res_circ_layer,    "Layer")
row_rawmat   <- make_composite_row(p_coia_rawmat,   res_len_rawmat,   res_circ_rawmat,   "Raw Material")
row_coretype <- make_composite_row(p_coia_coretype, res_len_coretype, res_circ_coretype, "Core Type")

rows_valid <- Filter(Negate(is.null), list(row_layer, row_rawmat, row_coretype))

if (length(rows_valid) > 0) {
  n_rows     <- length(rows_valid)
  n_subplots <- n_rows * 3  # 3 columns per row
  
  # ---- Step 1: main composite (no plot_annotation; manual tags below) ----
  # Tag each subplot A-I manually
  all_tags <- letters[seq_len(n_subplots)]
  tag_idx  <- 1L
  
  rows_tagged <- lapply(rows_valid, function(row) {
    # a row is a patchwork object and cannot be indexed cell-by-cell,
    # so tagging during rebuild in make_composite_row is more reliable.
    # wrap_elements freezes the whole row; tags are applied at the subplot level (see below)
    row
  })
  
  # Safer: tag the three subplots before make_composite_row returns
  # Rebuild tagged rows
  make_tagged_row <- function(p_coia, p_len, res_circ, tags) {
    p_rose <- make_rose_for_composite(res_circ)
    if (is.null(p_coia) || is.null(p_len) || is.null(p_rose)) return(NULL)
    
    p1 <- strip_margin(p_coia)    + labs(tag = tags[1]) + theme(plot.tag = element_text(size = 9, face = "bold"))
    p2 <- get_len_plot(p_len)     + labs(tag = tags[2]) + theme(plot.tag = element_text(size = 9, face = "bold"))
    p3 <- p_rose                  + labs(tag = tags[3]) + theme(plot.tag = element_text(size = 9, face = "bold"))
    
    (p1 | p2 | p3) + plot_layout(widths = c(5, 2, 2))
  }
  
  # Rebuild each tagged row
  inputs <- list(
    list(p_coia_layer,    res_len_layer,    res_circ_layer),
    list(p_coia_rawmat,   res_len_rawmat,   res_circ_rawmat),
    list(p_coia_coretype, res_len_coretype, res_circ_coretype)
  )
  valid_mask <- !sapply(rows_valid, is.null)  # rows_valid already filtered; length == n_rows
  
  tagged_rows <- vector("list", n_rows)
  for (i in seq_len(n_rows)) {
    tags_i <- letters[((i - 1) * 3 + 1):(i * 3)]
    inp    <- inputs[[i]]
    tagged_rows[[i]] <- make_tagged_row(inp[[1]], inp[[2]], inp[[3]], tags_i)
  }
  
  # ---- Step 2: assemble body (frozen so the external image does not break layout) ----
  p_main <- Reduce(`/`, tagged_rows) +
    plot_layout(heights = rep(1, n_rows))
  
  # ---- Step 3: external image, manual final tag ----
  next_tag    <- letters[n_subplots + 1]
  external_img <- png::readPNG(here("analysis/figures/source_panels/Axis_trajectory_SDG.png"))
  grob_img     <- grid::rasterGrob(external_img, interpolate = TRUE,
                                   x     = grid::unit(0.5055, "npc"),
                                   width = grid::unit(0.949,  "npc"))
  p_external   <- wrap_elements(full = grob_img) +
    labs(tag = next_tag) +
    theme(
      plot.tag    = element_text(size = 9, face = "bold"),
      plot.margin = margin(0, 0, 0, 0)
    )
  
  # ---- Step 4: final assembly ----
  p_final <- wrap_elements(full = p_main) / p_external +
    plot_layout(heights = c(n_rows, 1))

  cat(sprintf("Figure built: L_CoIA_composite.png (%d rows x 3 cols + external-image row)\n", n_rows))
  
} else {
  cat("  [skip] no valid subplots in any row; composite not built\n")
}
# ==============================================================================
# ========== Level 3: PERMANOVA ==========
# ==============================================================================

cat("\n\n")
cat("##  Level 3: PERMANOVA — group-structure analysis             ##\n")

run_permanova <- function(dist_mat, group_vec, group_name, domain_name,
                          n_perm = 9999) {
  group_vec <- as.character(group_vec)
  if (length(unique(group_vec)) < 2) {
    cat(sprintf("  [skip] %s ~ %s: too few group levels\n", domain_name, group_name))
    return(NULL)
  }
  df_tmp <- data.frame(group = factor(group_vec))
  res    <- adonis2(dist_mat ~ group, data = df_tmp,
                    permutations = n_perm, add = "lingoes")
  cat(sprintf("\n----- %s ~ %s -----\n", domain_name, group_name))
  print(res)
  tibble(domain   = domain_name,
         grouping = group_name,
         R2       = res$R2[1],
         F_value  = res$F[1],
         p_value  = res$`Pr(>F)`[1],
         df_group = res$Df[1],
         df_resid = res$Df[2])
}

run_permdisp <- function(dist_mat, group_vec, group_name, domain_name) {
  group_vec <- as.character(group_vec)
  if (length(unique(group_vec)) < 2) return(NULL)
  eig_vals    <- cmdscale(dist_mat, eig = TRUE)$eig
  correction  <- abs(min(c(eig_vals[eig_vals < 0], 0)))
  D_corrected <- as.dist(as.matrix(dist_mat) + correction)
  bd  <- betadisper(D_corrected, factor(group_vec))
  res <- permutest(bd, permutations = 9999)
  cat(sprintf("  PERMDISP — %s ~ %s: F = %.3f, p = %.3f\n",
              domain_name, group_name,
              res$tab$F[1], res$tab$`Pr(>F)`[1]))
  tibble(domain = domain_name, grouping = group_name,
         F_disp = res$tab$F[1], p_disp = res$tab$`Pr(>F)`[1])
}

pairwise_permanova <- function(dist_mat, group_vec, group_name,
                               domain_name, n_perm = 9999) {
  group_vec <- as.character(group_vec)
  groups    <- unique(group_vec)
  if (length(groups) < 2) return(NULL)
  map_dfr(combn(groups, 2, simplify = FALSE), function(pair) {
    idx    <- group_vec %in% pair
    d_sub  <- as.dist(as.matrix(dist_mat)[idx, idx])
    df_tmp <- data.frame(group = factor(group_vec[idx]))
    res    <- adonis2(d_sub ~ group, data = df_tmp,
                      permutations = n_perm, add = "lingoes")
    tibble(domain = domain_name, grouping = group_name,
           group1 = pair[1], group2 = pair[2],
           R2 = res$R2[1], F_value = res$F[1],
           p_raw = res$`Pr(>F)`[1])
  })
}

if (!is.null(meta_layer)) {
  D_morph_layer <- extract_subdist(D_morph_arch, meta_layer$ID)
  D_scar_layer  <- extract_subdist(D_scar_arch,  meta_layer$ID)
}
if (!is.null(meta_rawmat)) {
  D_morph_rawmat <- extract_subdist(D_morph_arch, meta_rawmat$ID)
  D_scar_rawmat  <- extract_subdist(D_scar_arch,  meta_rawmat$ID)
}
if (!is.null(meta_coretype)) {
  D_morph_coretype <- extract_subdist(D_morph_arch, meta_coretype$ID)
  D_scar_coretype  <- extract_subdist(D_scar_arch,  meta_coretype$ID)
}

cat("\n==== L3-1: global PERMANOVA ====\n")

permanova_results <- bind_rows(
  if (!is.null(meta_layer)) bind_rows(
    run_permanova(D_morph_layer,    meta_layer$layer,         "Layer",        "Morphology"),
    run_permanova(D_scar_layer,     meta_layer$layer,         "Layer",        "Scar Direction")
  ),
  if (!is.null(meta_rawmat)) bind_rows(
    run_permanova(D_morph_rawmat,   meta_rawmat$raw_material, "Raw Material", "Morphology"),
    run_permanova(D_scar_rawmat,    meta_rawmat$raw_material, "Raw Material", "Scar Direction")
  ),
  if (!is.null(meta_coretype)) bind_rows(
    run_permanova(D_morph_coretype, meta_coretype$core_type,  "Core Type",    "Morphology"),
    run_permanova(D_scar_coretype,  meta_coretype$core_type,  "Core Type",    "Scar Direction")
  )
)

cat("\n==== L3-2：PERMDISP ====\n")

permdisp_results <- bind_rows(
  if (!is.null(meta_layer)) bind_rows(
    run_permdisp(D_morph_layer,    meta_layer$layer,         "Layer",        "Morphology"),
    run_permdisp(D_scar_layer,     meta_layer$layer,         "Layer",        "Scar Direction")
  ),
  if (!is.null(meta_rawmat)) bind_rows(
    run_permdisp(D_morph_rawmat,   meta_rawmat$raw_material, "Raw Material", "Morphology"),
    run_permdisp(D_scar_rawmat,    meta_rawmat$raw_material, "Raw Material", "Scar Direction")
  ),
  if (!is.null(meta_coretype)) bind_rows(
    run_permdisp(D_morph_coretype, meta_coretype$core_type,  "Core Type",    "Morphology"),
    run_permdisp(D_scar_coretype,  meta_coretype$core_type,  "Core Type",    "Scar Direction")
  )
)

cat("\n==== L3-3: pairwise PERMANOVA ====\n")

pairwise_results <- bind_rows(
  if (!is.null(meta_layer)) bind_rows(
    pairwise_permanova(D_morph_layer,    meta_layer$layer,         "Layer",        "Morphology"),
    pairwise_permanova(D_scar_layer,     meta_layer$layer,         "Layer",        "Scar Direction")
  ),
  if (!is.null(meta_rawmat)) bind_rows(
    pairwise_permanova(D_morph_rawmat,   meta_rawmat$raw_material, "Raw Material", "Morphology"),
    pairwise_permanova(D_scar_rawmat,    meta_rawmat$raw_material, "Raw Material", "Scar Direction")
  ),
  if (!is.null(meta_coretype)) bind_rows(
    pairwise_permanova(D_morph_coretype, meta_coretype$core_type,  "Core Type",    "Morphology"),
    pairwise_permanova(D_scar_coretype,  meta_coretype$core_type,  "Core Type",    "Scar Direction")
  )
) %>%
  mutate(
    p_holm = p.adjust(p_raw, method = "holm"),
    significance = case_when(
      p_holm < 0.001 ~ "***", p_holm < 0.01 ~ "**",
      p_holm < 0.05  ~ "*",   TRUE           ~ "ns"
    )
  )

cat("\nPairwise PERMANOVA (Holm):\n")
pairwise_results %>%
  mutate(across(c(R2, F_value, p_raw, p_holm), ~ round(.x, 4))) %>%
  print(n = Inf)

# ==============================================================================
# ---- Save numeric results ----
# ==============================================================================

bind_rows(
  permanova_results %>%
    mutate(test = "PERMANOVA", statistic = F_value) %>%
    select(test, domain, grouping, R2, statistic, p_value),
  permdisp_results %>%
    mutate(test = "PERMDISP", R2 = NA,
           statistic = F_disp, p_value = p_disp) %>%
    select(test, domain, grouping, R2, statistic, p_value),
  pairwise_results %>%
    mutate(test = "Pairwise_PERMANOVA",
           statistic = F_value, p_value = p_holm) %>%
    select(test, domain, grouping, R2, statistic, p_value)
) %>%
  write_csv(here("analysis/data/derived_data/L3_permanova.csv"))
cat("Saved: L3_permanova.csv\n")

# ==============================================================================
# ---- Summary print ----
# ==============================================================================

cat("\n\n")
cat("##  Three-level analysis summary (SDG)                         ##\n")

cat("\n[Level 1: baseline]\n")
cat(sprintf("  Mantel r = %.4f, p = %.3f  ->  %s\n",
            mantel_global$statistic, mantel_global$signif,
            ifelse(mantel_global$signif < 0.05, "correlated", "independent (n.s.)")))
cat(sprintf("  RV       = %.4f, p = %.3f  ->  %s\n",
            coin_arch$RV, rv_test$pvalue,
            ifelse(rv_test$pvalue < 0.05, "covarying", "independent (n.s.)")))

cat("\n[Level 2A: grouped Mantel (Holm)]\n")
if (nrow(l2_mantel) > 0) {
  n_sig <- sum(l2_mantel$p_holm < 0.05, na.rm = TRUE)
  cat(sprintf("  %d subgroups tested, significant after Holm: %d\n",
              nrow(l2_mantel), n_sig))
  if (n_sig == 0) {
    cat("  -> all subgroups non-significant; morphology-technology independence is robust\n")
  } else {
    cat("  -> the following subgroups are significant (potential findings):\n")
    l2_mantel %>% filter(p_holm < 0.05) %>%
      select(group_var, group, n, mantel_r, p_holm, significance) %>%
      print()
  }
}

cat("\n[Level 2B] see L2_Arrow_Length_*.png and descriptive stats\n")
cat("[Level 2C] see L2_Arrow_Direction_rose_*.png and L2_circular_stats.csv\n")

cat("\n[Level 3: PERMANOVA]\n")
permanova_results %>%
  mutate(
    sig = case_when(
      p_value < 0.001 ~ "***", p_value < 0.01 ~ "**",
      p_value < 0.05  ~ "*",   TRUE            ~ "ns"
    )
  ) %>%
  select(domain, grouping, R2, F_value, p_value, sig) %>%
  mutate(across(c(R2, F_value, p_value), ~ round(.x, 4))) %>%
  arrange(domain, grouping) %>%
  print()

cat("\n[CoIA biplots] L1_CoIA_Biplot_layer.png / raw_material.png / core_type.png\n")
cat("[Composite] L_CoIA_composite.png (CoIA | arrow length | direction, 3 rows x 3 cols)\n")
cat("[Sankey] L1_CoIA_Sankey.png\n")

cat("\n\n========== SDG analysis complete ==========\n")
cat("Main output files:\n")
cat("  L1_results.csv\n")
cat("  L1_CoIA_Sankey.png\n")
cat("  L_CoIA_composite.png\n")
cat("  L2_grouped_mantel.csv\n")
cat("  L2_arrow_stats.csv\n")
cat("  L2_circular_stats.csv\n")
cat("  L3_permanova.csv\n")
cat("  CoIA_scores_full.csv\n")
cat("  CoIA_coords_full.csv\n")
cat("  PCA_CoIA_contribution.csv\n")
cat("  SDG_morph_ILR_scores.csv\n")
cat("  SDG_scar_ILR_scores.csv\n")
