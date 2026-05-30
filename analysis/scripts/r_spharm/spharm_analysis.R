# spharm_analysis.R
# SPHARM feature analysis: LDA visualisation + statistical tests.
# Sourced by the `spharm_analysis` target.
#
# Pipeline:
#   1. Load SPHARM power spectra (direction + morphology) and specimen metadata.
#   2. Filter power-spectrum features; split into EXP+IM / SDG+IM subsets.
#   3. EXP specimens: z-score standardisation + LDA visualisation.
#   4. Direction metrics (SPI, fabric E+I): KW + Dunn + PERMANOVA + PERMDISP.
#   5. Pairwise post-hoc summary bubble plot.
#   6. Save filtered data for downstream scripts.
#   7. Combined panel: three left plots + bubble plot on the right.
#
# Input:
#   - analysis/data/derived_data/SPHARM_direction.csv
#   - analysis/data/derived_data/SPHARM_morphology.csv
#   - analysis/data/derived_data/directions_aligned_svd.csv
#   - analysis/data/raw_data/SDG_core_metric.xlsx
#
# Returns (object, no figure written to disk): exp_method_compare_combined

library(here)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(readxl)
library(MASS)
library(vegan)
library(FSA)
library(RVAideMemoire)
library(ggtern)

conflicted::conflicts_prefer(ggplot2::aes)
conflicted::conflicts_prefer(ggplot2::theme_bw)
conflicted::conflicts_prefer(ggplot2::ggsave)
conflicted::conflicts_prefer(ggplot2::annotate)
conflicted::conflicts_prefer(dplyr::select)
conflicted::conflicts_prefer(dplyr::filter)
conflicted::conflicts_prefer(stats::sd)
set.seed(42)

# ==============================================================================
# Global parameters
# ==============================================================================

POWER_COLS_DIR   <- paste0("power_l", 1:6)
POWER_COLS_MORPH <- paste0("power_l", 1:8)
EXCLUDE_TYPES    <- c("Biface")
LEVALLOIS_MERGE  <- c("Levallois convergent", "Levallois laminar",
                      "Levallois preferential", "Levallois recurrent")

TYPOLOGY_COLORS <- c(
  "Levallois"      = "#4A6E8A",
  "Discoid"        = "#802520",
  "Unidirectional" = "#BA8530",
  "Multiplatform"  = "#8A7A68",
  "Bidirectional"  = "#788C4A"
)

TYPOLOGY_ORDER <- c(
  "Unidirectional", "Bidirectional", "Levallois", "Discoid", "Multiplatform"
)

# ==============================================================================
# 1. Load data
# ==============================================================================

SPHARM_direction  <- read_csv(here("analysis/data/derived_data/SPHARM_direction.csv"))
SPHARM_morphology <- read_csv(here("analysis/data/derived_data/SPHARM_morphology.csv"))

metric_data <- read_xlsx(here("analysis/data/raw_data/SDG_core_metric.xlsx")) %>%
  select(ID, Layer, Core_type_Li_merged, Raw_mat) %>%
  mutate(Layer = as.factor(Layer), Raw_mat = as.factor(Raw_mat))

SPHARM_morphology <- SPHARM_morphology %>%
  left_join(SPHARM_direction %>% select(ID, Typology), by = "ID")

# ==============================================================================
# 2. Filter features + split subsets
# ==============================================================================

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

# ==============================================================================
# 3. EXP specimens: z-score standardisation + LDA visualisation
# ==============================================================================

df_exp_dir   <- dir_splits$exp_im
df_exp_morph <- morph_splits$exp_im

scale_features <- function(df_target, cols) {
  ref_mat  <- df_target %>%
    filter(!str_starts(ID, "IM_")) %>%
    select(all_of(cols)) %>%
    as.matrix()
  col_mean <- colMeans(ref_mat)
  col_sd   <- apply(ref_mat, 2, sd)
  mat      <- df_target %>% select(all_of(cols)) %>% as.matrix()
  scale(mat, center = col_mean, scale = col_sd)
}

z_dir   <- scale_features(df_exp_dir,   POWER_COLS_DIR)
z_morph <- scale_features(df_exp_morph, POWER_COLS_MORPH)
colnames(z_dir)   <- paste0("dir_",   POWER_COLS_DIR)
colnames(z_morph) <- paste0("morph_", POWER_COLS_MORPH)

df_exp_only <- df_exp_dir %>%
  filter(!str_starts(ID, "IM_"), !Typology %in% EXCLUDE_TYPES) %>%
  mutate(
    Typology = case_when(
      Typology %in% LEVALLOIS_MERGE ~ "Levallois",
      TRUE ~ Typology
    ),
    Typology = droplevels(as.factor(Typology))
  )

non_im_idx <- !str_starts(df_exp_dir$ID, "IM_") &
  !df_exp_dir$Typology %in% EXCLUDE_TYPES

y_typology <- df_exp_only$Typology

cat(sprintf("EXP specimens retained: %d, classes: %d\n", nrow(df_exp_only), nlevels(y_typology)))
print(table(y_typology))

# --- LDA helper ---
run_lda_plot <- function(X, y, ids, filename) {
  fit      <- lda(X, grouping = y)
  prop_var <- fit$svd^2 / sum(fit$svd^2)
  
  scores <- predict(fit)$x %>%
    as.data.frame() %>%
    mutate(Typology = factor(y, levels = TYPOLOGY_ORDER), ID = ids)
  
  hull_df <- scores %>%
    group_by(Typology) %>%
    slice(chull(LD1, LD2)) %>%
    ungroup()
  
  p <- ggplot(scores, aes(x = LD1, y = LD2, color = Typology)) +
    geom_hline(yintercept = 0, color = "grey50", linewidth = 0.35, linetype = "dashed") +
    geom_vline(xintercept = 0, color = "grey50", linewidth = 0.35, linetype = "dashed") +
    geom_polygon(data = hull_df, aes(fill = Typology, group = Typology),
                 alpha = 0.25, color = NA) +
    geom_polygon(data = hull_df, aes(color = Typology, group = Typology),
                 fill = NA, linewidth = 0.01) +
    geom_point(size = 2.0, alpha = 0.88, stroke = 0.3, shape = 16) +
    scale_color_manual(values = TYPOLOGY_COLORS) +
    scale_fill_manual(values  = TYPOLOGY_COLORS) +
    scale_x_continuous(limits = c(-3, 4), expand = expansion(mult = 0.08),
                       breaks = seq(-3, 4, by = 1)) +
    scale_y_continuous(limits = c(-5.5, 3.5), expand = expansion(mult = 0.08),
                       breaks = seq(-5.5, 3.5, by = 1)) +
    labs(x = sprintf("LD1 (%.1f%%)", prop_var[1] * 100),
         y = sprintf("LD2 (%.1f%%)", prop_var[2] * 100)) +
    theme_bw() +
    theme(
      panel.grid        = element_blank(),
      legend.position   = c(0.9, 0.2),
      legend.justification  = c(0.5, 0.5),
      legend.background     = element_rect(fill = "transparent", colour = NA),
      legend.box.background = element_rect(fill = "transparent", colour = NA)
    )
  
  
  invisible(list(fit = fit, scores = scores, prop_var = prop_var))
}

res_morph <- run_lda_plot(z_morph[non_im_idx, ], y_typology, df_exp_only$ID,
                          "LDA_morph_by_Typology.png")
res_dir   <- run_lda_plot(z_dir[non_im_idx, ],   y_typology, df_exp_only$ID,
                          "LDA_dir_by_Typology.png")

# ==============================================================================
# 4. Direction metrics: SPI + fabric (Benn) + PERMANOVA + PERMDISP
# ==============================================================================

compute_SPI <- function(ux, uy, uz) {
  sqrt(sum(ux)^2 + sum(uy)^2 + sum(uz)^2) / length(ux)
}

compute_EI <- function(ux, uy, uz) {
  n      <- length(ux)
  U      <- cbind(ux, uy, uz)
  T_mat  <- (t(U) %*% U) / n
  eig    <- eigen(T_mat, symmetric = TRUE)
  lambda <- sort(eig$values, decreasing = TRUE)
  lambda <- pmax(lambda, 0)
  list(
    E       = ifelse(lambda[1] > 1e-10, 1 - lambda[2] / lambda[1], NA_real_),
    I       = ifelse(lambda[1] > 1e-10,     lambda[3] / lambda[1], NA_real_),
    lambda1 = lambda[1], lambda2 = lambda[2], lambda3 = lambda[3]
  )
}

raw_dirs <- read_csv(here("analysis/data/derived_data/directions_aligned_svd.csv")) %>%
  filter(str_starts(ID, "EXP"), !is.na(ux), !is.na(uy), !is.na(uz))

results <- raw_dirs %>%
  group_by(ID) %>%
  summarise(
    n_scars = n(),
    SPI     = compute_SPI(ux, uy, uz),
    EI      = list(compute_EI(ux, uy, uz)),
    .groups = "drop"
  ) %>%
  mutate(
    E       = map_dbl(EI, "E"),
    I       = map_dbl(EI, "I"),
    lambda1 = map_dbl(EI, "lambda1"),
    lambda2 = map_dbl(EI, "lambda2"),
    lambda3 = map_dbl(EI, "lambda3")
  ) %>%
  select(-EI) %>%
  arrange(ID)

cat("SPI/E/I computed, specimens:", nrow(results), "\n")
print(results %>%
        select(ID, n_scars, SPI, E, I) %>%
        mutate(across(c(SPI, E, I), \(x) round(x, 4))),
      n = Inf)

results_typed <- results %>%
  left_join(SPHARM_direction %>% select(ID, Typology), by = "ID") %>%
  filter(str_starts(ID, "EXP"), !Typology %in% EXCLUDE_TYPES) %>%
  mutate(
    Typology = case_when(
      Typology %in% LEVALLOIS_MERGE ~ "Levallois",
      TRUE ~ Typology
    ),
    Typology = droplevels(as.factor(Typology))
  ) %>%
  filter(complete.cases(SPI, E, I))

y_rei <- results_typed$Typology

cat(sprintf("Direction-metric specimens available: %d, classes: %d\n", nrow(results_typed), nlevels(y_rei)))
print(table(y_rei))

# --- SPI：Kruskal-Wallis + Dunn（Holm）---
kw_SPI <- kruskal.test(SPI ~ Typology, data = results_typed)
cat("\n--- Kruskal-Wallis test (SPI) ---\n")
print(kw_SPI)

cat("\n--- Dunn post-hoc (Holm) : SPI ---\n")
dunn_SPI <- FSA::dunnTest(SPI ~ Typology, data = results_typed, method = "holm")
print(dunn_SPI)

# --- SPI boxplot ---
p_SPI_box <- ggplot(results_typed,
                    aes(x = factor(Typology, levels = TYPOLOGY_ORDER),
                        y = SPI, fill = Typology, color = Typology)) +
  geom_boxplot(outlier.shape = 21, outlier.size = 2.5,
               alpha = 0.25, linewidth = 0.5) +
  geom_jitter(width = 0.15, size = 2.5, alpha = 0.7, shape = 16) +
  stat_summary(fun = mean, geom = "point",
               shape = 16, size = 4, color = "white") +
  annotate("text", x = Inf, y = Inf,
           label = "SPI: Kruskal-Wallis\nP < 0.001",
           hjust = 1.05, vjust = 1.2, size = 4, color = "grey40") +
  scale_fill_manual(values  = TYPOLOGY_COLORS) +
  scale_color_manual(values = TYPOLOGY_COLORS) +
  scale_x_discrete(
    limits = TYPOLOGY_ORDER,
    labels = c(
      "Unidirectional" = "Uni.",
      "Bidirectional"  = "Bi.",
      "Levallois"      = "Lev.",
      "Discoid"        = "Dis.",
      "Multiplatform"  = "Multi."
    ),
    expand = expansion(add = 0.6)
  ) +
  scale_y_continuous(limits = c(0.08, 1.0), breaks = seq(0.0, 1.0, by = 0.1)) +
  theme_bw() +
  theme(
    panel.grid     = element_blank(),
    axis.text.x    = element_text(size = 9.5),
    axis.text.y    = element_text(size = 9.5),
    legend.position = "none"
  ) +
  labs(x = NULL, y = "SPI")

# --- PERMANOVA + PERMDISP helper ---
run_permanova <- function(X, group_vec, label) {
  df_grp <- data.frame(Typology = group_vec)
  d      <- dist(X, method = "euclidean")
  
  set.seed(42)
  perm_global <- adonis2(X ~ Typology, data = df_grp,
                         method = "euclidean", permutations = 999)
  cat(sprintf("\n--- PERMANOVA : %s ---\n", label)); print(perm_global)
  
  set.seed(42)
  disp      <- betadisper(d, group_vec)
  disp_test <- permutest(disp, permutations = 999)
  cat(sprintf("\n--- PERMDISP : %s ---\n", label)); print(disp_test)
  
  disp_tukey <- TukeyHSD(disp)
  cat(sprintf("\n--- PERMDISP TukeyHSD : %s ---\n", label)); print(disp_tukey)
  
  set.seed(42)
  pairwise_res <- pairwise.perm.manova(d, group_vec, nperm = 999, p.method = "holm")
  cat(sprintf("\n--- Pairwise PERMANOVA (Holm) : %s ---\n", label))
  print(pairwise_res$p.value)
  
  invisible(list(global    = perm_global,
                 disp      = disp,
                 disp_test = disp_test,
                 disp_tukey = disp_tukey,
                 pairwise  = pairwise_res))
}

# --- Fabric: Benn ternary + PERMANOVA ---
X_EI    <- results_typed %>% select(E, I) %>% as.matrix()
perm_EI <- run_permanova(X_EI, y_rei, "Fabric (E+I)")

# Benn ternary coordinates (top = IS, left = PL, right = EL)
df_ternary <- results_typed %>%
  mutate(
    EL = 1 - (lambda2 / lambda1), 
    IS = lambda3 / lambda1,        
    PL = 1 - EL - IS,              
    Typology = factor(Typology, levels = TYPOLOGY_ORDER)
  )

p_benn <- ggtern(
  data = df_ternary,
  aes(x = PL, z = EL, y = IS, color = Typology, fill = Typology)
) +
  geom_polygon(
    data = df_ternary %>%
      group_by(Typology) %>%
      slice(chull(PL, EL)) %>%
      ungroup(),
    aes(group = Typology),
    alpha = 0.25, color = NA
  ) +
  geom_point(size = 1.8, alpha = 0.85, shape = 16) +
  scale_color_manual(values = TYPOLOGY_COLORS) +
  scale_fill_manual(values  = TYPOLOGY_COLORS) +
  labs(
    x     = "Planar",
    z     = "Linear",
    y     = "Isotropic",
    color = "Fabric: PERMANOVA\nP = 0.001",
    fill  = "Fabric: PERMANOVA\nP = 0.001"
  ) +
  ggtern::theme_bw(base_size = 4) +
  theme(
    tern.panel.grid.major = element_line(color = "grey85", linewidth = 0.15, linetype = "dashed"),
    tern.panel.grid.minor = element_line(color = "grey85", linewidth = 0.15, linetype = "dashed"),
    tern.axis.title.T     = element_text(size = 10, color = "grey20"),
    tern.axis.title.L     = element_text(size = 10, color = "grey20"),
    tern.axis.title.R     = element_text(size = 10, color = "grey20"),
    tern.axis.text.T = element_text(size = 9),
    tern.axis.text.L = element_text(size = 9),
    tern.axis.text.R = element_text(size = 9),
    legend.position       = "right",
    legend.title          = element_text(size = 11, color = "grey40", 
                                         margin = margin(b = 10)),
    legend.text           = element_text(size = 10, color = "grey40"),
    legend.key.size       = unit(1, "lines"),
    legend.background     = element_rect(fill = "transparent", colour = NA),
    legend.box.background = element_rect(fill = "transparent", colour = NA),
    plot.margin           = ggplot2::margin(2, 2, 2, 2)
  )
p_benn

# --- SPHARM PERMANOVA ---
perm_morph <- run_permanova(z_morph[non_im_idx, ], y_typology, "SPHARM morphology spectrum")
perm_dir   <- run_permanova(z_dir[non_im_idx, ],   y_typology, "SPHARM direction spectrum")

cat("\n========== PERMANOVA + PERMDISP summary ==========\n")
cat(sprintf("SPHARM morphology : R2 = %.3f, p = %.3f | PERMDISP p = %.3f\n",
            perm_morph$global$R2[1],
            perm_morph$global$`Pr(>F)`[1],
            perm_morph$disp_test$tab$`Pr(>F)`[1]))
cat(sprintf("SPHARM direction  : R2 = %.3f, p = %.3f | PERMDISP p = %.3f\n",
            perm_dir$global$R2[1],
            perm_dir$global$`Pr(>F)`[1],
            perm_dir$disp_test$tab$`Pr(>F)`[1]))
cat(sprintf("Fabric (E+I): R² = %.3f, p = %.3f | PERMDISP p = %.3f\n",
            perm_EI$global$R2[1],
            perm_EI$global$`Pr(>F)`[1],
            perm_EI$disp_test$tab$`Pr(>F)`[1]))
cat(sprintf("SPI KW      : H  = %.2f, df = %d, p < 0.001\n",
            kw_SPI$statistic, kw_SPI$parameter))

# ==============================================================================
# 5. Bubble plot: pairwise post-hoc summary
# ==============================================================================

pairs_full <- c(
  "Bidirectional –\nDiscoidal",    "Bidirectional –\nLevallois",
  "Bidirectional –\nMultiplatform","Bidirectional –\nUnidirectional",
  "Discoidal –\nLevallois",        "Discoidal –\nMultiplatform",
  "Discoidal –\nUnidirectional",   "Levallois –\nMultiplatform",
  "Levallois –\nUnidirectional",   "Multiplatform –\nUnidirectional"
)

method_levels <- c("SPI", "Fabric", "SP-SPHARM")

p_adj_vec <- c(
  # SPI Dunn（Holm）
  1.0000, 0.6534, 0.9665, 0.0058,
  0.1455, 0.1572, 0.0001,
  0.8816, 0.0040, 0.0421,
  # Fabric E+I（Holm）
  0.010, 0.010, 0.357, 0.072,
  0.055, 0.357, 0.010,
  0.804, 0.010, 0.012,
  # Direction SPHARM（Holm）
  0.458, 0.010, 0.045, 0.010,
  0.012, 0.458, 0.010,
  0.012, 0.028, 0.010
)

plot_data <- tibble(
  Pair   = factor(rep(pairs_full, 3), levels = rev(pairs_full)),
  Method = factor(rep(method_levels, each = length(pairs_full)),
                  levels = method_levels),
  p_adj  = p_adj_vec
) %>%
  mutate(
    Sig = case_when(
      p_adj <= 0.001 ~ "p ≤ 0.001",
      p_adj <= 0.01  ~ "p ≤ 0.01",
      p_adj <= 0.05  ~ "p ≤ 0.05",
      TRUE           ~ "ns"
    ),
    Sig = factor(Sig, levels = c("p ≤ 0.001", "p ≤ 0.01", "p ≤ 0.05", "ns"))
  )

p_bubble <- ggplot(
  plot_data %>% filter(Sig != "ns"),
  aes(x = Method, y = Pair)
) +
  geom_point(aes(fill = Sig),
             shape = 21, size = 10, color = "white", stroke = 0.1, alpha = 0.7) +
  geom_text(aes(label = case_when(
    Sig == "p ≤ 0.001" ~ "***",
    Sig == "p ≤ 0.01"  ~ "**",
    TRUE               ~ "*"
  ), color = Sig),
  size = 5, fontface = "bold") +
  scale_fill_manual(values  = c("p ≤ 0.001" = "#802520",
                                "p ≤ 0.01"  = "#B26538",
                                "p ≤ 0.05"  = "#BA8530"),
                    guide = "none") +
  scale_color_manual(values = c("p ≤ 0.001" = "#F5EDDC",
                                "p ≤ 0.01"  = "#F5EDDC",
                                "p ≤ 0.05"  = "#F5EDDC"),
                     guide = "none") +
  scale_x_discrete(position = "top", limits = method_levels) +
  scale_y_discrete(position = "right", limits = levels(plot_data$Pair),
                   expand = expansion(add = 0.6)) +
  theme_bw() +
  theme(
    panel.grid.major = element_line(color = "gray50", linewidth = 0.35,
                                    linetype = "dashed"),
    panel.grid.minor = element_blank(),
    axis.text.x      = element_text(face = "bold", size = 10, hjust = 0.5,
                                    margin = ggplot2::margin(b = 4)),
    axis.text.y      = element_text(size = 10, hjust = 0,
                                    margin = ggplot2::margin(r = -1)),
    axis.title       = element_blank(),
    legend.position  = "none",
    plot.margin      = ggplot2::margin(6, 6, 6, 6)
  )

# ==============================================================================
# 6. Save filtered data
# ==============================================================================

saveRDS(SPHARM_direction_filter,
        here("analysis/data/derived_data/SPHARM_direction_filter.rds"))
saveRDS(SPHARM_morphology_filter,
        here("analysis/data/derived_data/SPHARM_morphology_filter.rds"))

cat("Saved: SPHARM_direction_filter.rds\n")
cat("Saved: SPHARM_morphology_filter.rds\n")

# ==============================================================================
# 7. Combined panel: three left plots + bubble plot
# ==============================================================================
benn_grob <- ggplotGrob(p_benn)
p_benn_wrap <- wrap_elements(full = benn_grob)

hull_dir <- res_dir$scores %>%
  group_by(Typology) %>%
  slice(chull(LD1, LD2)) %>%
  ungroup()

p_dir_plot <- ggplot(res_dir$scores,
                     aes(x = LD1, y = LD2, color = Typology)) +
  geom_hline(yintercept = 0, color = "grey50", linewidth = 0.35, linetype = "dashed") +
  geom_vline(xintercept = 0, color = "grey50", linewidth = 0.35, linetype = "dashed") +
  annotate("text", x = Inf, y = Inf,
           label = "SP-SPHARM: PERMANOVA\nP = 0.001",
           hjust = 1.05, vjust = 1.2, size = 4, color = "grey40") +
  geom_polygon(data = hull_dir, aes(fill = Typology, group = Typology),
               alpha = 0.25, color = NA) +
  geom_polygon(data = hull_dir, aes(color = Typology, group = Typology),
               fill = NA, linewidth = 0.01) +
  geom_point(size = 2.0, alpha = 0.88, stroke = 0.3, shape = 16) +
  scale_color_manual(values = TYPOLOGY_COLORS) +
  scale_fill_manual(values  = TYPOLOGY_COLORS) +
  scale_x_continuous(limits = c(-2.5, 4.5), expand = expansion(mult = 0.08),
                     breaks = seq(-2.5, 4.5, by = 1)) +
  scale_y_continuous(limits = c(-5, 3), expand = expansion(mult = 0.08),
                     breaks = seq(-5, 3, by = 1)) +
  labs(x = sprintf("LD1 (%.1f%%)", res_dir$prop_var[1] * 100),
       y = sprintf("LD2 (%.1f%%)", res_dir$prop_var[2] * 100)) +
  theme_bw() +
  theme(panel.grid = element_blank(), legend.position = "none")

left_col <- (p_SPI_box / p_benn_wrap / p_dir_plot) +
  plot_layout(ncol = 1, heights = c(1, 1, 1))

exp_method_compare_combined <- (left_col | p_bubble) +
  plot_layout(ncol = 2, widths = c(1, 0.5)) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(plot.tag = element_text(face = "bold"))
  )
