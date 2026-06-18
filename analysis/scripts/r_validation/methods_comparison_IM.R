# methods_comparison_IM.R
# Method-validity check on ideal models (IM). Sourced by the `im_comparison` target.
#
# Pipeline:
#   Part A  Direction metrics — R (mean resultant length), E (elongation), I (isotropy).
#   Part B  SPHARM power spectra (IM only) — per-degree variance and per-type spectra.
#   Part C  Discriminative-power comparison — merge R/E/I with SPHARM power, then
#           pairwise standardised Euclidean distance heatmap (SPI vs fabric vs SPHARM).
#
# Input:
#   - analysis/data/raw_data/Scar_orientation_data.xlsx (sheet 1: IM specimens)
#   - analysis/data/derived_data/SPHARM_direction.csv
#
# Returns (objects, no figures written to disk): avg_dist, p_heatmap

library(here)
library(tidyverse)
library(readxl)
library(ggrepel)
library(grid)
library(patchwork)
conflicted::conflicts_prefer(dplyr::select)
conflicted::conflicts_prefer(dplyr::filter)
conflicted::conflicts_prefer(base::`%*%`)

# ==============================================================================
# Helpers
# ==============================================================================

compute_R <- function(dx, dy, dz) {
  resultant_magnitude <- sqrt(sum(dx)^2 + sum(dy)^2 + sum(dz)^2)
  total_length        <- sum(sqrt(dx^2 + dy^2 + dz^2))
  resultant_magnitude / total_length
}

compute_EI <- function(ux, uy, uz) {
  n      <- length(ux)
  U      <- cbind(ux, uy, uz)
  T_mat  <- (t(U) %*% U) / n
  eig    <- eigen(T_mat, symmetric = TRUE)
  lambda <- sort(eig$values, decreasing = TRUE)
  lambda <- pmax(lambda, 0)
  list(
    E = ifelse(lambda[1] > 1e-10, 1 - lambda[2] / lambda[1], NA_real_),
    I = ifelse(lambda[1] > 1e-10,     lambda[3] / lambda[1], NA_real_)
  )
}

GeomRoundTile <- ggplot2::ggproto(
  "GeomRoundTile", ggplot2::GeomTile,
  draw_panel = function(self, data, panel_params, coord,
                        radius = grid::unit(2, "pt")) {
    coords <- coord$transform(data, panel_params)
    grobs <- lapply(seq_len(nrow(coords)), function(i) {
      a <- coords$alpha[i];     if (is.null(a) || is.na(a)) a <- 1
      lw <- coords$linewidth[i]; if (is.null(lw))           lw <- 0.1
      grid::roundrectGrob(
        x = coords$xmin[i], y = coords$ymin[i],
        width  = coords$xmax[i] - coords$xmin[i],
        height = coords$ymax[i] - coords$ymin[i],
        just = c("left", "bottom"), r = radius,
        gp = grid::gpar(
          col  = coords$colour[i],
          fill = scales::alpha(coords$fill[i], a),
          lwd  = lw * ggplot2::.pt,
          lty  = coords$linetype[i] %||% 1
        )
      )
    })
    grid::gTree(children = do.call(grid::gList, grobs))
  }
)

geom_round_tile <- function(mapping = NULL, data = NULL, stat = "identity",
                            position = "identity", ..., radius = grid::unit(2, "pt"),
                            na.rm = FALSE, show.legend = NA, inherit.aes = TRUE) {
  ggplot2::layer(
    geom = GeomRoundTile, mapping = mapping, data = data, stat = stat,
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(radius = radius, na.rm = na.rm, ...)
  )
}

# ==============================================================================
# 1. Load data
# ==============================================================================

raw <- read_excel(
  here("analysis/data/raw_data/Scar_orientation_data.xlsx"), sheet = 1)

SPHARM_direction <- read_csv(
  here("analysis/data/derived_data/SPHARM_direction.csv"),
  show_col_types = FALSE)

# ==============================================================================
# Part A: direction metrics (R, E, I) — all specimens
# ==============================================================================

raw_dirs <- raw %>%
  mutate(
    dx     = End_X - Start_X,
    dy     = End_Y - Start_Y,
    dz     = End_Z - Start_Z,
    length = sqrt(dx^2 + dy^2 + dz^2)
  ) %>%
  filter(length > 1e-10) %>%
  mutate(
    ux = dx / length,
    uy = dy / length,
    uz = dz / length
  )

results <- raw_dirs %>%
  group_by(ID) %>%
  summarise(
    n_scars = n(),
    R       = compute_R(ux, uy, uz),
    E       = compute_EI(ux, uy, uz)$E,
    I       = compute_EI(ux, uy, uz)$I,
    .groups = "drop"
  ) %>%
  arrange(ID)

results %>%
  mutate(across(c(R, E, I), \(x) round(x, 4))) %>%
  print(n = Inf)

# ==============================================================================
# Part B: SPHARM power spectra — IM specimens only
# ==============================================================================

SPHARM_IM <- SPHARM_direction %>%
  filter(str_starts(ID, "IM_")) %>%
  select(ID, Typology, power_l1:power_l20)

cat("Ideal-model specimens:", nrow(SPHARM_IM), "\n")
print(SPHARM_IM$ID)

# --- B-1: per-degree variance ---
variance_IM <- SPHARM_IM %>%
  select(starts_with("power_l")) %>%
  summarise(across(everything(), var)) %>%
  pivot_longer(cols = everything(),
               names_to = "degree_label", values_to = "variance") %>%
  mutate(degree = as.integer(str_remove(degree_label, "power_l"))) %>%
  arrange(degree) %>%
  mutate(
    var_pct    = variance / sum(variance) * 100,
    var_cumsum = cumsum(var_pct)
  )

variance_IM %>%
  mutate(across(c(var_pct, var_cumsum), \(x) round(x, 2))) %>%
  print(n = Inf)

# --- B-2: per-type power spectra ---
df_long <- SPHARM_IM %>%
  pivot_longer(
    cols      = starts_with("power_l"),
    names_to  = "degree_label",
    values_to = "power"
  ) %>%
  mutate(
    degree   = as.integer(str_remove(degree_label, "power_l")),
    ID_label = ID %>%
      str_remove("^IM_") %>%
      str_replace_all("_", " ") %>%
      str_to_sentence()
  )

y_min <- min(df_long$power, na.rm = TRUE)
y_max <- max(df_long$power, na.rm = TRUE)

df_long %>%
  group_by(ID_label) %>%
  group_split() %>%
  walk(function(df_sub) {
    
    type_name <- as.character(df_sub$ID_label[1])
    
    p <- ggplot(df_sub, aes(x = degree, y = power)) +
      geom_line(color = "#B26538", linewidth = 0.5, alpha = 0.7) +
      geom_point(color = "#B26538", size = 2.5, shape = 16, alpha = 0.7) +
      scale_x_continuous(breaks = 1:20) +
      scale_y_continuous(
        labels = scales::label_scientific(digits = 2),
        limits = c(y_min, y_max)
      ) +
      theme_classic() +
      labs(
        title = type_name,
        x     = "Spherical Harmonic Degree (l)",
        y     = "Normalised Power"
      ) +
      theme(
        plot.title         = element_text(face = "bold", size = 13, hjust = 0.5),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor   = element_blank(),
        axis.text.x        = element_text(size = 9),
        axis.text.y        = element_text(size = 9)
      )
    
    file_name <- type_name %>%
      str_replace_all(" ", "_") %>%
      str_to_lower()
    
    
    cat("Built:", type_name, "\n")
  })

# ==============================================================================
# Part C: discriminative-power comparison
# ==============================================================================

# --- C-1: IM fabric (R/E/I) ---
fabric_IM <- raw %>%
  filter(str_starts(ID, "IM_")) %>%
  mutate(
    dx  = End_X - Start_X,
    dy  = End_Y - Start_Y,
    dz  = End_Z - Start_Z,
    len = sqrt(dx^2 + dy^2 + dz^2)
  ) %>%
  filter(len > 1e-10) %>%
  mutate(ux = dx / len, uy = dy / len, uz = dz / len) %>%
  group_by(ID) %>%
  summarise(
    R = compute_R(ux, uy, uz),
    E = compute_EI(ux, uy, uz)$E,
    I = compute_EI(ux, uy, uz)$I,
    .groups = "drop"
  )

# --- C-2: merge datasets ---
id_order <- c(
  "Cylindrical unipolar cortical", "Cylindrical unipolar scarred",
  "Cylindrical bipolar",
  "Conical unipolar cortical",     "Conical unipolar scarred",
  "Discoid",                       "Discoid unifacial",
  "Levallois preferential",        "Levallois convergent",
  "Levallois laminar",
  "Biface",                        "Multiplatform"
)

df_im <- SPHARM_IM %>%
  left_join(fabric_IM, by = "ID") %>%
  mutate(
    label = ID %>%
      str_remove("^IM_") %>%
      str_replace_all("_", " ") %>%
      str_to_sentence()
  )

cat("===== Merged dataset (n =", nrow(df_im), ") =====\n")
print(df_im %>% select(label, R, E, I, power_l1:power_l4), n = Inf)

# --- C-3: pairwise Euclidean distance heatmap ---
make_dist_df <- function(X, labels, method_name) {
  d <- as.matrix(dist(scale(X), method = "euclidean"))
  rownames(d) <- colnames(d) <- labels
  as.data.frame(d) %>%
    rownames_to_column("From") %>%
    pivot_longer(-From, names_to = "To", values_to = "distance") %>%
    mutate(method = method_name)
}

labs <- df_im$label

dist_all <- bind_rows(
  make_dist_df(df_im %>% select(R)                 %>% as.matrix(), labs, "SPI"),
  make_dist_df(df_im %>% select(E, I)              %>% as.matrix(), labs, "Fabric"),
  make_dist_df(df_im %>% select(power_l1:power_l4) %>% as.matrix(), labs, "SPHARM")
) %>%
  mutate(
    From   = factor(From,   levels = id_order),
    To     = factor(To,     levels = id_order),
    method = factor(method, levels = c("SPI", "Fabric", "SPHARM"))
  )

avg_dist <- dist_all %>%
  filter(From != To) %>%
  group_by(method) %>%
  summarise(mean_dist = round(mean(distance), 2), .groups = "drop")

cat("\n===== Standardised mean inter-class distance per method =====\n")
cat("(larger = stronger overall discrimination)\n")
print(avg_dist)

dist_all_upper <- dist_all %>%
  dplyr::filter(as.numeric(From) < as.numeric(To))

p_heatmap <- ggplot(dist_all_upper, aes(x = To, y = From, fill = distance)) +
  geom_round_tile(color = "white", linewidth = 0.1, radius = unit(3, "pt")) +
  geom_text(aes(label = sprintf("%.1f", distance)), size = 1.7) +
  facet_wrap(~ method, ncol = 3,
             labeller = as_labeller(c(
               SPI    = "Scar Pattern Index",
               Fabric = "Fabric metrics",
               SPHARM = "SP-SPHARM"
             ))) +
  scale_fill_gradient2(
    low      = "#5C7F71",
    mid      = "#F5EDDC",
    high     = "#802520",
    midpoint = 2,
    name     = "Euclidean\ndistance\n(standardised)"
  ) +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  scale_y_discrete(limits = rev) +
  theme_bw(base_size = 6) +
  labs(x = NULL, y = NULL) +
  theme(
    strip.text       = element_text(face = "bold", size = 6.5),
    strip.background = element_rect(fill = "#EBEBEB", color = "#EBEBEB"),
    axis.text        = element_text(size = 6),
    legend.title     = element_text(size = 6),
    legend.text      = element_text(size = 6),
    legend.key.size  = unit(0.32, "cm")
  )
