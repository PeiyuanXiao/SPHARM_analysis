# sdg_projection_explore.R

suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(readxl)
  library(compositions)
  library(patchwork)
  library(conflicted)
})
suppressMessages({
  conflicts_prefer(dplyr::filter, dplyr::select, dplyr::lag,
                   stats::sd, stats::var, stats::dist, stats::cor, stats::cov,
                   base::scale, base::norm, base::`%*%`, .quiet = TRUE)
})

set.seed(42)
SEED <- 42

OUT_DIR <- here("analysis/robustness/sdg_projection_explore")
FIG_DIR <- file.path(OUT_DIR, "figures")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

POWER_COLS_DIR   <- paste0("power_l", 1:6)
POWER_COLS_MORPH <- paste0("power_l", 1:8)

SPEC_LMAX <- 12

MESH_TARGET_FACES <- 20000

COUNT_RATIO_TRIGGER  <- 1.5
LENGTH_RATIO_TRIGGER <- 1.25

# =============================================================================
# Helper import — verbatim from the joint-MFA scripts, without running them
# =============================================================================
import_helpers <- function(path, wanted) {
  if (!file.exists(path))
    stop("helper source not found: ", path,
         " — stopping rather than reimplementing its helpers here.")
  got <- character(0)
  for (e in parse(path)) {
    if (!is.call(e)) next
    if (!as.character(e[[1]])[1] %in% c("<-", "=")) next
    if (!is.name(e[[2]])) next
    nm <- as.character(e[[2]])
    if (!nm %in% wanted) next
    rhs <- e[[3]]
    if (!(is.call(rhs) && identical(as.character(rhs[[1]])[1], "function"))) next
    eval(e, envir = globalenv())
    got <- c(got, nm)
  }
  miss <- setdiff(wanted, got)
  if (length(miss))
    stop("could not import from ", basename(path), ": ",
         paste(miss, collapse = ", "),
         " — the source script's structure changed; stopping.")
  got
}

EXP_SCRIPT <- here("analysis/robustness/joint_mfa_discrimination",
                   "joint_mfa_discrimination_stats.R")
SDG_SCRIPT <- here("analysis/robustness/joint_mfa_discrimination",
                   "joint_mfa_discrimination_SDG.R")

h_exp <- import_helpers(EXP_SCRIPT,
  c("replace_zeros", "make_ilr", "ilr_parts", "mfa_s1", "block_inertia",
    "filter_spharm", "split_by_group", "safe_filter_groups"))
cat(sprintf("Imported %d helpers verbatim from %s:\n  %s\n",
            length(h_exp), basename(EXP_SCRIPT), paste(sort(h_exp), collapse = ", ")))

# =============================================================================
# Inputs
# =============================================================================
cat("\nLoading committed main outputs...\n")

SPHARM_direction  <- read_csv(
  here("analysis/data/derived_data/SPHARM_direction.csv"),  show_col_types = FALSE)
SPHARM_morphology <- read_csv(
  here("analysis/data/derived_data/SPHARM_morphology.csv"), show_col_types = FALSE)
SPHARM_morphology <- SPHARM_morphology %>%
  left_join(SPHARM_direction %>% select(ID, Typology), by = "ID")

metric_data <- read_xlsx(here("analysis/data/raw_data/SDG_core_metric.xlsx")) %>%
  select(ID, Layer, Core_type_Li_merged, Raw_mat) %>%
  mutate(Layer = as.factor(Layer), Raw_mat = as.factor(Raw_mat))

assemblage_of <- function(id)
  case_when(str_starts(id, "EXP") ~ "EXP",
            str_starts(id, "SDG") ~ "SDG",
            str_starts(id, "IM_") ~ "IM",
            TRUE ~ NA_character_)

# =============================================================================
# STEP 0 — BATCH-EFFECT CHECK
# =============================================================================
cat("\n", strrep("=", 78), "\n", sep = "")
cat("STEP 0 — BATCH-EFFECT CHECK (equipment and recorder differ between groups)\n")
cat(strrep("=", 78), "\n", sep = "")

# ---- 0(a) scars per specimen and scar length --------------------------------
scar_xlsx <- here("analysis/data/raw_data/Scar_orientation_data.xlsx")
scar_raw  <- map_dfr(1:3, function(i)
  read_excel(scar_xlsx, sheet = i) %>% mutate(across(everything(), as.character)))

scars <- scar_raw %>%
  mutate(ID = str_trim(ID),
         group = assemblage_of(ID),
         across(c(Start_X, Start_Y, Start_Z, End_X, End_Y, End_Z), as.numeric),
         length_mm = sqrt((End_X - Start_X)^2 +
                          (End_Y - Start_Y)^2 +
                          (End_Z - Start_Z)^2)) %>%
  filter(group %in% c("EXP", "SDG"), !is.na(length_mm))

five_num <- function(x) {
  q <- stats::quantile(x, c(0, .25, .5, .75, 1), names = FALSE, na.rm = TRUE)
  tibble(min = q[1], Q1 = q[2], median = q[3], Q3 = q[4], max = q[5],
         mean = mean(x, na.rm = TRUE), n = length(x))
}

counts_per_spec <- scars %>% count(group, ID, name = "n_scars")
tab_counts <- counts_per_spec %>% group_by(group) %>%
  reframe(five_num(n_scars)) %>% rename(n_specimens = n)
tab_lengths <- scars %>% group_by(group) %>%
  reframe(five_num(length_mm)) %>% rename(n_scars = n)

cat("\n(a) SCARS PER SPECIMEN\n")
print(as.data.frame(tab_counts %>% mutate(across(where(is.numeric), ~ round(.x, 2)))))
cat("\n(a) SCAR LENGTH (mm)\n")
print(as.data.frame(tab_lengths %>% mutate(across(where(is.numeric), ~ round(.x, 2)))))

gv <- function(tab, g, col) tab[[col]][tab$group == g]
count_ratio  <- gv(tab_counts,  "EXP", "median") / gv(tab_counts,  "SDG", "median")
length_ratio <- gv(tab_lengths, "SDG", "median") / gv(tab_lengths, "EXP", "median")

ks_counts  <- suppressWarnings(stats::ks.test(
  counts_per_spec$n_scars[counts_per_spec$group == "EXP"],
  counts_per_spec$n_scars[counts_per_spec$group == "SDG"]))
ks_lengths <- suppressWarnings(stats::ks.test(
  scars$length_mm[scars$group == "EXP"],
  scars$length_mm[scars$group == "SDG"]))

cat(sprintf("\n  median scars/specimen  EXP/SDG ratio = %.2f  (trigger > %.2f)\n",
            count_ratio, COUNT_RATIO_TRIGGER))
cat(sprintf("  median scar length     SDG/EXP ratio = %.2f  (trigger > %.2f)\n",
            length_ratio, LENGTH_RATIO_TRIGGER))
cat(sprintf("  KS scars/specimen D = %.3f, p = %.3g | KS length D = %.3f, p = %.3g\n",
            ks_counts$statistic, ks_counts$p.value,
            ks_lengths$statistic, ks_lengths$p.value))

NEED_T10 <- (count_ratio > COUNT_RATIO_TRIGGER) ||
            (length_ratio > LENGTH_RATIO_TRIGGER) ||
            (1 / length_ratio > LENGTH_RATIO_TRIGGER)
cat(sprintf("\n  => >10 mm control required: %s\n", ifelse(NEED_T10, "YES", "no")))

p_counts <- ggplot(counts_per_spec, aes(group, n_scars, fill = group)) +
  geom_boxplot(width = .5, outlier.size = .8) +
  labs(title = "Step 0(a) scars per specimen", y = "scars", x = NULL) +
  theme_bw() + theme(legend.position = "none")
p_len <- ggplot(scars, aes(length_mm, colour = group)) +
  geom_density() + scale_x_log10() +
  geom_vline(xintercept = c(2, 10), linetype = "dashed", linewidth = .3) +
  labs(title = "Step 0(a) scar length (log scale; dashed = 2 and 10 mm)",
       x = "length (mm)") +
  theme_bw()

# ---- 0(b) mesh resolution after decimation ----------------------------------
mesh_csv <- file.path(OUT_DIR, "mesh_resolution.csv")
if (!file.exists(mesh_csv)) {
  cat("\n(b) SKIPPED — run mesh_resolution_check.py first (see header).\n")
  tab_mesh <- NULL
} else {
  mesh <- read_csv(mesh_csv, show_col_types = FALSE)
  tab_mesh <- mesh %>% group_by(group) %>%
    summarise(n = n(),
              faces_raw_med    = median(faces_raw),
              faces_raw_min    = min(faces_raw),
              faces_dec_med    = median(faces_decimated),
              faces_dec_min    = min(faces_decimated),
              faces_dec_max    = max(faces_decimated),
              verts_dec_med    = median(vertices_decimated),
              verts_dec_min    = min(vertices_decimated),
              verts_dec_max    = max(vertices_decimated),
              n_below_target   = sum(!decimation_active),
              .groups = "drop")
  cat("\n(b) MESH RESOLUTION (production decimation target = ",
      format(MESH_TARGET_FACES, big.mark = ","), " faces)\n", sep = "")
  print(as.data.frame(tab_mesh))
  cat(sprintf("\n  specimens whose raw mesh is already below the target: EXP %d, SDG %d\n",
              tab_mesh$n_below_target[tab_mesh$group == "EXP"],
              tab_mesh$n_below_target[tab_mesh$group == "SDG"]))
}

# ---- 0(c) mean power spectrum by degree -------------------------------------
spec_long <- function(df, cols_prefix, lmax, descriptor) {
  keep <- intersect(paste0(cols_prefix, 1:lmax), names(df))
  df %>%
    mutate(group = assemblage_of(ID)) %>%
    filter(group %in% c("EXP", "SDG")) %>%
    select(ID, group, all_of(keep)) %>%
    pivot_longer(all_of(keep), names_to = "l", values_to = "power") %>%
    mutate(l = as.integer(str_remove(l, cols_prefix)), descriptor = descriptor)
}

spec_all <- bind_rows(
  spec_long(SPHARM_morphology, "power_l", SPEC_LMAX, "M-SPHARM"),
  spec_long(SPHARM_direction,  "power_l", SPEC_LMAX, "SP-SPHARM"))

spec_summary <- spec_all %>%
  group_by(descriptor, group, l) %>%
  summarise(mean_power = mean(power), se = sd(power) / sqrt(n()),
            n = n(), .groups = "drop")

cat("\n(c) MEAN POWER BY DEGREE (l = 1-", SPEC_LMAX, ")\n", sep = "")
print(as.data.frame(spec_summary %>%
  select(descriptor, group, l, mean_power) %>%
  pivot_wider(names_from = group, values_from = mean_power) %>%
  mutate(ratio_SDG_EXP = SDG / EXP,
         across(c(EXP, SDG), ~ signif(.x, 4)),
         ratio_SDG_EXP = round(ratio_SDG_EXP, 3))))

hi <- spec_summary %>%
  filter(l > SPEC_LMAX / 2) %>%
  select(descriptor, group, l, mean_power) %>%
  pivot_wider(names_from = group, values_from = mean_power) %>%
  group_by(descriptor) %>%
  summarise(mean_log2_ratio_hi = mean(log2(SDG / EXP)), .groups = "drop")
cat("\n  high-degree (l > ", SPEC_LMAX / 2, ") mean log2(SDG/EXP) power ratio:\n", sep = "")
print(as.data.frame(hi %>% mutate(mean_log2_ratio_hi = round(mean_log2_ratio_hi, 3))))

inrange <- spec_summary %>%
  mutate(lmax_used = if_else(descriptor == "M-SPHARM", 8L, 6L)) %>%
  filter(l <= lmax_used) %>%
  select(descriptor, group, l, mean_power) %>%
  pivot_wider(names_from = group, values_from = mean_power) %>%
  group_by(descriptor) %>%
  summarise(mean_abs_log2_ratio = mean(abs(log2(SDG / EXP))),
            max_abs_log2_ratio  = max(abs(log2(SDG / EXP))),
            worst_l = l[which.max(abs(log2(SDG / EXP)))], .groups = "drop")
cat("\n  WITHIN the retained truncation (M l=1-8, SP l=1-6):\n")
print(as.data.frame(inrange %>% mutate(across(where(is.numeric), ~ round(.x, 3)))))

sp_by_spec <- SPHARM_direction %>%
  mutate(group = assemblage_of(ID)) %>%
  filter(group %in% c("EXP", "SDG")) %>%
  select(ID, group, all_of(paste0("power_l", 1:SPEC_LMAX))) %>%
  pivot_longer(-c(ID, group), names_to = "l", values_to = "power") %>%
  mutate(l = as.integer(str_remove(l, "power_l"))) %>%
  group_by(ID, group) %>%
  summarise(hi_share = sum(power[l > 6]) / sum(power),
            p3 = power[l == 3], p4 = power[l == 4], .groups = "drop") %>%
  left_join(counts_per_spec %>% select(ID, n_scars), by = "ID")

ct <- function(d, v) {
  r <- suppressWarnings(stats::cor.test(d$n_scars, d[[v]], method = "spearman"))
  sprintf("rho = %+.3f (p = %.3g)", r$estimate, r$p.value)
}
exp_only <- sp_by_spec %>% filter(group == "EXP", !is.na(n_scars))
sdg_only <- sp_by_spec %>% filter(group == "SDG", !is.na(n_scars))
cat("\n  SP-SPHARM vs scars per specimen, Spearman, WITHIN each batch:\n")
cat(sprintf("    EXP (n = %d)  high-degree share %s | power_l3 %s | power_l4 %s\n",
            nrow(exp_only), ct(exp_only, "hi_share"),
            ct(exp_only, "p3"), ct(exp_only, "p4")))
cat(sprintf("    SDG (n = %d)  high-degree share %s | power_l3 %s | power_l4 %s\n",
            nrow(sdg_only), ct(sdg_only, "hi_share"),
            ct(sdg_only, "p3"), ct(sdg_only, "p4")))
write_csv(sp_by_spec, file.path(OUT_DIR, "step0_sp_spectrum_vs_scarcount.csv"))

p_spec <- ggplot(spec_summary, aes(l, mean_power, colour = group)) +
  geom_line() + geom_point(size = 1) +
  facet_wrap(~ descriptor, scales = "free_y") +
  scale_y_log10() +
  labs(title = "Step 0(c) mean normalised power by degree (log scale)",
       y = "mean power", x = "degree l") +
  theme_bw()

ggsave(file.path(FIG_DIR, "step0_batch_check.png"),
       p_counts + p_len + p_spec + plot_layout(ncol = 1),
       width = 9, height = 11, dpi = 120)

write_csv(tab_counts,  file.path(OUT_DIR, "step0_scars_per_specimen.csv"))
write_csv(tab_lengths, file.path(OUT_DIR, "step0_scar_lengths.csv"))
write_csv(spec_summary, file.path(OUT_DIR, "step0_mean_power_by_degree.csv"))
if (!is.null(tab_mesh)) write_csv(tab_mesh, file.path(OUT_DIR, "step0_mesh_resolution_summary.csv"))

cat("\n", strrep("-", 78), "\n", sep = "")
cat("STEP 0 VERDICT\n")
cat(strrep("-", 78), "\n", sep = "")
cat(sprintf("  scars/specimen median   EXP %g vs SDG %g  (ratio %.2f)\n",
            gv(tab_counts, "EXP", "median"), gv(tab_counts, "SDG", "median"),
            count_ratio))
cat(sprintf("  scar length median (mm) EXP %.2f vs SDG %.2f  (ratio %.2f)\n",
            gv(tab_lengths, "EXP", "median"), gv(tab_lengths, "SDG", "median"),
            gv(tab_lengths, "SDG", "median") / gv(tab_lengths, "EXP", "median")))
cat(sprintf("  >10 mm control required : %s\n", ifelse(NEED_T10, "YES", "no")))
cat("\nSteps 1-3 are gated on this verdict; see the console report.\n")

# =============================================================================
# STEP 1 — EXP REFERENCE SPACE, SDG PROJECTED AS SUPPLEMENTARY INDIVIDUALS
# =============================================================================
cat("\n", strrep("=", 78), "\n", sep = "")
cat("STEP 1 — EXP REFERENCE SPACE + SDG SUPPLEMENTARY PROJECTION\n")
cat(strrep("=", 78), "\n", sep = "")

ASSEMBLAGE         <- "EXP"
EXCLUDE_TYPES      <- c("Biface")
LEVALLOIS_MERGE    <- c("Levallois convergent", "Levallois laminar",
                        "Levallois preferential", "Levallois recurrent")
TYPOLOGY_ORDER     <- c("Unidirectional", "Bidirectional", "Levallois",
                        "Discoid", "Multiplatform")
EXCLUDE_CORE_TYPES <- c("Handaxe", "Pick")
RESTRICT_LAYERS    <- NULL
TYPOLOGY_COLORS <- c("Levallois" = "#4A6E8A", "Discoid" = "#802520",
                     "Unidirectional" = "#BA8530", "Multiplatform" = "#8A7A68",
                     "Bidirectional" = "#788C4A")

core_meta <- read_excel(here("analysis/data/raw_data/SDG_core_metric.xlsx")) %>%
  select(ID = ID, raw_material = Raw_mat, core_type = Core_type_Li_merged) %>%
  mutate(across(everything(), ~ str_trim(as.character(.))))

invisible(import_helpers(EXP_SCRIPT, "build_blocks"))
invisible(import_helpers(SDG_SCRIPT, "build_blocks_sdg"))
cat("Imported the two data-preparation recipes verbatim ",
    "(build_blocks, build_blocks_sdg).\n", sep = "")

# ---- the projection machinery -----------------------------------------------
ilr_fixed <- function(power_df, parts) {
  X <- as.matrix(power_df)[, parts, drop = FALSE]
  as.matrix(ilr(replace_zeros(X)))
}

make_space <- function(ref, sup, space) {
  parts_M  <- ilr_parts(ref$morph)
  parts_SP <- ilr_parts(ref$scar)
  ZM_r  <- ilr_fixed(ref$morph, parts_M);  ZSP_r <- ilr_fixed(ref$scar, parts_SP)
  ZM_s  <- ilr_fixed(sup$morph, parts_M);  ZSP_s <- ilr_fixed(sup$scar, parts_SP)
  sM    <- mfa_s1(ZM_r)
  sSP   <- mfa_s1(ZSP_r)
  Zr <- switch(space, M = ZM_r, SP = ZSP_r, MFA = cbind(ZM_r / sM, ZSP_r / sSP))
  Zs <- switch(space, M = ZM_s, SP = ZSP_s, MFA = cbind(ZM_s / sM, ZSP_s / sSP))
  rownames(Zr) <- ref$ids; rownames(Zs) <- sup$ids
  pca <- prcomp(Zr, center = TRUE, scale. = FALSE)
  list(pca = pca, Zr = Zr, Zs = Zs,
       scores_ref = pca$x,
       scores_sup = stats::predict(pca, Zs),
       s1 = c(M = sM, SP = sSP))
}

typicality <- function(Zr, grp, Zs) {
  G <- levels(droplevels(grp)); d <- ncol(Zr)
  cent <- lapply(G, function(g) colMeans(Zr[grp == g, , drop = FALSE]))
  names(cent) <- G
  Sw <- Reduce(`+`, lapply(G, function(g) {
    Xg <- Zr[grp == g, , drop = FALSE]
    (nrow(Xg) - 1) * stats::cov(Xg)
  })) / (nrow(Zr) - length(G))

  invertible <- tryCatch({ solve(Sw); TRUE }, error = function(e) FALSE)
  if (!invertible) {
    warning("pooled within-group covariance is singular; falling back to ",
            "Euclidean distance (reported in the output as metric = 'euclidean')")
    Sw <- diag(d)
  }
  metric <- if (invertible) "mahalanobis" else "euclidean"

  d2 <- function(X) sapply(G, function(g) stats::mahalanobis(X, cent[[g]], Sw))
  D2_ref <- d2(Zr); D2_sup <- d2(Zs)
  if (is.null(dim(D2_sup))) D2_sup <- matrix(D2_sup, nrow = nrow(Zs),
                                             dimnames = list(NULL, G))
  own <- D2_ref[cbind(seq_len(nrow(Zr)), match(as.character(grp), G))]
  p_emp <- sapply(G, function(g) {
    ref_g <- own[grp == g]
    sapply(D2_sup[, g], function(x) mean(ref_g >= x))
  })
  if (is.null(dim(p_emp))) p_emp <- matrix(p_emp, nrow = nrow(Zs),
                                           dimnames = list(NULL, G))
  list(D2 = D2_sup, p_emp = p_emp,
       p_chisq = stats::pchisq(D2_sup, df = d, lower.tail = FALSE),
       metric = metric, d = d, own_ref = own, groups = G)
}

subset_blocks <- function(b, idx)
  list(ids = b$ids[idx], morph = b$morph[idx, , drop = FALSE],
       scar = b$scar[idx, , drop = FALSE], group = droplevels(b$group[idx]))

SPACES <- c(M = "M-SPHARM only", SP = "SP-SPHARM only", MFA = "MFA merged")

run_projection <- function(dir_frame, morph_frame, tag) {
  SPHARM_direction  <<- dir_frame
  SPHARM_morphology <<- morph_frame
  b_exp <- build_blocks("EXP")
  b_sdg <- build_blocks_sdg()

  cat(sprintf("\n[%s] EXP n = %d (%d groups: %s)\n", tag,
              length(b_exp$ids), nlevels(b_exp$group),
              paste(levels(b_exp$group), collapse = ", ")))
  cat(sprintf("[%s] SDG n = %d (%d groups: %s)\n", tag,
              length(b_sdg$ids), nlevels(b_sdg$group),
              paste(levels(b_sdg$group), collapse = ", ")))

  out <- map(names(SPACES), function(sp) {
    S  <- make_space(b_exp, b_sdg, sp)
    ty <- typicality(S$Zr, b_exp$group, S$Zs)
    list(space = sp, S = S, ty = ty, b_exp = b_exp, b_sdg = b_sdg)
  })
  names(out) <- names(SPACES)
  out
}

res <- run_projection(SPHARM_direction, SPHARM_morphology, "t>2mm")

for (sp in names(SPACES)) {
  S <- res[[sp]]$S
  ev <- S$pca$sdev^2; ev <- ev / sum(ev)
  cat(sprintf("  %-14s dims = %2d | PC1-2 explain %.1f%% | metric = %s\n",
              SPACES[[sp]], ncol(S$Zr), 100 * sum(ev[1:2]),
              res[[sp]]$ty$metric))
}

# =============================================================================
# STEP 2 — WHERE DO THE SDG CORES LAND
# =============================================================================
cat("\n", strrep("=", 78), "\n", sep = "")
cat("STEP 2 — TYPICALITY OF EACH SDG CORE RELATIVE TO THE EXP TYPES\n")
cat(strrep("=", 78), "\n", sep = "")

sdg_meta <- tibble(ID = res$MFA$b_sdg$ids,
                   sdg_type = as.character(res$MFA$b_sdg$group)) %>%
  left_join(core_meta %>% select(ID, raw_material), by = "ID") %>%
  left_join(metric_data %>% transmute(ID, layer = as.character(Layer)), by = "ID")

TYP_ALPHA <- 0.05

tidy_typ <- function(r, space) {
  G <- r$ty$groups
  as_tibble(r$ty$p_emp) %>% mutate(ID = r$b_sdg$ids, .before = 1) %>%
    pivot_longer(all_of(G), names_to = "exp_type", values_to = "p_emp") %>%
    left_join(as_tibble(r$ty$D2) %>% mutate(ID = r$b_sdg$ids, .before = 1) %>%
                pivot_longer(all_of(G), names_to = "exp_type", values_to = "D2"),
              by = c("ID", "exp_type")) %>%
    left_join(as_tibble(r$ty$p_chisq) %>% mutate(ID = r$b_sdg$ids, .before = 1) %>%
                pivot_longer(all_of(G), names_to = "exp_type", values_to = "p_chisq"),
              by = c("ID", "exp_type")) %>%
    mutate(space = space)
}
typ_long <- map_dfr(names(SPACES), ~ tidy_typ(res[[.x]], .x)) %>%
  left_join(sdg_meta, by = "ID") %>%
  mutate(space = factor(space, levels = names(SPACES)),
         exp_type = factor(exp_type, levels = TYPOLOGY_ORDER))

nearest <- typ_long %>% group_by(space, ID) %>%
  slice_max(p_emp, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(space, ID, sdg_type, raw_material, layer,
         nearest_exp_type = exp_type, nearest_p_emp = p_emp, nearest_D2 = D2)

sel <- typ_long %>% group_by(space, ID) %>%
  summarise(k = sum(p_emp > TYP_ALPHA),
            best = as.character(exp_type[which.max(p_emp)]), .groups = "drop")
cat("\n  SELECTIVITY — number of EXP types each SDG core is typical of:\n")
print(as.data.frame(sel %>% group_by(space) %>%
  summarise(median_k = median(k), mean_k = round(mean(k), 2),
            typical_of_none = sum(k == 0), typical_of_one = sum(k == 1),
            typical_of_all = sum(k == nlevels(typ_long$exp_type)),
            .groups = "drop")))
cat("\n  Levallois as the single best match, and how many of those are selective:\n")
print(as.data.frame(sel %>% group_by(space) %>%
  summarise(levallois_best = sum(best == "Levallois"),
            of_which_selective = sum(best == "Levallois" & k <= 2),
            .groups = "drop")))

cat(sprintf("\n  typicality threshold: p_emp > %.2f\n", TYP_ALPHA))
for (sp in names(SPACES)) {
  n_sdg <- length(res[[sp]]$b_sdg$ids)
  in_any <- typ_long %>% filter(space == sp) %>% group_by(ID) %>%
    summarise(any_typ = any(p_emp > TYP_ALPHA),
              lev_typ = any(p_emp > TYP_ALPHA & exp_type == "Levallois"),
              .groups = "drop")
  cat(sprintf("\n  %-14s : %d/%d SDG cores typical of at least one EXP type; %d of Levallois\n",
              SPACES[[sp]], sum(in_any$any_typ), n_sdg, sum(in_any$lev_typ)))
  cat("    nearest EXP type by SDG core type:\n")
  print(with(nearest %>% filter(space == sp),
             table(sdg_type, nearest_exp_type)))
}

# =============================================================================
# STEP 3 — EXP LEAVE-ONE-OUT BASELINE
# =============================================================================
cat("\n", strrep("=", 78), "\n", sep = "")
cat("STEP 3 — EXP LEAVE-ONE-OUT BASELINE (within-batch typicality)\n")
cat(strrep("=", 78), "\n", sep = "")

loo_baseline <- function(b_exp, space) {
  n <- length(b_exp$ids)
  map_dfr(seq_len(n), function(i) {
    ref <- subset_blocks(b_exp, -i)
    sup <- subset_blocks(b_exp,  i)
    S   <- make_space(ref, sup, space)
    ty  <- typicality(S$Zr, ref$group, S$Zs)
    own <- as.character(b_exp$group[i])
    tibble(ID = b_exp$ids[i], exp_type = own, space = space,
           p_own  = if (own %in% ty$groups) ty$p_emp[1, own] else NA_real_,
           p_max  = max(ty$p_emp[1, ]),
           nearest = ty$groups[which.max(ty$p_emp[1, ])])
  })
}

loo <- map_dfr(names(SPACES), ~ loo_baseline(res$MFA$b_exp, .x)) %>%
  mutate(space = factor(space, levels = names(SPACES)))

sdg_max <- typ_long %>% group_by(space, ID) %>%
  summarise(p_max = max(p_emp), .groups = "drop")

cat("\n  distribution of the maximum typicality probability over EXP types:\n")
cmp <- bind_rows(
  loo     %>% transmute(space, set = "EXP (leave-one-out)", p_max),
  sdg_max %>% transmute(space, set = "SDG (projected)",     p_max))
print(as.data.frame(cmp %>% group_by(space, set) %>%
  summarise(n = n(), median = median(p_max), Q1 = quantile(p_max, .25),
            Q3 = quantile(p_max, .75), frac_gt_alpha = mean(p_max > TYP_ALPHA),
            .groups = "drop") %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))))

cat("\n  Wilcoxon rank-sum, SDG vs EXP-LOO maximum typicality:\n")
for (sp in names(SPACES)) {
  a <- cmp$p_max[cmp$space == sp & cmp$set == "SDG (projected)"]
  b <- cmp$p_max[cmp$space == sp & cmp$set == "EXP (leave-one-out)"]
  w <- suppressWarnings(stats::wilcox.test(a, b))
  cat(sprintf("    %-14s W = %.0f, p = %.3g | median SDG %.3f vs EXP-LOO %.3f\n",
              SPACES[[sp]], w$statistic, w$p.value, median(a), median(b)))
}

# =============================================================================
# >10 mm CONTROL — same three spaces on the committed >10 mm scar subset
# =============================================================================
cat("\n", strrep("=", 78), "\n", sep = "")
cat(">10 mm SCAR CONTROL (committed Table S2 subset; not recomputed)\n")
cat(strrep("=", 78), "\n", sep = "")

t10_csv <- here("analysis/robustness/scar_threshold_sensitivity/spectra",
                "SPHARM_direction_t10.0.csv")
dir_t10 <- read_csv(t10_csv, show_col_types = FALSE) %>%
  slice(match(SPHARM_direction$ID, ID))
if (any(is.na(dir_t10$ID)))
  stop(">10 mm subset does not cover every specimen in SPHARM_direction.csv.")

dir_main   <- SPHARM_direction
morph_main <- SPHARM_morphology
res10 <- run_projection(dir_t10, morph_main, "t>10mm")
SPHARM_direction <<- dir_main

for (sp in names(SPACES)) {
  n_sdg <- length(res10[[sp]]$b_sdg$ids)
  p10 <- tidy_typ(res10[[sp]], sp) %>% group_by(ID) %>%
    summarise(any_typ = any(p_emp > TYP_ALPHA), .groups = "drop")
  p2 <- typ_long %>% filter(space == sp) %>% group_by(ID) %>%
    summarise(any_typ = any(p_emp > TYP_ALPHA), .groups = "drop")
  cat(sprintf("  %-14s typical of >=1 EXP type: >2 mm %d/%d | >10 mm %d/%d\n",
              SPACES[[sp]], sum(p2$any_typ), nrow(p2),
              sum(p10$any_typ), n_sdg))
}

# =============================================================================
# FIGURES (deliberately plain — this is an internal probe)
# =============================================================================
scatter_one <- function(r, sp) {
  ev <- r$S$pca$sdev^2; ev <- ev / sum(ev)
  ref <- as_tibble(r$S$scores_ref[, 1:2]) %>%
    setNames(c("PC1", "PC2")) %>% mutate(type = r$b_exp$group)
  sup <- as_tibble(r$S$scores_sup[, 1:2]) %>% setNames(c("PC1", "PC2"))
  ggplot(ref, aes(PC1, PC2)) +
    stat_ellipse(aes(colour = type), level = .68, linewidth = .3) +
    geom_point(aes(colour = type), size = 1.6) +
    geom_point(data = sup, shape = 4, colour = "grey20", size = 1.6,
               stroke = .6) +
    scale_colour_manual(values = TYPOLOGY_COLORS, name = "EXP type") +
    labs(title = SPACES[[sp]],
         subtitle = "coloured = EXP reference; x = SDG projected",
         x = sprintf("PC1 (%.1f%%)", 100 * ev[1]),
         y = sprintf("PC2 (%.1f%%)", 100 * ev[2])) +
    theme_bw(base_size = 9)
}
p_scatter <- wrap_plots(map(names(SPACES), ~ scatter_one(res[[.x]], .x)),
                        ncol = 3, guides = "collect")
ggsave(file.path(FIG_DIR, "step1_projection_scatter.png"), p_scatter,
       width = 13, height = 4.6, dpi = 130)

ord <- typ_long %>% filter(space == "MFA") %>% group_by(ID) %>%
  summarise(m = max(p_emp), .groups = "drop") %>% arrange(m) %>% pull(ID)
p_heat <- typ_long %>% mutate(ID = factor(ID, levels = ord)) %>%
  ggplot(aes(exp_type, ID, fill = p_emp)) +
  geom_tile() +
  facet_wrap(~ space, nrow = 1,
             labeller = labeller(space = SPACES)) +
  scale_fill_viridis_c(name = "typicality p", limits = c(0, 1)) +
  labs(title = "Typicality of each SDG core relative to the EXP types",
       subtitle = sprintf("empirical quantile in the EXP within-type distance distribution; ordered by MFA maximum (threshold p > %.2f)", TYP_ALPHA),
       x = NULL, y = NULL) +
  theme_bw(base_size = 8) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size = 4.5))
ggsave(file.path(FIG_DIR, "step2_typicality_heatmap.png"), p_heat,
       width = 11, height = 8, dpi = 130)

p_loo <- ggplot(cmp, aes(p_max, colour = set)) +
  stat_ecdf(linewidth = .7) +
  geom_vline(xintercept = TYP_ALPHA, linetype = "dashed", linewidth = .3) +
  facet_wrap(~ space, nrow = 1, labeller = labeller(space = SPACES)) +
  scale_colour_manual(values = c("EXP (leave-one-out)" = "#4A6E8A",
                                 "SDG (projected)"     = "#802520"),
                      name = NULL) +
  labs(title = "Step 3: maximum typicality, SDG projection vs EXP leave-one-out",
       subtitle = "dashed line = the p > 0.05 typicality threshold",
       x = "maximum typicality probability over EXP types", y = "ECDF") +
  theme_bw(base_size = 9) + theme(legend.position = "bottom")
ggsave(file.path(FIG_DIR, "step3_loo_baseline.png"), p_loo,
       width = 11, height = 4, dpi = 130)

# =============================================================================
# CSV OUTPUT
# =============================================================================
write_csv(
  typ_long %>%
    select(space, ID, sdg_type, raw_material, layer, exp_type, D2, p_emp, p_chisq) %>%
    arrange(space, ID, exp_type),
  file.path(OUT_DIR, "sdg_projection_typicality.csv"))
write_csv(nearest, file.path(OUT_DIR, "sdg_projection_nearest_exp_type.csv"))
write_csv(cmp,     file.path(OUT_DIR, "sdg_projection_loo_comparison.csv"))

cat("\nWrote 3 figures and 3 csv files to ", OUT_DIR, "\n", sep = "")

cat("\n"); print(sessionInfo())
