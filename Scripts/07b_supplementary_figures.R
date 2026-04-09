# =============================================================================
# 07b_supplementary_figures.R — Supplementary figures
# =============================================================================
#
# SFig 1: Temperature metric tournament (evolutionary vs acute)
# SFig 2: Acinetobacter species panel using collection-day temperature
# =============================================================================

library(tidyverse)
library(patchwork)
library(scales)

base_dir <- "/Users/bxc331/Desktop/AMR_open_project"
res_dir  <- file.path(base_dir, "Results")
out_dir  <- file.path(base_dir, "Datasets")
fig_dir  <- file.path(base_dir, "Figures")

theme_lph <- theme_bw(base_size = 8) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.25, colour = "grey92"),
    axis.text = element_text(size = 6.5),
    axis.title = element_text(size = 7.5),
    plot.tag = element_text(face = "bold", size = 10),
    plot.title = element_text(face = "bold", size = 9)
  )
theme_set(theme_lph)

analysis <- readRDS(file.path(out_dir, "analysis_ready.rds"))

# =============================================================================
# SFIG 0: FULL DATASET MAP (all 1,233 samples)
# =============================================================================
cat("--- SFig 0: Full dataset map ---\n")

samples_all <- readRDS(file.path(out_dir, "clean_samples.rds")) %>%
  left_join(analysis %>% select(genepid, T_30d), by = "genepid")

# Flag which samples are in the analysis subset (n=1,077)
# Must also intersect with resistome matrices to match PERMANOVA sample set
conf_vars <- c("T_30d", "gdp_pcap_ppp", "sanitation",
               "health_exp_gdp", "oop_health_exp", "immunization_dpt",
               "animal_amc_mgkg", "pop_density", "year", "Region")
resistome <- readRDS(file.path(out_dir, "resistome_matrices.rds"))
resistome_ids <- as.character(resistome$fg_cluster_clr$genepid)
analysis_cc <- analysis %>%
  mutate(genepid = as.character(genepid)) %>%
  filter(complete.cases(across(all_of(conf_vars))),
         genepid %in% resistome_ids)
samples_all <- samples_all %>%
  mutate(genepid = as.character(genepid),
         in_analysis = genepid %in% analysis_cc$genepid)

world <- map_data("world")

n_total <- nrow(samples_all)
n_analysis <- sum(samples_all$in_analysis)
n_excluded <- n_total - n_analysis

cat(sprintf("  Total: %d, Analysis: %d, Excluded: %d\n",
            n_total, n_analysis, n_excluded))

sfig0 <- ggplot() +
  geom_polygon(data = world, aes(x = long, y = lat, group = group),
               fill = "grey95", colour = "grey80", linewidth = 0.12) +
  geom_point(data = samples_all %>% filter(!in_analysis),
             aes(x = lon, y = lat), size = 1.0, alpha = 0.6,
             shape = 4, colour = "grey50") +
  geom_point(data = samples_all %>% filter(in_analysis, !is.na(T_30d)),
             aes(x = lon, y = lat, fill = T_30d),
             size = 1.2, alpha = 0.8, shape = 21, colour = "black", stroke = 0.15) +
  scale_fill_gradient2(low = "#4575B4", mid = "#FFFFBF", high = "#D73027",
                       midpoint = 15, name = "30-day mean\ntemperature (\u00B0C)") +
  coord_fixed(1.3, xlim = c(-170, 180), ylim = c(-55, 75)) +
  labs(title = sprintf("All %d sewage metagenomes from 112 countries", n_total),
       subtitle = sprintf("Coloured circles: %d samples in analysis subset; grey crosses: %d excluded (missing covariates)",
                           n_analysis, n_excluded)) +
  theme_void(base_size = 8) +
  theme(legend.position = "right",
        legend.key.height = unit(0.6, "cm"), legend.key.width = unit(0.3, "cm"),
        legend.text = element_text(size = 6), legend.title = element_text(size = 7),
        plot.title = element_text(face = "bold", size = 9),
        plot.subtitle = element_text(size = 7, colour = "grey30"))

ggsave(file.path(fig_dir, "SFig0_full_dataset_map.pdf"), sfig0,
       width = 8, height = 4.5)
cat("  Saved SFig0\n")

# =============================================================================
# SFIG 1: TEMPERATURE METRIC TOURNAMENT
# =============================================================================
cat("--- SFig 1: Metric tournament ---\n")

mt <- read_csv(file.path(res_dir, "layer1_metric_tournament.csv"),
               show_col_types = FALSE)

# Classify metrics
mt <- mt %>%
  mutate(
    timescale = case_when(
      str_detect(metric, "^BIO") ~ "Evolutionary\n(WorldClim)",
      str_detect(metric, "T_sampling|T_30d|T_days30") ~ "Acute\n(sampling period)",
      TRUE ~ "Intermediate\n(NASA POWER)"
    ),
    label = case_when(
      metric == "T_30d" ~ "Annual mean (BIO1)",
      metric == "BIO10_mean_temp_warmest_quarter" ~ "Warmest quarter mean (BIO10)",
      metric == "T_sampling" ~ "Collection-day temp",
      metric == "T_30d" ~ "30-day mean before sampling",
      metric == "T_90d" ~ "90-day mean before sampling",
      metric == "T_365d" ~ "365-day mean before sampling",
      metric == "T_annual_mean" ~ "Annual mean (NASA POWER)",
      metric == "T_max" ~ "Annual max temp",
      metric == "T_days30" ~ "Days >30\u00b0C",
      metric == "T_95pct" ~ "95th percentile temp",
      metric == "T_amplitude" ~ "Annual amplitude",
      metric == "T_variability" ~ "Temperature variability",
      TRUE ~ metric
    )
  )

sfig1 <- ggplot(mt, aes(x = R2, y = reorder(label, R2), fill = timescale)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = sprintf("%.4f", R2)), hjust = -0.1, size = 2.2) +
  scale_fill_manual(values = c("Evolutionary\n(WorldClim)" = "#E66101",
                                "Intermediate\n(NASA POWER)" = "#FDB863",
                                "Acute\n(sampling period)" = "#B2ABD2"),
                    name = "Timescale") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.25))) +
  labs(tag = "a",
       title = "Temperature metric tournament (functional resistome)",
       x = "PERMANOVA R\u00b2",
       y = NULL)

ggsave(file.path(fig_dir, "SFig1_metric_tournament.pdf"), sfig1,
       width = 7, height = 4)
cat("  Saved SFig1\n")

# =============================================================================
# SFIG 2: ACINETOBACTER PANEL WITH COLLECTION-DAY TEMPERATURE
# =============================================================================
cat("--- SFig 2: Acinetobacter with collection-day temp ---\n")

motu <- readRDS(file.path(out_dir, "motu_species_cache.rds"))
motu_ids <- as.character(as.integer(motu$genepid))
motu_mat <- as.matrix(motu[, -1])
total_reads <- rowSums(motu_mat)
motu_ra <- motu_mat / total_reads
all_cols <- names(motu)[-1]

merge_df <- tibble(genepid = motu_ids) %>%
  inner_join(analysis %>% mutate(genepid = as.character(genepid)) %>%
               select(genepid, T_30d, T_sampling), by = "genepid") %>%
  filter(!is.na(T_30d))

idx <- match(merge_df$genepid, motu_ids)
motu_ra_matched <- motu_ra[idx, ]

get_ra <- function(pattern) {
  cols <- all_cols[str_detect(all_cols, fixed(pattern))]
  if (length(cols) == 0) return(rep(0, nrow(motu_ra_matched)))
  rowSums(motu_ra_matched[, cols, drop = FALSE])
}

species <- c("A. baumannii" = "Acinetobacter baumannii",
             "A. johnsonii" = "Acinetobacter johnsonii",
             "A. lwoffii" = "Acinetobacter lwoffii")

acin_genus <- rowSums(motu_ra_matched[, all_cols[str_detect(all_cols, "^Acinetobacter")], drop = FALSE])

# Build data for both temp metrics
build_acin_data <- function(temp_col, temp_label) {
  temp_vals <- merge_df[[temp_col]]
  ok <- !is.na(temp_vals)

  plot_df <- tibble(taxon = "Acinetobacter (all spp.)",
                    temp = temp_vals[ok], ra = acin_genus[ok])
  for (i in seq_along(species)) {
    sp_ra <- get_ra(species[i])
    plot_df <- bind_rows(plot_df, tibble(
      taxon = names(species)[i], temp = temp_vals[ok], ra = sp_ra[ok]
    ))
  }

  # Stats
  stats <- plot_df %>%
    group_by(taxon) %>%
    summarise(
      rho = cor(ra, temp, method = "spearman"),
      p = suppressWarnings(cor.test(ra, temp, method = "spearman")$p.value),
      .groups = "drop"
    ) %>%
    mutate(
      p_label = ifelse(p < 0.001, "p<0.001", sprintf("p=%.3f", p)),
      legend_label = sprintf("%s (rho=%+.2f, %s)", taxon, rho, p_label)
    )

  plot_df <- plot_df %>%
    left_join(stats %>% select(taxon, legend_label), by = "taxon") %>%
    filter(ra > 0)

  list(data = plot_df, stats = stats, temp_label = temp_label)
}

bio1_data <- build_acin_data("T_30d", "30-day mean temperature (\u00b0C)")
samp_data <- build_acin_data("T_sampling", "Collection-day temperature (\u00b0C)")

sp_cols_3 <- c("#0072B2", "#009E73", "#D55E00")
col_map_bio1 <- setNames(c("grey55", sp_cols_3), unique(bio1_data$data$legend_label))
col_map_samp <- setNames(c("grey55", sp_cols_3), unique(samp_data$data$legend_label))

make_acin_plot <- function(d, col_map, tag, title) {
  ggplot(d$data, aes(x = temp, y = log10(ra), colour = legend_label)) +
    geom_point(alpha = 0.08, size = 0.3) +
    geom_smooth(method = "lm", se = TRUE, linewidth = 0.7, alpha = 0.15) +
    scale_colour_manual(values = col_map, name = NULL) +
    guides(colour = guide_legend(override.aes = list(alpha = 1, size = 0, linewidth = 1))) +
    labs(tag = tag, title = title, x = d$temp_label, y = "Relative abundance (log10)") +
    theme(legend.key = element_rect(fill = NA, colour = NA),
          legend.position.inside = c(0.98, 0.02),
          legend.justification = c(1, 0),
          legend.background = element_rect(fill = alpha("white", 0.85), colour = NA),
          legend.key.size = unit(0.25, "cm"),
          legend.text = element_text(size = 5.5),
          plot.title = element_text(face = "bold.italic", size = 8))
}

sfig2a <- make_acin_plot(bio1_data, col_map_bio1, "a", "Acinetobacter — 30-day mean temperature")
sfig2b <- make_acin_plot(samp_data, col_map_samp, "b", "Acinetobacter — Collection-day temperature")

sfig2 <- sfig2a | sfig2b

ggsave(file.path(fig_dir, "SFig2_acinetobacter_temp_comparison.pdf"), sfig2,
       width = 10, height = 4.5)
cat("  Saved SFig2\n")

cat("\nDone.\n")

# =============================================================================
# SFIG 3: ADJUSTED SPECIES CORRELATIONS (raw vs adjusted comparison)
# =============================================================================
cat("--- SFig 3: Adjusted species correlations ---\n")

adj_sp <- read_csv(file.path(res_dir, "adjusted_species_temp.csv"),
                   show_col_types = FALSE) %>%
  filter(!is.na(rho_adj))

adj_sp_long <- adj_sp %>%
  select(genus, display, rho_raw, rho_adj) %>%
  pivot_longer(c(rho_raw, rho_adj), names_to = "model", values_to = "rho") %>%
  mutate(model = recode(model,
                        rho_raw = "Unadjusted",
                        rho_adj = "Adjusted for\nconfounders"))

sfig3 <- ggplot(adj_sp_long, aes(x = rho, y = reorder(display, rho),
                                  shape = model)) +
  geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_line(aes(group = display), colour = "grey70", linewidth = 0.3) +
  geom_point(aes(fill = model), size = 2.5, colour = "black", stroke = 0.3) +
  scale_shape_manual(values = c("Unadjusted" = 21, "Adjusted for\nconfounders" = 24)) +
  scale_fill_manual(values = c("Unadjusted" = "grey70",
                                "Adjusted for\nconfounders" = "#2C3E50")) +
  labs(title = "Species temperature associations: unadjusted vs covariate-adjusted",
       x = "Spearman rho with 30-day mean temperature",
       y = NULL, shape = NULL, fill = NULL) +
  theme(legend.position = "bottom")

ggsave(file.path(fig_dir, "SFig3_adjusted_species.pdf"), sfig3,
       width = 7, height = 5)
cat("  Saved SFig3\n")

# =============================================================================
# SFIG 4: ADJUSTED ARG FAMILY CORRELATIONS
# =============================================================================
cat("--- SFig 4: Adjusted ARG correlations ---\n")

adj_arg <- read_csv(file.path(res_dir, "adjusted_arg_temp.csv"),
                    show_col_types = FALSE) %>%
  mutate(label = case_when(
    arg_family == "dfr" ~ "dfr (trimethoprim)",
    arg_family == "qnr" ~ "qnr (fluoroquinolone)",
    arg_family == "tet" ~ "tet (tetracycline)",
    arg_family == "blaTEM" ~ "blaTEM (beta-lactam)",
    arg_family == "aac" ~ "aac (aminoglycoside)",
    arg_family == "blaSHV" ~ "blaSHV (beta-lactam)",
    arg_family == "blaCTX-M" ~ "blaCTX-M (beta-lactam)",
    arg_family == "aph" ~ "aph (aminoglycoside)",
    arg_family == "ant" ~ "ant (aminoglycoside)",
    arg_family == "blaOXA" ~ "blaOXA (carbapenem)",
    arg_family == "erm" ~ "erm (macrolide)",
    TRUE ~ arg_family
  ))

adj_arg_long <- adj_arg %>%
  select(label, rho_raw, rho_adj) %>%
  pivot_longer(c(rho_raw, rho_adj), names_to = "model", values_to = "rho") %>%
  mutate(model = recode(model,
                        rho_raw = "Unadjusted",
                        rho_adj = "Adjusted for\nconfounders"))

sfig4 <- ggplot(adj_arg_long, aes(x = rho, y = reorder(label, rho),
                                    shape = model)) +
  geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_line(aes(group = label), colour = "grey70", linewidth = 0.3) +
  geom_point(aes(fill = model), size = 2.5, colour = "black", stroke = 0.3) +
  scale_shape_manual(values = c("Unadjusted" = 21, "Adjusted for\nconfounders" = 24)) +
  scale_fill_manual(values = c("Unadjusted" = "grey70",
                                "Adjusted for\nconfounders" = "#2C3E50")) +
  labs(title = "Resistance gene family temperature associations: unadjusted vs adjusted",
       x = "Spearman rho with 30-day mean temperature",
       y = "Resistance gene family", shape = NULL, fill = NULL) +
  theme(legend.position = "bottom")

ggsave(file.path(fig_dir, "SFig4_adjusted_args.pdf"), sfig4,
       width = 7, height = 4.5)
cat("  Saved SFig4\n")

cat("\nAll supplementary figures done.\n")

# =============================================================================
# SFIG 5: FIG 2 REPLICATED WITH COLLECTION-DAY TEMPERATURE
# =============================================================================
cat("--- SFig 5: Species panels with collection-day temp ---\n")

motu <- readRDS(file.path(out_dir, "motu_species_cache.rds"))
motu_ids <- as.character(as.integer(motu$genepid))
motu_mat <- as.matrix(motu[, -1])
total_reads <- rowSums(motu_mat)
motu_ra <- motu_mat / total_reads
all_cols <- names(motu)[-1]

# Require BOTH temperature metrics so n is comparable to Fig 2
merge_df_samp <- tibble(genepid = motu_ids) %>%
  inner_join(analysis %>% mutate(genepid = as.character(genepid)) %>%
               select(genepid, T_sampling, T_30d), by = "genepid") %>%
  filter(!is.na(T_sampling), !is.na(T_30d))
cat("  Samples with both temp metrics:", nrow(merge_df_samp), "\n")

idx_s <- match(merge_df_samp$genepid, motu_ids)
motu_ra_samp <- motu_ra[idx_s, ]

get_ra_s <- function(pattern, use_fixed = FALSE) {
  if (use_fixed) {
    cols <- all_cols[str_detect(all_cols, fixed(pattern))]
  } else {
    cols <- all_cols[str_detect(all_cols, pattern)]
  }
  if (length(cols) == 0) return(rep(0, nrow(motu_ra_samp)))
  rowSums(motu_ra_samp[, cols, drop = FALSE])
}

sp_cols <- c("#0072B2", "#009E73", "#D55E00", "#CC79A7")

make_genus_panel_samp <- function(genus_display, genus_pattern, species_list, tag) {
  genus_ra <- get_ra_s(genus_pattern)
  plot_df <- tibble(taxon = paste0(genus_display, " (all spp.)"),
                    temp = merge_df_samp$T_sampling, ra = genus_ra)
  for (i in seq_along(species_list)) {
    sp_ra <- get_ra_s(species_list[i], use_fixed = TRUE)
    plot_df <- bind_rows(plot_df, tibble(
      taxon = names(species_list)[i],
      temp = merge_df_samp$T_sampling, ra = sp_ra))
  }
  if (length(species_list) == 1) {
    genus_sum <- sum(genus_ra)
    sp_sum <- sum(get_ra_s(species_list[1], use_fixed = TRUE))
    if (abs(genus_sum - sp_sum) / max(genus_sum, 1e-10) < 0.01) {
      plot_df <- plot_df %>% filter(taxon != paste0(genus_display, " (all spp.)"))
    }
  }
  stats <- plot_df %>%
    group_by(taxon) %>%
    summarise(rho = cor(ra, temp, method = "spearman"),
              p = suppressWarnings(cor.test(ra, temp, method = "spearman")$p.value),
              .groups = "drop") %>%
    mutate(p_label = ifelse(p < 0.001, "p<0.001", sprintf("p=%.3f", p)),
           legend_label = sprintf("%s\nrho=%+.2f, %s", taxon, rho, p_label))
  plot_df <- plot_df %>%
    left_join(stats %>% select(taxon, legend_label), by = "taxon")
  plot_df_nz <- plot_df %>% filter(ra > 0)

  taxa <- unique(plot_df_nz$legend_label)
  is_genus <- str_detect(taxa, "all spp")
  colour_map <- setNames(rep("grey55", length(taxa)), taxa)
  sp_idx <- which(!is_genus)
  for (j in seq_along(sp_idx)) colour_map[taxa[sp_idx[j]]] <- sp_cols[min(j, length(sp_cols))]

  ggplot(plot_df_nz, aes(x = temp, y = log10(ra), colour = legend_label)) +
    geom_point(alpha = 0.08, size = 0.3) +
    geom_smooth(method = "lm", se = TRUE, linewidth = 0.7, alpha = 0.15) +
    scale_colour_manual(values = colour_map, name = NULL) +
    guides(colour = guide_legend(override.aes = list(alpha = 1, size = 0, linewidth = 1))) +
    labs(tag = tag, title = genus_display,
         x = "Collection-day temperature (\u00b0C)",
         y = "Relative abundance (log10)") +
    theme(legend.key = element_rect(fill = NA, colour = NA),
          legend.position = "inside",
          legend.position.inside = c(0.98, 0.02),
          legend.justification = c(1, 0),
          legend.background = element_rect(fill = alpha("white", 0.85), colour = NA),
          legend.key.size = unit(0.25, "cm"),
          legend.text = element_text(size = 5),
          plot.title = element_text(face = "bold.italic", size = 8))
}

# --- WHO BPPL 2024 Critical ---
s5a <- make_genus_panel_samp("Acinetobacter", "^Acinetobacter",
  c("A. baumannii" = "Acinetobacter baumannii",
    "A. johnsonii" = "Acinetobacter johnsonii",
    "A. lwoffii" = "Acinetobacter lwoffii"), "A")
s5b <- make_genus_panel_samp("Klebsiella", "^Klebsiella",
  c("K. pneumoniae" = "Klebsiella pneumoniae"), "B")
s5c <- make_genus_panel_samp("Escherichia", "^Escherichia",
  c("E. coli" = "Escherichia coli"), "C")
s5d <- make_genus_panel_samp("Enterobacter", "^Enterobacter",
  c("E. cloacae" = "Enterobacter cloacae"), "D")
s5e <- make_genus_panel_samp("Mycobacterium", "^Mycobacterium",
  c("M. tuberculosis" = "Mycobacterium tuberculosis"), "E")

# --- WHO BPPL 2024 High ---
s5f <- make_genus_panel_samp("Enterococcus", "^Enterococcus",
  c("E. faecium" = "Enterococcus faecium",
    "E. faecalis" = "Enterococcus faecalis"), "F")
s5g <- make_genus_panel_samp("Staphylococcus", "^Staphylococcus",
  c("S. aureus" = "Staphylococcus aureus"), "G")
s5h <- make_genus_panel_samp("Pseudomonas", "^Pseudomonas",
  c("P. aeruginosa" = "Pseudomonas aeruginosa"), "H")
s5i <- make_genus_panel_samp("Salmonella", "^Salmonella",
  c("S. enterica" = "Salmonella enterica"), "I")
s5j <- make_genus_panel_samp("Neisseria", "^Neisseria",
  c("N. gonorrhoeae" = "Neisseria gonorrhoeae"), "J")

# --- WHO BPPL 2024 Medium ---
s5k <- make_genus_panel_samp("Streptococcus", "^Streptococcus",
  c("S. pneumoniae" = "Streptococcus pneumoniae",
    "S. pyogenes" = "Streptococcus pyogenes"), "K")
s5l <- make_genus_panel_samp("Haemophilus", "^Haemophilus",
  c("H. influenzae" = "Haemophilus influenzae"), "L")

sfig5 <- (s5a | s5b | s5c | s5d) / (s5e | s5f | s5g | s5h) / (s5i | s5j | s5k | s5l)
ggsave(file.path(fig_dir, "SFig5_species_collection_day_temp.pdf"), sfig5,
       width = 14, height = 10.5)
cat("  Saved SFig5\n")

# =============================================================================
# SFIG 6: FIG 3A REPLICATED WITH COLLECTION-DAY TEMPERATURE
# =============================================================================
cat("--- SFig 6: ARG families with collection-day temp ---\n")

# Rebuild ARG family abundances and correlate with T_sampling
library(readxl)
counts <- readRDS(file.path(out_dir, "clean_ARG_counts.rds"))
sd2 <- read_excel(file.path(out_dir, "Martiny_et_al/supplementary/41467_2025_66070_MOESM5_ESM.xlsx"),
                  sheet = "SD2 ARG metadata")

gene_map <- sd2 %>%
  filter(str_detect(fa_name, "^resfinder")) %>%
  mutate(real_name = str_extract(fa_name, "resfinder\\|(.+)", group = 1) %>% str_to_lower()) %>%
  distinct(cluster_representative_98, real_name)

cluster_family <- gene_map %>% mutate(arg_family = NA_character_)
for (fam in c("dfr","qnr","tet","blaTEM","aac","blaSHV","blaCTX-M","aph","blaOXA","erm")) {
  pat <- switch(fam, dfr="^dfr", qnr="^qnr", tet="^tet", blaTEM="^blatem",
    aac="^aac", blaSHV="^blashv", `blaCTX-M`="^blactx-m", aph="^aph",
    blaOXA="^blaoxa", erm="^erm")
  cluster_family <- cluster_family %>%
    mutate(arg_family = ifelse(is.na(arg_family) & str_detect(real_name, pat), fam, arg_family))
}
cluster_family <- cluster_family %>%
  mutate(arg_family = ifelse(is.na(arg_family) & str_detect(real_name, "^ant\\("), "ant", arg_family)) %>%
  filter(!is.na(arg_family)) %>% distinct(cluster_representative_98, arg_family)

total_abund <- counts %>%
  group_by(genepid) %>%
  summarise(total_arg = sum(fragmentCountAln_adj, na.rm = TRUE), .groups = "drop") %>%
  mutate(genepid = as.character(genepid))

family_abund <- counts %>%
  inner_join(cluster_family, by = "cluster_representative_98", relationship = "many-to-many") %>%
  group_by(genepid, arg_family) %>%
  summarise(abundance = sum(fragmentCountAln_adj, na.rm = TRUE), .groups = "drop") %>%
  mutate(genepid = as.character(genepid))

# Require both temp metrics for fair comparison with main figures
analysis_samp <- analysis %>%
  mutate(genepid = as.character(genepid)) %>%
  filter(!is.na(T_sampling), !is.na(T_30d))
cat("  Samples with both temp metrics:", nrow(analysis_samp), "\n")

arg_cors_samp <- map_dfr(unique(cluster_family$arg_family), function(fam) {
  fam_data <- family_abund %>% filter(arg_family == fam) %>%
    right_join(total_abund, by = "genepid") %>%
    mutate(abundance = replace_na(abundance, 0), ra = abundance / total_arg) %>%
    inner_join(analysis_samp %>% select(genepid, T_sampling), by = "genepid")
  if (nrow(fam_data) < 50) return(NULL)
  ct <- cor.test(fam_data$ra, fam_data$T_sampling, method = "spearman")
  tibble(arg_family = fam, rho = ct$estimate, p = ct$p.value)
}) %>%
  mutate(
    sig = p.adjust(p, method = "BH") < 0.05,
    label = case_when(
      arg_family == "dfr" ~ "dfr (trimethoprim)",
      arg_family == "qnr" ~ "qnr (fluoroquinolone)",
      arg_family == "tet" ~ "tet (tetracycline)",
      arg_family == "blaTEM" ~ "blaTEM (beta-lactam)",
      arg_family == "aac" ~ "aac (aminoglycoside)",
      arg_family == "blaSHV" ~ "blaSHV (beta-lactam)",
      arg_family == "blaCTX-M" ~ "blaCTX-M (beta-lactam)",
      arg_family == "aph" ~ "aph (aminoglycoside)",
      arg_family == "ant" ~ "ant (aminoglycoside)",
      arg_family == "blaOXA" ~ "blaOXA (carbapenem)",
      arg_family == "erm" ~ "erm (macrolide)", TRUE ~ arg_family),
    p_label = ifelse(p < 0.001, "p<0.001",
                     ifelse(p < 0.01, sprintf("p=%.3f", p), sprintf("p=%.2f", p))),
    point_col = ifelse(sig, "#2C3E50", "grey65")
  )

sfig6 <- ggplot(arg_cors_samp, aes(x = rho, y = reorder(label, rho))) +
  geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_segment(aes(xend = 0, yend = reorder(label, rho)),
               colour = "grey75", linewidth = 0.4) +
  geom_point(aes(fill = point_col), size = 2.8,
             shape = 21, colour = "black", stroke = 0.3) +
  geom_text(aes(label = p_label, x = rho + ifelse(rho > 0, 0.03, -0.03)),
            hjust = ifelse(arg_cors_samp$rho > 0, 0, 1),
            size = 2, colour = "grey30") +
  scale_fill_identity() +
  labs(title = "Resistance gene families vs collection-day temperature (unadjusted)",
       x = "Spearman rho with collection-day temperature",
       y = "Resistance gene family")

ggsave(file.path(fig_dir, "SFig6_args_collection_day_temp.pdf"), sfig6,
       width = 7, height = 4.5)
cat("  Saved SFig6\n")

# =============================================================================
# SFIG 7: DIRECTED ACYCLIC GRAPH (DAG)
# =============================================================================
cat("--- SFig 7: DAG ---\n")

# Build DAG using ggplot — no extra packages needed
dag_nodes <- tibble(
  label = c("30-day mean\ntemperature", "Bacterial\ncommunity\ncomposition",
            "Resistance gene\ncomposition\n(resistome)", "Socioeconomic\nfactors\n(GDP, WASH, health)",
            "Antimicrobial\nconsumption\n(human + animal)", "WHO\nRegion"),
  x = c(0, 2, 4, 2, 4, 0),
  y = c(2, 2, 2, 4, 4, 4),
  fill = c("#D73027", "#FDB863", "#4575B4", "#999999", "#999999", "#999999")
)

dag_edges <- tibble(
  x = c(0, 2, 0, 2, 4, 0, 2, 0),
  y = c(2, 2, 2, 4, 4, 4, 4, 4),
  xend = c(2, 4, 4, 0, 2, 2, 4, 4),
  yend = c(2, 2, 2, 2, 2, 2, 2, 2),
  type = c("main", "main", "direct", "confounder", "confounder",
           "confounder", "confounder", "confounder"),
  ltype = c("solid", "solid", "dashed", "solid", "solid",
            "solid", "solid", "solid")
)

# Simplify: key paths only
dag_edges2 <- tibble(
  x =    c(0,   2,   0,   2,   4,   0),
  y =    c(2,   2,   2,   4,   4,   4),
  xend = c(2,   4,   4,   0,   2,   2),
  yend = c(2,   2,   2,   2,   2,   2),
  label = c("a: temp changes\ncommunity", "b: community\ncarries ARGs",
            "c': direct\n(residual)", "confounders", "confounders", "confounders"),
  colour = c("#D73027", "#E66101", "#D73027", "grey50", "grey50", "grey50"),
  ltype = c("solid", "solid", "dashed", "solid", "solid", "solid"),
  lwidth = c(1.2, 1.2, 0.6, 0.5, 0.5, 0.5)
)

sfig7 <- ggplot() +
  # Edges
  geom_segment(data = dag_edges2,
               aes(x = x, y = y, xend = xend, yend = yend,
                   colour = colour, linetype = ltype, linewidth = lwidth),
               arrow = arrow(length = unit(0.15, "cm"), type = "closed")) +
  scale_colour_identity() +
  scale_linetype_identity() +
  scale_linewidth_identity() +
  # Nodes
  geom_label(data = dag_nodes, aes(x = x, y = y, label = label, fill = fill),
             colour = "white", fontface = "bold", size = 2.5,
             label.padding = unit(0.4, "lines"), label.size = 0) +
  scale_fill_identity() +
  # Annotations for mediation
  annotate("text", x = 1, y = 1.7, label = "Path a", size = 2.3,
           fontface = "italic", colour = "#D73027") +
  annotate("text", x = 3, y = 1.7, label = "Path b", size = 2.3,
           fontface = "italic", colour = "#E66101") +
  annotate("text", x = 2, y = 1.3, label = "Path c' (direct)", size = 2.3,
           fontface = "italic", colour = "#D73027") +
  annotate("text", x = 2, y = 0.7,
           label = "Indirect (a x b): ~65% of total effect\nDirect (c'): ~35% of total effect",
           size = 2.5, fontface = "bold") +
  coord_cartesian(xlim = c(-0.8, 4.8), ylim = c(0.3, 4.7)) +
  theme_void() +
  labs(title = "Assumed causal structure (directed acyclic graph)")

ggsave(file.path(fig_dir, "SFig7_DAG.pdf"), sfig7, width = 7, height = 4.5)
cat("  Saved SFig7\n")

# =============================================================================
# SFIG 8: PROJECTED CHANGE IN KEY SPECIES UNDER WARMING
# =============================================================================
cat("--- SFig 8: Warming projections ---\n")

# Use the adjusted linear models to project change under +1, +1.5, +2 C warming
# For each species, we have the adjusted rho — but we need the actual slope (beta)
# from a linear model to project. Use the raw data + region-adjusted model.

analysis <- readRDS(file.path(out_dir, "analysis_ready.rds")) %>%
  mutate(genepid = as.character(genepid))

motu <- readRDS(file.path(out_dir, "motu_species_cache.rds"))
motu_ids <- as.character(as.integer(motu$genepid))
motu_mat <- as.matrix(motu[, -1])
total_reads <- rowSums(motu_mat)
motu_ra <- motu_mat / total_reads
all_cols <- names(motu)[-1]

merge_proj <- tibble(genepid = motu_ids) %>%
  inner_join(analysis %>% select(genepid, T_30d, Region),
             by = "genepid") %>%
  filter(!is.na(T_30d))

idx_p <- match(merge_proj$genepid, motu_ids)

species_proj <- list(
  "A. baumannii" = "Acinetobacter baumannii",
  "K. pneumoniae" = "Klebsiella pneumoniae",
  "E. coli" = "Escherichia coli",
  "Enterococcus (all spp.)" = "^Enterococcus"
)

proj_results <- map_dfr(names(species_proj), function(sp) {
  pat <- species_proj[[sp]]
  if (str_starts(pat, "\\^")) {
    cols <- all_cols[str_detect(all_cols, pat)]
  } else {
    cols <- all_cols[str_detect(all_cols, fixed(pat))]
  }
  ra <- rowSums(motu_ra[idx_p, cols, drop = FALSE])

  df <- tibble(ra = ra, temp = merge_proj$T_30d,
               region = merge_proj$Region) %>%
    filter(ra > 0)

  mod <- lm(log10(ra) ~ temp + region, data = df)
  beta <- coef(mod)["temp"]

  # Predict % change for warming scenarios
  # At mean RA, a beta of X per degree means:
  # change = 10^(beta * delta_T) - 1 (multiplicative change on original scale)
  tibble(
    species = sp,
    beta_per_degree = beta,
    warming = c(1, 1.5, 2),
    pct_change = (10^(beta * c(1, 1.5, 2)) - 1) * 100
  )
})

sfig8 <- ggplot(proj_results, aes(x = factor(warming), y = pct_change, fill = species)) +
  geom_col(position = "dodge", width = 0.7) +
  geom_hline(yintercept = 0, colour = "grey50") +
  scale_fill_manual(values = c("A. baumannii" = "#D55E00",
                                "K. pneumoniae" = "#0072B2",
                                "E. coli" = "#009E73",
                                "Enterococcus (all spp.)" = "#CC79A7"),
                    name = NULL) +
  labs(title = "Projected change in species relative abundance under warming scenarios",
       subtitle = "Based on region-adjusted linear models (cross-sectional, illustrative only)",
       x = "Temperature increase (\u00b0C)",
       y = "Projected change in relative abundance (%)") +
  theme(plot.subtitle = element_text(size = 7, face = "italic"))

ggsave(file.path(fig_dir, "SFig8_warming_projections.pdf"), sfig8, width = 7, height = 4)
cat("  Saved SFig8\n")

# =============================================================================
# SFIG 9: REGIONAL BREAKDOWN — SPECIES-TEMPERATURE BY WHO REGION
# =============================================================================
cat("--- SFig 9: Regional breakdown ---\n")

# For key species, compute rho per WHO Region
species_regional <- list(
  "A. baumannii" = "Acinetobacter baumannii",
  "K. pneumoniae" = "Klebsiella pneumoniae",
  "E. coli" = "Escherichia coli",
  "Enterococcus (all spp.)" = "^Enterococcus"
)

regional_cors <- map_dfr(names(species_regional), function(sp) {
  pat <- species_regional[[sp]]
  if (str_starts(pat, "\\^")) {
    cols <- all_cols[str_detect(all_cols, pat)]
  } else {
    cols <- all_cols[str_detect(all_cols, fixed(pat))]
  }
  ra <- rowSums(motu_ra[idx_p, cols, drop = FALSE])

  df <- tibble(ra = ra, temp = merge_proj$T_30d,
               region = merge_proj$Region)

  # Per region
  df %>%
    group_by(region) %>%
    filter(n() >= 20) %>%
    summarise(
      rho = cor(ra, temp, method = "spearman"),
      p = suppressWarnings(cor.test(ra, temp, method = "spearman")$p.value),
      n = n(),
      .groups = "drop"
    ) %>%
    mutate(species = sp)
})

# Add global estimate
global_cors <- map_dfr(names(species_regional), function(sp) {
  pat <- species_regional[[sp]]
  if (str_starts(pat, "\\^")) {
    cols <- all_cols[str_detect(all_cols, pat)]
  } else {
    cols <- all_cols[str_detect(all_cols, fixed(pat))]
  }
  ra <- rowSums(motu_ra[idx_p, cols, drop = FALSE])
  ct <- cor.test(ra, merge_proj$T_30d, method = "spearman")
  tibble(region = "Global", rho = ct$estimate, p = ct$p.value,
         n = length(ra), species = sp)
})

all_regional <- bind_rows(regional_cors, global_cors) %>%
  mutate(sig = p < 0.05,
         region = factor(region, levels = c(sort(unique(regional_cors$region)), "Global")))

sfig9 <- ggplot(all_regional, aes(x = rho, y = region, colour = species)) +
  geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_point(aes(shape = sig, size = ifelse(region == "Global", 3.5, 2.5)),
             position = position_dodge(width = 0.6)) +
  scale_colour_manual(values = c("A. baumannii" = "#D55E00",
                                  "K. pneumoniae" = "#0072B2",
                                  "E. coli" = "#009E73",
                                  "Enterococcus (all spp.)" = "#CC79A7"),
                      name = NULL) +
  scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 1),
                     labels = c("TRUE" = "p<0.05", "FALSE" = "NS"),
                     name = NULL) +
  scale_size_identity() +
  labs(title = "Species-temperature associations by WHO Region",
       x = "Spearman rho with 30-day mean temperature",
       y = NULL) +
  theme(legend.position = "bottom")

ggsave(file.path(fig_dir, "SFig9_regional_breakdown.pdf"), sfig9, width = 8, height = 5)
cat("  Saved SFig9\n")
