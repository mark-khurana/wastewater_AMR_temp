# =============================================================================
# 09c_city_species_figures.R — Per-city species panels (mirrors Fig 2)
# =============================================================================
#
# SFig15–19: One figure per city, each with 8 WHO priority genera showing
# genus total (grey) + key species (coloured) vs temperature.
# Spearman rho + p annotated for each taxon.
# =============================================================================

library(tidyverse)
library(patchwork)
library(scales)
library(lubridate)

base_dir <- "/Users/bxc331/Desktop/AMR_open_project"
data_dir <- file.path(base_dir, "Datasets/Becsei_et_al")
fig_dir  <- file.path(base_dir, "Figures")

cat("=== 09c_city_species_figures.R ===\n\n")

theme_lph <- theme_bw(base_size = 8) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.25, colour = "grey92"),
    axis.text = element_text(size = 6.5),
    axis.title = element_text(size = 7.5),
    plot.tag = element_text(face = "bold", size = 10),
    plot.title = element_text(face = "bold.italic", size = 8)
  )
theme_set(theme_lph)

sp_cols <- c("#0072B2", "#009E73", "#D55E00", "#CC79A7")

# =============================================================================
# LOAD DATA
# =============================================================================

meta <- read.csv(file.path(data_dir, "meta_location.csv"), stringsAsFactors = FALSE)
gen_abund <- read.csv(file.path(data_dir, "genomic_genus_abundance.csv"),
                       row.names = 1, check.names = FALSE)
sp_abund  <- read.csv(file.path(data_dir, "genomic_species_abundance.csv"),
                       row.names = 1, check.names = FALSE)

df <- meta %>%
  mutate(collection_date = as.Date(collection_date)) %>%
  filter(!is.na(collection_date))

agg_key <- df %>%
  group_by(city, collection_date) %>%
  summarise(sample_ids = list(complete_name), .groups = "drop")

agg_matrix <- function(mat, keys) {
  out <- matrix(NA, nrow(keys), ncol(mat)); colnames(out) <- colnames(mat)
  for (i in seq_len(nrow(keys)))
    out[i, ] <- colMeans(mat[keys$sample_ids[[i]], , drop = FALSE])
  out
}

gen_agg <- agg_matrix(gen_abund, agg_key)
sp_agg  <- agg_matrix(sp_abund, agg_key)

# Temperature (30-day mean)
daily_temps <- readRDS(file.path(data_dir, "becsei_climate.rds"))

ts_data <- agg_key %>%
  select(city, collection_date)

ts_data$T_air_30d <- NA_real_
for (i in seq_len(nrow(ts_data))) {
  ct <- daily_temps %>% filter(city == ts_data$city[i])
  w30 <- ct$t_mean[ct$date >= (ts_data$collection_date[i] - 30) &
                    ct$date <= ts_data$collection_date[i]]
  if (length(w30) > 0) ts_data$T_air_30d[i] <- mean(w30, na.rm = TRUE)
}

# Relative abundances
gen_total <- rowSums(gen_agg)
sp_total  <- rowSums(sp_agg)

# =============================================================================
# GENUS/SPECIES DEFINITIONS (same as Fig 2)
# =============================================================================

panel_defs <- list(
  # WHO BPPL 2024 Critical
  list(genus = "Acinetobacter", pattern = "^Acinetobacter", tag = "A",
       species = c("A. baumannii" = "Acinetobacter baumannii",
                   "A. johnsonii" = "Acinetobacter johnsonii",
                   "A. lwoffii" = "Acinetobacter lwoffii")),
  list(genus = "Klebsiella", pattern = "^Klebsiella", tag = "B",
       species = c("K. pneumoniae" = "Klebsiella pneumoniae")),
  list(genus = "Escherichia", pattern = "^Escherichia", tag = "C",
       species = c("E. coli" = "Escherichia coli")),
  list(genus = "Enterobacter", pattern = "^Enterobacter", tag = "D",
       species = c("E. cloacae" = "Enterobacter cloacae")),
  list(genus = "Mycobacterium", pattern = "^Mycobacterium", tag = "E",
       species = c("M. tuberculosis" = "Mycobacterium tuberculosis")),
  # WHO BPPL 2024 High
  list(genus = "Enterococcus", pattern = "^Enterococcus", tag = "F",
       species = c("E. faecium" = "Enterococcus faecium",
                   "E. faecalis" = "Enterococcus faecalis")),
  list(genus = "Staphylococcus", pattern = "^Staphylococcus", tag = "G",
       species = c("S. aureus" = "Staphylococcus aureus")),
  list(genus = "Pseudomonas", pattern = "^Pseudomonas", tag = "H",
       species = c("P. aeruginosa" = "Pseudomonas aeruginosa")),
  list(genus = "Salmonella", pattern = "^Salmonella", tag = "I",
       species = c("S. enterica" = "Salmonella enterica")),
  list(genus = "Neisseria", pattern = "^Neisseria", tag = "J",
       species = c("N. gonorrhoeae" = "Neisseria gonorrhoeae")),
  # WHO BPPL 2024 Medium
  list(genus = "Streptococcus", pattern = "^Streptococcus", tag = "K",
       species = c("S. pneumoniae" = "Streptococcus pneumoniae",
                   "S. pyogenes" = "Streptococcus pyogenes")),
  list(genus = "Haemophilus", pattern = "^Haemophilus", tag = "L",
       species = c("H. influenzae" = "Haemophilus influenzae"))
)

# =============================================================================
# HELPER: build one genus panel for a given city subset
# =============================================================================

make_city_panel <- function(def, city_idx, temp_vals) {
  # Genus RA
  gen_cols <- grep(def$pattern, colnames(gen_agg), value = TRUE)
  genus_ra <- rowSums(gen_agg[city_idx, gen_cols, drop = FALSE]) / gen_total[city_idx]

  plot_df <- tibble(
    taxon = paste0(def$genus, " (all spp.)"),
    temp = temp_vals,
    ra = genus_ra
  )

  for (i in seq_along(def$species)) {
    sp_col <- grep(paste0("^", def$species[i], "$"), colnames(sp_agg), value = TRUE)
    sp_ra <- if (length(sp_col) > 0) sp_agg[city_idx, sp_col[1]] / sp_total[city_idx] else rep(0, length(city_idx))
    plot_df <- bind_rows(plot_df, tibble(
      taxon = names(def$species)[i],
      temp = temp_vals,
      ra = sp_ra
    ))
  }

  # Spearman rho + p for each taxon
  stats <- plot_df %>%
    group_by(taxon) %>%
    summarise(
      rho = suppressWarnings(cor.test(ra, temp, method = "spearman")$estimate),
      p = suppressWarnings(cor.test(ra, temp, method = "spearman")$p.value),
      .groups = "drop"
    ) %>%
    mutate(
      p_label = ifelse(p < 0.001, "p<0.001",
                ifelse(p < 0.01, sprintf("p=%.3f", p),
                       sprintf("p=%.2f", p))),
      legend_label = sprintf("%s\nrho=%+.2f, %s", taxon, rho, p_label)
    )

  plot_df <- plot_df %>%
    left_join(stats %>% select(taxon, legend_label), by = "taxon")

  # Colours
  taxa <- unique(plot_df$legend_label)
  is_genus <- str_detect(taxa, "all spp")
  colour_map <- setNames(rep("grey55", length(taxa)), taxa)
  sp_idx <- which(!is_genus)
  for (j in seq_along(sp_idx))
    colour_map[taxa[sp_idx[j]]] <- sp_cols[min(j, length(sp_cols))]

  plot_nonzero <- plot_df %>% filter(ra > 0)

  ggplot(plot_nonzero, aes(x = temp, y = log10(ra), colour = legend_label)) +
    geom_point(alpha = 0.15, size = 0.5) +
    geom_smooth(method = "lm", se = TRUE, linewidth = 0.7, alpha = 0.15) +
    scale_colour_manual(values = colour_map, name = NULL) +
    guides(colour = guide_legend(override.aes = list(alpha = 1, size = 0, linewidth = 1))) +
    labs(tag = def$tag,
         title = def$genus,
         x = "30-day mean air temperature (\u00B0C)",
         y = "Relative abundance (log10)") +
    theme(legend.key = element_rect(fill = NA, colour = NA),
          legend.position = "inside",
          legend.position.inside = c(0.98, 0.02),
          legend.justification = c(1, 0),
          legend.background = element_rect(fill = alpha("white", 0.85), colour = NA),
          legend.key.size = unit(0.25, "cm"),
          legend.text = element_text(size = 4.5),
          legend.spacing.y = unit(0.05, "cm"))
}

# =============================================================================
# GENERATE ONE FIGURE PER CITY
# =============================================================================

cities <- c("Copenhagen", "Rotterdam", "Budapest", "Bologna", "Rome")
sfig_nums <- c("SFig15", "SFig16", "SFig17", "SFig18", "SFig19")

for (ci in seq_along(cities)) {
  city_name <- cities[ci]
  cat(sprintf("--- %s: %s ---\n", sfig_nums[ci], city_name))

  city_rows <- which(ts_data$city == city_name & !is.na(ts_data$T_air_30d))
  temp_city <- ts_data$T_air_30d[city_rows]

  if (length(city_rows) < 10) {
    cat("  Skipping (too few samples)\n")
    next
  }

  panels <- lapply(panel_defs, make_city_panel,
                   city_idx = city_rows, temp_vals = temp_city)

  fig <- (panels[[1]] | panels[[2]] | panels[[3]] | panels[[4]]) /
         (panels[[5]] | panels[[6]] | panels[[7]] | panels[[8]]) /
         (panels[[9]] | panels[[10]] | panels[[11]] | panels[[12]]) +
    plot_annotation(
      title = sprintf("%s (n=%d samples)", city_name, length(city_rows)),
      theme = theme(plot.title = element_text(face = "bold", size = 10))
    )

  fname <- sprintf("%s_%s_species.pdf", sfig_nums[ci],
                   tolower(gsub(" ", "_", city_name)))
  ggsave(file.path(fig_dir, fname), fig, width = 14, height = 10.5)
  cat(sprintf("  Saved %s\n", fname))
}

cat("\nDone.\n")
