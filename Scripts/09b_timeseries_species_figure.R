# =============================================================================
# 09b_timeseries_species_figure.R — Fig 5: Species-level temporal patterns
# =============================================================================
#
# Mirrors Fig 2 from Part 1: 8 WHO priority genera, each panel showing
# genus total (grey) + key species (coloured) vs temperature.
# Uses Becsei et al. 5-city biweekly data (n=162 city-dates).
#
# Species correlations are city-adjusted (partial Spearman).
# =============================================================================

library(tidyverse)
library(patchwork)
library(scales)
library(lubridate)
library(compositions)
library(ggrepel)
library(showtext)
showtext_auto()
font_add_google("Noto Sans", "notosans")

base_dir <- "/Users/bxc331/Desktop/AMR_open_project"
data_dir <- file.path(base_dir, "Datasets/Becsei_et_al")
fig_dir  <- file.path(base_dir, "Figures")

cat("=== 09b_timeseries_species_figure.R ===\n\n")

theme_lph <- theme_bw(base_size = 8, base_family = "notosans") +
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

# Aggregate replicates
agg_matrix <- function(mat, keys) {
  out <- matrix(NA, nrow(keys), ncol(mat)); colnames(out) <- colnames(mat)
  for (i in seq_len(nrow(keys)))
    out[i, ] <- colMeans(mat[keys$sample_ids[[i]], , drop = FALSE])
  out
}

gen_agg <- agg_matrix(gen_abund, agg_key)
sp_agg  <- agg_matrix(sp_abund, agg_key)

# Temperature
daily_temps <- readRDS(file.path(data_dir, "becsei_climate.rds"))
city_coords <- tibble(
  city = c("Copenhagen", "Rotterdam", "Bologna", "Budapest", "Rome"),
  lat  = c(55.6761, 51.9225, 44.4949, 47.4979, 41.9028),
  lon  = c(12.5683, 4.4792, 11.3426, 19.0402, 12.4964)
)

ts_data <- agg_key %>%
  select(city, collection_date) %>%
  left_join(city_coords, by = "city")

ts_data$T_air_30d <- NA_real_
for (i in seq_len(nrow(ts_data))) {
  ct <- daily_temps %>% filter(city == ts_data$city[i])
  w30 <- ct$t_mean[ct$date >= (ts_data$collection_date[i] - 30) &
                    ct$date <= ts_data$collection_date[i]]
  if (length(w30) > 0) ts_data$T_air_30d[i] <- mean(w30, na.rm = TRUE)
}

idx <- which(!is.na(ts_data$T_air_30d))
temp <- ts_data$T_air_30d[idx]
city_f <- factor(ts_data$city[idx])

# Relative abundance (proportion of total reads)
gen_total <- rowSums(gen_agg[idx, ])
sp_total  <- rowSums(sp_agg[idx, ])

# Helper: get relative abundance for a taxon (genus or species)
get_gen_ra <- function(pattern) {
  cols <- grep(pattern, colnames(gen_agg), value = TRUE)
  if (length(cols) == 0) return(rep(0, length(idx)))
  rowSums(gen_agg[idx, cols, drop = FALSE]) / gen_total
}

get_sp_ra <- function(name) {
  cols <- grep(paste0("^", name, "$"), colnames(sp_agg), value = TRUE)
  if (length(cols) == 0) return(rep(0, length(idx)))
  sp_agg[idx, cols[1]] / sp_total
}

# City-adjusted partial Spearman
partial_rho <- function(y) {
  y_r <- residuals(lm(y ~ city_f))
  t_r <- residuals(lm(temp ~ city_f))
  ct <- suppressWarnings(cor.test(t_r, y_r, method = "spearman"))
  list(rho = ct$estimate, p = ct$p.value)
}

# =============================================================================
# BUILD PANELS — same structure as Fig 2
# =============================================================================

make_panel <- function(genus_display, genus_pattern, species_list, tag) {
  # genus_pattern: for genus-level columns
  # species_list: named vector c("display" = "species name in data")

  genus_ra <- get_gen_ra(genus_pattern)

  plot_df <- tibble(
    taxon = paste0(genus_display, " (all spp.)"),
    temp = temp,
    ra = genus_ra
  )

  for (i in seq_along(species_list)) {
    sp_ra <- get_sp_ra(species_list[i])
    plot_df <- bind_rows(plot_df, tibble(
      taxon = names(species_list)[i],
      temp = temp,
      ra = sp_ra
    ))
  }

  # Compute stats for each taxon
  stats <- plot_df %>%
    group_by(taxon) %>%
    summarise(
      rho_raw = partial_rho(ra)$rho,
      p_raw = partial_rho(ra)$p,
      .groups = "drop"
    ) %>%
    mutate(
      p_label = ifelse(p_raw < 0.001, "p<0.001", sprintf("p=%.3f", p_raw)),
      legend_label = paste0(taxon, "\n\u03C1=", sprintf("%+.2f", rho_raw), ", ", p_label, " (adj.)")
    )

  plot_df <- plot_df %>%
    left_join(stats %>% select(taxon, legend_label), by = "taxon")

  # Colour map: genus in grey, species in colour
  taxa <- unique(plot_df$legend_label)
  is_genus <- str_detect(taxa, "all spp")
  colour_map <- setNames(rep("grey55", length(taxa)), taxa)
  sp_idx <- which(!is_genus)
  for (j in seq_along(sp_idx))
    colour_map[taxa[sp_idx[j]]] <- sp_cols[min(j, length(sp_cols))]

  # Filter zeros for log plotting
  plot_nonzero <- plot_df %>% filter(ra > 0)

  ggplot(plot_nonzero, aes(x = temp, y = log10(ra), colour = legend_label)) +
    geom_point(alpha = 0.08, size = 0.3) +
    geom_smooth(method = "lm", se = TRUE, linewidth = 0.7, alpha = 0.15) +
    scale_colour_manual(values = colour_map, name = NULL) +
    guides(colour = guide_legend(override.aes = list(alpha = 1, size = 0, linewidth = 1))) +
    labs(tag = tag,
         title = genus_display,
         x = "30-day mean air temperature (\u00B0C)",
         y = "Relative abundance (log10)") +
    theme(legend.key = element_rect(fill = NA, colour = NA),
          legend.position = "inside",
          legend.position.inside = c(0.98, 0.02),
          legend.justification = c(1, 0),
          legend.background = element_rect(fill = alpha("white", 0.85), colour = NA),
          legend.key.size = unit(0.25, "cm"),
          legend.text = element_text(size = 5),
          legend.spacing.y = unit(0.1, "cm"))
}

cat("Building panels...\n")

# --- WHO BPPL 2024 Critical ---
p_a <- make_panel("Acinetobacter", "^Acinetobacter",
  c("A. baumannii" = "Acinetobacter baumannii",
    "A. johnsonii" = "Acinetobacter johnsonii",
    "A. lwoffii" = "Acinetobacter lwoffii"), "A")

p_b <- make_panel("Klebsiella", "^Klebsiella",
  c("K. pneumoniae" = "Klebsiella pneumoniae"), "B")

p_c <- make_panel("Escherichia", "^Escherichia",
  c("E. coli" = "Escherichia coli"), "C")

p_d <- make_panel("Enterobacter", "^Enterobacter",
  c("E. cloacae" = "Enterobacter cloacae"), "D")

p_e <- make_panel("Mycobacterium", "^Mycobacterium",
  c("M. tuberculosis" = "Mycobacterium tuberculosis"), "E")

# --- WHO BPPL 2024 High ---
p_f <- make_panel("Enterococcus", "^Enterococcus",
  c("E. faecium" = "Enterococcus faecium",
    "E. faecalis" = "Enterococcus faecalis"), "F")

p_g <- make_panel("Staphylococcus", "^Staphylococcus",
  c("S. aureus" = "Staphylococcus aureus"), "G")

p_h <- make_panel("Pseudomonas", "^Pseudomonas",
  c("P. aeruginosa" = "Pseudomonas aeruginosa"), "H")

p_i <- make_panel("Salmonella", "^Salmonella",
  c("S. enterica" = "Salmonella enterica"), "I")

p_j <- make_panel("Neisseria", "^Neisseria",
  c("N. gonorrhoeae" = "Neisseria gonorrhoeae"), "J")

# --- WHO BPPL 2024 Medium ---
p_k <- make_panel("Streptococcus", "^Streptococcus",
  c("S. pneumoniae" = "Streptococcus pneumoniae",
    "S. pyogenes" = "Streptococcus pyogenes"), "K")

p_l <- make_panel("Haemophilus", "^Haemophilus",
  c("H. influenzae" = "Haemophilus influenzae"), "L")

fig5 <- (p_a | p_b | p_c | p_d) / (p_e | p_f | p_g | p_h) / (p_i | p_j | p_k | p_l)

ggsave(file.path(fig_dir, "SFig13_species_temporal.pdf"), fig5,
       width = 14, height = 10.5)
cat("Saved SFig13\n")
