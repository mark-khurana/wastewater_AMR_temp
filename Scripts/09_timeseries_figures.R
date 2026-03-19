# Temporal replication figures (Figure 3, Figures S3–S4).

library(tidyverse)
library(patchwork)
library(scales)
library(lubridate)
library(compositions)
library(showtext)
library(here)
library(vegan)
library(ggrepel)
showtext_auto()
font_add_google("Noto Sans", "notosans")

theme_lph <- theme_bw(base_size = 10, base_family = "notosans") +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.25, colour = "grey92"),
    legend.position = "bottom",
    legend.key.size = unit(0.35, "cm"),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 9, face = "bold"),
    plot.title = element_text(face = "bold", size = 10),
    axis.text = element_text(size = 8.5),
    axis.title = element_text(size = 9.5),
    plot.tag = element_text(face = "bold", size = 12)
  )
theme_set(theme_lph)

city_cols <- c("Copenhagen" = "#0072B2", "Rotterdam" = "#009E73",
               "Budapest" = "#E69F00", "Bologna" = "#D55E00", "Rome" = "#CC79A7")

concordance <- read.csv(here("Results", "timeseries_concordance.csv"))

meta <- read.csv(here("Datasets", "Becsei_et_al", "meta_location.csv"), stringsAsFactors = FALSE)
gen_abund_raw <- read.csv(here("Datasets", "Becsei_et_al", "genomic_genus_abundance.csv"),
                           row.names = 1, check.names = FALSE)
res_class_raw <- read.csv(here("Datasets", "Becsei_et_al", "resfinder_class_abundance.csv"),
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

gen_agg <- agg_matrix(gen_abund_raw, agg_key)
res_agg <- agg_matrix(res_class_raw, agg_key)

do_clr <- function(mat) {
  mat <- as.data.frame(mat)
  for (j in seq_len(ncol(mat))) {
    nz <- mat[[j]][mat[[j]] > 0]
    if (length(nz) > 0) mat[[j]][mat[[j]] == 0] <- min(nz) / 2
  }
  out <- as.data.frame(clr(as.matrix(mat)))
  colnames(out) <- colnames(mat)
  out
}

gen_mat <- as.data.frame(gen_agg)
keep_gen <- (colSums(gen_mat > 0) / nrow(gen_mat)) > 0.2 & colSums(gen_mat) > 100
gen_mat <- gen_mat[, keep_gen]
gen_clr <- do_clr(gen_mat)
res_clr <- do_clr(res_agg)

gen_ra <- gen_mat / rowSums(gen_mat)
gen_bray <- vegdist(gen_ra, method = "bray")
gen_pcoa <- cmdscale(gen_bray, k = 5, eig = TRUE)

daily_temps <- readRDS(here("Datasets", "Becsei_et_al", "becsei_climate.rds"))
city_coords <- tibble(
  city = c("Copenhagen", "Rotterdam", "Bologna", "Budapest", "Rome"),
  lat  = c(55.6761, 51.9225, 44.4949, 47.4979, 41.9028),
  lon  = c(12.5683, 4.4792, 11.3426, 19.0402, 12.4964)
)

ts_data <- agg_key %>%
  select(city, collection_date) %>%
  mutate(month = month(collection_date)) %>%
  left_join(city_coords, by = "city")

ts_data$T_air_30d <- NA_real_
for (i in seq_len(nrow(ts_data))) {
  ct <- daily_temps %>% filter(city == ts_data$city[i])
  w30 <- ct$t_mean[ct$date >= (ts_data$collection_date[i] - 30) &
                    ct$date <= ts_data$collection_date[i]]
  if (length(w30) > 0) ts_data$T_air_30d[i] <- mean(w30, na.rm = TRUE)
}

ts_data$city_f <- factor(ts_data$city,
  levels = c("Copenhagen", "Rotterdam", "Budapest", "Bologna", "Rome"))
ts_data$bac_PC1 <- gen_pcoa$points[, 1]
ts_data$bac_PC2 <- gen_pcoa$points[, 2]
ts_data$bac_PC3 <- gen_pcoa$points[, 3]

part1_genera <- c("Acinetobacter", "Klebsiella", "Escherichia", "Enterococcus",
                  "Pseudomonas", "Staphylococcus", "Salmonella", "Clostridioides")

adj_cross <- read.csv(here("Results", "adjusted_species_temp.csv")) %>%
  filter(grepl("all spp\\.|E\\. coli", display)) %>%
  mutate(genus = genus) %>%
  select(genus, rho_cross_adj = rho_adj, p_cross_adj = p_adj)

conc_plot <- concordance %>%
  mutate(is_part1 = genus %in% part1_genera) %>%
  left_join(adj_cross, by = "genus")

part1_conc <- conc_plot %>% filter(is_part1)
n_conc_part1 <- sum(part1_conc$concordant)
n_part1 <- nrow(part1_conc)

fig3a <- ggplot(conc_plot %>% filter(is_part1),
                aes(x = rho_cross_adj, y = rho_temporal)) +
  geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = "grey70", linewidth = 0.3) +
  geom_abline(slope = 1, intercept = 0, colour = "grey85", linetype = "dashed", linewidth = 0.3) +
  geom_point(size = 3, colour = "#2C3E50") +
  geom_text_repel(
    aes(label = genus), size = 3.2, fontface = "italic",
    max.overlaps = 25, segment.size = 0.2, segment.colour = "grey50",
    min.segment.length = 0.1
  ) +
  scale_x_continuous(limits = c(-0.7, 0.7)) +
  scale_y_continuous(limits = c(-0.7, 0.7)) +
  coord_equal() +
  labs(tag = "A",
       x = expression(paste("Cross-sectional ", rho, " (adjusted)")),
       y = expression(paste("Within-city ", rho))) +
  theme(plot.margin = margin(5.5, 5.5, 5.5, 5.5))

mod_idx <- which(!is.na(ts_data$T_air_30d))
res_dist <- dist(res_clr[mod_idx, ])

perm_df <- ts_data[mod_idx, ] %>%
  mutate(city_f = factor(city), month_f = factor(month))
perm_df$bac_PC1 <- gen_pcoa$points[mod_idx, 1]
perm_df$bac_PC2 <- gen_pcoa$points[mod_idx, 2]
perm_df$bac_PC3 <- gen_pcoa$points[mod_idx, 3]
perm_df$bac_PC4 <- gen_pcoa$points[mod_idx, 4]
perm_df$bac_PC5 <- gen_pcoa$points[mod_idx, 5]

set.seed(42)
perm_total <- adonis2(res_dist ~ T_air_30d + city_f + month_f,
                       data = perm_df, permutations = 999)
r2_total <- perm_total["T_air_30d", "R2"]
p_total  <- perm_total["T_air_30d", "Pr(>F)"]

set.seed(42)
perm_direct <- adonis2(res_dist ~ bac_PC1 + bac_PC2 + bac_PC3 + bac_PC4 + bac_PC5 + T_air_30d + city_f + month_f,
                        data = perm_df, permutations = 999)
r2_direct <- perm_direct["T_air_30d", "R2"]
p_direct  <- perm_direct["T_air_30d", "Pr(>F)"]

r2_indirect <- r2_total - r2_direct
pct_mediated <- round(100 * r2_indirect / r2_total, 0)

med_data <- tibble(
  path = factor(c("Indirect\n(via bacteriome)", "Direct"),
                levels = c("Indirect\n(via bacteriome)", "Direct")),
  R2 = c(r2_indirect, r2_direct)
)

fig3b <- ggplot(med_data, aes(x = "Resistome", y = R2, fill = path)) +
  geom_col(position = "stack", width = 0.45) +
  scale_fill_manual(values = c("Direct" = "#CC6677",
                                "Indirect\n(via bacteriome)" = "#88CCEE")) +
  scale_y_continuous(labels = percent_format(accuracy = 0.01),
                     expand = expansion(mult = c(0, 0.25))) +
  geom_text(aes(label = sprintf("%.2f%%", 100 * R2)),
            position = position_stack(vjust = 0.5), size = 3.5) +
  annotate("text", x = 1, y = -0.003,
           label = sprintf("PERMANOVA p=%s\nadjusted for city + month",
                           format.pval(p_total, digits = 2)),
           vjust = 1, size = 2.8, colour = "grey40") +
  labs(tag = "B",
       title = "Temporal resistome mediation\nvia bacteriome (5 cities)",
       y = "Temperature R\u00B2",
       x = NULL, fill = NULL) +
  theme(legend.position = "bottom",
        legend.key.size = unit(0.35, "cm"),
        legend.text = element_text(size = 8),
        legend.margin = margin(t = -2),
        plot.title = element_text(size = 9.5))

fig2_genera <- c("Acinetobacter", "Klebsiella", "Escherichia", "Enterococcus",
                 "Pseudomonas", "Staphylococcus", "Salmonella", "Clostridioides")

gen_stats <- tibble()
gen_city_stats <- tibble()

for (g in fig2_genera) {
  plot_idx <- which(!is.na(ts_data$T_air_30d))
  temp_vals <- ts_data$T_air_30d[plot_idx]
  clr_vals <- gen_clr[[g]][plot_idx]
  city_vals <- factor(ts_data$city[plot_idx])

  temp_resid <- residuals(lm(temp_vals ~ city_vals))
  clr_resid <- residuals(lm(clr_vals ~ city_vals))
  ct <- cor.test(temp_resid, clr_resid, method = "spearman")

  gen_stats <- bind_rows(gen_stats, tibble(
    genus = g, rho = ct$estimate, p_value = ct$p.value
  ))

  for (cc in levels(city_vals)) {
    city_idx <- which(city_vals == cc)
    if (length(city_idx) < 8) next
    ct_c <- suppressWarnings(cor.test(temp_vals[city_idx], clr_vals[city_idx],
                                       method = "spearman"))
    gen_city_stats <- bind_rows(gen_city_stats, tibble(
      genus = g, city = cc,
      rho_city = ct_c$estimate, p_city = ct_c$p.value
    ))
  }
}

gen_stats <- gen_stats %>%
  mutate(
    genus = factor(genus, levels = fig2_genera),
    p_label = ifelse(p_value < 0.001, "p<0.001",
              ifelse(p_value < 0.01, sprintf("p=%.3f", p_value),
                     sprintf("p=%.2f", p_value))),
    overall_anno = paste0("Overall: \u03C1=", sprintf("%+.2f", rho), ", ", p_label)
  )

genera_plot <- ts_data %>%
  filter(!is.na(T_air_30d)) %>%
  select(city, T_air_30d, city_f)

for (g in fig2_genera) {
  genera_plot[[g]] <- gen_clr[[g]][!is.na(ts_data$T_air_30d)]
}

genera_long <- genera_plot %>%
  pivot_longer(cols = all_of(fig2_genera), names_to = "genus", values_to = "clr") %>%
  mutate(genus = factor(genus, levels = fig2_genera))

city_short <- c("Copenhagen" = "CPH", "Rotterdam" = "ROT",
                "Budapest" = "BUD", "Bologna" = "BOL", "Rome" = "ROM")

city_rho_df <- gen_city_stats %>%
  mutate(
    short = city_short[city],
    city_f = factor(city, levels = names(city_short)),
    sig = ifelse(p_city < 0.001, "***",
          ifelse(p_city < 0.01, "**",
          ifelse(p_city < 0.05, "*", "")))
  ) %>%
  arrange(genus, city_f)

city_anno <- city_rho_df %>%
  group_by(genus) %>%
  summarise(
    line2 = paste(sprintf("%s %+.2f%s", short[1:3], rho_city[1:3], sig[1:3]), collapse = " "),
    line3 = paste(sprintf("%s %+.2f%s", short[4:5], rho_city[4:5], sig[4:5]), collapse = " "),
    .groups = "drop"
  )

anno_lines <- gen_stats %>%
  left_join(city_anno, by = "genus") %>%
  mutate(full_anno = sprintf("%s\n%s\n%s", overall_anno, line2, line3))

anno_df <- anno_lines %>%
  left_join(
    genera_long %>% group_by(genus) %>%
      summarise(x_pos = -Inf, y_pos = Inf, .groups = "drop"),
    by = "genus"
  )

fig3c <- ggplot(genera_long, aes(x = T_air_30d, y = clr)) +
  geom_point(aes(colour = city_f), size = 0.8, alpha = 0.3) +
  geom_smooth(aes(colour = city_f), method = "lm", se = TRUE,
              linewidth = 0.7, alpha = 0.15) +
  geom_label(data = anno_df,
             aes(x = x_pos, y = y_pos, label = full_anno),
             inherit.aes = FALSE, hjust = 0, vjust = 1, size = 2.8,
             colour = "grey25", lineheight = 1.1,
             fill = alpha("white", 0.8), label.size = 0,
             label.padding = unit(0.2, "lines")) +
  facet_wrap(~genus, scales = "free_y", ncol = 2) +
  scale_colour_manual(values = city_cols, name = "City") +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.08))) +
  labs(tag = "C",
       x = "30-day mean air temperature (\u00B0C)",
       y = "Abundance (CLR-transformed)") +
  theme(legend.position = "bottom",
        legend.margin = margin(t = -2),
        strip.text = element_text(face = "bold.italic", size = 9))

left_col <- (fig3a / fig3b) + plot_layout(heights = c(1.2, 1))
fig3 <- (left_col | fig3c) + plot_layout(widths = c(1, 1.8))

ggsave(here("Figures", "Fig3_temporal_replication.pdf"), fig3,
       width = 12, height = 10)
ggsave(here("Figures", "Fig3_temporal_replication.png"), fig3,
       width = 12, height = 10, dpi = 300)

ts_dm <- ts_data %>%
  filter(!is.na(T_air_30d)) %>%
  mutate(BPC1 = gen_pcoa$points[, 1]) %>%
  group_by(city) %>%
  mutate(bac_PC1_dm = BPC1 - mean(BPC1)) %>%
  ungroup()

dm_range <- range(ts_dm$bac_PC1_dm, na.rm = TRUE)
t_range <- range(ts_dm$T_air_30d, na.rm = TRUE)
dm_to_t <- function(x) (x - dm_range[1]) / diff(dm_range) * diff(t_range) + t_range[1]
t_to_dm <- function(x) (x - t_range[1]) / diff(t_range) * diff(dm_range) + dm_range[1]

sfig10 <- ggplot(ts_dm, aes(x = collection_date)) +
  geom_line(aes(y = T_air_30d, colour = "Temperature"),
            linewidth = 0.6, alpha = 0.8) +
  geom_point(aes(y = T_air_30d, colour = "Temperature"), size = 0.8, alpha = 0.6) +
  geom_line(aes(y = dm_to_t(bac_PC1_dm), colour = "Bacteriome PC1\n(city-demeaned)"),
            linewidth = 0.6, alpha = 0.8, linetype = "dashed") +
  geom_point(aes(y = dm_to_t(bac_PC1_dm), colour = "Bacteriome PC1\n(city-demeaned)"),
             size = 0.8, alpha = 0.6) +
  scale_colour_manual(values = c("Temperature" = "#D73027",
                                  "Bacteriome PC1\n(city-demeaned)" = "#4575B4"),
                      name = NULL) +
  scale_y_continuous(
    name = "30-day mean air temperature (\u00B0C)",
    sec.axis = sec_axis(~ t_to_dm(.), name = "Bacteriome PC1 (demeaned)")
  ) +
  facet_wrap(~city_f, nrow = 1, scales = "free_x") +
  labs(x = NULL) +
  theme(legend.position = "bottom",
        strip.text = element_text(face = "bold", size = 7),
        axis.title.y.right = element_text(colour = "#4575B4", size = 7),
        axis.text.y.right = element_text(colour = "#4575B4"),
        axis.title.y.left = element_text(colour = "#D73027", size = 7),
        axis.text.y.left = element_text(colour = "#D73027"),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 5))

ggsave(here("Figures", "SFig10_city_timeseries.pdf"), sfig10,
       width = 11, height = 3.5)

drug_results <- read.csv(here("Results", "timeseries_drug_models.csv"))
drug_plot <- drug_results %>%
  mutate(
    sig = p_adj < 0.05,
    point_col = ifelse(sig, "#2C3E50", "grey55"),
    p_label = ifelse(p_value < 0.001, "***",
              ifelse(p_value < 0.01, "**",
              ifelse(p_value < 0.05, "*", "")))
  )

sfig11 <- ggplot(drug_plot, aes(x = beta_temp, y = reorder(drug_class, beta_temp))) +
  geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_segment(aes(xend = 0, yend = reorder(drug_class, beta_temp)),
               colour = "grey75", linewidth = 0.4) +
  geom_point(aes(fill = point_col), size = 2.8,
             shape = 21, colour = "black", stroke = 0.3) +
  geom_text(aes(label = p_label, x = beta_temp + sign(beta_temp) * 0.005),
            hjust = ifelse(drug_plot$beta_temp > 0, 0, 1), size = 3) +
  scale_fill_identity() +
  labs(title = "Resistome drug-class associations with temperature (within-city)",
       x = "Abundance change per \u00B0C (adjusted for city + month)",
       y = NULL)

ggsave(here("Figures", "SFig11_drug_class.pdf"), sfig11,
       width = 7, height = 4)

top_drugs <- c("Beta-Lactam", "Fosfomycin", "Folate Pathway Antagonist",
               "Quinolone", "Tetracycline", "Aminoglycoside")

ts_drugs <- ts_data %>%
  filter(!is.na(T_air_30d)) %>%
  select(city, T_air_30d, city_f)
for (dc in top_drugs) ts_drugs[[dc]] <- res_clr[[dc]][!is.na(ts_data$T_air_30d)]

ts_drugs_long <- ts_drugs %>%
  pivot_longer(cols = all_of(top_drugs), names_to = "drug_class", values_to = "clr") %>%
  mutate(drug_class = factor(drug_class, levels = top_drugs))

sfig12 <- ggplot(ts_drugs_long, aes(x = T_air_30d, y = clr, colour = city_f)) +
  geom_point(size = 0.8, alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.5, alpha = 0.12) +
  facet_wrap(~drug_class, scales = "free_y", nrow = 2) +
  scale_colour_manual(values = city_cols, name = "City") +
  labs(x = "30-day mean air temperature (\u00B0C)",
       y = "Abundance (CLR-transformed)") +
  theme(legend.position = "bottom",
        strip.text = element_text(face = "bold", size = 7))

ggsave(here("Figures", "SFig12_city_drug_class.pdf"), sfig12,
       width = 9, height = 5)

ts_data$T_air_day <- NA_real_
for (i in seq_len(nrow(ts_data))) {
  ct <- daily_temps %>% filter(city == ts_data$city[i])
  hit <- which(ct$date == ts_data$collection_date[i])
  if (length(hit) > 0) ts_data$T_air_day[i] <- ct$t_mean[hit[1]]
}

genera_plot_day <- ts_data %>%
  filter(!is.na(T_air_day)) %>%
  select(city, T_air_day, city_f)

for (g in fig2_genera) {
  genera_plot_day[[g]] <- gen_clr[[g]][!is.na(ts_data$T_air_day)]
}

genera_long_day <- genera_plot_day %>%
  pivot_longer(cols = all_of(fig2_genera), names_to = "genus", values_to = "clr") %>%
  mutate(genus = factor(genus, levels = fig2_genera))

idx_day <- which(!is.na(ts_data$T_air_day))
temp_day <- ts_data$T_air_day[idx_day]
city_day <- factor(ts_data$city[idx_day])

gen_stats_day <- tibble()
gen_city_day <- tibble()

for (g in fig2_genera) {
  clr_vals <- gen_clr[[g]][idx_day]

  t_r <- residuals(lm(temp_day ~ city_day))
  c_r <- residuals(lm(clr_vals ~ city_day))
  ct <- cor.test(t_r, c_r, method = "spearman")
  gen_stats_day <- bind_rows(gen_stats_day, tibble(
    genus = g, rho = ct$estimate, p_value = ct$p.value
  ))

  for (cc in levels(city_day)) {
    ci <- which(city_day == cc)
    if (length(ci) < 8) next
    ct_c <- suppressWarnings(cor.test(temp_day[ci], clr_vals[ci], method = "spearman"))
    gen_city_day <- bind_rows(gen_city_day, tibble(
      genus = g, city = cc, rho_city = ct_c$estimate, p_city = ct_c$p.value
    ))
  }
}

gen_stats_day <- gen_stats_day %>%
  mutate(
    genus = factor(genus, levels = fig2_genera),
    p_label = ifelse(p_value < 0.001, "p<0.001",
              ifelse(p_value < 0.01, sprintf("p=%.3f", p_value),
                     sprintf("p=%.2f", p_value))),
    overall_anno = paste0("Overall: \u03C1=", sprintf("%+.2f", rho), ", ", p_label)
  )

city_rho_day <- gen_city_day %>%
  mutate(
    short = city_short[city],
    city_f = factor(city, levels = names(city_short)),
    sig = ifelse(p_city < 0.001, "***",
          ifelse(p_city < 0.01, "**",
          ifelse(p_city < 0.05, "*", "")))
  ) %>%
  arrange(genus, city_f)

city_anno_day <- city_rho_day %>%
  group_by(genus) %>%
  summarise(
    line2 = paste(sprintf("%s %+.2f%s", short[1:3], rho_city[1:3], sig[1:3]), collapse = " "),
    line3 = paste(sprintf("%s %+.2f%s", short[4:5], rho_city[4:5], sig[4:5]), collapse = " "),
    .groups = "drop"
  )

anno_day <- gen_stats_day %>%
  left_join(city_anno_day, by = "genus") %>%
  mutate(
    full_anno = sprintf("%s\n%s\n%s", overall_anno, line2, line3),
    x_pos = -Inf, y_pos = Inf
  )

sfig14 <- ggplot(genera_long_day, aes(x = T_air_day, y = clr)) +
  geom_point(aes(colour = city_f), size = 0.8, alpha = 0.5) +
  geom_smooth(aes(colour = city_f), method = "lm", se = TRUE,
              linewidth = 0.6, alpha = 0.12) +
  geom_label(data = anno_day,
             aes(x = x_pos, y = y_pos, label = full_anno),
             inherit.aes = FALSE, hjust = 0, vjust = 1, size = 1.65,
             colour = "grey25", lineheight = 1.1,
             fill = alpha("white", 0.8), label.size = 0,
             label.padding = unit(0.2, "lines")) +
  facet_wrap(~genus, scales = "free_y", nrow = 2) +
  scale_colour_manual(values = city_cols, name = "City") +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.08))) +
  labs(title = "Same-day air temperature (sensitivity analysis)",
       x = "Collection-day air temperature (\u00B0C)",
       y = "Abundance (CLR-transformed)") +
  theme(legend.position = "bottom",
        strip.text = element_text(face = "bold.italic", size = 7))

ggsave(here("Figures", "SFig14_sameday_temp.pdf"), sfig14,
       width = 11, height = 6)

mod_temp <- ts_data$T_air_30d[mod_idx]
mod_city <- factor(ts_data$city[mod_idx])
mod_month <- factor(ts_data$month[mod_idx])

gen_stats_madj <- tibble()
gen_city_madj <- tibble()

for (g in fig2_genera) {
  clr_vals <- gen_clr[[g]][mod_idx]

  temp_resid <- residuals(lm(mod_temp ~ mod_city + mod_month))
  clr_resid <- residuals(lm(clr_vals ~ mod_city + mod_month))
  ct <- cor.test(temp_resid, clr_resid, method = "spearman")
  gen_stats_madj <- bind_rows(gen_stats_madj, tibble(
    genus = g, rho = ct$estimate, p_value = ct$p.value
  ))

  for (cc in levels(mod_city)) {
    ci <- which(mod_city == cc)
    if (length(ci) < 10) next
    m_c <- mod_month[ci]
    if (length(unique(m_c)) > 2) {
      t_r <- residuals(lm(mod_temp[ci] ~ m_c))
      c_r <- residuals(lm(clr_vals[ci] ~ m_c))
    } else {
      t_r <- mod_temp[ci]
      c_r <- clr_vals[ci]
    }
    ct_c <- suppressWarnings(cor.test(t_r, c_r, method = "spearman"))
    gen_city_madj <- bind_rows(gen_city_madj, tibble(
      genus = g, city = cc, rho_city = ct_c$estimate, p_city = ct_c$p.value
    ))
  }
}

gen_stats_madj <- gen_stats_madj %>%
  mutate(
    genus = factor(genus, levels = fig2_genera),
    p_label = ifelse(p_value < 0.001, "p<0.001",
              ifelse(p_value < 0.01, sprintf("p=%.3f", p_value),
                     sprintf("p=%.2f", p_value))),
    overall_anno = paste0("Overall: \u03C1=", sprintf("%+.2f", rho), ", ", p_label)
  )

city_rho_madj <- gen_city_madj %>%
  mutate(
    short = city_short[city],
    city_f = factor(city, levels = names(city_short)),
    sig = ifelse(p_city < 0.001, "***",
          ifelse(p_city < 0.01, "**",
          ifelse(p_city < 0.05, "*", "")))
  ) %>%
  arrange(genus, city_f)

city_anno_madj <- city_rho_madj %>%
  group_by(genus) %>%
  summarise(
    line2 = paste(sprintf("%s %+.2f%s", short[1:3], rho_city[1:3], sig[1:3]), collapse = " "),
    line3 = paste(sprintf("%s %+.2f%s", short[4:5], rho_city[4:5], sig[4:5]), collapse = " "),
    .groups = "drop"
  )

anno_madj <- gen_stats_madj %>%
  left_join(city_anno_madj, by = "genus") %>%
  mutate(full_anno = sprintf("%s\n%s\n%s", overall_anno, line2, line3),
         x_pos = -Inf, y_pos = Inf)

sfig20 <- ggplot(genera_long, aes(x = T_air_30d, y = clr)) +
  geom_point(aes(colour = city_f), size = 0.8, alpha = 0.5) +
  geom_smooth(aes(colour = city_f), method = "lm", se = TRUE,
              linewidth = 0.6, alpha = 0.12) +
  geom_label(data = anno_madj,
             aes(x = x_pos, y = y_pos, label = full_anno),
             inherit.aes = FALSE, hjust = 0, vjust = 1, size = 1.65,
             colour = "grey25", lineheight = 1.1,
             fill = alpha("white", 0.8), label.size = 0,
             label.padding = unit(0.2, "lines")) +
  facet_wrap(~genus, scales = "free_y", nrow = 2) +
  scale_colour_manual(values = city_cols, name = "City") +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.08))) +
  labs(title = "Month-adjusted sensitivity (city + month residualization)",
       x = "30-day mean air temperature (\u00B0C)",
       y = "Abundance (CLR-transformed)") +
  theme(legend.position = "bottom",
        strip.text = element_text(face = "bold.italic", size = 7))

ggsave(here("Figures", "SFig20_month_adjusted.pdf"), sfig20,
       width = 11, height = 6)
