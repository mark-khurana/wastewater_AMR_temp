# Main figures (Figures 1–3) and supplementary PCoA ordination (Figure S2).

library(tidyverse)
library(patchwork)
library(scales)
library(showtext)
library(here)
showtext_auto()
font_add_google("Noto Sans", "notosans")

theme_lph <- theme_bw(base_size = 8, base_family = "notosans") +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.25, colour = "grey92"),
    legend.position = "bottom",
    legend.key.size = unit(0.3, "cm"),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 7, face = "bold"),
    plot.title = element_text(face = "bold.italic", size = 8),
    axis.text = element_text(size = 6.5),
    axis.title = element_text(size = 7.5),
    plot.tag = element_text(face = "bold", size = 10)
  )
theme_set(theme_lph)

pal_resistome <- c("Functional" = "#CC6677", "Acquired" = "#6699CC")
analysis <- readRDS(here("Datasets", "analysis_ready.rds"))

conf_vars <- c("T_30d", "gdp_pcap_ppp", "sanitation",
               "health_exp_gdp", "oop_health_exp", "immunization_dpt",
               "animal_amc_mgkg", "pop_density", "year", "Region")
resistome <- readRDS(here("Datasets", "resistome_matrices.rds"))
resistome_ids <- as.character(resistome$fg_cluster_clr$genepid)
analysis_cc <- analysis %>%
  mutate(genepid = as.character(genepid)) %>%
  filter(complete.cases(across(all_of(conf_vars))),
         genepid %in% resistome_ids)

samples_all <- readRDS(here("Datasets", "clean_samples.rds")) %>%
  mutate(genepid = as.character(genepid)) %>%
  select(-any_of("Region")) %>%
  left_join(analysis %>% mutate(genepid = as.character(genepid)) %>%
              select(genepid, T_30d, WHO_Region = Region), by = "genepid")

samples_main <- samples_all %>%
  filter(genepid %in% analysis_cc$genepid)

world <- map_data("world")

who_cols <- c("AFRO" = "#66C2A5", "AMRO" = "#FC8D62", "EMRO" = "#E78AC3",
              "EURO" = "#8DA0CB", "SEARO" = "#A6D854", "WPRO" = "#FFD92F")

fig1a <- ggplot() +
  geom_polygon(data = world, aes(x = long, y = lat, group = group),
               fill = "grey95", colour = "grey80", linewidth = 0.12) +
  geom_point(data = samples_main %>% filter(!is.na(WHO_Region)),
             aes(x = lon, y = lat, fill = WHO_Region),
             size = 1.2, alpha = 0.8, shape = 21, colour = "black", stroke = 0.15) +
  scale_fill_manual(values = who_cols, name = "WHO Region") +
  coord_fixed(1.3, xlim = c(-170, 180), ylim = c(-55, 75)) +
  labs(tag = "A") +
  theme_void(base_size = 8) +
  theme(legend.position = "right",
        legend.key.height = unit(0.3, "cm"), legend.key.width = unit(0.3, "cm"),
        legend.text = element_text(size = 6), legend.title = element_text(size = 7),
        plot.tag = element_text(face = "bold", size = 10))

layer2_perm <- readRDS(here("Results", "layer2_partial_permanova.rds"))
layer1_perm <- readRDS(here("Results", "layer1_permanova.rds"))

n_layer2 <- layer2_perm$perm_fg$n
step_data <- tibble(
  Resistome = rep(c("Functional", "Acquired"), 2),
  Adjustment = factor(rep(c("+ WHO Region\n+ year",
                             "+ All covariates"), each = 2),
                      levels = c("+ WHO Region\n+ year",
                                 "+ All covariates")),
  R2 = c(layer2_perm$perm_fg$naive["T_30d", "R2"],
         layer2_perm$perm_acq$naive["T_30d", "R2"],
         layer2_perm$perm_fg$full["T_30d", "R2"],
         layer2_perm$perm_acq$full["T_30d", "R2"])
)

fig1b <- ggplot(step_data, aes(x = Adjustment, y = R2, fill = Resistome)) +
  geom_col(position = "dodge", width = 0.65) +
  geom_text(aes(label = sprintf("%.2f%%", 100*R2)),
            position = position_dodge(width = 0.65), vjust = -0.4, size = 2.2) +
  scale_fill_manual(values = pal_resistome) +
  scale_y_continuous(labels = percent_format(accuracy = 0.01),
                     expand = expansion(mult = c(0, 0.15))) +
  labs(tag = "B", y = "Temperature R\u00B2\n(PERMANOVA)", x = NULL, fill = NULL) +
  theme(legend.position.inside = c(0.8, 0.85),
        legend.background = element_rect(fill = "white", colour = NA))

layer3_med <- readRDS(here("Results", "layer3_mediation.rds"))
med_adj <- tibble(
  resistome = c("Functional", "Acquired"),
  r2_direct = c(layer3_med$med_fg$r2_cprime, layer3_med$med_acq$r2_cprime),
  r2_indirect = c(layer3_med$med_fg$indirect, layer3_med$med_acq$indirect)
)

med_long <- med_adj %>%
  select(resistome, r2_direct, r2_indirect) %>%
  pivot_longer(-resistome, names_to = "path", values_to = "R2") %>%
  mutate(path = factor(path,
                       levels = c("r2_indirect", "r2_direct"),
                       labels = c("Indirect\n(via bacteriome)", "Direct")))

pct_labels <- med_adj %>%
  mutate(pct_med = 100 * r2_indirect / (r2_direct + r2_indirect),
         pct = sprintf("%.0f%% mediated", pct_med),
         total = r2_direct + r2_indirect)

fig1c <- ggplot(med_long, aes(x = resistome, y = R2, fill = path)) +
  geom_col(position = "stack", width = 0.45) +
  scale_fill_manual(values = c("Direct" = "#CC6677",
                                "Indirect\n(via bacteriome)" = "#88CCEE")) +
  scale_y_continuous(labels = percent_format(accuracy = 0.01),
                     expand = expansion(mult = c(0, 0.25))) +
  geom_text(aes(label = sprintf("%.2f%%", 100*R2)),
            position = position_stack(vjust = 0.5), size = 2.2) +
  labs(tag = "C", y = "Temperature R\u00B2\n(adjusted, mediation)",
       x = NULL, fill = NULL) +
  theme(legend.position.inside = c(0.65, 0.8),
        legend.background = element_rect(fill = "white", colour = NA))

motu <- readRDS(here("Datasets", "motu_species_cache.rds"))
motu_ids <- as.character(as.integer(motu$genepid))
motu_mat <- as.matrix(motu[, -1])
total_reads <- rowSums(motu_mat)
motu_ra <- motu_mat / total_reads
all_cols <- names(motu)[-1]

merge_df <- tibble(genepid = motu_ids) %>%
  inner_join(analysis %>% mutate(genepid = as.character(genepid)) %>%
               select(genepid, T_30d), by = "genepid") %>%
  filter(!is.na(T_30d))

idx <- match(merge_df$genepid, motu_ids)
motu_ra_matched <- motu_ra[idx, ]

get_ra <- function(pattern, use_fixed = FALSE) {
  if (use_fixed) {
    cols <- all_cols[str_detect(all_cols, fixed(pattern))]
  } else {
    cols <- all_cols[str_detect(all_cols, pattern)]
  }
  if (length(cols) == 0) return(rep(0, nrow(motu_ra_matched)))
  rowSums(motu_ra_matched[, cols, drop = FALSE])
}

sp_cols <- c("#0072B2", "#009E73", "#D55E00", "#CC79A7")

make_genus_panel <- function(genus_display, genus_pattern, species_list, tag) {
  genus_ra <- get_ra(genus_pattern)

  plot_df <- tibble(
    taxon = paste0(genus_display, " (all spp.)"),
    temp = merge_df$T_30d,
    ra = genus_ra
  )

  for (i in seq_along(species_list)) {
    sp_ra <- get_ra(species_list[i], use_fixed = TRUE)
    plot_df <- bind_rows(plot_df, tibble(
      taxon = names(species_list)[i],
      temp = merge_df$T_30d,
      ra = sp_ra
    ))
  }

  if (length(species_list) == 1) {
    genus_ra_sum <- sum(genus_ra)
    sp_ra_sum <- sum(get_ra(species_list[1], use_fixed = TRUE))
    if (abs(genus_ra_sum - sp_ra_sum) / max(genus_ra_sum, 1e-10) < 0.01) {
      plot_df <- plot_df %>% filter(taxon != paste0(genus_display, " (all spp.)"))
    }
  }

  adj_sp <- tryCatch(
    read_csv(here("Results", "adjusted_species_temp.csv"), show_col_types = FALSE),
    error = function(e) NULL
  )

  stats <- plot_df %>%
    group_by(taxon) %>%
    summarise(
      rho_raw = cor(ra, temp, method = "spearman"),
      p_raw = suppressWarnings(cor.test(ra, temp, method = "spearman")$p.value),
      .groups = "drop"
    )

  if (!is.null(adj_sp)) {
    stats <- stats %>%
      left_join(
        adj_sp %>% select(display, rho_adj, p_adj),
        by = c("taxon" = "display")
      ) %>%
      left_join(
        adj_sp %>% filter(str_detect(display, "all spp")) %>%
          select(display, rho_adj2 = rho_adj, p_adj2 = p_adj),
        by = c("taxon" = "display")
      ) %>%
      mutate(
        rho = coalesce(rho_adj, rho_adj2, rho_raw),
        p = coalesce(p_adj, p_adj2, p_raw)
      ) %>%
      select(taxon, rho, p)
  } else {
    stats <- stats %>% rename(rho = rho_raw, p = p_raw)
  }

  stats <- stats %>%
    mutate(
      p_label = ifelse(p < 0.001, "p<0.001", sprintf("p=%.3f", p)),
      legend_label = paste0(taxon, "\n\u03C1=", sprintf("%+.2f", rho), ", ", p_label)
    )

  plot_df <- plot_df %>%
    left_join(stats %>% select(taxon, legend_label), by = "taxon")

  taxa <- unique(plot_df$legend_label)
  is_genus <- str_detect(taxa, "all spp")
  colour_map <- setNames(rep("grey55", length(taxa)), taxa)
  sp_idx <- which(!is_genus)
  for (j in seq_along(sp_idx)) {
    colour_map[taxa[sp_idx[j]]] <- sp_cols[min(j, length(sp_cols))]
  }

  plot_df_nonzero <- plot_df %>% filter(ra > 0)

  ggplot(plot_df_nonzero, aes(x = temp, y = log10(ra), colour = legend_label)) +
    geom_point(alpha = 0.03, size = 0.2) +
    geom_smooth(method = "lm", se = TRUE, linewidth = 0.8, alpha = 0.2) +
    scale_colour_manual(values = colour_map, name = NULL) +
    guides(colour = guide_legend(override.aes = list(alpha = 1, size = 0, linewidth = 1.2))) +
    labs(tag = tag,
         title = genus_display,
         x = "30-day mean temperature (\u00B0C)",
         y = "Relative abundance (log10)") +
    theme(legend.key = element_rect(fill = NA, colour = NA),
          legend.position = "inside",
          legend.position.inside = c(0.98, 0.02),
          legend.justification = c(1, 0),
          legend.background = element_rect(fill = alpha("white", 0.92), colour = NA),
          legend.key.size = unit(0.3, "cm"),
          legend.text = element_text(size = 6.5),
          legend.spacing.y = unit(0.1, "cm"),
          plot.title = element_text(face = "bold.italic", size = 8))
}

p_a <- make_genus_panel("Acinetobacter", "^Acinetobacter",
  c("A. baumannii" = "Acinetobacter baumannii",
    "A. johnsonii" = "Acinetobacter johnsonii",
    "A. lwoffii" = "Acinetobacter lwoffii"), "A")

p_b <- make_genus_panel("Klebsiella", "^Klebsiella",
  c("K. pneumoniae" = "Klebsiella pneumoniae"), "B")

p_c <- make_genus_panel("Escherichia", "^Escherichia",
  c("E. coli" = "Escherichia coli"), "C")

p_d <- make_genus_panel("Enterococcus", "^Enterococcus",
  c("E. faecium" = "Enterococcus faecium",
    "E. faecalis" = "Enterococcus faecalis"), "D")

p_e <- make_genus_panel("Pseudomonas", "^Pseudomonas",
  c("P. aeruginosa" = "Pseudomonas aeruginosa"), "E")

p_f <- make_genus_panel("Staphylococcus", "^Staphylococcus",
  c("S. aureus" = "Staphylococcus aureus"), "F")

p_g <- make_genus_panel("Salmonella", "^Salmonella",
  c("S. enterica" = "Salmonella enterica"), "G")

p_h <- make_genus_panel("Clostridioides", "^Clostridioides",
  c("C. difficile" = "Clostridioides difficile"), "H")

strip_y <- theme(axis.title.y = element_blank())
strip_x <- theme(axis.title.x = element_blank())

fig2 <- ((p_a + strip_x) | (p_b + strip_x + strip_y) | (p_c + strip_x + strip_y) | (p_d + strip_x + strip_y)) /
        (p_e | (p_f + strip_y) | (p_g + strip_y) | (p_h + strip_y))

ggsave(here("Figures", "Fig2_species_pathogens.pdf"), fig2,
       width = 12, height = 7)

adj_arg <- read_csv(here("Results", "adjusted_arg_temp.csv"),
                    show_col_types = FALSE)

arg_data <- adj_arg %>%
  mutate(
    rho = rho_adj,
    p = p_adj,
    sig = p < 0.05,
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
      arg_family == "erm" ~ "erm (macrolide)",
      TRUE ~ arg_family
    ),
    p_label = ifelse(p < 0.001, "p<0.001",
                     ifelse(p < 0.01, sprintf("p=%.3f", p),
                            sprintf("p=%.2f", p))),
    point_col = ifelse(sig, "#2C3E50", "grey65")
  )

fig3a <- ggplot(arg_data, aes(x = rho, y = reorder(label, rho))) +
  geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_segment(aes(xend = 0, yend = reorder(label, rho)),
               colour = "grey75", linewidth = 0.4) +
  geom_point(aes(fill = point_col), size = 2.8,
             shape = 21, colour = "black", stroke = 0.3) +
  geom_text(aes(label = p_label, x = rho + ifelse(rho > 0, 0.03, -0.03)),
            hjust = ifelse(arg_data$rho > 0, 0, 1),
            size = 2, colour = "grey30") +
  scale_fill_identity() +
  scale_x_continuous(limits = c(-0.25, 0.35), breaks = seq(-0.2, 0.3, 0.1)) +
  labs(tag = "A",
       x = expression(paste("Spearman ", rho, " with 30-day mean temperature")),
       y = "Resistance gene family")

layer3_data <- readRDS(here("Results", "layer3_mediation.rds"))

pcoa_df <- layer3_data$bact_pcs %>%
  mutate(genepid = as.character(genepid)) %>%
  left_join(analysis %>% mutate(genepid = as.character(genepid)) %>%
              select(genepid, T_30d), by = "genepid") %>%
  filter(!is.na(T_30d)) %>%
  mutate(temp_group = cut(T_30d,
                          breaks = c(-Inf, 10, 20, Inf),
                          labels = c("<10\u00B0C", "10-20\u00b0C", ">20\u00B0C")))

eig_pct <- round(100 * layer3_data$pcoa_eig[1:2] /
                    sum(layer3_data$pcoa_eig[layer3_data$pcoa_eig > 0]), 1)

fig3b <- ggplot(pcoa_df, aes(x = BPC1, y = BPC2)) +
  geom_point(aes(fill = T_30d),
             alpha = 0.5, size = 1.0, shape = 21, stroke = 0.1, colour = "grey60") +
  stat_density_2d(aes(colour = temp_group), linewidth = 0.6, bins = 4) +
  scale_fill_gradient2(low = "#4575B4", mid = "#FFFFBF", high = "#D73027",
                       midpoint = 15, name = "Temp (\u00B0C)",
                       guide = guide_colourbar(order = 1)) +
  scale_colour_manual(values = c("<10\u00B0C" = "#4575B4",
                                  "10-20\u00b0C" = "#D4A017",
                                  ">20\u00B0C" = "#D73027"),
                      name = "Density",
                      guide = guide_legend(order = 2)) +
  labs(tag = "A", title = "Bacteriome",
       x = sprintf("PCoA1 (%.1f%%)", eig_pct[1]),
       y = sprintf("PCoA2 (%.1f%%)", eig_pct[2])) +
  theme(legend.position = "right",
        legend.key.height = unit(0.4, "cm"), legend.key.width = unit(0.25, "cm"),
        plot.title = element_text(face = "bold", size = 8))

resistome <- readRDS(here("Datasets", "resistome_matrices.rds"))
fg_clr <- resistome$fg_cluster_clr %>%
  mutate(genepid = as.character(genepid))

shared_ids <- intersect(fg_clr$genepid, pcoa_df$genepid)
fg_sub <- fg_clr %>% filter(genepid %in% shared_ids) %>% arrange(genepid)
clim_sub <- pcoa_df %>% filter(genepid %in% shared_ids) %>% arrange(genepid)

fg_mat <- fg_sub %>% select(-genepid) %>% as.matrix()
fg_dist <- dist(fg_mat)
fg_pcoa <- cmdscale(fg_dist, k = 2, eig = TRUE)
fg_eig_pct <- round(100 * fg_pcoa$eig[1:2] /
                       sum(fg_pcoa$eig[fg_pcoa$eig > 0]), 1)

res_pcoa_df <- tibble(
  RPC1 = fg_pcoa$points[, 1],
  RPC2 = fg_pcoa$points[, 2],
  BIO1 = clim_sub$T_30d,
  temp_group = clim_sub$temp_group
)

fig3c <- ggplot(res_pcoa_df, aes(x = RPC1, y = RPC2)) +
  geom_point(aes(fill = BIO1),
             alpha = 0.5, size = 1.0, shape = 21, stroke = 0.1, colour = "grey60") +
  stat_density_2d(aes(colour = temp_group), linewidth = 0.6, bins = 4) +
  scale_fill_gradient2(low = "#4575B4", mid = "#FFFFBF", high = "#D73027",
                       midpoint = 15, name = "Temp (\u00B0C)",
                       guide = guide_colourbar(order = 1)) +
  scale_colour_manual(values = c("<10\u00B0C" = "#4575B4",
                                  "10-20\u00b0C" = "#D4A017",
                                  ">20\u00B0C" = "#D73027"),
                      name = "Density",
                      guide = guide_legend(order = 2)) +
  labs(tag = "B", title = "Resistome",
       x = sprintf("PCoA1 (%.1f%%)", fg_eig_pct[1]),
       y = sprintf("PCoA2 (%.1f%%)", fg_eig_pct[2])) +
  theme(legend.position = "right",
        legend.key.height = unit(0.4, "cm"), legend.key.width = unit(0.25, "cm"),
        plot.title = element_text(face = "bold", size = 8))

ggsave(here("Figures", "SFig_adjusted_args_lollipop.pdf"), fig3a,
       width = 5, height = 4)

sfig_pcoa <- fig3b | fig3c
ggsave(here("Figures", "SFig_pcoa_bacteriome_resistome.pdf"), sfig_pcoa,
       width = 7.5, height = 4)

fig1 <- fig1a / (fig1b | fig1c) +
  plot_layout(heights = c(1, 0.8))
ggsave(here("Figures", "Fig1_temperature_resistome.pdf"), fig1,
       width = 7.5, height = 6)
ggsave(here("Figures", "Fig1_temperature_resistome.png"), fig1,
       width = 7.5, height = 6, dpi = 300)

ggsave(here("Figures", "Fig2_species_pathogens.png"), fig2,
       width = 12, height = 7, dpi = 300)
ggsave(here("Figures", "SFig_adjusted_args_lollipop.png"), fig3a,
       width = 5, height = 4, dpi = 300)
ggsave(here("Figures", "SFig_pcoa_bacteriome_resistome.png"), sfig_pcoa,
       width = 7.5, height = 4, dpi = 300)
