# =============================================================================
# 12_reviewer_supp_outputs.R
# Builds the NEW supplementary items for the reviewer follow-up, numbered to
# continue the current supplementary (existing: Figs S1-S4, Tables S1-S7):
#   Figure S5 : Species-level within-city temperature associations (longitudinal)
#   Figure S6 : Bacteriome composition varies by location (Q4)
#   Table  S8 : Cross-sectional vs longitudinal species-level correlations
# Depends on: Results/reviewer_q3_species_long.rds, Results/reviewer_q4_location.rds
# =============================================================================

library(tidyverse)
library(patchwork)
library(scales)
library(vegan)
library(showtext)
showtext_auto()
font_add_google("Noto Sans", "notosans")
showtext_opts(dpi = 300)

base_dir <- "/Users/bxc331/Desktop/AMR_open_project"
data_dir <- file.path(base_dir, "Datasets/Becsei_et_al")
fig_dir  <- file.path(base_dir, "supplementary_figures_export")
dir.create(fig_dir, showWarnings = FALSE)

theme_lph <- theme_bw(base_size = 8, base_family = "notosans") +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linewidth = 0.25, colour = "grey92"),
        axis.text = element_text(size = 6.5),
        axis.title = element_text(size = 7.5),
        plot.tag = element_text(face = "bold", size = 10),
        plot.title = element_text(face = "bold.italic", size = 8),
        legend.position = "bottom")
theme_set(theme_lph)

sp_cols <- c("#0072B2", "#009E73", "#D55E00", "#CC79A7")

q3 <- readRDS(file.path(base_dir, "Results", "reviewer_q3_species_long.rds"))
q4 <- readRDS(file.path(base_dir, "Results", "reviewer_q4_location.rds"))

# =============================================================================
# TABLE S8 : cross-sectional (adjusted) vs longitudinal (within-city) species rho
# =============================================================================
cs <- q3$cs %>% transmute(display, cs_rho_adj = rho_adj, cs_p_adj = p_adj)

# map long species (full names) to display abbreviations used cross-sectionally
disp_map <- c(
  "Acinetobacter baumannii"    = "A. baumannii",
  "Klebsiella pneumoniae"      = "K. pneumoniae",
  "Escherichia coli"           = "E. coli",
  "Enterobacter cloacae"       = "E. cloacae",
  "Mycobacterium tuberculosis" = "M. tuberculosis",
  "Enterococcus faecium"       = "E. faecium",
  "Enterococcus faecalis"      = "E. faecalis",
  "Staphylococcus aureus"      = "S. aureus",
  "Pseudomonas aeruginosa"     = "P. aeruginosa",
  "Salmonella enterica"        = "S. enterica",
  "Neisseria gonorrhoeae"      = "N. gonorrhoeae",
  "Streptococcus pneumoniae"   = "S. pneumoniae",
  "Streptococcus pyogenes"     = "S. pyogenes",
  "Haemophilus influenzae"     = "H. influenzae")

tableS8 <- q3$q3 %>%
  mutate(display = disp_map[species]) %>%
  left_join(cs, by = "display") %>%
  transmute(Species = display,
            cs_rho = cs_rho_adj, cs_p = cs_p_adj,
            lo_rho = long_rho, lo_p = long_p,
            `Cross-sectional adjusted rho` =
              ifelse(is.na(cs_rho_adj), "—", sprintf("%+.2f", cs_rho_adj)),
            `Cross-sectional p` =
              ifelse(is.na(cs_p_adj), "—", signif(cs_p_adj, 2)),
            `Longitudinal within-city rho` =
              ifelse(is.na(long_rho), "n.d.", sprintf("%+.2f", long_rho)),
            `Longitudinal p` =
              ifelse(is.na(long_p), "n.d.", signif(long_p, 2)),
            Concordant = case_when(
              is.na(cs_rho_adj) | is.na(long_rho) ~ "—",
              sign(cs_rho_adj) == sign(long_rho) ~ "Yes",
              TRUE ~ "No"))

write.csv(tableS8 %>% select(-cs_rho, -cs_p, -lo_rho, -lo_p),
          file.path(base_dir, "Results", "tableS8_species_longitudinal.csv"),
          row.names = FALSE)
cat("=== Table S8 ===\n")
print(as.data.frame(tableS8 %>% select(-cs_rho, -cs_p, -lo_rho, -lo_p)))

# =============================================================================
# FIGURE S5 : species-level within-city temperature associations
# Compact lollipop: cross-sectional adjusted rho vs longitudinal within-city rho
# =============================================================================
plot_df <- tableS8 %>%
  transmute(Species, cs = cs_rho, lo = lo_rho) %>%
  filter(!is.na(lo), !is.na(cs)) %>%   # paired species only
  arrange(cs) %>%
  mutate(Species = factor(Species, levels = Species))

figS5 <- plot_df %>%
  pivot_longer(c(cs, lo), names_to = "dataset", values_to = "rho") %>%
  mutate(dataset = recode(dataset,
                          cs = "Cross-sectional (adjusted)",
                          lo = "Longitudinal (within-city)")) %>%
  ggplot(aes(x = rho, y = Species, colour = dataset)) +
  geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_line(aes(group = Species), colour = "grey75", linewidth = 0.4) +
  geom_point(size = 2) +
  scale_colour_manual(values = c("Cross-sectional (adjusted)" = "#4A7298",
                                 "Longitudinal (within-city)" = "#D55E00"),
                      name = NULL) +
  labs(x = expression("Spearman "*rho*" with 30-day mean temperature"),
       y = NULL) +
  theme(axis.text.y = element_text(face = "italic"))

ggsave(file.path(fig_dir, "FigureS5_species_longitudinal.pdf"), figS5,
       width = 6.5, height = 4.5)
ggsave(file.path(fig_dir, "FigureS5_species_longitudinal.png"), figS5,
       width = 6.5, height = 4.5, dpi = 300)
cat("Saved Figure S5\n")

# =============================================================================
# FIGURE S6 : bacteriome composition varies by location
# Panel A: cross-sectional PCoA coloured by WHO region
# Panel B: 5-city PCoA coloured by city (shows between-city separation)
# =============================================================================
# --- Panel B from 5-city data (already computed in q4) ---
bray5 <- q4$bray5
cf    <- q4$cf
pco5  <- cmdscale(as.dist(bray5), k = 2, eig = TRUE)
eig5  <- round(100 * pco5$eig[1:2] / sum(pco5$eig[pco5$eig > 0]), 1)
pcoB <- tibble(PCo1 = pco5$points[, 1], PCo2 = pco5$points[, 2], City = cf)

city_cols <- c(Copenhagen = "#0072B2", Rotterdam = "#009E73",
               Budapest = "#D55E00", Bologna = "#CC79A7", Rome = "#E69F00")

pB <- ggplot(pcoB, aes(PCo1, PCo2, colour = City)) +
  geom_point(size = 1.4, alpha = 0.8) +
  stat_ellipse(level = 0.68, linewidth = 0.4) +
  scale_colour_manual(values = city_cols, name = NULL) +
  labs(tag = "B", title = "Five-city longitudinal bacteriome",
       x = sprintf("PCo1 (%.1f%%)", eig5[1]),
       y = sprintf("PCo2 (%.1f%%)", eig5[2]))

# --- Panel A cross-sectional PCoA coloured by WHO region ---
# recompute lightweight from q4-stored cross-sectional objects if present,
# else load analysis-ready + genus cache
analysis <- readRDS(file.path(base_dir, "Datasets", "analysis_ready.rds"))
genus_raw <- readRDS(file.path(base_dir, "Datasets", "genus_counts_cache.rds")) %>%
  mutate(genepid = as.character(as.integer(genepid)))
resistome <- readRDS(file.path(base_dir, "Datasets", "resistome_matrices.rds"))
genera <- setdiff(names(genus_raw), "genepid")
prev <- colSums(genus_raw[, genera] > 0) / nrow(genus_raw)
keep_genera <- names(prev[prev >= 0.05])
genus_filt <- genus_raw %>% select(genepid, all_of(keep_genera))
shared_ids <- Reduce(intersect, list(genus_filt$genepid,
                                     resistome$fg_cluster_clr$genepid, analysis$genepid))
genus_sub <- genus_filt %>% filter(genepid %in% shared_ids) %>% arrange(genepid)
clim_sub  <- analysis %>% filter(genepid %in% shared_ids) %>% arrange(genepid)
genus_ra  <- genus_sub %>% mutate(across(-genepid, ~ . / rowSums(across(-genepid))))
keep <- complete.cases(clim_sub[, c("T_30d","Region","year","gdp_pcap_ppp","sanitation",
                                     "health_exp_gdp","oop_health_exp","immunization_dpt",
                                     "animal_amc_mgkg","pop_density")])
genus_ra <- genus_ra[keep, ]; clim_sub <- clim_sub[keep, ]
brayCS <- vegdist(genus_ra %>% select(-genepid), method = "bray")
pcoCS <- cmdscale(brayCS, k = 2, eig = TRUE)
eigCS <- round(100 * pcoCS$eig[1:2] / sum(pcoCS$eig[pcoCS$eig > 0]), 1)
pcoA <- tibble(PCo1 = pcoCS$points[, 1], PCo2 = pcoCS$points[, 2],
               Region = clim_sub$WHO_region)

pA <- ggplot(pcoA, aes(PCo1, PCo2, colour = Region)) +
  geom_point(size = 0.8, alpha = 0.6) +
  scale_colour_brewer(palette = "Dark2", name = NULL) +
  labs(tag = "A", title = "Global cross-sectional bacteriome",
       x = sprintf("PCo1 (%.1f%%)", eigCS[1]),
       y = sprintf("PCo2 (%.1f%%)", eigCS[2])) +
  guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1)))

figS6 <- pA + pB + plot_layout(widths = c(1, 1))
ggsave(file.path(fig_dir, "FigureS6_bacteriome_by_location.pdf"), figS6,
       width = 9, height = 4.2)
ggsave(file.path(fig_dir, "FigureS6_bacteriome_by_location.png"), figS6,
       width = 9, height = 4.2, dpi = 300)
cat("Saved Figure S6\n")

cat("\nDONE building supp outputs.\n")
