# PCoA variance scree plot and mediation sensitivity to axis count (k = 3, 5, 10)

library(tidyverse)
library(vegan)
library(patchwork)
library(showtext)
library(here)
showtext_auto()
font_add_google("Noto Sans", "notosans")

theme_lph <- theme_bw(base_size = 8, base_family = "notosans") +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.25, colour = "grey92"),
    legend.position  = "bottom",
    legend.key.size  = unit(0.3, "cm"),
    legend.text      = element_text(size = 6),
    legend.title     = element_text(size = 7, face = "bold"),
    plot.title       = element_text(face = "bold.italic", size = 8),
    axis.text        = element_text(size = 6.5),
    axis.title       = element_text(size = 7.5),
    plot.tag         = element_text(face = "bold", size = 10)
  )
theme_set(theme_lph)

med_data <- readRDS(here("Results", "layer3_mediation.rds"))
eig <- med_data$pcoa_eig
pos_eig <- eig[eig > 0]
pct_var <- 100 * pos_eig / sum(pos_eig)

n_show <- 10
scree_df <- tibble(
  axis     = 1:n_show,
  pct      = pct_var[1:n_show],
  cum_pct  = cumsum(pct_var[1:n_show])
)

cum_col <- "#2E6E4E"

p_scree <- ggplot(scree_df, aes(x = axis)) +
  geom_col(aes(y = pct), fill = "#4A7298", width = 0.7) +
  geom_line(aes(y = cum_pct / max(cum_pct) * max(pct)),
            colour = cum_col, linewidth = 0.7) +
  geom_point(aes(y = cum_pct / max(cum_pct) * max(pct)),
             colour = cum_col, size = 1.8) +
  geom_vline(xintercept = 5.5, linetype = "dashed", colour = "grey40",
             linewidth = 0.4) +
  annotate("text", x = 5.7, y = max(scree_df$pct) * 0.95,
           label = "k = 5", hjust = 0, size = 3, family = "notosans",
           colour = "grey30") +
  scale_x_continuous(breaks = 1:n_show, expand = expansion(add = 0.5)) +
  scale_y_continuous(
    name = "Variance explained (%)",
    sec.axis = sec_axis(
      ~ . / max(scree_df$pct) * max(scree_df$cum_pct),
      name = "Cumulative variance (%)"
    ),
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(x = "PCoA axis") +
  theme(axis.title       = element_text(size = 10),
        axis.text        = element_text(size = 9),
        axis.title.y.right = element_text(colour = cum_col, size = 10),
        axis.text.y.right  = element_text(colour = cum_col, size = 9))

ggsave(here("Figures", "eFigure_pcoa_variance.png"), p_scree,
       width = 5.5, height = 3.5, dpi = 300, bg = "white")
ggsave(here("Figures", "eFigure_pcoa_variance.pdf"), p_scree,
       width = 5.5, height = 3.5, device = cairo_pdf, bg = "white")

analysis  <- readRDS(here("Datasets", "analysis_ready.rds"))
resistome <- readRDS(here("Datasets", "resistome_matrices.rds"))
genus_raw <- readRDS(here("Datasets", "genus_counts_cache.rds"))

genus_raw <- genus_raw %>% mutate(genepid = as.character(as.integer(genepid)))

genera <- setdiff(names(genus_raw), "genepid")
prev <- colSums(genus_raw[, genera] > 0) / nrow(genus_raw)
keep_genera <- names(prev[prev >= 0.05])
genus_filt <- genus_raw %>% select(genepid, all_of(keep_genera))

shared_ids <- Reduce(intersect, list(
  genus_filt$genepid,
  resistome$fg_cluster_clr$genepid,
  analysis$genepid
))

genus_sub <- genus_filt %>% filter(genepid %in% shared_ids) %>% arrange(genepid)
clim_sub  <- analysis   %>% filter(genepid %in% shared_ids) %>% arrange(genepid)

genus_ra <- genus_sub %>%
  mutate(across(-genepid, ~ . / rowSums(across(-genepid))))

keep <- complete.cases(clim_sub[, c("T_30d", "Region", "year",
                                     "gdp_pcap_ppp", "sanitation", "health_exp_gdp",
                                     "oop_health_exp", "immunization_dpt",
                                     "animal_amc_mgkg", "pop_density")])
genus_ra <- genus_ra[keep, ]
clim_sub <- clim_sub[keep, ]

bray_dist <- vegdist(genus_ra %>% select(-genepid), method = "bray")
pcoa_full <- cmdscale(bray_dist, k = 10, eig = TRUE)
bact_pcs_all <- as.data.frame(pcoa_full$points)
names(bact_pcs_all) <- paste0("BPC", 1:10)
bact_pcs_all$genepid <- as.character(genus_ra$genepid)

prep_resistome <- function(clr_mat, bact_pcs_df, clim_df, label) {
  clr_mat <- clr_mat %>% mutate(genepid = as.character(genepid))
  clim_df <- clim_df %>% mutate(genepid = as.character(genepid))
  bact_pcs_df <- bact_pcs_df %>% mutate(genepid = as.character(genepid))

  shared <- Reduce(intersect, list(clr_mat$genepid, bact_pcs_df$genepid, clim_df$genepid))
  clr_s  <- clr_mat     %>% filter(genepid %in% shared) %>% arrange(genepid)
  bpc_s  <- bact_pcs_df %>% filter(genepid %in% shared) %>% arrange(genepid)
  clim_s <- clim_df     %>% filter(genepid %in% shared) %>% arrange(genepid)

  cc <- complete.cases(clim_s[, c("T_30d", "Region", "year",
                                    "gdp_pcap_ppp", "sanitation", "health_exp_gdp",
                                    "oop_health_exp", "immunization_dpt",
                                    "animal_amc_mgkg", "pop_density")])
  clr_s <- clr_s[cc, ]; bpc_s <- bpc_s[cc, ]; clim_s <- clim_s[cc, ]

  dist_mat <- dist(clr_s %>% select(-genepid))

  med_df <- clim_s %>%
    left_join(bpc_s, by = "genepid") %>%
    mutate(year_f = factor(year))

  list(dist_mat = dist_mat, med_df = med_df, n = nrow(clr_s), label = label)
}

prep_fg  <- prep_resistome(resistome$fg_cluster_clr, bact_pcs_all, analysis, "Functional")
prep_acq <- prep_resistome(resistome$acq_cluster_clr, bact_pcs_all, analysis, "Acquired")

run_mediation_k <- function(prep, k) {
  bpc_cols <- paste0("BPC", 1:k)

  perm_c <- adonis2(prep$dist_mat ~ T_30d + Region + year_f,
                    data = prep$med_df, permutations = 99, by = "margin")

  bpc_formula <- as.formula(
    paste("prep$dist_mat ~ T_30d +", paste(bpc_cols, collapse = " + "), "+ Region + year_f")
  )
  perm_cp <- adonis2(bpc_formula, data = prep$med_df, permutations = 99, by = "margin")

  r2_c  <- perm_c["T_30d", "R2"]
  r2_cp <- perm_cp["T_30d", "R2"]
  indirect <- r2_c - r2_cp
  prop_med <- ifelse(r2_c > 0, indirect / r2_c * 100, NA)

  tibble(resistome = prep$label, k = k,
         r2_total = r2_c, r2_direct = r2_cp,
         r2_indirect = indirect, pct_mediated = prop_med)
}

k_values <- c(3, 5, 10)
results <- list()

for (kk in k_values) {
  results[[length(results) + 1]] <- run_mediation_k(prep_fg, kk)
  results[[length(results) + 1]] <- run_mediation_k(prep_acq, kk)
}

sens_df <- bind_rows(results)

saveRDS(sens_df, here("Results", "pcoa_mediation_sensitivity.rds"))

library(flextable)
library(officer)

tab_df <- sens_df %>%
  mutate(
    `PCoA axes (k)` = k,
    `Cumulative variance (%)` = sapply(k, function(kk) {
      sprintf("%.1f", sum(pct_var[1:kk]))
    }),
    `Total R2 (%)` = sprintf("%.2f", r2_total * 100),
    `Direct R2 (%)` = sprintf("%.2f", r2_direct * 100),
    `Mediated R2 (%)` = sprintf("%.2f", r2_indirect * 100),
    `Proportion mediated (%)` = sprintf("%.1f", pct_mediated),
    Resistome = resistome
  ) %>%
  select(Resistome, `PCoA axes (k)`, `Cumulative variance (%)`,
         `Total R2 (%)`, `Direct R2 (%)`, `Mediated R2 (%)`,
         `Proportion mediated (%)`)

ft <- flextable(tab_df) %>%
  merge_v(j = "Resistome") %>%
  valign(j = "Resistome", valign = "top") %>%
  bold(j = "Resistome") %>%
  bold(i = ~ `PCoA axes (k)` == 5) %>%
  add_footer_lines(paste0(
    "Sensitivity of the mediation analysis to the number of bacteriome ",
    "principal coordinates retained. Total R\u00B2: marginal PERMANOVA R\u00B2 ",
    "for 30-day mean temperature adjusting for WHO region and sampling year. ",
    "Direct R\u00B2: temperature R\u00B2 after additionally conditioning on k ",
    "bacteriome principal coordinates. Proportion mediated = (Total \u2212 Direct) / Total \u00D7 100. ",
    "The primary analysis (k=5, bold) is shown alongside 3 and 10 axes. ",
    "Cumulative variance refers to the proportion of bacteriome Bray\u2013Curtis ",
    "variation captured by the retained axes. 99 permutations; R\u00B2 is deterministic."
  )) %>%
  fontsize(size = 9, part = "footer") %>%
  autofit()

doc <- read_docx() %>%
  body_add_par(
    "Supplementary Table. Sensitivity of mediation results to the number of bacteriome principal coordinates.",
    style = "heading 2") %>%
  body_add_par("") %>%
  body_add_flextable(ft)

print(doc, target = here("Tables", "eTable_pcoa_mediation_sensitivity.docx"))
