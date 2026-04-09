# Layer 3: Species-level mediation analysis (temperature -> bacteriome -> resistome).
# Outputs: Results/layer3_genus_temp_correlations.csv, Results/layer3_mediation_summary.csv, Results/layer3_mediation.rds

library(here)
library(tidyverse)
library(vegan)
library(compositions)
library(broom)

analysis  <- readRDS(here("Datasets", "analysis_ready.rds"))
resistome <- readRDS(here("Datasets", "resistome_matrices.rds"))

genus_file <- here("Datasets", "genus_counts_cache.rds")

if (file.exists(genus_file)) {
  genus_raw <- readRDS(genus_file)
} else {
  zip_path <- here("Datasets", "Martiny_et_al", "motus_counts.zip")
  tmp_dir <- tempdir()
  unzip(zip_path, files = "motus_counts/motus_agg_pad/genus_motus_pad_agg_counts.csv",
        exdir = tmp_dir)
  genus_raw <- read_csv(file.path(tmp_dir, "motus_counts/motus_agg_pad/genus_motus_pad_agg_counts.csv"),
                        show_col_types = FALSE)
  saveRDS(genus_raw, genus_file)
}

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
clim_sub <- analysis %>% filter(genepid %in% shared_ids) %>% arrange(genepid)

genus_ra <- genus_sub %>%
  mutate(across(-genepid, ~ . / rowSums(across(-genepid))))

temp_cors <- map_dfr(keep_genera, function(g) {
  x <- genus_ra[[g]]
  y <- clim_sub$T_30d
  ok <- !is.na(x) & !is.na(y)
  if (sum(ok) < 30) return(NULL)
  ct <- cor.test(x[ok], y[ok], method = "spearman")
  tibble(genus = g, rho = ct$estimate, p_value = ct$p.value, n = sum(ok))
}) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH"),
         category = case_when(
           rho > 0.2 & p_adj < 0.05 ~ "thermophilic",
           rho < -0.2 & p_adj < 0.05 ~ "psychrophilic",
           TRUE ~ "neutral"
         )) %>%
  arrange(desc(abs(rho)))

genus_ra_mat <- genus_ra %>% filter(genepid %in% shared_ids) %>% arrange(genepid)
clim_match <- analysis %>% filter(genepid %in% shared_ids) %>% arrange(genepid)

keep <- complete.cases(clim_match[, c("T_30d", "Region", "year",
                                       "gdp_pcap_ppp", "sanitation", "health_exp_gdp",
                                       "oop_health_exp", "immunization_dpt",
                                       "animal_amc_mgkg", "pop_density")])
genus_ra_mat <- genus_ra_mat[keep, ]
clim_match <- clim_match[keep, ]

bray_dist <- vegdist(genus_ra_mat %>% select(-genepid), method = "bray")

clim_match <- clim_match %>% mutate(year_f = factor(year))

perm_bact <- adonis2(bray_dist ~ T_30d + Region + year_f,
                     data = clim_match, permutations = 999, by = "margin")

pcoa <- cmdscale(bray_dist, k = 5, eig = TRUE)
bact_pcs <- as.data.frame(pcoa$points)
names(bact_pcs) <- paste0("BPC", 1:5)
bact_pcs$genepid <- as.character(genus_ra_mat$genepid)

eig_pct <- round(100 * pcoa$eig[1:5] / sum(pcoa$eig[pcoa$eig > 0]), 1)

run_mediation <- function(clr_mat, bact_pcs_df, clim_df, label) {
  clr_mat <- clr_mat %>% mutate(genepid = as.character(genepid))
  clim_df <- clim_df %>% mutate(genepid = as.character(genepid))
  bact_pcs_df <- bact_pcs_df %>% mutate(genepid = as.character(genepid))
  shared <- Reduce(intersect, list(clr_mat$genepid, bact_pcs_df$genepid, clim_df$genepid))
  clr_sub <- clr_mat %>% filter(genepid %in% shared) %>% arrange(genepid)
  bpc_sub <- bact_pcs_df %>% filter(genepid %in% shared) %>% arrange(genepid)
  clim_sub <- clim_df %>% filter(genepid %in% shared) %>% arrange(genepid)

  keep <- complete.cases(clim_sub[, c("T_30d", "Region", "year",
                                       "gdp_pcap_ppp", "sanitation", "health_exp_gdp",
                                       "oop_health_exp", "immunization_dpt",
                                       "animal_amc_mgkg", "pop_density")])
  clr_sub <- clr_sub[keep, ]; bpc_sub <- bpc_sub[keep, ]; clim_sub <- clim_sub[keep, ]

  dist_mat <- dist(clr_sub %>% select(-genepid))

  med_df <- clim_sub %>%
    left_join(bpc_sub, by = "genepid") %>%
    mutate(year_f = factor(year))

  perm_c <- adonis2(dist_mat ~ T_30d + Region + year_f,
                    data = med_df, permutations = 999, by = "margin")

  perm_cprime <- adonis2(dist_mat ~ T_30d + BPC1 + BPC2 + BPC3 +
                           BPC4 + BPC5 + Region + year_f,
                         data = med_df, permutations = 999, by = "margin")

  r2_c <- perm_c["T_30d", "R2"]
  r2_cprime <- perm_cprime["T_30d", "R2"]
  indirect <- r2_c - r2_cprime
  prop_mediated <- ifelse(r2_c > 0, indirect / r2_c, NA)

  path_a <- map_dfr(paste0("BPC", 1:5), function(pc) {
    mod <- lm(reformulate("T_30d", pc), data = med_df)
    tid <- tidy(mod) %>% filter(term == "T_30d")
    tibble(pc = pc, beta = tid$estimate, p = tid$p.value,
           r2 = summary(mod)$r.squared)
  })

  list(
    perm_c = perm_c, perm_cprime = perm_cprime,
    r2_c = r2_c, r2_cprime = r2_cprime,
    indirect = indirect, prop_mediated = prop_mediated,
    path_a = path_a, n = nrow(clr_sub), label = label
  )
}

med_fg  <- run_mediation(resistome$fg_cluster_clr, bact_pcs, analysis, "FG")
med_acq <- run_mediation(resistome$acq_cluster_clr, bact_pcs, analysis, "Acquired")

mediation_summary <- tibble(
  resistome = c("FG", "Acquired"),
  r2_total = c(med_fg$r2_c, med_acq$r2_c),
  r2_direct = c(med_fg$r2_cprime, med_acq$r2_cprime),
  r2_indirect = c(med_fg$indirect, med_acq$indirect),
  pct_mediated = c(med_fg$prop_mediated, med_acq$prop_mediated) * 100
)

focus_genera <- c(
  # WHO BPPL 2024 — Critical
  "Acinetobacter", "Escherichia", "Klebsiella", "Enterobacter", "Mycobacterium",
  # WHO BPPL 2024 — High
  "Enterococcus", "Staphylococcus", "Pseudomonas", "Salmonella", "Neisseria",
  # WHO BPPL 2024 — Medium
  "Streptococcus", "Haemophilus",
  # Additional (foodborne / 2017 WHO)
  "Campylobacter", "Vibrio", "Clostridioides"
)

focus_results <- temp_cors %>%
  filter(genus %in% focus_genera | str_detect(genus, paste(focus_genera, collapse = "|")))

write_csv(temp_cors, here("Results", "layer3_genus_temp_correlations.csv"))
write_csv(mediation_summary, here("Results", "layer3_mediation_summary.csv"))
saveRDS(list(med_fg = med_fg, med_acq = med_acq,
             perm_bact = perm_bact, bact_pcs = bact_pcs,
             temp_cors = temp_cors, pcoa_eig = pcoa$eig),
        here("Results", "layer3_mediation.rds"))
