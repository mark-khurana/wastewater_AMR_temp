# =============================================================================
# 07c_adjusted_correlations.R — Covariate-adjusted species + ARG correlations
# =============================================================================
#
# For supplementary: repeat Fig 2 and Fig 3a correlations after residualizing
# on confounders (partial correlation approach).
#
# Method: For each species/ARG family, fit:
#   log(abundance) ~ log_gdp + sanitation + health_exp + oop_exp +
#                     immunization + log_animal_amc + Region
# Take residuals. Also residualize temperature on same confounders.
# Then Spearman correlation of residuals = partial correlation.
#
# Output: Results/adjusted_species_temp.csv, Results/adjusted_arg_temp.csv
#         Figures/SFig3_adjusted_species.pdf, Figures/SFig4_adjusted_args.pdf
# =============================================================================

library(tidyverse)
library(patchwork)

base_dir <- "/Users/bxc331/Desktop/AMR_open_project"
out_dir  <- file.path(base_dir, "Datasets")
res_dir  <- file.path(base_dir, "Results")
fig_dir  <- file.path(base_dir, "Figures")

cat("=== 07c_adjusted_correlations.R ===\n\n")

analysis <- readRDS(file.path(out_dir, "analysis_ready.rds")) %>%
  mutate(genepid = as.character(genepid),
         log_gdp = log(gdp_pcap_ppp),
         log_animal_amc = log(animal_amc_mgkg + 1),
         log_pop_density = log(pop_density + 1))

# Confounders to adjust for (including collection year)
conf_formula <- ~ log_gdp + sanitation + health_exp_gdp + oop_health_exp +
  immunization_dpt + log_animal_amc + Region + year_f

# Subset to complete cases, matching the main PERMANOVA sample
conf_vars <- c("T_30d", "log_gdp", "sanitation", "health_exp_gdp",
               "oop_health_exp", "immunization_dpt", "log_animal_amc", "log_pop_density",
               "Region", "year")
resistome <- readRDS(file.path(out_dir, "resistome_matrices.rds"))
resistome_ids <- as.character(resistome$fg_cluster_clr$genepid)
analysis_cc <- analysis %>%
  filter(complete.cases(across(all_of(conf_vars))),
         genepid %in% resistome_ids) %>%
  mutate(year_f = factor(year))
cat("Complete cases for adjustment:", nrow(analysis_cc), "\n\n")

# Residualize temperature on confounders (including year)
temp_mod <- lm(T_30d ~ log_gdp + sanitation + health_exp_gdp +
                 oop_health_exp + immunization_dpt + log_pop_density +
                 log_animal_amc + Region + year_f,
               data = analysis_cc)
analysis_cc$temp_resid <- residuals(temp_mod)

# =============================================================================
# A. ADJUSTED SPECIES CORRELATIONS
# =============================================================================
cat("=== A. Adjusted species correlations ===\n")

motu <- readRDS(file.path(out_dir, "motu_species_cache.rds"))
motu_ids <- as.character(as.integer(motu$genepid))
motu_mat <- as.matrix(motu[, -1])
total_reads <- rowSums(motu_mat)
motu_ra <- motu_mat / total_reads
all_cols <- names(motu)[-1]

# Match to analysis_cc
shared <- intersect(motu_ids, analysis_cc$genepid)
idx_motu <- match(shared, motu_ids)
idx_clim <- match(shared, analysis_cc$genepid)
motu_ra_sub <- motu_ra[idx_motu, ]
clim_sub <- analysis_cc[idx_clim, ]

cat("Matched samples:", nrow(clim_sub), "\n")

get_ra <- function(pattern, use_fixed = FALSE) {
  if (use_fixed) {
    cols <- all_cols[str_detect(all_cols, fixed(pattern))]
  } else {
    cols <- all_cols[str_detect(all_cols, pattern)]
  }
  if (length(cols) == 0) return(rep(0, nrow(motu_ra_sub)))
  rowSums(motu_ra_sub[, cols, drop = FALSE])
}

# Species to test
species_list <- list(
  list(genus = "Acinetobacter", display = "A. baumannii", pattern = "Acinetobacter baumannii"),
  list(genus = "Acinetobacter", display = "A. johnsonii", pattern = "Acinetobacter johnsonii"),
  list(genus = "Acinetobacter", display = "Acinetobacter (all spp.)", pattern = "^Acinetobacter"),
  list(genus = "Klebsiella", display = "K. pneumoniae", pattern = "Klebsiella pneumoniae"),
  list(genus = "Klebsiella", display = "Klebsiella (all spp.)", pattern = "^Klebsiella"),
  list(genus = "Escherichia", display = "E. coli", pattern = "Escherichia coli"),
  list(genus = "Enterococcus", display = "E. faecium", pattern = "Enterococcus faecium"),
  list(genus = "Enterococcus", display = "Enterococcus (all spp.)", pattern = "^Enterococcus"),
  list(genus = "Pseudomonas", display = "P. aeruginosa", pattern = "Pseudomonas aeruginosa"),
  list(genus = "Pseudomonas", display = "Pseudomonas (all spp.)", pattern = "^Pseudomonas"),
  list(genus = "Staphylococcus", display = "S. aureus", pattern = "Staphylococcus aureus"),
  list(genus = "Staphylococcus", display = "Staphylococcus (all spp.)", pattern = "^Staphylococcus"),
  list(genus = "Salmonella", display = "S. enterica", pattern = "Salmonella enterica"),
  list(genus = "Salmonella", display = "Salmonella (all spp.)", pattern = "^Salmonella"),
  list(genus = "Clostridioides", display = "C. difficile", pattern = "Clostridioides difficile"),
  list(genus = "Clostridioides", display = "Clostridioides (all spp.)", pattern = "^Clostridioides"),
  # WHO BPPL 2024 additions
  list(genus = "Mycobacterium", display = "M. tuberculosis", pattern = "Mycobacterium tuberculosis"),
  list(genus = "Mycobacterium", display = "Mycobacterium (all spp.)", pattern = "^Mycobacterium"),
  list(genus = "Neisseria", display = "N. gonorrhoeae", pattern = "Neisseria gonorrhoeae"),
  list(genus = "Neisseria", display = "Neisseria (all spp.)", pattern = "^Neisseria"),
  list(genus = "Haemophilus", display = "H. influenzae", pattern = "Haemophilus influenzae"),
  list(genus = "Haemophilus", display = "Haemophilus (all spp.)", pattern = "^Haemophilus"),
  list(genus = "Enterobacter", display = "E. cloacae", pattern = "Enterobacter cloacae"),
  list(genus = "Enterobacter", display = "Enterobacter (all spp.)", pattern = "^Enterobacter"),
  list(genus = "Streptococcus", display = "S. pneumoniae", pattern = "Streptococcus pneumoniae"),
  list(genus = "Streptococcus", display = "S. pyogenes", pattern = "Streptococcus pyogenes"),
  list(genus = "Streptococcus", display = "Streptococcus (all spp.)", pattern = "^Streptococcus")
)

adj_species <- map_dfr(species_list, function(sp) {
  use_fixed <- !str_starts(sp$pattern, "\\^")
  ra <- get_ra(sp$pattern, use_fixed = use_fixed)

  # Unadjusted
  ct_raw <- cor.test(ra, clim_sub$T_30d, method = "spearman")

  # Adjusted: residualize abundance on confounders, then correlate with temp residuals
  df <- tibble(ra = ra, temp_resid = clim_sub$temp_resid) %>%
    bind_cols(clim_sub %>% select(log_gdp, sanitation, health_exp_gdp,
                                   oop_health_exp, immunization_dpt,
                                   log_animal_amc, Region, year_f))
  # Only use non-zero for the abundance residualization
  if (sum(ra > 0) < 50) {
    return(tibble(genus = sp$genus, display = sp$display,
                  rho_raw = ct_raw$estimate, p_raw = ct_raw$p.value,
                  rho_adj = NA_real_, p_adj = NA_real_))
  }

  abund_mod <- lm(log(ra + 1e-8) ~ log_gdp + sanitation + health_exp_gdp +
                    oop_health_exp + immunization_dpt + log_animal_amc +
                    Region + year_f,
                  data = df)
  df$abund_resid <- residuals(abund_mod)

  ct_adj <- cor.test(df$abund_resid, df$temp_resid, method = "spearman")

  tibble(genus = sp$genus, display = sp$display,
         rho_raw = ct_raw$estimate, p_raw = ct_raw$p.value,
         rho_adj = ct_adj$estimate, p_adj = ct_adj$p.value)
})

cat("\n--- Species: raw vs adjusted correlations ---\n")
adj_species %>%
  mutate(across(c(rho_raw, rho_adj), ~ round(.x, 3)),
         across(c(p_raw, p_adj), ~ signif(.x, 3))) %>%
  print(n = 20)

write_csv(adj_species, file.path(res_dir, "adjusted_species_temp.csv"))

# =============================================================================
# B. ADJUSTED ARG FAMILY CORRELATIONS
# =============================================================================
cat("\n=== B. Adjusted ARG family correlations ===\n")

arg_raw <- read_csv(file.path(res_dir, "layer3c_arg_family_temp.csv"),
                    show_col_types = FALSE)

# Rebuild ARG family abundances for the cc subset
library(readxl)
counts <- readRDS(file.path(out_dir, "clean_ARG_counts.rds"))
sd2 <- read_excel(file.path(out_dir, "Martiny_et_al/supplementary/41467_2025_66070_MOESM5_ESM.xlsx"),
                  sheet = "SD2 ARG metadata")

gene_map <- sd2 %>%
  filter(str_detect(fa_name, "^resfinder")) %>%
  mutate(real_name = str_extract(fa_name, "resfinder\\|(.+)", group = 1) %>% str_to_lower()) %>%
  distinct(cluster_representative_98, real_name)

arg_families <- tribble(
  ~family, ~pattern,
  "dfr", "^dfr", "qnr", "^qnr", "tet", "^tet",
  "blaTEM", "^blatem", "aac", "^aac", "blaSHV", "^blashv",
  "blaCTX-M", "^blactx-m", "aph", "^aph", "ant", "^ant\\\\(",
  "blaOXA", "^blaoxa", "erm", "^erm"
)

# Simpler approach: use the family abundances from the raw results
# and just compute partial correlations
# Load the family abundances by reconstructing from counts
cluster_family <- gene_map %>%
  mutate(arg_family = NA_character_)

for (fam in c("dfr", "qnr", "tet", "blaTEM", "aac", "blaSHV", "blaCTX-M", "aph", "blaOXA", "erm")) {
  pat <- switch(fam,
    dfr = "^dfr", qnr = "^qnr", tet = "^tet", blaTEM = "^blatem",
    aac = "^aac", blaSHV = "^blashv", `blaCTX-M` = "^blactx-m",
    aph = "^aph", blaOXA = "^blaoxa", erm = "^erm")
  cluster_family <- cluster_family %>%
    mutate(arg_family = ifelse(is.na(arg_family) & str_detect(real_name, pat),
                                fam, arg_family))
}

# ant needs special pattern
cluster_family <- cluster_family %>%
  mutate(arg_family = ifelse(is.na(arg_family) & str_detect(real_name, "^ant\\("),
                              "ant", arg_family))

cluster_family <- cluster_family %>% filter(!is.na(arg_family)) %>%
  distinct(cluster_representative_98, arg_family)

total_abund <- counts %>%
  group_by(genepid) %>%
  summarise(total_arg = sum(fragmentCountAln_adj, na.rm = TRUE), .groups = "drop") %>%
  mutate(genepid = as.character(genepid))

family_abund <- counts %>%
  inner_join(cluster_family, by = "cluster_representative_98") %>%
  group_by(genepid, arg_family) %>%
  summarise(abundance = sum(fragmentCountAln_adj, na.rm = TRUE), .groups = "drop") %>%
  mutate(genepid = as.character(genepid))

# For each family, compute raw and adjusted correlation
adj_args <- map_dfr(unique(cluster_family$arg_family), function(fam) {
  fam_data <- family_abund %>% filter(arg_family == fam) %>%
    right_join(total_abund, by = "genepid") %>%
    mutate(abundance = replace_na(abundance, 0),
           ra = abundance / total_arg) %>%
    inner_join(analysis_cc %>% select(genepid, T_30d, temp_resid,
                                       log_gdp, sanitation, health_exp_gdp,
                                       oop_health_exp, immunization_dpt,
                                       log_animal_amc, Region, year_f),
               by = "genepid")

  if (nrow(fam_data) < 50) return(NULL)

  ct_raw <- cor.test(fam_data$ra, fam_data$T_30d, method = "spearman")

  abund_mod <- lm(ra ~ log_gdp + sanitation + health_exp_gdp + oop_health_exp +
                    immunization_dpt + log_animal_amc + Region + year_f,
                  data = fam_data)
  fam_data$abund_resid <- residuals(abund_mod)

  ct_adj <- cor.test(fam_data$abund_resid, fam_data$temp_resid, method = "spearman")

  tibble(arg_family = fam,
         rho_raw = ct_raw$estimate, p_raw = ct_raw$p.value,
         rho_adj = ct_adj$estimate, p_adj = ct_adj$p.value)
})

cat("\n--- ARG families: raw vs adjusted correlations ---\n")
adj_args %>%
  mutate(across(c(rho_raw, rho_adj), ~ round(.x, 3)),
         across(c(p_raw, p_adj), ~ signif(.x, 3))) %>%
  print(n = 15)

write_csv(adj_args, file.path(res_dir, "adjusted_arg_temp.csv"))

cat("\nSaved adjusted correlations to Results/\n")
