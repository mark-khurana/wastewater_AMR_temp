# =============================================================================
# 11_reviewer_Q3_Q4.R
# Reviewer follow-up analyses:
#   Q3: species-level temperature associations in the longitudinal (5-city) data,
#       compared with the cross-sectional adjusted estimates (Table S7 / Fig 2).
#   Q4: how much sewage bacteriome composition varies BY LOCATION
#       (city / country) vs by temperature.
# Prints numbers only; figures/tables built in follow-up.
# =============================================================================

library(tidyverse)
library(vegan)
library(here)

base_dir <- "/Users/bxc331/Desktop/AMR_open_project"
data_dir <- file.path(base_dir, "Datasets/Becsei_et_al")

cat("################  Q3: SPECIES-LEVEL LONGITUDINAL  ################\n\n")

meta <- read.csv(file.path(data_dir, "meta_location.csv"), stringsAsFactors = FALSE)
gen_abund <- read.csv(file.path(data_dir, "genomic_genus_abundance.csv"),
                      row.names = 1, check.names = FALSE)
sp_abund  <- read.csv(file.path(data_dir, "genomic_species_abundance.csv"),
                      row.names = 1, check.names = FALSE)

df <- meta %>% mutate(collection_date = as.Date(collection_date)) %>%
  filter(!is.na(collection_date))

agg_key <- df %>% group_by(city, collection_date) %>%
  summarise(sample_ids = list(complete_name), .groups = "drop")

agg_matrix <- function(mat, keys) {
  out <- matrix(NA, nrow(keys), ncol(mat)); colnames(out) <- colnames(mat)
  for (i in seq_len(nrow(keys)))
    out[i, ] <- colMeans(mat[keys$sample_ids[[i]], , drop = FALSE])
  out
}
gen_agg <- agg_matrix(gen_abund, agg_key)
sp_agg  <- agg_matrix(sp_abund, agg_key)

daily_temps <- readRDS(file.path(data_dir, "becsei_climate.rds"))
ts_data <- agg_key %>% select(city, collection_date)
ts_data$T_air_30d <- NA_real_
for (i in seq_len(nrow(ts_data))) {
  ct <- daily_temps %>% filter(city == ts_data$city[i])
  w30 <- ct$t_mean[ct$date >= (ts_data$collection_date[i] - 30) &
                     ct$date <= ts_data$collection_date[i]]
  if (length(w30) > 0) ts_data$T_air_30d[i] <- mean(w30, na.rm = TRUE)
}
idx    <- which(!is.na(ts_data$T_air_30d))
temp   <- ts_data$T_air_30d[idx]
city_f <- factor(ts_data$city[idx])
cat("n city-date obs with temp:", length(idx), "\n")
cat("cities:", paste(levels(city_f), collapse = ", "), "\n\n")

gen_total <- rowSums(gen_agg[idx, ])
sp_total  <- rowSums(sp_agg[idx, ])

get_sp_ra <- function(name) {
  cols <- grep(paste0("^", name, "$"), colnames(sp_agg), value = TRUE)
  if (length(cols) == 0) return(rep(NA_real_, length(idx)))
  sp_agg[idx, cols[1]] / sp_total
}

# City-adjusted partial Spearman (within-city association with temperature)
partial_rho <- function(y) {
  ok <- !is.na(y)
  if (sum(ok) < 20) return(c(rho = NA, p = NA))
  y_r <- residuals(lm(y[ok] ~ city_f[ok]))
  t_r <- residuals(lm(temp[ok] ~ city_f[ok]))
  ct <- suppressWarnings(cor.test(t_r, y_r, method = "spearman"))
  c(rho = unname(ct$estimate), p = ct$p.value)
}

who_species <- c(
  "Acinetobacter baumannii", "Klebsiella pneumoniae", "Escherichia coli",
  "Enterobacter cloacae", "Mycobacterium tuberculosis",
  "Enterococcus faecium", "Enterococcus faecalis", "Staphylococcus aureus",
  "Pseudomonas aeruginosa", "Salmonella enterica", "Neisseria gonorrhoeae",
  "Streptococcus pneumoniae", "Streptococcus pyogenes", "Haemophilus influenzae"
)

q3 <- map_dfr(who_species, function(sp) {
  r <- partial_rho(get_sp_ra(sp))
  tibble(species = sp, long_rho = r["rho"], long_p = r["p"])
})
print(as.data.frame(q3), digits = 3)

# Cross-sectional adjusted species estimates for side-by-side
cs <- read.csv(file.path(base_dir, "Results", "adjusted_species_temp.csv"),
               stringsAsFactors = FALSE)
cat("\nCross-sectional species file columns:\n"); print(names(cs))
saveRDS(list(q3 = q3, cs = cs), file.path(base_dir, "Results", "reviewer_q3_species_long.rds"))

cat("\n\n################  Q4: BACTERIOME COMPOSITION BY LOCATION  ################\n\n")

# ---- (a) Cross-sectional: variance in bacteriome explained by country vs temp ----
analysis <- readRDS(file.path(base_dir, "Datasets", "analysis_ready.rds"))
genus_raw <- readRDS(file.path(base_dir, "Datasets", "genus_counts_cache.rds")) %>%
  mutate(genepid = as.character(as.integer(genepid)))
resistome <- readRDS(file.path(base_dir, "Datasets", "resistome_matrices.rds"))

genera <- setdiff(names(genus_raw), "genepid")
prev <- colSums(genus_raw[, genera] > 0) / nrow(genus_raw)
keep_genera <- names(prev[prev >= 0.05])
genus_filt <- genus_raw %>% select(genepid, all_of(keep_genera))

shared_ids <- Reduce(intersect, list(
  genus_filt$genepid, resistome$fg_cluster_clr$genepid, analysis$genepid))

genus_sub <- genus_filt %>% filter(genepid %in% shared_ids) %>% arrange(genepid)
clim_sub  <- analysis %>% filter(genepid %in% shared_ids) %>% arrange(genepid)
genus_ra  <- genus_sub %>% mutate(across(-genepid, ~ . / rowSums(across(-genepid))))

keep <- complete.cases(clim_sub[, c("T_30d", "Region", "year",
                                     "gdp_pcap_ppp", "sanitation", "health_exp_gdp",
                                     "oop_health_exp", "immunization_dpt",
                                     "animal_amc_mgkg", "pop_density")])
genus_ra <- genus_ra[keep, ]
clim_sub <- clim_sub[keep, ]
cat("Cross-sectional samples:", nrow(clim_sub),
    "| cities:", n_distinct(clim_sub$city),
    "| countries:", n_distinct(clim_sub$country), "\n\n")

bray <- vegdist(genus_ra %>% select(-genepid), method = "bray")
clim_sub <- clim_sub %>% mutate(year_f = factor(year),
                                country_f = factor(country))

# NB: R2 is exact from observed data; permutations only affect the p-value.
# Country factor has ~100 levels, so use a modest permutation count for speed.
set.seed(1)
cat("--- PERMANOVA: bacteriome ~ country + T_30d (marginal) ---\n"); flush.console()
pm1 <- adonis2(bray ~ country_f + T_30d, data = clim_sub,
               permutations = 99, by = "margin")
print(pm1); flush.console()

cat("\n--- PERMANOVA: bacteriome ~ T_30d alone ---\n"); flush.console()
pm2 <- adonis2(bray ~ T_30d, data = clim_sub, permutations = 999, by = "margin")
print(pm2); flush.console()

# ---- (b) 5-city: between-city vs within-city compositional dissimilarity ----
cat("\n--- 5-city genus composition: between- vs within-city Bray-Curtis ---\n")
gen_ra_5 <- gen_agg[idx, ] / rowSums(gen_agg[idx, ])
bray5 <- as.matrix(vegdist(gen_ra_5, method = "bray"))
cf <- as.character(city_f)
ut <- upper.tri(bray5)
same <- outer(cf, cf, `==`)
within_bc  <- mean(bray5[ut & same])
between_bc <- mean(bray5[ut & !same])
cat(sprintf("Mean within-city Bray-Curtis:  %.3f\n", within_bc))
cat(sprintf("Mean between-city Bray-Curtis: %.3f\n", between_bc))

set.seed(1)
pm5 <- adonis2(as.dist(bray5) ~ city_f, permutations = 999)
cat("\n5-city PERMANOVA bacteriome ~ city:\n"); print(pm5)

# Copenhagen vs Budapest specifically (reviewer's example)
cph <- which(cf == "Copenhagen"); bud <- which(cf == "Budapest")
cat(sprintf("\nCopenhagen vs Budapest mean between-Bray-Curtis: %.3f\n",
            mean(bray5[cph, bud])))
cat(sprintf("Copenhagen mean temp: %.1f C | Budapest mean temp: %.1f C\n",
            mean(temp[cph]), mean(temp[bud])))

saveRDS(list(pm_country = pm1, pm_temp = pm2, within_bc = within_bc,
             between_bc = between_bc, pm5 = pm5,
             bray5 = bray5, cf = cf, temp5 = temp,
             gen_ra_5 = gen_ra_5),
        file.path(base_dir, "Results", "reviewer_q4_location.rds"))

cat("\nDONE.\n")
