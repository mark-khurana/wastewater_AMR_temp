# Within-city temporal analysis using Becsei et al. longitudinal data (5 European cities)

library(tidyverse)
library(lubridate)
library(vegan)
library(compositions)
library(nasapower)
library(here)

meta <- read.csv(here("Datasets", "Becsei_et_al", "meta_location.csv"), stringsAsFactors = FALSE)
res_class <- read.csv(here("Datasets", "Becsei_et_al", "resfinder_class_abundance.csv"),
                       row.names = 1, check.names = FALSE)
res_gene <- read.csv(here("Datasets", "Becsei_et_al", "resfinder_gene_abundance.csv"),
                      row.names = 1, check.names = FALSE)
gen_abund <- read.csv(here("Datasets", "Becsei_et_al", "genomic_genus_abundance.csv"),
                       row.names = 1, check.names = FALSE)

df <- meta %>%
  mutate(collection_date = as.Date(collection_date)) %>%
  filter(!is.na(collection_date))

agg_key <- df %>%
  group_by(city, collection_date) %>%
  summarise(sample_ids = list(complete_name), .groups = "drop")

agg_matrix <- function(mat, keys) {
  out <- matrix(NA, nrow = nrow(keys), ncol = ncol(mat))
  colnames(out) <- colnames(mat)
  for (i in seq_len(nrow(keys))) {
    ids <- keys$sample_ids[[i]]
    out[i, ] <- colMeans(mat[ids, , drop = FALSE])
  }
  out
}

res_class_agg <- agg_matrix(res_class, agg_key)
res_gene_agg  <- agg_matrix(res_gene, agg_key)
gen_agg       <- agg_matrix(gen_abund, agg_key)

ts_data <- agg_key %>%
  select(city, collection_date) %>%
  mutate(month = month(collection_date), year = year(collection_date))

city_coords <- tibble(
  city = c("Copenhagen", "Rotterdam", "Bologna", "Budapest", "Rome"),
  lat  = c(55.6761, 51.9225, 44.4949, 47.4979, 41.9028),
  lon  = c(12.5683, 4.4792, 11.3426, 19.0402, 12.4964)
)

climate_cache <- here("Datasets", "Becsei_et_al", "becsei_climate.rds")
if (file.exists(climate_cache)) {
  daily_temps <- readRDS(climate_cache)
} else {
  daily_list <- list()
  for (i in seq_len(nrow(city_coords))) {
    cc <- city_coords[i, ]
    yearly_dfs <- list()
    for (yr in 2018:2021) {
      d <- tryCatch(
        get_power(community = "ag", lonlat = c(cc$lon, cc$lat),
                  pars = c("T2M", "T2M_MAX", "T2M_MIN"),
                  dates = c(paste0(yr, "-01-01"), paste0(yr, "-12-31")),
                  temporal_api = "daily"),
        error = function(e) { warning(e$message); NULL })
      if (!is.null(d))
        yearly_dfs[[as.character(yr)]] <- d %>%
          transmute(date = YYYYMMDD, t_mean = T2M, t_max = T2M_MAX, t_min = T2M_MIN)
      Sys.sleep(0.5)
    }
    city_daily <- bind_rows(yearly_dfs) %>% mutate(city = cc$city)
    daily_list[[cc$city]] <- city_daily
  }
  daily_temps <- bind_rows(daily_list)
  saveRDS(daily_temps, climate_cache)
}

ts_data <- ts_data %>% left_join(city_coords, by = "city")
ts_data$T_air_30d <- NA_real_
for (i in seq_len(nrow(ts_data))) {
  ct <- daily_temps %>% filter(city == ts_data$city[i])
  w30 <- ct$t_mean[ct$date >= (ts_data$collection_date[i] - 30) &
                    ct$date <= ts_data$collection_date[i]]
  if (length(w30) > 0) ts_data$T_air_30d[i] <- mean(w30, na.rm = TRUE)
}

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

res_clr <- do_clr(res_class_agg)

gen_mat <- as.data.frame(gen_agg)
keep_gen <- (colSums(gen_mat > 0) / nrow(gen_mat)) > 0.2 & colSums(gen_mat) > 100
gen_mat <- gen_mat[, keep_gen]
gen_clr <- do_clr(gen_mat)

gen_ra <- gen_mat / rowSums(gen_mat)
gen_bray <- vegdist(gen_ra, method = "bray")
gen_pcoa <- cmdscale(gen_bray, k = 5, eig = TRUE)
gen_eig_pct <- round(100 * gen_pcoa$eig[1:5] / sum(gen_pcoa$eig[gen_pcoa$eig > 0]), 1)
ts_data$bac_PC1 <- gen_pcoa$points[, 1]
ts_data$bac_PC2 <- gen_pcoa$points[, 2]
ts_data$bac_PC3 <- gen_pcoa$points[, 3]
ts_data$bac_PC4 <- gen_pcoa$points[, 4]
ts_data$bac_PC5 <- gen_pcoa$points[, 5]

gene_names <- colnames(res_gene_agg)
family_map <- list(
  dfr  = grep("^dfr", gene_names, value = TRUE),
  aac  = grep("^aac", gene_names, value = TRUE),
  aph  = grep("^aph", gene_names, value = TRUE),
  ant  = grep("^ant", gene_names, value = TRUE),
  tet  = grep("^tet", gene_names, value = TRUE),
  qnr  = grep("^qnr", gene_names, value = TRUE),
  erm  = grep("^erm", gene_names, value = TRUE),
  blaTEM  = grep("^blaTEM", gene_names, value = TRUE),
  blaSHV  = grep("^blaSHV", gene_names, value = TRUE),
  `blaCTX-M` = grep("^blaCTX", gene_names, value = TRUE),
  blaOXA  = grep("^blaOXA", gene_names, value = TRUE)
)

gene_family_mat <- matrix(NA, nrow = nrow(res_gene_agg), ncol = length(family_map))
colnames(gene_family_mat) <- names(family_map)
for (fam in names(family_map)) {
  cols <- family_map[[fam]]
  if (length(cols) == 1) {
    gene_family_mat[, fam] <- res_gene_agg[, cols]
  } else if (length(cols) > 1) {
    gene_family_mat[, fam] <- rowSums(res_gene_agg[, cols])
  }
}

gene_fam_clr <- do_clr(gene_family_mat)

fit_fe <- function(y, data) {
  data$y <- y
  m <- lm(y ~ T_air_30d + city_f + month_f, data = data)
  ms <- summary(m)
  tibble(
    beta_temp = coef(m)["T_air_30d"],
    se = ms$coefficients["T_air_30d", "Std. Error"],
    t_stat = ms$coefficients["T_air_30d", "t value"],
    p_value = ms$coefficients["T_air_30d", "Pr(>|t|)"],
    r_squared = ms$r.squared
  )
}

mod_data <- ts_data %>%
  filter(!is.na(T_air_30d)) %>%
  mutate(city_f = factor(city), month_f = factor(month))
n_mod <- nrow(mod_data)

genus_results <- map_dfr(colnames(gen_clr), function(g) {
  y <- gen_clr[[g]][!is.na(ts_data$T_air_30d)]
  fe_result <- fit_fe(y, mod_data)

  temp_resid <- residuals(lm(T_air_30d ~ city_f, data = mod_data))
  abund_resid <- residuals(lm(y ~ city_f, data = mod_data))
  rho_city_adj <- suppressWarnings(cor.test(abund_resid, temp_resid, method = "spearman"))

  fe_result %>%
    mutate(genus = g,
           rho_temporal = rho_city_adj$estimate,
           rho_temporal_p = rho_city_adj$p.value)
}) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH"),
         direction = ifelse(beta_temp > 0, "warm", "cool")) %>%
  arrange(p_value)

martiny_genus <- read.csv(here("Results", "layer3_genus_temp_correlations.csv")) %>%
  filter(!grepl("mOTU|ext_", genus)) %>%
  select(genus, rho_cross = rho, p_cross = p_value, category)

concordance <- inner_join(
  genus_results %>% select(genus, beta_temp, rho_temporal, rho_temporal_p, p_value, p_adj, direction),
  martiny_genus, by = "genus"
) %>%
  mutate(cross_dir = ifelse(rho_cross > 0, "warm", "cool"),
         concordant = direction == cross_dir)

rho_all <- cor.test(concordance$rho_temporal, concordance$rho_cross, method = "spearman")

sig_mart <- concordance %>% filter(p_cross < 0.05 / nrow(martiny_genus))
if (nrow(sig_mart) > 10) {
  rho_sig <- cor.test(sig_mart$rho_temporal, sig_mart$rho_cross, method = "spearman")
}

sig_both <- concordance %>% filter(p_value < 0.05, p_cross < 0.05)

drug_results <- map_dfr(colnames(res_clr), function(dc) {
  fit_fe(res_clr[[dc]][!is.na(ts_data$T_air_30d)], mod_data) %>%
    mutate(drug_class = dc)
}) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH"),
         direction = ifelse(beta_temp > 0, "warm", "cool")) %>%
  arrange(p_value)

martiny_drug <- read.csv(here("Results", "layer1_drug_class_models.csv")) %>%
  filter(resistome == "Acquired") %>%
  select(drug_class, beta_cross = beta_temp, p_cross = p_value)

drug_name_map <- c(
  "Aminoglycoside" = "aminoglycoside",
  "Beta-Lactam" = "beta_lactam",
  "Folate Pathway Antagonist" = "folate_pathway_antagonist",
  "Fosfomycin" = "fosfomycin",
  "Glycopeptide" = "glycopeptide",
  "Macrolide" = "macrolide",
  "Quinolone" = "quinolone",
  "Tetracycline" = "tetracycline",
  "Polymyxin" = "polymyxin",
  "Steroid Antibacterial" = "steroid_antibacterial",
  "Lincosamide" = "lincosamide"
)

drug_conc <- drug_results %>%
  mutate(martiny_name = drug_name_map[drug_class]) %>%
  filter(!is.na(martiny_name)) %>%
  inner_join(martiny_drug, by = c("martiny_name" = "drug_class")) %>%
  mutate(cross_dir = ifelse(beta_cross > 0, "warm", "cool"),
         concordant = direction == cross_dir)

cruz_loya <- tibble(
  arg_family = c("dfr","aac","aph","ant","tet","qnr","erm","blaTEM","blaSHV","blaCTX-M","blaOXA"),
  prediction = c("heat","heat","heat","heat","cold","cold","cold","neutral","neutral","neutral","heat"),
  group = c("Folate antagonist","Aminoglycoside","Aminoglycoside","Aminoglycoside",
            "Tetracycline","Quinolone","Macrolide","Beta-lactam","Beta-lactam","Beta-lactam","Beta-lactam")
)

gene_fam_results <- map_dfr(colnames(gene_fam_clr), function(fam) {
  fit_fe(gene_fam_clr[[fam]][!is.na(ts_data$T_air_30d)], mod_data) %>%
    mutate(arg_family = fam)
}) %>%
  left_join(cruz_loya, by = "arg_family") %>%
  mutate(
    direction = ifelse(beta_temp > 0, "warm", "cool"),
    prediction_met = case_when(
      prediction == "heat" & beta_temp > 0 ~ TRUE,
      prediction == "cold" & beta_temp < 0 ~ TRUE,
      prediction == "neutral" ~ NA,
      TRUE ~ FALSE
    )
  ) %>%
  arrange(p_value)

martiny_arg <- read.csv(here("Results", "layer3c_arg_family_temp.csv"))

gene_fam_results <- gene_fam_results %>%
  left_join(
    martiny_arg %>% select(arg_family, rho_cross = rho, rho_p_cross = rho_p),
    by = "arg_family"
  ) %>%
  mutate(
    cross_dir = ifelse(rho_cross > 0, "warm", "cool"),
    concordant_cross = direction == cross_dir
  )

mediation_results <- map_dfr(colnames(res_clr), function(dc) {
  y <- res_clr[[dc]][!is.na(ts_data$T_air_30d)]
  mod_data$y <- y

  m_base <- lm(y ~ T_air_30d + city_f + month_f, data = mod_data)
  m_med  <- lm(y ~ T_air_30d + bac_PC1 + bac_PC2 + bac_PC3 + bac_PC4 + bac_PC5 + city_f + month_f,
               data = mod_data)

  b_base <- coef(m_base)["T_air_30d"]
  b_med  <- coef(m_med)["T_air_30d"]

  tibble(
    drug_class = dc,
    beta_base = b_base,
    p_base = summary(m_base)$coefficients["T_air_30d", 4],
    beta_mediated = b_med,
    p_mediated = summary(m_med)$coefficients["T_air_30d", 4],
    bac_PC1_p = summary(m_med)$coefficients["bac_PC1", 4],
    attenuation_pct = ifelse(abs(b_base) > 1e-6, 100 * (1 - b_med / b_base), NA)
  )
}) %>% arrange(p_base)

sig_med <- mediation_results %>% filter(p_base < 0.05)

write.csv(genus_results, here("Results", "timeseries_genus_models.csv"), row.names = FALSE)
write.csv(concordance, here("Results", "timeseries_concordance.csv"), row.names = FALSE)
write.csv(drug_results, here("Results", "timeseries_drug_models.csv"), row.names = FALSE)
write.csv(drug_conc, here("Results", "timeseries_drug_concordance.csv"), row.names = FALSE)
write.csv(gene_fam_results, here("Results", "timeseries_gene_family_models.csv"), row.names = FALSE)
write.csv(mediation_results, here("Results", "timeseries_mediation.csv"), row.names = FALSE)
