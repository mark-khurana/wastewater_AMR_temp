# Layer 2: Partial PERMANOVA with socioeconomic confounders, variance partitioning, driver comparison, AMC sensitivity.
# Outputs: Results/layer2_partial_permanova.rds, Results/layer2_variance_partitioning.csv, Results/layer2_driver_comparison.csv, Results/layer2_dissociation.csv, Results/layer2_amc_sensitivity.csv

library(here)
library(tidyverse)
library(vegan)
library(broom)

analysis  <- readRDS(here("Datasets", "analysis_ready.rds"))
resistome <- readRDS(here("Datasets", "resistome_matrices.rds"))

conf_vars <- c("gdp_pcap_ppp", "sanitation", "pop_density", "health_exp_gdp",
               "oop_health_exp", "tourism_arrivals", "immunization_dpt",
               "human_amc_ddd", "animal_amc_mgkg")

run_partial_permanova <- function(clr_mat, climate_df, label) {
  shared <- intersect(clr_mat$genepid, climate_df$genepid)
  clr_sub <- clr_mat %>% filter(genepid %in% shared) %>% arrange(genepid)
  clim_sub <- climate_df %>% filter(genepid %in% shared) %>% arrange(genepid)

  keep <- complete.cases(clim_sub[, c("T_30d", "Region",
                                       "gdp_pcap_ppp", "sanitation",
                                       "health_exp_gdp", "oop_health_exp",
                                       "immunization_dpt", "animal_amc_mgkg",
                                       "pop_density", "year")])
  clr_sub <- clr_sub[keep, ]
  clim_sub <- clim_sub[keep, ]

  dist_mat <- dist(clr_sub %>% select(-genepid))

  clim_sub <- clim_sub %>%
    mutate(log_gdp = log(gdp_pcap_ppp),
           log_pop_density = log(pop_density + 1),
           log_animal_amc = log(animal_amc_mgkg + 1),
           year_f = factor(year))

  naive <- adonis2(dist_mat ~ T_30d + Region + year_f,
                   data = clim_sub, permutations = 999, by = "margin")

  full <- adonis2(dist_mat ~ T_30d + log_gdp + sanitation +
                    health_exp_gdp + oop_health_exp + immunization_dpt +
                    log_pop_density + log_animal_amc + Region + year_f,
                  data = clim_sub, permutations = 999, by = "margin")

  list(naive = naive, full = full, n = nrow(clr_sub), label = label)
}

perm_fg  <- run_partial_permanova(resistome$fg_cluster_clr, analysis, "FG resistome")
perm_acq <- run_partial_permanova(resistome$acq_cluster_clr, analysis, "Acquired resistome")

run_varpart <- function(clr_mat, climate_df, label) {
  shared <- intersect(clr_mat$genepid, climate_df$genepid)
  clr_sub <- clr_mat %>% filter(genepid %in% shared) %>% arrange(genepid)
  clim_sub <- climate_df %>% filter(genepid %in% shared) %>% arrange(genepid)

  keep <- complete.cases(clim_sub[, c("T_30d", "Region",
                                       "gdp_pcap_ppp", "sanitation",
                                       "health_exp_gdp", "oop_health_exp",
                                       "immunization_dpt", "animal_amc_mgkg",
                                       "pop_density", "year")])
  clr_sub <- clr_sub[keep, ]
  clim_sub <- clim_sub[keep, ]

  Y <- clr_sub %>% select(-genepid) %>% as.matrix()

  temp_vars <- clim_sub %>%
    transmute(BIO1 = T_30d)

  year_dummies <- model.matrix(~ factor(year), data = clim_sub)[, -1, drop = FALSE]

  socio_vars <- clim_sub %>%
    transmute(log_gdp = log(gdp_pcap_ppp),
              sanitation = sanitation,
              health_exp = health_exp_gdp,
              oop_exp = oop_health_exp,
              immunization = immunization_dpt,
              log_pop_density = log(pop_density + 1),
              log_animal_amc = log(animal_amc_mgkg + 1)) %>%
    bind_cols(as_tibble(year_dummies))

  vp <- varpart(Y, temp_vars, socio_vars)

  fractions <- tibble(
    resistome = label,
    temp_pure = vp$part$indfract$Adj.R.squared[1],
    shared = vp$part$indfract$Adj.R.squared[2],
    socio_pure = vp$part$indfract$Adj.R.squared[3],
    residual = vp$part$indfract$Adj.R.squared[4]
  )

  list(varpart = vp, fractions = fractions)
}

vp_fg  <- run_varpart(resistome$fg_cluster_clr, analysis, "FG")
vp_acq <- run_varpart(resistome$acq_cluster_clr, analysis, "Acquired")

varpart_summary <- bind_rows(vp_fg$fractions, vp_acq$fractions)

cruz_loya <- resistome$cruz_loya

run_driver_comparison <- function(class_clr, climate_df, resistome_type) {
  shared <- intersect(class_clr$genepid, climate_df$genepid)
  clr_sub <- class_clr %>% filter(genepid %in% shared) %>% arrange(genepid)
  clim_sub <- climate_df %>% filter(genepid %in% shared) %>% arrange(genepid)

  drug_classes <- setdiff(names(clr_sub), "genepid")

  map_dfr(drug_classes, function(dc) {
    df <- tibble(genepid = clr_sub$genepid, clr_abundance = clr_sub[[dc]]) %>%
      left_join(clim_sub %>%
                  select(genepid, T_30d, Region,
                         gdp_pcap_ppp, sanitation, health_exp_gdp,
                         oop_health_exp, immunization_dpt, animal_amc_mgkg,
                         pop_density, year),
                by = "genepid") %>%
      filter(complete.cases(.)) %>%
      mutate(log_gdp = log(gdp_pcap_ppp),
             log_pop_density = log(pop_density + 1),
             log_animal_amc = log(animal_amc_mgkg + 1),
             year_f = factor(year))

    if (nrow(df) < 30) return(NULL)

    m_naive <- lm(clr_abundance ~ T_30d + Region + year_f, data = df)
    t_naive <- tidy(m_naive) %>% filter(term == "T_30d")

    m_full <- lm(clr_abundance ~ T_30d + log_gdp + sanitation +
                   health_exp_gdp + oop_health_exp + immunization_dpt +
                   log_pop_density + log_animal_amc + Region + year_f, data = df)
    t_full <- tidy(m_full) %>% filter(term == "T_30d")

    tibble(
      resistome = resistome_type, drug_class = dc, n = nrow(df),
      beta_naive = t_naive$estimate, p_naive = t_naive$p.value,
      beta_full = t_full$estimate, p_full = t_full$p.value,
      beta_change_pct = 100 * (t_full$estimate - t_naive$estimate) / abs(t_naive$estimate),
      r2_naive = summary(m_naive)$r.squared,
      r2_full = summary(m_full)$r.squared,
      r2_gain = summary(m_full)$r.squared - summary(m_naive)$r.squared
    )
  }) %>%
    left_join(cruz_loya, by = c("drug_class" = "class")) %>%
    mutate(
      stressor_group = replace_na(stressor_group, "unclassified"),
      p_adj_naive = p.adjust(p_naive, method = "BH"),
      p_adj_full = p.adjust(p_full, method = "BH")
    ) %>%
    arrange(p_adj_full)
}

driver_fg  <- run_driver_comparison(resistome$fg_class_clr, analysis, "FG")
driver_acq <- run_driver_comparison(resistome$acq_class_clr, analysis, "Acquired")
all_drivers <- bind_rows(driver_fg, driver_acq)

extract_dissociation <- function(perm_result) {
  full <- perm_result$full
  predictors <- rownames(full)[!rownames(full) %in% c("Residual", "Total")]

  map_dfr(predictors, function(pred) {
    tibble(
      predictor = pred,
      R2 = full[pred, "R2"],
      p_value = full[pred, "Pr(>F)"]
    )
  })
}

dissoc_fg <- extract_dissociation(perm_fg) %>% mutate(resistome = "FG")
dissoc_acq <- extract_dissociation(perm_acq) %>% mutate(resistome = "Acquired")
dissociation <- bind_rows(dissoc_fg, dissoc_acq)

run_amc_sensitivity <- function(clr_mat, climate_df, label) {
  shared <- intersect(clr_mat$genepid, climate_df$genepid)
  clr_sub <- clr_mat %>% filter(genepid %in% shared) %>% arrange(genepid)
  clim_sub <- climate_df %>% filter(genepid %in% shared) %>% arrange(genepid)

  keep <- complete.cases(clim_sub[, c("T_30d", "Region",
                                       "gdp_pcap_ppp", "sanitation",
                                       "health_exp_gdp", "oop_health_exp",
                                       "immunization_dpt", "human_amc_ddd",
                                       "animal_amc_mgkg", "pop_density",
                                       "year")])
  clr_sub <- clr_sub[keep, ]
  clim_sub <- clim_sub[keep, ]

  if (nrow(clr_sub) < 50) {
    return(NULL)
  }

  dist_mat <- dist(clr_sub %>% select(-genepid))

  clim_sub <- clim_sub %>%
    mutate(log_gdp = log(gdp_pcap_ppp),
           log_pop_density = log(pop_density + 1),
           log_human_amc = log(human_amc_ddd),
           log_animal_amc = log(animal_amc_mgkg + 1),
           year_f = factor(year))

  no_amc <- adonis2(dist_mat ~ T_30d + log_gdp + sanitation +
                      health_exp_gdp + oop_health_exp + immunization_dpt +
                      log_pop_density + log_animal_amc + Region + year_f,
                    data = clim_sub, permutations = 999, by = "margin")

  with_amc <- adonis2(dist_mat ~ T_30d + log_gdp + sanitation +
                        health_exp_gdp + oop_health_exp + immunization_dpt +
                        log_pop_density + log_animal_amc + log_human_amc +
                        Region + year_f,
                      data = clim_sub, permutations = 999, by = "margin")

  list(no_amc = no_amc, with_amc = with_amc, n = nrow(clr_sub), label = label)
}

amc_fg  <- run_amc_sensitivity(resistome$fg_cluster_clr, analysis, "FG")
amc_acq <- run_amc_sensitivity(resistome$acq_cluster_clr, analysis, "Acquired")

if (!is.null(amc_fg) && !is.null(amc_acq)) {
  amc_sensitivity <- tibble(
    resistome = c("FG", "Acquired"),
    n_samples = c(amc_fg$n, amc_acq$n),
    r2_bio1_no_amc = c(amc_fg$no_amc["T_30d", "R2"],
                        amc_acq$no_amc["T_30d", "R2"]),
    r2_bio1_with_amc = c(amc_fg$with_amc["T_30d", "R2"],
                          amc_acq$with_amc["T_30d", "R2"]),
    r2_human_amc = c(amc_fg$with_amc["log_human_amc", "R2"],
                      amc_acq$with_amc["log_human_amc", "R2"]),
    p_human_amc = c(amc_fg$with_amc["log_human_amc", "Pr(>F)"],
                     amc_acq$with_amc["log_human_amc", "Pr(>F)"])
  )
}

saveRDS(list(perm_fg = perm_fg, perm_acq = perm_acq,
             vp_fg = vp_fg$varpart, vp_acq = vp_acq$varpart,
             amc_fg = amc_fg, amc_acq = amc_acq),
        here("Results", "layer2_partial_permanova.rds"))
write_csv(varpart_summary, here("Results", "layer2_variance_partitioning.csv"))
write_csv(all_drivers, here("Results", "layer2_driver_comparison.csv"))
write_csv(dissociation, here("Results", "layer2_dissociation.csv"))
if (exists("amc_sensitivity")) {
  write_csv(amc_sensitivity, here("Results", "layer2_amc_sensitivity.csv"))
}
