# Layer 1: Naive temperature-resistome associations (PERMANOVA, drug-class models, GAMs, metric tournament).
# Outputs: Results/layer1_permanova.rds, Results/layer1_drug_class_models.csv, Results/layer1_gam_results.csv, Results/layer1_metric_tournament.csv

library(here)
library(tidyverse)
library(vegan)
library(mgcv)
library(broom)

analysis  <- readRDS(here("Datasets", "analysis_ready.rds"))
resistome <- readRDS(here("Datasets", "resistome_matrices.rds"))

run_permanova <- function(clr_mat, climate_df, label) {
  shared <- intersect(clr_mat$genepid, climate_df$genepid)
  clr_sub <- clr_mat %>% filter(genepid %in% shared) %>% arrange(genepid)
  clim_sub <- climate_df %>% filter(genepid %in% shared) %>% arrange(genepid)

  keep <- !is.na(clim_sub$T_30d) & !is.na(clim_sub$Region)
  clr_sub <- clr_sub[keep, ]
  clim_sub <- clim_sub[keep, ]

  dist_mat <- dist(clr_sub %>% select(-genepid))

  result <- adonis2(
    dist_mat ~ T_30d + Region,
    data = clim_sub, permutations = 999, by = "margin"
  )

  result
}

perm_fg  <- run_permanova(resistome$fg_cluster_clr, analysis, "FG resistome")
perm_acq <- run_permanova(resistome$acq_cluster_clr, analysis, "Acquired resistome")

r2_fg  <- perm_fg["T_30d", "R2"]
r2_acq <- perm_acq["T_30d", "R2"]

cruz_loya <- resistome$cruz_loya

run_class_models <- function(class_clr, climate_df, resistome_type) {
  shared <- intersect(class_clr$genepid, climate_df$genepid)
  clr_sub <- class_clr %>% filter(genepid %in% shared) %>% arrange(genepid)
  clim_sub <- climate_df %>% filter(genepid %in% shared) %>% arrange(genepid)

  drug_classes <- setdiff(names(clr_sub), "genepid")

  map_dfr(drug_classes, function(dc) {
    df <- tibble(genepid = clr_sub$genepid, clr_abundance = clr_sub[[dc]]) %>%
      left_join(clim_sub %>% select(genepid, T_30d, Region),
                by = "genepid") %>%
      filter(!is.na(T_30d), !is.na(Region))

    if (nrow(df) < 30) return(NULL)

    mod <- lm(clr_abundance ~ T_30d + Region, data = df)
    tid <- tidy(mod) %>% filter(term == "T_30d")

    tibble(
      resistome = resistome_type, drug_class = dc, n = nrow(df),
      beta_temp = tid$estimate, se = tid$std.error,
      t_value = tid$statistic, p_value = tid$p.value,
      r2_full = summary(mod)$r.squared
    )
  }) %>%
    left_join(cruz_loya, by = c("drug_class" = "class")) %>%
    mutate(
      stressor_group = replace_na(stressor_group, "unclassified"),
      p_adj = p.adjust(p_value, method = "BH")
    ) %>%
    arrange(p_adj)
}

fg_class_results  <- run_class_models(resistome$fg_class_clr, analysis, "FG")
acq_class_results <- run_class_models(resistome$acq_class_clr, analysis, "Acquired")
all_class_results <- bind_rows(fg_class_results, acq_class_results)

run_gams <- function(class_clr, climate_df, resistome_type, target_classes = NULL) {
  shared <- intersect(class_clr$genepid, climate_df$genepid)
  clr_sub <- class_clr %>% filter(genepid %in% shared) %>% arrange(genepid)
  clim_sub <- climate_df %>% filter(genepid %in% shared) %>% arrange(genepid)

  if (is.null(target_classes)) {
    target_classes <- intersect(setdiff(names(clr_sub), "genepid"), cruz_loya$class)
  }

  map_dfr(target_classes, function(dc) {
    df <- tibble(clr_abundance = clr_sub[[dc]],
                 BIO1 = clim_sub$T_30d,
                 Region = clim_sub$Region) %>%
      filter(!is.na(BIO1), !is.na(Region))
    if (nrow(df) < 50) return(NULL)

    gam_mod <- gam(clr_abundance ~ s(BIO1, k = 5) + Region, data = df, method = "REML")
    lm_mod  <- gam(clr_abundance ~ BIO1 + Region, data = df, method = "REML")
    gam_sum <- summary(gam_mod)

    tibble(
      resistome = resistome_type, drug_class = dc, n = nrow(df),
      edf = gam_sum$s.table[1, "edf"],
      p_smooth = gam_sum$s.table[1, "p-value"],
      aic_linear = AIC(lm_mod), aic_gam = AIC(gam_mod),
      delta_aic = AIC(lm_mod) - AIC(gam_mod),
      deviance_explained = gam_sum$dev.expl,
      nonlinear = gam_sum$s.table[1, "edf"] > 1.5 & gam_sum$s.table[1, "p-value"] < 0.05
    )
  })
}

gam_results <- bind_rows(
  run_gams(resistome$fg_class_clr, analysis, "FG"),
  run_gams(resistome$acq_class_clr, analysis, "Acquired")
)

temp_metrics <- c(
  "T_30d", "BIO10_mean_temp_warmest_quarter",
  "T_sampling", "T_30d", "T_90d", "T_365d",
  "T_annual_mean", "T_max", "T_days30", "T_95pct",
  "T_amplitude", "T_variability"
)

metric_comparison <- map_dfr(temp_metrics, function(metric) {
  shared <- intersect(resistome$fg_cluster_clr$genepid, analysis$genepid)
  clr_sub <- resistome$fg_cluster_clr %>% filter(genepid %in% shared) %>% arrange(genepid)
  clim_sub <- analysis %>% filter(genepid %in% shared) %>% arrange(genepid)

  keep <- !is.na(clim_sub[[metric]]) & !is.na(clim_sub$Region)
  clr_sub <- clr_sub[keep, ]; clim_sub <- clim_sub[keep, ]

  dist_mat <- dist(clr_sub %>% select(-genepid))
  clim_sub$temp_var <- clim_sub[[metric]]

  result <- adonis2(dist_mat ~ temp_var + Region, data = clim_sub,
                    permutations = 999, by = "margin")

  tibble(metric = metric, n = nrow(clr_sub),
         R2 = result["temp_var", "R2"],
         F_value = result["temp_var", "F"],
         p_value = result["temp_var", "Pr(>F)"])
})

best <- metric_comparison %>% slice_max(R2, n = 1)

saveRDS(list(perm_fg = perm_fg, perm_acq = perm_acq), here("Results", "layer1_permanova.rds"))
write_csv(all_class_results, here("Results", "layer1_drug_class_models.csv"))
write_csv(gam_results, here("Results", "layer1_gam_results.csv"))
write_csv(metric_comparison, here("Results", "layer1_metric_tournament.csv"))
