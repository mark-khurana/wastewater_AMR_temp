# Re-runs only the dispersion test from 13_reviewer_sensitivity.R at 999 permutations, so the
# reported p-value is not sitting on the 99-permutation floor. betadisper is cheap; the
# PERMANOVA arms in script 13 are not, and keep the 99 permutations already documented in the
# supplementary methods for sensitivity analyses.
# Output: Results/rev_dispersion.csv (overwrites the 99-permutation version)

library(here); library(tidyverse); library(vegan); library(compositions)

set.seed(20260923)
N_PERM <- 999

analysis <- readRDS(here("Datasets", "analysis_ready.rds"))
counts   <- readRDS(here("Datasets", "clean_ARG_counts.rds"))

model_vars <- c("T_30d", "Region", "gdp_pcap_ppp", "sanitation", "health_exp_gdp",
                "oop_health_exp", "immunization_dpt", "animal_amc_mgkg",
                "pop_density", "year")

build_matrix <- function(df) {
  df %>%
    filter(!is.na(cluster_representative_98)) %>%
    group_by(genepid, cluster_representative_98) %>%
    summarise(abundance = sum(fragmentCountAln_adj, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = cluster_representative_98, values_from = abundance, values_fill = 0)
}

mats <- list(
  Functional = build_matrix(counts %>% filter(functional_amr == "True")),
  Acquired   = build_matrix(counts %>% filter(resfinder == "True"))
)

out <- list()

for (res_type in names(mats)) {
  m      <- mats[[res_type]]
  # Primary zero replacement (03_merge_analysis_ready.R): 65% of the per-cluster minimum
  # non-zero value, taken over all samples.
  dl     <- apply(as.matrix(m[, -1]), 2, function(x) min(x[x > 0]))
  shared <- intersect(m$genepid, analysis$genepid)
  m      <- m        %>% filter(genepid %in% shared) %>% arrange(genepid)
  clim   <- analysis %>% filter(genepid %in% shared) %>% arrange(genepid)
  keep   <- complete.cases(clim[, model_vars])
  m      <- as.matrix(m[keep, -1, drop = FALSE])
  clim   <- clim[keep, ]

  for (j in seq_len(ncol(m))) m[m[, j] == 0, j] <- 0.65 * dl[j]
  d <- dist(as.matrix(clr(m)))

  groupings <- list(
    "Temperature quartile" = cut(clim$T_30d, breaks = quantile(clim$T_30d, 0:4 / 4),
                                 include.lowest = TRUE, labels = paste0("Q", 1:4)),
    "WHO region"           = factor(clim$Region)
  )

  for (g in names(groupings)) {
    message(sprintf("[disp999] %s | %s", res_type, g))
    bd <- betadisper(d, groupings[[g]], type = "centroid")
    pt <- permutest(bd, permutations = N_PERM)
    mn <- tapply(bd$distances, bd$group, mean)
    out[[length(out) + 1]] <- tibble(
      resistome = res_type, grouping = g, n_groups = nlevels(groupings[[g]]),
      permutations = N_PERM,
      F_value = pt$tab$F[1], p_value = pt$tab$`Pr(>F)`[1],
      min_group_mean_dist = min(mn), max_group_mean_dist = max(mn),
      spread_ratio = max(mn) / min(mn),
      group_means = paste(sprintf("%s=%.2f", names(mn), mn), collapse = "; ")
    )
  }
}

bind_rows(out) %>% write_csv(here("Results", "rev_dispersion.csv"))
message("done")
