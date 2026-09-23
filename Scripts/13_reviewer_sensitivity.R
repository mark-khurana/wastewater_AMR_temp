# Reviewer-requested sensitivity analyses for the JAC-AMR-2026-454 revision.
#   (1) Homogeneity of multivariate dispersion - the PERMANOVA assumption Reviewer 1 asks about.
#   (2) Zero-handling sensitivity for the CLR transformation.
#   (3) Aitchison vs Bray-Curtis on the same resistome.
# Model specification is lifted from 05_layer2_covariates.R so the arms are directly comparable
# to the primary adjusted result.
# Outputs: Results/rev_dispersion.csv, Results/rev_zero_sensitivity.csv, Results/rev_braycurtis.csv

library(here)
library(tidyverse)
library(vegan)
library(compositions)
library(zCompositions)

set.seed(20260923)
N_PERM <- as.integer(Sys.getenv("N_PERM", "999"))

analysis <- readRDS(here("Datasets", "analysis_ready.rds"))
counts   <- readRDS(here("Datasets", "clean_ARG_counts.rds"))

model_vars <- c("T_30d", "Region", "gdp_pcap_ppp", "sanitation", "health_exp_gdp",
                "oop_health_exp", "immunization_dpt", "animal_amc_mgkg",
                "pop_density", "year")

build_matrix <- function(df, group_col, value_col) {
  df %>%
    filter(!is.na(.data[[group_col]])) %>%
    group_by(genepid, .data[[group_col]]) %>%
    summarise(abundance = sum(.data[[value_col]], na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = all_of(group_col), values_from = abundance, values_fill = 0)
}

fg_counts  <- counts %>% filter(functional_amr == "True")
acq_counts <- counts %>% filter(resfinder == "True")

mats_adj <- list(
  Functional = build_matrix(fg_counts,  "cluster_representative_98", "fragmentCountAln_adj"),
  Acquired   = build_matrix(acq_counts, "cluster_representative_98", "fragmentCountAln_adj")
)
mats_raw <- list(
  Functional = build_matrix(fg_counts,  "cluster_representative_98", "fragmentCountAln"),
  Acquired   = build_matrix(acq_counts, "cluster_representative_98", "fragmentCountAln")
)

# Restrict a genepid-keyed matrix to samples with complete covariate data, and return the
# aligned covariate frame with the model transforms applied.
align <- function(mat_with_id) {
  shared  <- intersect(mat_with_id$genepid, analysis$genepid)
  m       <- mat_with_id %>% filter(genepid %in% shared) %>% arrange(genepid)
  clim    <- analysis   %>% filter(genepid %in% shared) %>% arrange(genepid)
  keep    <- complete.cases(clim[, model_vars])
  clim    <- clim[keep, ]
  list(
    mat  = as.matrix(m[keep, -1, drop = FALSE]),
    clim = clim %>% mutate(log_gdp          = log(gdp_pcap_ppp),
                           log_pop_density  = log(pop_density + 1),
                           log_animal_amc   = log(animal_amc_mgkg + 1),
                           year_f           = factor(year))
  )
}

# The fully adjusted PERMANOVA from 05_layer2_covariates.R. by = "margin" is set explicitly
# because vegan >= 2.7 defaults to by = NULL, which returns NA per-term R2.
adjusted_permanova <- function(dist_mat, clim, perms = N_PERM) {
  fit <- adonis2(dist_mat ~ T_30d + log_gdp + sanitation + health_exp_gdp +
                   oop_health_exp + immunization_dpt + log_pop_density +
                   log_animal_amc + Region + year_f,
                 data = clim, permutations = perms, by = "margin")
  tibble(n = nrow(clim), R2 = fit["T_30d", "R2"], p = fit["T_30d", "Pr(>F)"])
}

clr_mat <- function(m) as.matrix(clr(m))

# --------------------------------------------------------------------------------------
# (2) Zero-handling sensitivity
# --------------------------------------------------------------------------------------
# compositions::zeroreplace(x, d, a) imputes every zero with a * d, where d is the detection
# limit in the units of x and a defaults to 2/3. The published pipeline passes d = 0.65, which
# imputes a flat 0.433 - see Results/rev_zero_diagnostics.csv for how that compares with the
# observed non-zero values.
zero_arms <- function(mat_adj, mat_raw) {
  list(
    "Multiplicative replacement, d = 0.65 (published pipeline)" = function() clr_mat(zeroreplace(mat_adj, d = 0.65)),
    "Multiplicative replacement, 65% of per-cluster minimum detected value" = function() {
      dl <- apply(mat_adj, 2, function(x) { p <- x[x > 0]; if (length(p)) min(p) else NA_real_ })
      dl[is.na(dl)] <- min(mat_adj[mat_adj > 0])
      out <- mat_adj
      for (j in seq_len(ncol(out))) out[out[, j] == 0, j] <- 0.65 * dl[j]
      clr_mat(out)
    },
    "Multiplicative replacement, 65% of global minimum detected value" = function() {
      out <- mat_adj
      out[out == 0] <- 0.65 * min(mat_adj[mat_adj > 0])
      clr_mat(out)
    },
    "Bayesian-multiplicative replacement (CZM) on raw read counts" = function() {
      keep <- colSums(mat_raw) > 0
      clr_mat(as.matrix(cmultRepl(mat_raw[, keep, drop = FALSE],
                                  label = 0, method = "CZM", output = "prop",
                                  z.warning = 1, suppress.print = TRUE)))
    },
    "Pseudocount of 1 read added to raw read counts" = function() clr_mat(mat_raw + 1),
    "Clusters present in >= 50% of samples only (d = 0.65)" = function() {
      keep <- colMeans(mat_adj > 0) >= 0.5
      clr_mat(zeroreplace(mat_adj[, keep, drop = FALSE], d = 0.65))
    }
  )
}

zero_results <- list()
zero_diag    <- list()

for (res_type in names(mats_adj)) {
  a_adj <- align(mats_adj[[res_type]])
  a_raw <- align(mats_raw[[res_type]])
  stopifnot(nrow(a_adj$mat) == nrow(a_raw$mat))

  pos <- a_adj$mat[a_adj$mat > 0]
  zero_diag[[res_type]] <- tibble(
    resistome                    = res_type,
    n_samples                    = nrow(a_adj$mat),
    n_clusters                   = ncol(a_adj$mat),
    zero_fraction                = mean(a_adj$mat == 0),
    imputed_value_d065           = (2 / 3) * 0.65,
    median_nonzero               = median(pos),
    pct_nonzero_below_imputed    = 100 * mean(pos < (2 / 3) * 0.65),
    n_clusters_prevalent_50pct   = sum(colMeans(a_adj$mat > 0) >= 0.5)
  )

  for (arm in names(zero_arms(a_adj$mat, a_raw$mat))) {
    message(sprintf("[zero] %-12s | %s", res_type, arm))
    X <- zero_arms(a_adj$mat, a_raw$mat)[[arm]]()
    out <- adjusted_permanova(dist(X), a_adj$clim)
    zero_results[[length(zero_results) + 1]] <-
      out %>% mutate(resistome = res_type, zero_method = arm, n_clusters = ncol(X), .before = 1)
  }
}

bind_rows(zero_diag)    %>% write_csv(here("Results", "rev_zero_diagnostics.csv"))
bind_rows(zero_results) %>% write_csv(here("Results", "rev_zero_sensitivity.csv"))

# --------------------------------------------------------------------------------------
# (1) Homogeneity of multivariate dispersion
# --------------------------------------------------------------------------------------
# Tested on the Aitchison distance used for the primary analysis. Temperature is continuous,
# so it is binned into quartiles purely to make the dispersion test definable.
disp_results <- list()

for (res_type in names(mats_adj)) {
  a <- align(mats_adj[[res_type]])
  d <- dist(clr_mat(zeroreplace(a$mat, d = 0.65)))

  groupings <- list(
    "Temperature quartile" = cut(a$clim$T_30d, breaks = quantile(a$clim$T_30d, 0:4 / 4),
                                 include.lowest = TRUE, labels = paste0("Q", 1:4)),
    "WHO region"           = factor(a$clim$Region)
  )

  for (g in names(groupings)) {
    message(sprintf("[disp] %-12s | %s", res_type, g))
    bd <- betadisper(d, groupings[[g]], type = "centroid")
    pt <- permutest(bd, permutations = N_PERM)
    disp_results[[length(disp_results) + 1]] <- tibble(
      resistome    = res_type,
      grouping     = g,
      n_groups     = nlevels(groupings[[g]]),
      F_value      = pt$tab$F[1],
      p_value      = pt$tab$`Pr(>F)`[1],
      min_group_mean_dist = min(tapply(bd$distances, bd$group, mean)),
      max_group_mean_dist = max(tapply(bd$distances, bd$group, mean)),
      spread_ratio = max(tapply(bd$distances, bd$group, mean)) /
                     min(tapply(bd$distances, bd$group, mean)),
      group_means  = paste(sprintf("%s=%.2f", levels(bd$group),
                                   tapply(bd$distances, bd$group, mean)), collapse = "; ")
    )
  }
}

bind_rows(disp_results) %>% write_csv(here("Results", "rev_dispersion.csv"))

# --------------------------------------------------------------------------------------
# (3) Aitchison vs Bray-Curtis on the same resistome
# --------------------------------------------------------------------------------------
bc_results <- list()

for (res_type in names(mats_adj)) {
  a <- align(mats_adj[[res_type]])
  rel <- a$mat / rowSums(a$mat)

  for (metric in c("Aitchison", "Bray-Curtis")) {
    message(sprintf("[dist] %-12s | %s", res_type, metric))
    d <- if (metric == "Aitchison") dist(clr_mat(zeroreplace(a$mat, d = 0.65)))
         else vegdist(rel, method = "bray")
    bc_results[[length(bc_results) + 1]] <-
      adjusted_permanova(d, a$clim) %>%
      mutate(resistome = res_type, metric = metric, .before = 1)
  }
}

bind_rows(bc_results) %>% write_csv(here("Results", "rev_braycurtis.csv"))

message("done. R ", getRversion(), " | vegan ", packageVersion("vegan"),
        " | ", N_PERM, " permutations")
