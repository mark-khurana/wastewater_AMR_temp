# Merge climate metrics, confounders, and ARG counts; build CLR-transformed resistome matrices.
# Outputs: Datasets/analysis_ready.rds, Datasets/resistome_matrices.rds

library(tidyverse)
library(compositions)
library(here)

climate     <- readRDS(here("Datasets", "climate_metrics.rds"))
confounders <- readRDS(here("Datasets", "country_confounders.rds"))
counts      <- readRDS(here("Datasets", "clean_ARG_counts.rds"))

analysis <- climate %>%
  left_join(confounders, by = "genepid")

fg_counts  <- counts %>% filter(functional_amr == "True")
acq_counts <- counts %>% filter(resfinder == "True")

build_matrix <- function(df, group_col, value_col = "fragmentCountAln_adj") {
  df %>%
    filter(!is.na(.data[[group_col]])) %>%
    group_by(genepid, .data[[group_col]]) %>%
    summarise(abundance = sum(.data[[value_col]], na.rm = TRUE),
              .groups = "drop") %>%
    pivot_wider(names_from = all_of(group_col),
                values_from = abundance, values_fill = 0)
}

fg_class_mat  <- build_matrix(fg_counts, "class")
acq_class_mat <- build_matrix(acq_counts, "class")

fg_cluster_mat  <- build_matrix(fg_counts, "cluster_representative_98")
acq_cluster_mat <- build_matrix(acq_counts, "cluster_representative_98")

clr_transform <- function(mat_with_id) {
  ids <- mat_with_id$genepid
  mat <- mat_with_id %>% select(-genepid) %>% as.matrix()
  mat_replaced <- zeroreplace(mat, d = 0.65)
  clr_vals <- as.data.frame(clr(mat_replaced))
  clr_vals$genepid <- ids
  clr_vals
}

fg_cluster_clr  <- clr_transform(fg_cluster_mat)
acq_cluster_clr <- clr_transform(acq_cluster_mat)
fg_class_clr    <- clr_transform(fg_class_mat)
acq_class_clr   <- clr_transform(acq_class_mat)

cruz_loya <- tribble(
  ~class,                        ~stressor_group,
  "aminoglycoside",              "heat-shock",
  "folate_pathway_antagonist",   "heat-shock",
  "fluoroquinolone",             "cold-shock",
  "tetracycline",                "cold-shock",
  "macrolide",                   "cold-shock",
  "quinolone",                   "cold-shock",
  "beta_lactam",                 "neutral",
  "polymyxin",                   "neutral",
  "phenicol",                    "neutral",
  "glycopeptide",                "neutral"
)

saveRDS(analysis, here("Datasets", "analysis_ready.rds"))

resistome <- list(
  fg_cluster_clr  = fg_cluster_clr,
  acq_cluster_clr = acq_cluster_clr,
  fg_class_clr    = fg_class_clr,
  acq_class_clr   = acq_class_clr,
  fg_class_mat    = fg_class_mat,
  acq_class_mat   = acq_class_mat,
  cruz_loya       = cruz_loya
)
saveRDS(resistome, here("Datasets", "resistome_matrices.rds"))
