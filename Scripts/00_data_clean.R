# Clean Martiny et al. sample metadata and annotate ARG counts with drug classes.
# Outputs: Datasets/clean_samples.rds, Datasets/clean_ARG_counts.rds, Datasets/exclusion_log.csv

library(tidyverse)
library(readxl)
library(here)

sd1 <- read_excel(
  here("Datasets", "Martiny_et_al", "supplementary", "41467_2025_66070_MOESM4_ESM.xlsx"),
  sheet = "SD1 Sample metadata"
)

sd2 <- read_excel(
  here("Datasets", "Martiny_et_al", "supplementary", "41467_2025_66070_MOESM5_ESM.xlsx"),
  sheet = "SD2 ARG metadata"
)

counts <- read.csv(here("Datasets", "Martiny_et_al", "counts_clusters.csv"))

sd1 <- sd1 %>%
  mutate(
    collection_date_raw = as.numeric(collection_date),
    collection_date = case_when(
      is.na(collection_date_raw) ~ NA_Date_,
      collection_date_raw < 10000 ~ NA_Date_,
      TRUE ~ as.Date(collection_date_raw, origin = "1899-12-30")
    )
  ) %>%
  select(-collection_date_raw)

samples <- sd1 %>%
  arrange(genepid, complete_name) %>%
  distinct(genepid, .keep_all = TRUE) %>%
  select(genepid, complete_name, country, city, lat, lon,
         Region, year, collection_date)

exclusion_log <- tibble(
  genepid = numeric(), country = character(), city = character(),
  reason = character()
)

missing_coords <- samples %>%
  filter(is.na(lat) | is.na(lon)) %>%
  mutate(reason = "Missing GPS coordinates (lat or lon is NA)")

exclusion_log <- bind_rows(exclusion_log,
  missing_coords %>% select(genepid, country, city, reason))

samples <- samples %>% filter(!is.na(lat), !is.na(lon))

cluster_class <- sd2 %>%
  distinct(cluster_representative_98, class) %>%
  group_by(cluster_representative_98) %>%
  slice(1) %>%
  ungroup()

counts_annotated <- counts %>%
  left_join(cluster_class, by = "cluster_representative_98")

saveRDS(samples, here("Datasets", "clean_samples.rds"))
saveRDS(counts_annotated, here("Datasets", "clean_ARG_counts.rds"))
write_csv(exclusion_log, here("Datasets", "exclusion_log.csv"))
