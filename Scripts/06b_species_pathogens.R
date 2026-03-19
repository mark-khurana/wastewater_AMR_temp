# WHO priority pathogen species-level temperature associations (genus vs species contrast).
# Outputs: Results/layer3b_species_temp.csv, Results/layer3b_genus_vs_species.csv

library(here)
library(tidyverse)

analysis <- readRDS(here("Datasets", "analysis_ready.rds"))

motu_cache <- here("Datasets", "motu_species_cache.rds")

if (file.exists(motu_cache)) {
  motu <- readRDS(motu_cache)
} else {
  tmp_dir <- tempdir()
  unzip(here("Datasets", "Martiny_et_al", "motus_counts.zip"),
        files = "motus_counts/motus_agg_pad/mOTU_motus_pad_agg_counts.csv",
        exdir = tmp_dir, overwrite = TRUE)
  motu <- read_csv(file.path(tmp_dir, "motus_counts/motus_agg_pad/mOTU_motus_pad_agg_counts.csv"),
                    show_col_types = FALSE)
  saveRDS(motu, motu_cache)
}

targets <- tribble(
  ~display_name,           ~pattern,                      ~priority,    ~genus_name,
  "A. baumannii",          "Acinetobacter baumannii",     "Critical",   "Acinetobacter",
  "K. pneumoniae",         "Klebsiella pneumoniae",       "Critical",   "Klebsiella",
  "P. aeruginosa",         "Pseudomonas aeruginosa",      "Critical",   "Pseudomonas",
  "E. coli",               "Escherichia coli",            "Critical",   "Escherichia",
  "E. faecium",            "Enterococcus faecium",        "High",       "Enterococcus",
  "E. faecalis",           "Enterococcus faecalis",       "High",       "Enterococcus",
  "S. aureus",             "Staphylococcus aureus",       "High",       "Staphylococcus",
  "Salmonella enterica",   "Salmonella enterica",         "High",       "Salmonella",
  "N. gonorrhoeae",        "Neisseria gonorrhoeae",       "High",       "Neisseria",
  "C. jejuni",             "Campylobacter jejuni",        "High",       "Campylobacter"
)

all_cols <- names(motu)[-1]

species_matches <- map(targets$pattern, function(pat) {
  all_cols[str_detect(all_cols, fixed(pat))]
})
names(species_matches) <- targets$display_name

genus_matches <- map(unique(targets$genus_name), function(g) {
  all_cols[str_detect(all_cols, paste0("^", g))]
})
names(genus_matches) <- unique(targets$genus_name)

motu_ids <- as.character(as.integer(motu$genepid))
motu_mat <- as.matrix(motu[, -1])
total_reads <- rowSums(motu_mat)

motu_ra <- motu_mat / total_reads

merge_df <- tibble(genepid = motu_ids, total_reads = total_reads) %>%
  inner_join(analysis %>% mutate(genepid = as.character(genepid)) %>%
               select(genepid, T_30d, Region),
             by = "genepid") %>%
  filter(!is.na(T_30d))

idx <- match(merge_df$genepid, motu_ids)
motu_ra_matched <- motu_ra[idx, ]

species_cors <- map_dfr(seq_len(nrow(targets)), function(i) {
  cols <- species_matches[[targets$display_name[i]]]
  if (length(cols) == 0) return(NULL)

  if (length(cols) == 1) {
    ra <- motu_ra_matched[, cols]
  } else {
    ra <- rowSums(motu_ra_matched[, cols, drop = FALSE])
  }

  prev <- mean(ra > 0)
  if (sum(ra > 0) < 30) return(NULL)

  ct <- cor.test(ra, merge_df$T_30d, method = "spearman")

  tibble(
    name = targets$display_name[i],
    level = "species",
    priority = targets$priority[i],
    genus = targets$genus_name[i],
    rho = ct$estimate,
    p_value = ct$p.value,
    prevalence = prev,
    n = length(ra)
  )
})

genus_cors <- map_dfr(names(genus_matches), function(g) {
  cols <- genus_matches[[g]]
  if (length(cols) == 0) return(NULL)

  ra <- rowSums(motu_ra_matched[, cols, drop = FALSE])
  prev <- mean(ra > 0)

  ct <- cor.test(ra, merge_df$T_30d, method = "spearman")

  tibble(
    name = paste(g, "(genus)"),
    level = "genus",
    priority = NA_character_,
    genus = g,
    rho = ct$estimate,
    p_value = ct$p.value,
    prevalence = prev,
    n = length(ra)
  )
})

all_cors <- bind_rows(species_cors, genus_cors) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  arrange(genus, level)

contrast <- all_cors %>%
  select(name, level, genus, rho, p_adj) %>%
  pivot_wider(names_from = level, values_from = c(rho, p_adj),
              names_glue = "{.value}_{level}") %>%
  filter(!is.na(rho_species)) %>%
  mutate(direction_match = sign(rho_species) == sign(rho_genus),
         rho_diff = rho_species - rho_genus)

genus_species_compare <- all_cors %>%
  filter(genus %in% (all_cors %>% filter(level == "species") %>% pull(genus))) %>%
  select(name, level, genus, rho, p_adj, prevalence)

write_csv(all_cors, here("Results", "layer3b_species_temp.csv"))
write_csv(genus_species_compare, here("Results", "layer3b_genus_vs_species.csv"))
