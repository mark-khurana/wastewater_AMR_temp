# Sensitivity analysis: robustness of temperature-resistome association across averaging windows

library(tidyverse)
library(vegan)
library(flextable)
library(officer)
library(here)

analysis  <- readRDS(here("Datasets", "analysis_ready.rds"))
resistome <- readRDS(here("Datasets", "resistome_matrices.rds"))

windows <- c("T_sampling", "T_30d", "T_90d", "T_365d", "T_annual_mean")
window_labels <- c(
  T_sampling    = "Collection-day",
  T_30d         = "30-day mean",
  T_90d         = "90-day mean",
  T_365d        = "365-day mean",
  T_annual_mean = "Annual mean (year-matched)"
)

cor_data <- analysis %>% select(all_of(windows)) %>% drop_na()
cor_mat  <- cor(cor_data, use = "complete.obs")

run_window_permanova <- function(clr_mat, climate_df, temp_var, resistome_label) {
  shared <- intersect(clr_mat$genepid, climate_df$genepid)
  clr_sub  <- clr_mat %>% filter(genepid %in% shared) %>% arrange(genepid)
  clim_sub <- climate_df %>% filter(genepid %in% shared) %>% arrange(genepid)

  keep <- complete.cases(clim_sub[, c(temp_var, "Region",
                                       "gdp_pcap_ppp", "sanitation",
                                       "health_exp_gdp", "oop_health_exp",
                                       "immunization_dpt", "animal_amc_mgkg",
                                       "pop_density", "year")])
  clr_sub  <- clr_sub[keep, ]
  clim_sub <- clim_sub[keep, ]

  dist_mat <- dist(clr_sub %>% select(-genepid))

  clim_sub <- clim_sub %>%
    mutate(temp = .data[[temp_var]],
           log_gdp = log(gdp_pcap_ppp),
           log_pop_density = log(pop_density + 1),
           log_animal_amc = log(animal_amc_mgkg + 1),
           year_f = factor(year))

  result <- adonis2(dist_mat ~ temp + log_gdp + sanitation +
                      health_exp_gdp + oop_health_exp + immunization_dpt +
                      log_pop_density + log_animal_amc + Region + year_f,
                    data = clim_sub, permutations = 99, by = "margin")

  tibble(
    Resistome = resistome_label,
    Window = window_labels[temp_var],
    n = nrow(clr_sub),
    R2 = result["temp", "R2"],
    R2_pct = sprintf("%.3f", result["temp", "R2"] * 100),
    F_val = sprintf("%.2f", result["temp", "F"]),
    p = result["temp", "Pr(>F)"]
  )
}

results <- map_dfr(windows, function(w) {
  bind_rows(
    run_window_permanova(resistome$fg_cluster_clr, analysis, w, "Functional"),
    run_window_permanova(resistome$acq_cluster_clr, analysis, w, "Acquired")
  )
})

results <- results %>%
  mutate(p_fmt = ifelse(p < 0.001, "\u22640\u00B7001",
                   ifelse(p < 0.01, sprintf("%.3f", p),
                          sprintf("%.2f", p))))

saveRDS(list(permanova = results, cor_matrix = cor_mat, cor_n = nrow(cor_data)),
        here("Results", "sensitivity_temp_windows.rds"))

tab_perm <- results %>%
  select(Resistome, Window, n, R2_pct, F_val, p_fmt)

ft_perm <- flextable(tab_perm) %>%
  set_header_labels(R2_pct = "R\u00B2 (%)", F_val = "F", p_fmt = "p") %>%
  merge_v(j = "Resistome") %>%
  valign(j = "Resistome", valign = "top") %>%
  bold(j = "Resistome")

cor_long <- as.data.frame(cor_mat) %>%
  rownames_to_column("Window_1") %>%
  pivot_longer(-Window_1, names_to = "Window_2", values_to = "r") %>%
  filter(as.integer(factor(Window_1, levels = windows)) <
         as.integer(factor(Window_2, levels = windows))) %>%
  mutate(Window_1 = window_labels[Window_1],
         Window_2 = window_labels[Window_2],
         r = sprintf("%.3f", r))

ft_cor <- flextable(cor_long) %>%
  set_header_labels(Window_1 = "Window A", Window_2 = "Window B",
                    r = "Pearson r")

doc <- read_docx() %>%
  body_add_par("Supplementary Table 9. Sensitivity of the temperature\u2013resistome association to the choice of temperature averaging window.", style = "heading 2") %>%
  body_add_par("") %>%
  body_add_par("Panel A. Adjusted PERMANOVA results by temperature window", style = "heading 3") %>%
  body_add_flextable(ft_perm %>%
    add_footer_lines(paste0(
      "Each row shows the marginal (Type III) PERMANOVA R\u00B2 for temperature, ",
      "adjusting for GDP per capita (log), basic sanitation, health expenditure (% GDP), ",
      "out-of-pocket health expenditure, DPT immunisation, population density (log), ",
      "animal antimicrobial consumption (log), WHO Region, and collection year. ",
      "999 permutations. Aitchison distance (Euclidean on CLR-transformed abundances)."
    )) %>%
    fontsize(size = 9, part = "footer") %>%
    autofit()
  ) %>%
  body_add_par("") %>%
  body_add_par("Panel B. Pairwise Pearson correlations between temperature windows", style = "heading 3") %>%
  body_add_par(sprintf("(n = %d samples with all windows available)", nrow(cor_data)),
               style = "Normal") %>%
  body_add_flextable(ft_cor %>% autofit())

print(doc, target = here("Tables", "eTable9_temp_windows.docx"))
