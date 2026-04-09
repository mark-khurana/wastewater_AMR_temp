# =============================================================================
# 10b_tables_docx.R — Generate all tables as Word (.docx) files
# =============================================================================

library(tidyverse)
library(flextable)
library(officer)

base_dir <- "/Users/bxc331/Desktop/AMR_open_project"
out_dir  <- file.path(base_dir, "Datasets")
res_dir  <- file.path(base_dir, "Results")
tab_dir  <- file.path(base_dir, "Tables")

cat("=== 10b_tables_docx.R ===\n\n")

# Load data
analysis <- readRDS(file.path(out_dir, "analysis_ready.rds"))
resistome <- readRDS(file.path(out_dir, "resistome_matrices.rds"))
resistome_ids <- as.character(resistome$fg_cluster_clr$genepid)

conf_vars <- c("T_30d", "gdp_pcap_ppp", "sanitation",
               "health_exp_gdp", "oop_health_exp", "immunization_dpt",
               "animal_amc_mgkg", "pop_density", "year", "Region")

cc <- analysis %>%
  mutate(genepid = as.character(genepid)) %>%
  filter(complete.cases(across(all_of(conf_vars))),
         genepid %in% resistome_ids)

fmt_range <- function(lo, hi, dig = 1) {
  sprintf("%.*f to %.*f", dig, lo, dig, hi)
}

# =============================================================================
# TABLE 1: Sampling design (MAIN)
# =============================================================================
cat("--- Table 1 ---\n")

regions <- sort(unique(cc$Region))

tab1 <- bind_rows(
  map_dfr(regions, function(reg) {
    d <- cc %>% filter(Region == reg)
    tibble(
      Region = reg,
      Samples = nrow(d),
      Cities = n_distinct(d$city),
      Countries = n_distinct(d$country),
      temp_range = fmt_range(min(d$T_30d),
                             max(d$T_30d))
    )
  }),
  tibble(
    Region = "Total",
    Samples = nrow(cc),
    Cities = n_distinct(cc$city),
    Countries = n_distinct(cc$country),
    temp_range = fmt_range(min(cc$T_30d),
                           max(cc$T_30d))
  )
)

year_dist <- cc %>% count(year) %>% mutate(label = sprintf("%d (n=%d)", year, n))
year_note <- paste(year_dist$label, collapse = "; ")

ft1 <- flextable(tab1) %>%
  set_header_labels(temp_range = "Temperature range (\u00B0C)") %>%
  bold(i = nrow(tab1), bold = TRUE) %>%
  hline(i = nrow(tab1) - 1, border = fp_border(color = "black")) %>%
  set_caption("Table 1. Characteristics of 1,066 sewage metagenomes from 104 countries included in the cross-sectional analysis, by WHO region.") %>%
  add_footer_lines(paste0(
    "Samples were collected between 2016 and 2021 (", year_note, "). ",
    "Temperature: 30-day mean air temperature preceding sample collection (NASA POWER T2M). ",
    "Country-level confounders adjusted for in the model: GDP per capita (PPP), basic sanitation coverage, ",
    "current health expenditure (% GDP), out-of-pocket health expenditure, DPT (diphtheria\u2013tetanus\u2013pertussis) ",
    "immunisation coverage (proxy for healthcare delivery capacity), population density, ",
    "and animal antimicrobial consumption. See Supplementary Table 5 for sources and rationale."
  )) %>%
  fontsize(size = 9, part = "footer") %>%
  autofit()

save_as_docx(ft1, path = file.path(tab_dir, "Table1_descriptive.docx"),
             pr_section = prop_section(page_size = page_size(orient = "landscape")))
cat("  Saved Table1_descriptive.docx\n")

# =============================================================================
# eTABLE 1: PERMANOVA results (SUPPLEMENTARY)
# =============================================================================
cat("--- eTable 1 ---\n")

perm <- readRDS(file.path(res_dir, "layer2_partial_permanova.rds"))

format_perm <- function(perm_obj, label) {
  df <- as.data.frame(perm_obj)
  df$Term <- rownames(df)
  df %>%
    filter(Term != "Total") %>%
    mutate(
      Resistome = label,
      R2_pct = round(R2 * 100, 3),
      F_val = round(`F`, 2),
      p = ifelse(is.na(`Pr(>F)`), "\u2014",
            ifelse(`Pr(>F)` < 0.001, "<0.001",
              ifelse(`Pr(>F)` < 0.01, sprintf("%.3f", `Pr(>F)`),
                     sprintf("%.2f", `Pr(>F)`))))
    ) %>%
    select(Resistome, Term, df = Df, R2_pct, F_val, p)
}

perm_all <- bind_rows(
  format_perm(perm$perm_fg$full, "Functional"),
  format_perm(perm$perm_acq$full, "Acquired")
)

# Clean term names
perm_all <- perm_all %>%
  mutate(Term = case_when(
    Term == "T_30d" ~ "30-day mean temperature",
    Term == "log_gdp" ~ "GDP per capita (log)",
    Term == "sanitation" ~ "Basic sanitation (%)",
    Term == "health_exp_gdp" ~ "Health expenditure (% GDP)",
    Term == "oop_health_exp" ~ "Out-of-pocket health exp. (%)",
    Term == "immunization_dpt" ~ "DPT immunisation (%)",
    Term == "log_pop_density" ~ "Population density (log)",
    Term == "log_animal_amc" ~ "Animal AMC (log mg/kg)",
    Term == "Region" ~ "WHO Region",
    Term == "year_f" ~ "Collection year",
    Term == "Residual" ~ "Residual",
    TRUE ~ Term
  ))

ft2 <- flextable(perm_all) %>%
  set_header_labels(R2_pct = "R\u00B2 (%)", F_val = "F") %>%
  merge_v(j = "Resistome") %>%
  valign(j = "Resistome", valign = "top") %>%
  bold(j = "Resistome") %>%
  set_caption("Supplementary Table 1. PERMANOVA results: variance in resistome composition explained by temperature and covariates (marginal, Type III). n=1,066.") %>%
  add_footer_lines(paste0(
    "PERMANOVA on Aitchison distance (Euclidean on CLR-transformed gene cluster abundances). ",
    "999 permutations, marginal (Type III) terms. Both resistome types use the same samples and covariates."
  )) %>%
  fontsize(size = 9, part = "footer") %>%
  autofit()

save_as_docx(ft2, path = file.path(tab_dir, "eTable1_permanova.docx"))
cat("  Saved eTable1_permanova.docx\n")

# =============================================================================
# eTABLE 2: Top genus correlations (SUPPLEMENTARY)
# =============================================================================
cat("--- eTable 2 ---\n")

genus_corr <- read.csv(file.path(res_dir, "layer3_genus_temp_correlations.csv"))

top_warm <- genus_corr %>% arrange(desc(rho)) %>% head(15)
top_cool <- genus_corr %>% arrange(rho) %>% head(15)

top_genera <- bind_rows(
  top_warm %>% mutate(Direction = "Warm-associated"),
  top_cool %>% mutate(Direction = "Cool-associated")
) %>%
  mutate(
    Genus = str_replace(genus, " ext_mOTU.*", " (unclassified)"),
    Genus = str_replace(Genus, "aceae gen\\. .*", "aceae (unclassified)"),
    Genus = str_replace(Genus, "ales gen\\. .*", "ales (unclassified)"),
    rho_val = round(rho, 3),
    p_adj_fmt = ifelse(p_adj < 0.001, "<0.001", sprintf("%.3f", p_adj)),
    n = n
  ) %>%
  select(Direction, Genus, rho_val, p_adj_fmt, n)

ft3 <- flextable(top_genera) %>%
  set_header_labels(rho_val = "\u03C1", p_adj_fmt = "p (adj.)") %>%
  merge_v(j = "Direction") %>%
  valign(j = "Direction", valign = "top") %>%
  italic(j = "Direction") %>%
  set_caption(paste0(
    "Supplementary Table 2. Top 15 warm-associated and 15 cool-associated genera ",
    "by Spearman correlation with 30-day mean temperature (unadjusted). ",
    "Full results for all ", nrow(genus_corr), " genera available as Supplementary Data."
  )) %>%
  add_footer_lines(paste0(
    "Unadjusted Spearman correlations between genus-level relative abundance and 30-day mean temperature. ",
    "p-values Bonferroni-adjusted for ", nrow(genus_corr), " tests."
  )) %>%
  fontsize(size = 9, part = "footer") %>%
  autofit()

save_as_docx(ft3, path = file.path(tab_dir, "eTable2_genus_correlations.docx"))
cat("  Saved eTable2_genus_correlations.docx\n")

# =============================================================================
# eTABLE 3: Species correlations (SUPPLEMENTARY)
# =============================================================================
cat("--- eTable 3 ---\n")

sp <- read_csv(file.path(res_dir, "adjusted_species_temp.csv"), show_col_types = FALSE)

fmt_p <- function(p) {
  ifelse(is.na(p), "\u2014",
    ifelse(p < 0.001, "<0.001",
      ifelse(p < 0.01, sprintf("%.3f", p),
             sprintf("%.2f", p))))
}

sp_tab <- sp %>%
  mutate(
    Genus = genus,
    Taxon = display,
    rho_raw_fmt = round(rho_raw, 3),
    p_raw_fmt = fmt_p(p_raw),
    rho_adj_fmt = ifelse(is.na(rho_adj), "\u2014", sprintf("%.3f", rho_adj)),
    p_adj_fmt = fmt_p(p_adj)
  ) %>%
  select(Genus, Taxon, rho_raw_fmt, p_raw_fmt, rho_adj_fmt, p_adj_fmt)

ft4 <- flextable(sp_tab) %>%
  set_header_labels(rho_raw_fmt = "\u03C1 (raw)", p_raw_fmt = "p (raw)",
                    rho_adj_fmt = "\u03C1 (adj.)", p_adj_fmt = "p (adj.)") %>%
  merge_v(j = "Genus") %>%
  valign(j = "Genus", valign = "top") %>%
  italic(j = "Genus") %>%
  set_caption("Supplementary Table 3. Species-level temperature associations: unadjusted and covariate-adjusted Spearman correlations (n=1,066).") %>%
  add_footer_lines(paste0(
    "Adjustment: GDP (log), basic sanitation, health expenditure (% GDP), out-of-pocket health expenditure, ",
    "DPT immunisation, population density (log), animal antimicrobial consumption (log mg/kg), WHO Region, ",
    "and collection year. Partial correlations computed by residualising both abundance and temperature on all covariates."
  )) %>%
  fontsize(size = 9, part = "footer") %>%
  autofit()

save_as_docx(ft4, path = file.path(tab_dir, "eTable3_species_correlations.docx"))
cat("  Saved eTable3_species_correlations.docx\n")

# =============================================================================
# eTABLE 4: Sensitivity analyses (SUPPLEMENTARY)
# =============================================================================
cat("--- eTable 4 ---\n")

sens <- tibble(
  Analysis = c(
    "Main model (30-day mean temp.)",
    "Without collection year",
    "AMC subset only",
    "Collection-day temperature",
    "Within-city temporal (30-day)",
    "Within-city (same-day temp.)",
    "Month-adjusted temporal",
    "Humidity adjustment"
  ),
  n = c("1,066", "1,066", "522", "1,066", "162", "162", "162", "1,064"),
  temp_r2 = c(
    "FG: 0.51%, Acq: 0.56%",
    "FG: 0.72%, Acq: 0.88%",
    "FG: 0.53%, Acq: 0.61%",
    "\u2014",
    "6.3% (total)",
    "\u2014",
    "\u2014",
    "FG: 0.52%, Acq: 0.58%"
  ),
  key_finding = c(
    "Both p=0.001; ~68% mediated via bacteriome",
    "Year explains ~2.7% independently",
    "Survives with similar magnitude",
    "Stronger for A. baumannii (\u03C1=+0.38)",
    "73% mediated via bacteriome",
    "Consistent direction and magnitude",
    "Only Enterococcus retains significance",
    "Temperature R\u00B2 unchanged after RH adjustment (Suppl. Table 7)"
  )
)

ft5 <- flextable(sens) %>%
  set_header_labels(temp_r2 = "Temperature R\u00B2",
                    key_finding = "Key finding") %>%
  set_caption("Supplementary Table 4. Summary of sensitivity analyses for the temperature\u2013resistome association.") %>%
  add_footer_lines(paste0(
    "FG: functional genomic resistome. Acq: acquired resistome. ",
    "R\u00B2: marginal PERMANOVA R\u00B2 for temperature, adjusting for all covariates. ",
    "AMC subset: restricted to samples with animal antimicrobial consumption data. ",
    "Collection-day temperature: NASA POWER air temperature on sampling date. ",
    "Within-city temporal: biweekly sampling in 5 European cities (2019\u20132021), ",
    "with city and calendar month as fixed effects."
  )) %>%
  fontsize(size = 9, part = "footer") %>%
  autofit()

save_as_docx(ft5, path = file.path(tab_dir, "eTable4_sensitivity.docx"))
cat("  Saved eTable4_sensitivity.docx\n")

# =============================================================================
# eTABLE 5: Data sources (SUPPLEMENTARY)
# =============================================================================
cat("--- eTable 5 ---\n")

sources <- tibble(
  Category = c(
    "Outcome",
    "Exposure",
    "Exposure",
    "Confounder",
    "Confounder",
    "Confounder",
    "Confounder",
    "Confounder",
    "Confounder",
    "Confounder",
    "Confounder",
    "Replication",
    "Reference"
  ),
  Variable = c(
    "Sewage resistome and bacteriome",
    "30-day mean air temperature (primary)",
    "Collection-day air temperature (sensitivity)",
    "GDP per capita, PPP",
    "Basic sanitation services",
    "Current health expenditure",
    "Out-of-pocket health expenditure",
    "DPT immunisation coverage",
    "Population density",
    "Animal antimicrobial consumption",
    "WHO Region",
    "5-city biweekly sewage metagenomes",
    "WHO Bacterial Priority Pathogen List"
  ),
  Rationale = c(
    "Metagenomic profiling of ARGs and bacterial taxa from urban sewage",
    "Thermal environment preceding sample collection; hypothesised driver of bacterial community composition",
    "Acute thermal conditions at time of sampling",
    "Overall economic development; confounds temperature via latitude\u2013wealth gradient",
    "Water, sanitation, and hygiene infrastructure; affects faecal pathogen transmission",
    "Healthcare system investment; associated with AMR stewardship capacity",
    "Proxy for healthcare access barriers; higher OOP associated with self-medication and over-the-counter antibiotic use",
    "Proxy for healthcare delivery capacity (cold chain, community health infrastructure)",
    "Urbanisation and crowding; affects transmission dynamics",
    "Veterinary antibiotic use in livestock; contributes to environmental AMR burden",
    "Broad geographic and cultural confounding",
    "Independent within-city temporal replication of cross-sectional findings",
    "Species selection for pathogen-level analyses"
  ),
  Source = c(
    "Martiny et al. 2025 (Nat Commun)",
    "NASA POWER v2.0 (T2M)",
    "NASA POWER v2.0 (T2M)",
    "World Bank WDI",
    "World Bank WDI",
    "World Bank WDI",
    "World Bank WDI",
    "World Bank WDI",
    "World Bank WDI",
    "Mulchandani et al. 2023 via OWID",
    "WHO",
    "Becsei et al. 2026",
    "WHO 2024"
  ),
  year_range = c(
    "2016\u20132021", "2016\u20132021", "2016\u20132021",
    "2018", "2017\u20132020", "2018\u20132019", "2018\u20132019",
    "2018\u20132020", "2018\u20132020", "2017\u20132020",
    "\u2014", "2019\u20132021", "\u2014"
  ),
  Coverage = c(
    "112 countries", "98%", "98%",
    "98%", "97%", "97%", "97%", "97%", "98%",
    "109 countries", "100%", "5 cities", "\u2014"
  )
)

ft6 <- flextable(sources) %>%
  set_header_labels(year_range = "Year(s)") %>%
  merge_v(j = "Category") %>%
  valign(j = "Category", valign = "top") %>%
  bold(j = "Category") %>%
  set_caption("Supplementary Table 5. Data sources, covariates, and rationale for inclusion in the analytical model.") %>%
  add_footer_lines(paste0(
    "Coverage: percentage of the 1,066-sample analysis subset with non-missing data, or number of countries/cities. ",
    "All country-level covariates were matched to sample country of origin and are time-invariant within the analysis ",
    "(single cross-sectional value per country, not year-matched). ",
    "WDI: World Development Indicators. OWID: Our World in Data. PPP: purchasing power parity. ",
    "DPT: diphtheria\u2013tetanus\u2013pertussis. OOP: out-of-pocket."
  )) %>%
  fontsize(size = 8, part = "body") %>%
  fontsize(size = 8, part = "footer") %>%
  width(j = "Rationale", width = 3) %>%
  autofit()

save_as_docx(ft6, path = file.path(tab_dir, "eTable5_datasources.docx"),
             pr_section = prop_section(page_size = page_size(orient = "landscape")))
cat("  Saved eTable5_datasources.docx\n")

# =============================================================================
# eTABLE 7: Humidity sensitivity — adjusted PERMANOVA (SUPPLEMENTARY)
# =============================================================================
cat("--- eTable 7 ---\n")

hum <- readRDS(file.path(res_dir, "sensitivity_humidity.rds"))

format_perm_hum <- function(perm_obj, model_label, resistome_label) {
  df <- as.data.frame(perm_obj)
  df$Term <- rownames(df)
  df %>%
    filter(Term != "Total") %>%
    mutate(
      Resistome = resistome_label,
      Model = model_label,
      R2_pct = sprintf("%.3f", R2 * 100),
      F_val = ifelse(is.na(`F`), "\u2014", sprintf("%.2f", `F`)),
      p = ifelse(is.na(`Pr(>F)`), "\u2014",
            ifelse(`Pr(>F)` < 0.001, "\u22640\u00B7001",
              ifelse(`Pr(>F)` < 0.01, sprintf("%.3f", `Pr(>F)`),
                     sprintf("%.2f", `Pr(>F)`))))
    ) %>%
    select(Resistome, Model, Term, Df, R2_pct, F_val, p)
}

clean_hum_terms <- function(df) {
  df %>% mutate(Term = case_when(
    Term == "T_30d" ~ "30-day mean temperature",
    Term == "RH2M_30d" ~ "30-day mean relative humidity",
    Term == "log_gdp" ~ "GDP per capita (log)",
    Term == "sanitation" ~ "Basic sanitation (%)",
    Term == "health_exp_gdp" ~ "Health expenditure (% GDP)",
    Term == "oop_health_exp" ~ "Out-of-pocket health exp. (%)",
    Term == "immunization_dpt" ~ "DPT immunisation (%)",
    Term == "log_pop_density" ~ "Population density (log)",
    Term == "log_animal_amc" ~ "Animal AMC (log mg/kg)",
    Term == "Region" ~ "WHO Region",
    Term == "year_f" ~ "Collection year",
    Term == "Residual" ~ "Residual",
    TRUE ~ Term
  ))
}

hum_perm_all <- bind_rows(
  format_perm_hum(hum$fg$perm_ref,   "Without humidity", "Functional"),
  format_perm_hum(hum$fg$perm_main,  "With humidity",    "Functional"),
  format_perm_hum(hum$acq$perm_ref,  "Without humidity", "Acquired"),
  format_perm_hum(hum$acq$perm_main, "With humidity",    "Acquired")
) %>% clean_hum_terms()

ft7 <- flextable(hum_perm_all) %>%
  set_header_labels(R2_pct = "R\u00B2 (%)", F_val = "F", Df = "df") %>%
  merge_v(j = c("Resistome", "Model")) %>%
  valign(j = c("Resistome", "Model"), valign = "top") %>%
  bold(j = "Resistome") %>%
  italic(j = "Model") %>%
  set_caption(paste0(
    "Supplementary Table 7. Sensitivity analysis: PERMANOVA results with and without ",
    "adjustment for relative humidity (n=", hum$fg$n, ")."
  )) %>%
  add_footer_lines(paste0(
    "Without humidity: 30-day mean temperature + 7 socioeconomic confounders + WHO Region + collection year. ",
    "With humidity: same model plus 30-day mean relative humidity (NASA POWER RH2M). ",
    "PERMANOVA on Aitchison distance (Euclidean on CLR-transformed gene cluster abundances), ",
    "999 permutations, marginal (Type III) tests. ",
    "Temperature and relative humidity are moderately correlated in this dataset ",
    "(Pearson r = \u22120\u00B771, variance inflation factor = 2\u00B70). ",
    "The temperature R\u00B2 is essentially unchanged after humidity adjustment ",
    "(functional: 0\u00B7512% \u2192 0\u00B7515%; acquired: 0\u00B7564% \u2192 0\u00B7579%), ",
    "confirming that the temperature\u2013resistome association is not confounded by moisture availability."
  )) %>%
  fontsize(size = 8, part = "body") %>%
  fontsize(size = 8, part = "footer") %>%
  autofit()

save_as_docx(ft7, path = file.path(tab_dir, "eTable7_humidity_sensitivity.docx"),
             pr_section = prop_section(page_size = page_size(orient = "landscape")))
cat("  Saved eTable7_humidity_sensitivity.docx\n")

# =============================================================================
# eTABLE 8: Grid distances (already generated separately)
# =============================================================================
# eTable8_grid_distances.docx is generated by Scripts/eTable8_grid_distances.R

# =============================================================================
# eTABLE 9: Temperature window sensitivity (SUPPLEMENTARY)
# =============================================================================
cat("--- eTable 9 ---\n")

tw <- readRDS(file.path(res_dir, "sensitivity_temp_windows.rds"))

tw_tab <- tw$permanova %>%
  select(Resistome, Window, n, R2_pct, F_val, p_fmt)

ft9_perm <- flextable(tw_tab) %>%
  set_header_labels(R2_pct = "R\u00B2 (%)", F_val = "F", p_fmt = "p") %>%
  merge_v(j = "Resistome") %>%
  valign(j = "Resistome", valign = "top") %>%
  bold(j = "Resistome") %>%
  bold(i = ~ Window == "30-day mean") %>%
  fontsize(size = 9, part = "all") %>%
  font(fontname = "Times New Roman", part = "all") %>%
  padding(padding.top = 1, padding.bottom = 1, part = "body") %>%
  border_remove() %>%
  hline_top(border = fp_border(color = "black", width = 1), part = "header") %>%
  hline_bottom(border = fp_border(color = "black", width = 1), part = "header") %>%
  hline_bottom(border = fp_border(color = "black", width = 1), part = "body") %>%
  align(j = 3:6, align = "center", part = "all") %>%
  autofit()

# Correlation panel
cor_long <- as.data.frame(tw$cor_matrix) %>%
  rownames_to_column("Window_1") %>%
  pivot_longer(-Window_1, names_to = "Window_2", values_to = "r") %>%
  filter(as.integer(factor(Window_1, levels = c("T_sampling","T_30d","T_90d","T_365d","T_annual_mean"))) <
         as.integer(factor(Window_2, levels = c("T_sampling","T_30d","T_90d","T_365d","T_annual_mean"))))

wl <- c(T_sampling = "Collection-day", T_30d = "30-day mean",
        T_90d = "90-day mean", T_365d = "365-day mean",
        T_annual_mean = "Annual mean")
cor_long <- cor_long %>%
  mutate(Window_1 = wl[Window_1], Window_2 = wl[Window_2],
         r = sprintf("%.3f", r))

ft9_cor <- flextable(cor_long) %>%
  set_header_labels(Window_1 = "Window A", Window_2 = "Window B", r = "Pearson r") %>%
  fontsize(size = 9, part = "all") %>%
  font(fontname = "Times New Roman", part = "all") %>%
  border_remove() %>%
  hline_top(border = fp_border(color = "black", width = 1), part = "header") %>%
  hline_bottom(border = fp_border(color = "black", width = 1), part = "header") %>%
  hline_bottom(border = fp_border(color = "black", width = 1), part = "body") %>%
  align(j = 3, align = "center", part = "all") %>%
  autofit()

doc9 <- read_docx() %>%
  body_add_par("Supplementary Table 9. Sensitivity of the temperature\u2013resistome association to the choice of temperature averaging window.", style = "heading 2") %>%
  body_add_par("") %>%
  body_add_par("Panel A. Adjusted PERMANOVA results by temperature window", style = "heading 3") %>%
  body_add_flextable(ft9_perm %>%
    add_footer_lines(paste0(
      "Each row shows the marginal (Type III) PERMANOVA R\u00B2 for the temperature variable, ",
      "adjusting for GDP per capita (log), basic sanitation, health expenditure (% GDP), ",
      "out-of-pocket health expenditure, DPT immunisation, population density (log), ",
      "animal antimicrobial consumption (log), WHO Region, and collection year. ",
      "99 permutations; R\u00B2 is deterministic. Primary analysis (30-day mean, bold) shown for reference."
    )) %>%
    fontsize(size = 8, part = "footer") %>%
    font(fontname = "Times New Roman", part = "footer") %>%
    autofit()
  ) %>%
  body_add_par("") %>%
  body_add_par("Panel B. Pairwise Pearson correlations between temperature windows", style = "heading 3") %>%
  body_add_par(sprintf("(n=%d samples with all windows available)", tw$cor_n), style = "Normal") %>%
  body_add_flextable(ft9_cor)

print(doc9, target = file.path(tab_dir, "eTable9_temp_windows.docx"))
cat("  Saved eTable9_temp_windows.docx\n")

# =============================================================================
# eTABLE 10: PCoA mediation sensitivity (SUPPLEMENTARY)
# =============================================================================
cat("--- eTable 10 ---\n")

sens_pcoa <- readRDS(file.path(res_dir, "pcoa_mediation_sensitivity.rds"))
med_data  <- readRDS(file.path(res_dir, "layer3_mediation.rds"))
eig <- med_data$pcoa_eig
pos_eig <- eig[eig > 0]
pct_var <- 100 * pos_eig / sum(pos_eig)

pcoa_tab <- sens_pcoa %>%
  mutate(
    Resistome = resistome,
    k = as.integer(k),
    cum_var = sapply(k, function(kk) sprintf("%.1f", sum(pct_var[1:kk]))),
    total_r2 = sprintf("%.2f", r2_total * 100),
    direct_r2 = sprintf("%.2f", r2_direct * 100),
    mediated_r2 = sprintf("%.2f", r2_indirect * 100),
    prop_med = sprintf("%.1f", pct_mediated)
  ) %>%
  select(Resistome, k, cum_var, total_r2, direct_r2, mediated_r2, prop_med)

ft10 <- flextable(pcoa_tab) %>%
  set_header_labels(k = "PCoA axes (k)",
                    cum_var = "Cumulative\nvariance (%)",
                    total_r2 = "Total R\u00B2 (%)",
                    direct_r2 = "Direct R\u00B2 (%)",
                    mediated_r2 = "Mediated R\u00B2 (%)",
                    prop_med = "Proportion\nmediated (%)") %>%
  merge_v(j = "Resistome") %>%
  valign(j = "Resistome", valign = "top") %>%
  bold(j = "Resistome") %>%
  bold(i = ~ k == 5) %>%
  fontsize(size = 9, part = "all") %>%
  font(fontname = "Times New Roman", part = "all") %>%
  padding(padding.top = 1, padding.bottom = 1, part = "body") %>%
  border_remove() %>%
  hline_top(border = fp_border(color = "black", width = 1), part = "header") %>%
  hline_bottom(border = fp_border(color = "black", width = 1), part = "header") %>%
  hline_bottom(border = fp_border(color = "black", width = 1), part = "body") %>%
  align(j = 2:7, align = "center", part = "all") %>%
  add_footer_lines(paste0(
    "Sensitivity of mediation results to the number of bacteriome principal coordinates retained. ",
    "Total R\u00B2: marginal PERMANOVA R\u00B2 for 30-day mean temperature adjusting for WHO region and sampling year. ",
    "Direct R\u00B2: temperature R\u00B2 after additionally conditioning on k bacteriome principal coordinates. ",
    "Proportion mediated = (Total \u2212 Direct) / Total \u00D7 100. ",
    "Primary analysis (k=5, bold) shown alongside 3 and 10 axes. ",
    "Cumulative variance: proportion of bacteriome Bray\u2013Curtis variation captured by retained axes. ",
    "99 permutations; R\u00B2 is deterministic."
  )) %>%
  fontsize(size = 8, part = "footer") %>%
  font(fontname = "Times New Roman", part = "footer") %>%
  autofit()

doc10 <- read_docx() %>%
  body_add_par("Supplementary Table 10. Sensitivity of mediation results to the number of bacteriome principal coordinates.", style = "heading 2") %>%
  body_add_par("") %>%
  body_add_flextable(ft10)

print(doc10, target = file.path(tab_dir, "eTable10_pcoa_mediation_sensitivity.docx"))
cat("  Saved eTable10_pcoa_mediation_sensitivity.docx\n")

# =============================================================================
# eTABLE 11: Complete-case analysis (SUPPLEMENTARY)
# =============================================================================
cat("--- eTable 11 ---\n")

resistome_ids_all <- as.character(resistome$fg_cluster_clr$genepid)
analysis$has_resistome <- as.character(analysis$genepid) %in% resistome_ids_all
covs_cc <- c("T_30d", "gdp_pcap_ppp", "sanitation", "health_exp_gdp",
             "oop_health_exp", "immunization_dpt", "animal_amc_mgkg", "pop_density")
a_pool <- analysis %>% filter(has_resistome)
a_pool$included <- complete.cases(a_pool[, covs_cc])

# Panel A: Flow
flow <- tibble(
  col1 = c("Total samples in dataset",
           "With matched resistome data",
           "With complete covariate data (final)"),
  col2 = c(as.character(nrow(analysis)),
           sprintf("%d (%d excluded)", nrow(a_pool), nrow(analysis) - nrow(a_pool)),
           sprintf("%d (%d excluded)", sum(a_pool$included), sum(!a_pool$included)))
)

# Panel B: By WHO region
region_tab <- a_pool %>%
  group_by(Region) %>%
  summarise(total = n(), retained = sum(included), excluded = total - retained, .groups = "drop") %>%
  mutate(col1 = Region,
         col2 = sprintf("%d/%d (%.0f%%)", retained, total, 100 * excluded / total)) %>%
  select(col1, col2)

# Fisher test
region_mat <- a_pool %>% group_by(Region) %>%
  summarise(incl = sum(included), excl = sum(!included), .groups = "drop")
fisher_region <- fisher.test(region_mat %>% select(incl, excl) %>% as.matrix(),
                             simulate.p.value = TRUE, B = 10000)

# Panel C: By income
income_tab <- a_pool %>%
  filter(!is.na(income)) %>%
  group_by(income) %>%
  summarise(total = n(), retained = sum(included), excluded = total - retained, .groups = "drop") %>%
  mutate(col1 = income,
         col2 = sprintf("%d/%d (%.0f%%)", retained, total, 100 * excluded / total)) %>%
  select(col1, col2)

income_mat <- a_pool %>% filter(!is.na(income)) %>% group_by(income) %>%
  summarise(incl = sum(included), excl = sum(!included), .groups = "drop")
fisher_income <- fisher.test(income_mat %>% select(incl, excl) %>% as.matrix(),
                             simulate.p.value = TRUE, B = 10000)

# Panel D: Covariates
cov_labels <- c(
  T_30d = "30-day mean temperature",
  gdp_pcap_ppp = "GDP per capita (PPP)",
  sanitation = "Basic sanitation access",
  health_exp_gdp = "Health expenditure (% GDP)",
  oop_health_exp = "Out-of-pocket health expenditure",
  immunization_dpt = "DPT immunisation coverage",
  animal_amc_mgkg = "Animal antimicrobial consumption",
  pop_density = "Population density"
)

cov_tab <- tibble(
  col1 = cov_labels[covs_cc],
  col2 = sprintf("%d (%.1f%%)",
    sapply(covs_cc, function(v) sum(is.na(a_pool[[v]]))),
    100 * sapply(covs_cc, function(v) sum(is.na(a_pool[[v]]))) / nrow(a_pool))
)

# Assemble
blank <- tibble(col1 = "", col2 = "")
combined <- bind_rows(
  tibble(col1 = "Sample derivation", col2 = ""),
  flow,
  blank,
  tibble(col1 = "Exclusions by WHO region", col2 = ""),
  region_tab,
  blank,
  tibble(col1 = "Exclusions by World Bank income group", col2 = ""),
  income_tab,
  blank,
  tibble(col1 = "Missing data by covariate", col2 = ""),
  cov_tab
)

header_rows <- which(combined$col1 %in% c("Sample derivation",
  "Exclusions by WHO region", "Exclusions by World Bank income group",
  "Missing data by covariate"))

ft11 <- flextable(combined) %>%
  set_header_labels(col1 = "", col2 = "n (% excluded)") %>%
  bold(i = header_rows, part = "body") %>%
  italic(i = header_rows, part = "body") %>%
  fontsize(size = 9, part = "all") %>%
  font(fontname = "Times New Roman", part = "all") %>%
  padding(padding.top = 1, padding.bottom = 1, part = "body") %>%
  border_remove() %>%
  hline_top(border = fp_border(color = "black", width = 1), part = "header") %>%
  hline_bottom(border = fp_border(color = "black", width = 1), part = "header") %>%
  hline_bottom(border = fp_border(color = "black", width = 1), part = "body") %>%
  width(j = 1, width = 3.2) %>%
  width(j = 2, width = 1.8) %>%
  align(j = 2, align = "center", part = "all") %>%
  add_footer_lines(paste0(
    "Of ", nrow(analysis), " samples, ", nrow(analysis) - nrow(a_pool),
    " lacked matched resistome data and ", sum(!a_pool$included),
    " had missing covariates, yielding ", sum(a_pool$included),
    " in the final analysis (", n_distinct(a_pool$country[a_pool$included]),
    " of ", n_distinct(analysis$country), " countries). ",
    "Data for exclusions by region and income are retained/total (% excluded). ",
    "Exclusion rates differed by WHO region (Fisher\u2019s exact p",
    ifelse(fisher_region$p.value < 0.001, "<0\u00B7001",
           sprintf("=%.3f", fisher_region$p.value)),
    ") and income group (p",
    ifelse(fisher_income$p.value < 0.001, "<0\u00B7001",
           sprintf("=%.3f", fisher_income$p.value)),
    "), driven primarily by C\u00f4te d\u2019Ivoire (n=18) and Taiwan (n=8), ",
    "which lacked World Bank data entirely. ",
    "For covariates, data are n (%) missing out of ", nrow(a_pool), " resistome-matched samples."
  )) %>%
  fontsize(size = 8, part = "footer") %>%
  font(fontname = "Times New Roman", part = "footer")

doc11 <- read_docx() %>%
  body_add_par("Supplementary Table 11. Sample derivation and completeness of covariate data.", style = "heading 2") %>%
  body_add_par("") %>%
  body_add_flextable(ft11)

print(doc11, target = file.path(tab_dir, "eTable11_complete_case.docx"))
cat("  Saved eTable11_complete_case.docx\n")

cat("\nAll .docx tables generated.\n")
