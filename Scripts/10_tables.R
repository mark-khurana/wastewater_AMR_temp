# =============================================================================
# 10_tables.R — Generate all LaTeX tables for Lancet Planetary Health submission
# =============================================================================

library(tidyverse)

base_dir <- "/Users/bxc331/Desktop/AMR_open_project"
out_dir  <- file.path(base_dir, "Datasets")
res_dir  <- file.path(base_dir, "Results")
tab_dir  <- file.path(base_dir, "Tables")
if (!dir.exists(tab_dir)) dir.create(tab_dir)

cat("=== 10_tables.R ===\n\n")

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

cat("Analysis subset: n =", nrow(cc), "\n\n")

# Helper
fmt_iqr <- function(med, lo, hi, dig = 1) {
  sprintf("%.*f (%.*f--%.*f)", dig, med, dig, lo, dig, hi)
}

fmt_range <- function(lo, hi, dig = 1) {
  sprintf("%.*f to %.*f", dig, lo, dig, hi)
}

# =============================================================================
# TABLE 1: Sampling design by World Bank Region (MAIN)
# =============================================================================
cat("--- Table 1: Sampling design ---\n")

regions <- sort(unique(cc$Region))

tab1 <- bind_rows(
  map_dfr(regions, function(reg) {
    d <- cc %>% filter(Region == reg)
    tibble(
      Region = reg,
      n = nrow(d),
      n_cities = n_distinct(d$city),
      n_countries = n_distinct(d$country),
      years = paste(sort(unique(d$year)), collapse = ", "),
      temp_range = fmt_range(min(d$T_30d), max(d$T_30d)),
      temp_med = round(median(d$T_30d), 1)
    )
  }),
  {
    tibble(
      Region = "Total",
      n = nrow(cc),
      n_cities = n_distinct(cc$city),
      n_countries = n_distinct(cc$country),
      years = paste(sort(unique(cc$year)), collapse = ", "),
      temp_range = fmt_range(min(cc$T_30d), max(cc$T_30d)),
      temp_med = round(median(cc$T_30d), 1)
    )
  }
)

# Year distribution as a separate summary
year_tab <- cc %>% count(year) %>% mutate(pct = round(100 * n / sum(n), 1))

tex1 <- c(
  "\\begin{table}[ht]",
  "\\centering",
  "\\caption{Sampling characteristics of 1,077 sewage metagenomes from 104 countries, by World Bank region.}",
  "\\label{tab:descriptive}",
  "\\small",
  "\\begin{tabular}{l rrrr rr}",
  "\\hline",
  paste0("\\textbf{Region} & \\textbf{Samples} & \\textbf{Cities} & ",
         "\\textbf{Countries} & \\textbf{Collection years} & ",
         "\\textbf{Temperature range} & \\textbf{Median} \\\\"),
  paste0(" & & & & & ($^\\circ$C) & ($^\\circ$C) \\\\"),
  "\\hline"
)

for (i in seq_len(nrow(tab1))) {
  r <- tab1[i, ]
  reg_escaped <- gsub("&", "\\\\&", r$Region)
  lab <- if (r$Region == "Total") paste0("\\hline\n\\textbf{Total}") else reg_escaped
  # Compact year range
  yr_compact <- paste0(min(cc$year[cc$Region == r$Region | r$Region == "Total"]),
                       "--", max(cc$year[cc$Region == r$Region | r$Region == "Total"]))
  if (r$Region == "Total") yr_compact <- paste0(min(cc$year), "--", max(cc$year))
  tex1 <- c(tex1, sprintf("%s & %d & %d & %d & %s & %s & %.1f \\\\",
    lab, r$n, r$n_cities, r$n_countries, yr_compact, r$temp_range, r$temp_med))
}

# Year distribution row
year_str <- paste(sprintf("%d (n=%d)", year_tab$year, year_tab$n), collapse = "; ")

tex1 <- c(tex1,
  "\\hline",
  "\\end{tabular}",
  "\\vspace{3mm}",
  "\\raggedright",
  "\\scriptsize",
  paste0("Temperature: 30-day mean air temperature preceding sample collection (NASA POWER T2M). ",
         "Year distribution: ", year_str, ". ",
         "Country-level socioeconomic covariates (GDP, sanitation, health expenditure, immunisation, ",
         "population density, animal antimicrobial consumption) are listed in Supplementary Table 5. ",
         "All covariates are matched to sample country of origin and do not vary between samples within a country."),
  "\\end{table}"
)

writeLines(tex1, file.path(tab_dir, "Table1_descriptive.tex"))
cat("  Saved Table1_descriptive.tex\n")

# =============================================================================
# eTABLE 1: Full PERMANOVA results (SUPPLEMENTARY)
# =============================================================================
cat("--- eTable 1: PERMANOVA results ---\n")

perm <- readRDS(file.path(res_dir, "layer2_partial_permanova.rds"))

format_perm <- function(perm_obj, label) {
  df <- as.data.frame(perm_obj)
  df$term <- rownames(df)
  df <- df %>%
    filter(term != "Total") %>%
    mutate(
      R2_pct = round(.data$R2 * 100, 3),
      F_val = round(.data$`F`, 2),
      p = ifelse(is.na(.data$`Pr(>F)`), NA_character_,
            ifelse(.data$`Pr(>F)` < 0.001, "$<$0.001",
              ifelse(.data$`Pr(>F)` < 0.01, sprintf("%.3f", .data$`Pr(>F)`),
                     sprintf("%.2f", .data$`Pr(>F)`)))),
      resistome = label
    ) %>%
    select(resistome, term, Df, R2_pct, F_val, p)
  df
}

# Extract the full model results
perm_fg_full  <- format_perm(perm$perm_fg$full, "Functional")
perm_acq_full <- format_perm(perm$perm_acq$full, "Acquired")

perm_all <- bind_rows(perm_fg_full, perm_acq_full)

# Clean term names
perm_all <- perm_all %>%
  mutate(term = case_when(
    term == "T_30d" ~ "30-day mean temperature",
    term == "log_gdp" ~ "GDP per capita (log)",
    term == "sanitation" ~ "Basic sanitation (\\%)",
    term == "health_exp_gdp" ~ "Health expenditure (\\% GDP)",
    term == "oop_health_exp" ~ "Out-of-pocket health exp. (\\%)",
    term == "immunization_dpt" ~ "DPT immunisation (\\%)",
    term == "log_pop_density" ~ "Population density (log)",
    term == "log_animal_amc" ~ "Animal AMC (log mg/kg)",
    term == "Region" ~ "WHO Region",
    term == "year_f" ~ "Collection year",
    term == "Residual" ~ "Residual",
    TRUE ~ term
  ))

tex2 <- c(
  "\\begin{table}[ht]",
  "\\centering",
  "\\caption{PERMANOVA results: variance in resistome composition explained by temperature and covariates (marginal, Type III). n=1,077.}",
  "\\label{tab:permanova}",
  "\\footnotesize",
  "\\begin{tabular}{ll rrrr}",
  "\\hline",
  "\\textbf{Resistome} & \\textbf{Term} & \\textbf{df} & \\textbf{R$^2$ (\\%)} & \\textbf{F} & \\textbf{p} \\\\",
  "\\hline"
)

current_res <- ""
for (i in seq_len(nrow(perm_all))) {
  r <- perm_all[i, ]
  res_lab <- if (r$resistome != current_res) {
    current_res <<- r$resistome
    if (i > 1) tex2 <- c(tex2, "\\hline")
    paste0("\\textbf{", r$resistome, "}")
  } else ""

  p_val <- if (is.na(r$p)) "---" else r$p
  f_val <- if (is.na(r$F_val)) "---" else sprintf("%.2f", r$F_val)
  tex2 <- c(tex2, sprintf("%s & %s & %d & %.3f & %s & %s \\\\",
    res_lab, r$term, r$Df, r$R2_pct, f_val, p_val))
}

tex2 <- c(tex2,
  "\\hline",
  "\\end{tabular}",
  "\\vspace{2mm}",
  "\\raggedright",
  "\\scriptsize",
  "PERMANOVA on Aitchison distance (Euclidean on CLR-transformed gene cluster abundances). 999 permutations, marginal (Type III) terms. Both resistome types use the same samples and covariates.",
  "\\end{table}"
)

writeLines(tex2, file.path(tab_dir, "eTable1_permanova.tex"))
cat("  Saved eTable1_permanova.tex\n")

# =============================================================================
# eTABLE 2: Full genus-level correlations (SUPPLEMENTARY)
# =============================================================================
cat("--- eTable 2: Full genus correlations ---\n")

genus_corr <- read.csv(file.path(res_dir, "layer3_genus_temp_correlations.csv"))
cat("  Total genera:", nrow(genus_corr), "\n")

# This is too large for a LaTeX table — save as CSV for supplementary data file
# and create a summary LaTeX table of top/bottom genera
top_warm <- genus_corr %>% arrange(desc(rho)) %>% head(15)
top_cool <- genus_corr %>% arrange(rho) %>% head(15)
top_genera <- bind_rows(
  top_warm %>% mutate(direction = "Warm-associated"),
  top_cool %>% mutate(direction = "Cool-associated")
)

# Clean genus names
top_genera <- top_genera %>%
  mutate(
    genus_clean = str_replace(genus, " ext_mOTU.*", " (unclassified)"),
    genus_clean = str_replace(genus_clean, "aceae gen\\. .*", "aceae (unclassified)"),
    genus_clean = str_replace(genus_clean, "ales gen\\. .*", "ales (unclassified)")
  )

tex3 <- c(
  "\\begin{table}[ht]",
  "\\centering",
  "\\caption{Top 15 warm-associated and 15 cool-associated genera by Spearman correlation with 30-day mean temperature (unadjusted). Full results for all genera available as Supplementary Data.}",
  "\\label{tab:genus_corr}",
  "\\footnotesize",
  "\\begin{tabular}{llrrr}",
  "\\hline",
  "\\textbf{Direction} & \\textbf{Genus} & \\textbf{$\\rho$} & \\textbf{p (adj.)} & \\textbf{n} \\\\",
  "\\hline"
)

current_dir <- ""
for (i in seq_len(nrow(top_genera))) {
  r <- top_genera[i, ]
  dir_lab <- if (r$direction != current_dir) {
    current_dir <<- r$direction
    if (i > 1) tex3 <- c(tex3, "\\hline")
    paste0("\\textit{", r$direction, "}")
  } else ""

  p_str <- if (r$p_adj < 0.001) "$<$0.001" else sprintf("%.3f", r$p_adj)
  tex3 <- c(tex3, sprintf("%s & %s & %.3f & %s & %d \\\\",
    dir_lab, r$genus_clean, r$rho, p_str, r$n))
}

tex3 <- c(tex3,
  "\\hline",
  "\\end{tabular}",
  "\\vspace{2mm}",
  "\\raggedright",
  "\\scriptsize",
  paste0("Unadjusted Spearman correlations between genus-level relative abundance and 30-day mean temperature. ",
         "p-values Bonferroni-adjusted for ", nrow(genus_corr), " tests. ",
         "Full table of all genera provided as Supplementary Data file."),
  "\\end{table}"
)

writeLines(tex3, file.path(tab_dir, "eTable2_genus_correlations.tex"))
# Also save full CSV
write_csv(genus_corr, file.path(tab_dir, "eTable2_full_genus_correlations.csv"))
cat("  Saved eTable2_genus_correlations.tex + CSV\n")

# =============================================================================
# eTABLE 3: Species-level raw + adjusted correlations (SUPPLEMENTARY)
# =============================================================================
cat("--- eTable 3: Species correlations ---\n")

sp <- read_csv(file.path(res_dir, "adjusted_species_temp.csv"), show_col_types = FALSE)

tex4 <- c(
  "\\begin{table}[ht]",
  "\\centering",
  "\\caption{Species-level temperature associations: unadjusted and covariate-adjusted Spearman correlations for WHO priority pathogen species (n=1,077).}",
  "\\label{tab:species_corr}",
  "\\footnotesize",
  "\\begin{tabular}{ll rr rr}",
  "\\hline",
  " & & \\multicolumn{2}{c}{\\textbf{Unadjusted}} & \\multicolumn{2}{c}{\\textbf{Adjusted}} \\\\",
  "\\textbf{Genus} & \\textbf{Taxon} & \\textbf{$\\rho$} & \\textbf{p} & \\textbf{$\\rho$} & \\textbf{p} \\\\",
  "\\hline"
)

fmt_p <- function(p) {
  if (is.na(p)) return("---")
  if (p < 0.001) return("$<$0.001")
  if (p < 0.01) return(sprintf("%.3f", p))
  sprintf("%.2f", p)
}

current_genus <- ""
for (i in seq_len(nrow(sp))) {
  r <- sp[i, ]
  g_lab <- if (r$genus != current_genus) {
    current_genus <<- r$genus
    paste0("\\textit{", r$genus, "}")
  } else ""

  rho_adj <- if (is.na(r$rho_adj)) "---" else sprintf("%.3f", r$rho_adj)
  tex4 <- c(tex4, sprintf("%s & %s & %.3f & %s & %s & %s \\\\",
    g_lab, r$display, r$rho_raw, fmt_p(r$p_raw), rho_adj, fmt_p(r$p_adj)))
}

tex4 <- c(tex4,
  "\\hline",
  "\\end{tabular}",
  "\\vspace{2mm}",
  "\\raggedright",
  "\\scriptsize",
  paste0("Adjustment: GDP (log), basic sanitation, health expenditure (\\% GDP), out-of-pocket health expenditure, ",
         "DPT immunisation, population density (log), animal antimicrobial consumption (log mg/kg), WHO Region, ",
         "and collection year. Partial correlations computed by residualising both abundance and temperature on all covariates."),
  "\\end{table}"
)

writeLines(tex4, file.path(tab_dir, "eTable3_species_correlations.tex"))
cat("  Saved eTable3_species_correlations.tex\n")

# =============================================================================
# eTABLE 4: Sensitivity analysis summary (SUPPLEMENTARY)
# =============================================================================
cat("--- eTable 4: Sensitivity analyses ---\n")

tex5 <- c(
  "\\begin{table}[ht]",
  "\\centering",
  "\\caption{Summary of sensitivity analyses for the temperature--resistome association.}",
  "\\label{tab:sensitivity}",
  "\\footnotesize",
  "\\begin{tabular}{p{4.5cm} p{2cm} p{3cm} p{4cm}}",
  "\\hline",
  "\\textbf{Analysis} & \\textbf{n} & \\textbf{Temperature R$^2$} & \\textbf{Key finding} \\\\",
  "\\hline",
  "Main model (30-day mean) & 1,066 & FG: 0.51\\%, Acq: 0.56\\% & Both p=0.001; $\\sim$68\\% mediated via bacteriome \\\\",
  "Without collection year & 1,066 & FG: 0.72\\%, Acq: 0.88\\% & Year explains $\\sim$2.7\\% independently \\\\",
  "AMC subset only & 522 & FG: 0.53\\%, Acq: 0.61\\% & Survives with similar magnitude \\\\",
  "Collection-day temperature & 1,066 & --- & Stronger for \\textit{A. baumannii} ($\\rho$=+0.38) \\\\",
  "Within-city temporal (30-day) & 162 & 6.3\\% (total) & 73\\% mediated via bacteriome \\\\",
  "Within-city (same-day temp) & 162 & --- & Consistent direction and magnitude \\\\",
  "Month-adjusted temporal & 162 & --- & Only Enterococcus retains significance \\\\",
  "\\hline",
  "\\end{tabular}",
  "\\vspace{2mm}",
  "\\raggedright",
  "\\scriptsize",
  paste0("FG: functional genomic resistome. Acq: acquired resistome. R$^2$: marginal PERMANOVA R$^2$ for temperature, ",
         "adjusting for all covariates. AMC subset: restricted to samples with animal antimicrobial consumption data. ",
         "Collection-day temperature: NASA POWER air temperature on sampling date. ",
         "Within-city temporal: biweekly sampling in 5 European cities (2019--2021), adjusted for city and calendar month."),
  "\\end{table}"
)

writeLines(tex5, file.path(tab_dir, "eTable4_sensitivity.tex"))
cat("  Saved eTable4_sensitivity.tex\n")

# =============================================================================
# eTABLE 5: Data sources (SUPPLEMENTARY)
# =============================================================================
cat("--- eTable 5: Data sources ---\n")

tex6 <- c(
  "\\begin{table}[ht]",
  "\\centering",
  "\\caption{Data sources for covariates and exposures.}",
  "\\label{tab:datasources}",
  "\\footnotesize",
  "\\begin{tabular}{p{3.5cm} p{3.5cm} p{2cm} p{2cm} p{2cm}}",
  "\\hline",
  "\\textbf{Variable} & \\textbf{Source} & \\textbf{Year(s)} & \\textbf{Coverage} & \\textbf{Resolution} \\\\",
  "\\hline",
  "Sewage metagenomes & Martiny et al. 2025 (Nat Commun) & 2016--2021 & 112 countries & Sample \\\\",
  "30-day mean temperature & NASA POWER v2.0 (T2M) & 2016--2021 & Global & 0.5$^\\circ$ \\\\",
  "Collection-day temperature & NASA POWER v2.0 & Daily & Global & 0.5$^\\circ$ \\\\",
  "GDP per capita (PPP) & World Bank WDI & 2018 & 98\\% & Country \\\\",
  "Basic sanitation (\\%) & World Bank WDI & 2017--2020 & 97\\% & Country \\\\",
  "Health expenditure (\\% GDP) & World Bank WDI & 2018--2019 & 97\\% & Country \\\\",
  "Out-of-pocket health exp. & World Bank WDI & 2018--2019 & 97\\% & Country \\\\",
  "DPT immunisation (\\%) & World Bank WDI & 2018--2020 & 97\\% & Country \\\\",
  "Population density & World Bank WDI & 2018--2020 & 98\\% & Country \\\\",
  "Animal AMC (mg/kg) & Mulchandani et al. 2023 (OWID) & 2017--2020 & 42 countries & Country \\\\",
  "5-city temporal data & Becsei et al. 2026 & 2019--2021 & 5 cities & Biweekly \\\\",
  "WHO Priority Pathogens & WHO BPPL 2024 & --- & --- & --- \\\\",
  "\\hline",
  "\\end{tabular}",
  "\\vspace{2mm}",
  "\\raggedright",
  "\\scriptsize",
  paste0("Coverage: percentage of the 1,077-sample analysis subset with non-missing data for that variable. ",
         "WDI: World Development Indicators. OWID: Our World in Data. ",
         "AMC: antimicrobial consumption. PPP: purchasing power parity. ",
         "Country-level covariates matched to sample country of origin."),
  "\\end{table}"
)

writeLines(tex6, file.path(tab_dir, "eTable5_datasources.tex"))
cat("  Saved eTable5_datasources.tex\n")

cat("\nAll tables generated.\n")
