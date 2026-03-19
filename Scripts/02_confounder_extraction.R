# Download country-level socioeconomic confounders (WDI) and antimicrobial consumption data.
# Outputs: Datasets/country_confounders.rds

library(tidyverse)
library(WDI)
library(here)

samples <- readRDS(here("Datasets", "clean_samples.rds"))
our_countries <- unique(samples$country)
our_years <- range(samples$year)

indicators <- c(
  gdp_pcap_ppp     = "NY.GDP.PCAP.PP.CD",
  sanitation       = "SH.STA.BASS.ZS",
  pop_density      = "EN.POP.DNST",
  health_exp_gdp   = "SH.XPD.CHEX.GD.ZS",
  oop_health_exp   = "SH.XPD.OOPC.CH.ZS",
  tourism_arrivals = "ST.INT.ARVL",
  immunization_dpt = "SH.IMM.IDPT",
  pm25             = "EN.ATM.PM25.MC.M3",
  school_enroll    = "SE.SEC.ENRR"
)

cache_file <- here("Datasets", "wdi_cache_v3.rds")

if (file.exists(cache_file)) {
  wdi_list <- readRDS(cache_file)
} else {
  wdi_list <- list()
  for (nm in names(indicators)) {
    result <- tryCatch({
      WDI(country = "all",
          indicator = setNames(indicators[nm], nm),
          start = our_years[1], end = our_years[2],
          extra = TRUE)
    }, error = function(e) { NULL })

    if (!is.null(result)) {
      wdi_list[[nm]] <- result %>%
        filter(region != "Aggregates") %>%
        select(country, iso2c, iso3c, year, income, all_of(nm))
    }
    Sys.sleep(1)
  }
  saveRDS(wdi_list, cache_file)
}

wdi <- wdi_list[[1]] %>% select(country, iso2c, iso3c, year, income)
for (nm in names(wdi_list)) {
  wdi <- wdi %>%
    left_join(wdi_list[[nm]] %>% select(country, year, all_of(nm)),
              by = c("country", "year"))
}
wdi <- distinct(wdi)

human_amc_file <- here("Datasets", "human_amc_owid.csv")
if (file.exists(human_amc_file)) {
  human_amc_clean <- read_csv(human_amc_file, show_col_types = FALSE) %>%
    rename(country_amc = country)
} else {
  human_amc_clean <- NULL
}

animal_amc_file <- here("Datasets", "animal_amc_owid.csv")
if (file.exists(animal_amc_file)) {
  animal_amc_clean <- read_csv(animal_amc_file, show_col_types = FALSE) %>%
    rename(country_amc = country)
} else {
  animal_amc_clean <- NULL
}

name_map <- tribble(
  ~martiny_name,        ~wdi_name,
  "USA",                "United States",
  "UK",                 "United Kingdom",
  "Czech Republic",     "Czechia",
  "South Korea",        "Korea, Rep.",
  "Iran",               "Iran, Islamic Rep.",
  "Russia",             "Russian Federation",
  "Venezuela",          "Venezuela, RB",
  "Laos",               "Lao PDR",
  "Vietnam",            "Viet Nam",
  "Syria",              "Syrian Arab Republic",
  "Egypt",              "Egypt, Arab Rep.",
  "Congo",              "Congo, Rep.",
  "DR Congo",           "Congo, Dem. Rep.",
  "Ivory Coast",        "Cote d'Ivoire",
  "C\u221a\u00a5te d'Ivoire", "Cote d'Ivoire",
  "Gambia",             "Gambia, The",
  "Brunei",             "Brunei Darussalam",
  "Slovakia",           "Slovak Republic",
  "Turkey",             "Turkiye",
  "Kyrgyzstan",         "Kyrgyz Republic",
  "Yemen",              "Yemen, Rep.",
  "Palestine",          "West Bank and Gaza",
  "Hong Kong",          "Hong Kong SAR, China",
  "United States of America", "United States",
  "Bosnia and Herz.",   "Bosnia and Herzegovina",
  "Saint Lucia",        "St. Lucia",
  "Dem. Rep. Congo",    "Congo, Dem. Rep."
)

samples_mapped <- samples %>%
  left_join(name_map, by = c("country" = "martiny_name")) %>%
  mutate(wdi_country = coalesce(wdi_name, country)) %>%
  select(-wdi_name)

confounders <- samples_mapped %>%
  select(genepid, country, wdi_country, year) %>%
  left_join(wdi, by = c("wdi_country" = "country", "year" = "year")) %>%
  select(genepid, country, wdi_country, year, iso3c, income,
         gdp_pcap_ppp, sanitation, pop_density, health_exp_gdp,
         oop_health_exp, tourism_arrivals, immunization_dpt,
         pm25, school_enroll)

if (!is.null(human_amc_clean)) {
  hamc_lookup <- human_amc_clean %>% select(iso3c, year, human_amc_ddd)
  confounders <- confounders %>%
    left_join(hamc_lookup, by = c("iso3c", "year"))

  country_mean_hamc <- human_amc_clean %>%
    group_by(iso3c) %>%
    summarise(human_amc_ddd_mean = mean(human_amc_ddd, na.rm = TRUE), .groups = "drop")
  confounders <- confounders %>%
    left_join(country_mean_hamc, by = "iso3c") %>%
    mutate(human_amc_ddd = coalesce(human_amc_ddd, human_amc_ddd_mean)) %>%
    select(-human_amc_ddd_mean)
} else {
  confounders$human_amc_ddd <- NA_real_
}

if (!is.null(animal_amc_clean)) {
  aamc_lookup <- animal_amc_clean %>%
    group_by(iso3c) %>%
    summarise(animal_amc_mgkg = mean(animal_amc_mgkg, na.rm = TRUE), .groups = "drop")
  confounders <- confounders %>%
    left_join(aamc_lookup, by = "iso3c")
} else {
  confounders$animal_amc_mgkg <- NA_real_
}

confounders_final <- confounders %>%
  select(genepid, income,
         gdp_pcap_ppp, sanitation, pop_density, health_exp_gdp,
         oop_health_exp, tourism_arrivals, immunization_dpt,
         pm25, school_enroll,
         human_amc_ddd, animal_amc_mgkg)

saveRDS(confounders_final, here("Datasets", "country_confounders.rds"))
