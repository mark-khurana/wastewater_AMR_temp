# Temperature and antimicrobial resistance in global wastewater

Code for: *Association between temperature and antimicrobial resistance in global wastewater: a metagenomic analysis of 104 countries*

## Data

The analysis uses publicly available data. Place raw files in `Datasets/` before running.

- **Sewage metagenomes** (cross-sectional): Martiny et al. 2025 ([doi:10.5281/zenodo.14652832](https://doi.org/10.5281/zenodo.14652832)). Download into `Datasets/Martiny_et_al/`.
- **Sewage metagenomes** (longitudinal, 5 European cities): Becsei et al. 2026 ([doi:10.1093/database/baaf089](https://doi.org/10.1093/database/baaf089)). Download into `Datasets/Becsei_et_al/`.
- **Country-level covariates**: World Bank World Development Indicators ([databank.worldbank.org](https://datatopics.worldbank.org/world-development-indicators/)). Downloaded automatically by `02_confounder_extraction.R` via the `WDI` package.
- **Temperature**: NASA POWER ([power.larc.nasa.gov](https://power.larc.nasa.gov)). Downloaded automatically by `01_climate_extraction.R` and `01b_year_matched_temperature.R` via the `nasapower` package.
- **Animal antimicrobial consumption**: Mulchandani et al. 2023 ([doi:10.1371/journal.pgph.0001305](https://doi.org/10.1371/journal.pgph.0001305)).
- **Human antimicrobial consumption**: WHO GLASS ([who.int/initiatives/glass](https://www.who.int/initiatives/glass)).

## Requirements

R (>= 4.3). Key packages: `tidyverse`, `vegan`, `compositions`, `nasapower`, `WDI`, `geodata`, `terra`, `patchwork`, `showtext`, `flextable`, `officer`, `here`.

Install all at once:

```r
install.packages(c("tidyverse", "vegan", "compositions", "nasapower", "WDI",
                   "geodata", "terra", "patchwork", "showtext", "flextable",
                   "officer", "here", "readxl", "scales", "broom", "mgcv"))
```

## Scripts

Run in order. Scripts that call external APIs (01, 01b, 02) cache results locally so re-runs are fast.

| Script | Description |
|--------|-------------|
| `00_data_clean.R` | Load and clean Martiny et al. sample metadata and ARG counts |
| `01_climate_extraction.R` | Extract WorldClim bioclimatic variables and NASA POWER temperatures |
| `01b_year_matched_temperature.R` | Year-matched 30-day mean temperatures from NASA POWER |
| `02_confounder_extraction.R` | Download country-level socioeconomic covariates from World Bank |
| `03_merge_analysis_ready.R` | Merge all data, apply CLR transformation, produce analysis-ready dataset |
| `04_layer1_naive_models.R` | PERMANOVA: temperature + WHO region + year (naive models) |
| `05_layer2_covariates.R` | PERMANOVA with full covariate adjustment and variance partitioning |
| `06_layer3_species_mediation.R` | Bacteriome mediation analysis and genus-level temperature correlations |
| `06b_species_pathogens.R` | Species-level correlations for WHO priority pathogens |
| `07_figures.R` | Figures 1, 2, and supplementary PCoA ordination (Figure S2) |
| `08_timeseries_5cities.R` | Within-city temporal replication using 5-city longitudinal data |
| `09_timeseries_figures.R` | Figure 3 and supplementary temporal figures (Figures S3, S4) |
| `sensitivity_humidity.R` | Humidity sensitivity analysis (Table S5) |
| `sensitivity_temp_windows.R` | Temperature averaging window sensitivity (Table S6) |
| `eFigure_pcoa_sensitivity.R` | PCoA scree plot and axis sensitivity (Figure S1) |

## Output

- `Figures/` — PDF and PNG figures
- `Results/` — CSV and RDS intermediate results
