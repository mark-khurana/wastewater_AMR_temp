# Extract WorldClim bioclimatic variables and NASA POWER daily temperature metrics.
# Outputs: Datasets/climate_metrics.rds, Datasets/nasapower_daily_cache.rds

library(tidyverse)
library(lubridate)
library(terra)
library(geodata)
library(nasapower)
library(here)

samples <- readRDS(here("Datasets", "clean_samples.rds"))

locations <- samples %>%
  distinct(lat, lon, .keep_all = TRUE) %>%
  select(city, country, lat, lon)

wc_path <- here("Datasets", "worldclim_cache")

bio <- worldclim_global(var = "bio", res = 10, path = wc_path)

pts <- vect(locations, geom = c("lon", "lat"), crs = "EPSG:4326")
bio_vals <- extract(bio, pts) %>% as_tibble() %>% select(-ID)

bio_names <- c(
  "BIO1_annual_mean_temp", "BIO2_mean_diurnal_range", "BIO3_isothermality",
  "BIO4_temp_seasonality", "BIO5_max_temp_warmest_month",
  "BIO6_min_temp_coldest_month", "BIO7_temp_annual_range",
  "BIO8_mean_temp_wettest_quarter", "BIO9_mean_temp_driest_quarter",
  "BIO10_mean_temp_warmest_quarter", "BIO11_mean_temp_coldest_quarter",
  "BIO12_annual_precip", "BIO13_precip_wettest_month",
  "BIO14_precip_driest_month", "BIO15_precip_seasonality",
  "BIO16_precip_wettest_quarter", "BIO17_precip_driest_quarter",
  "BIO18_precip_warmest_quarter", "BIO19_precip_coldest_quarter"
)
names(bio_vals) <- bio_names

locations_bio <- bind_cols(locations, bio_vals)

loc_date_ranges <- samples %>%
  group_by(lat, lon) %>%
  summarise(
    earliest = min(collection_date, na.rm = TRUE) - days(365),
    latest = pmax(
      max(collection_date, na.rm = TRUE),
      as.Date(paste0(max(year, na.rm = TRUE), "-12-31"))
    ),
    .groups = "drop"
  ) %>%
  mutate(earliest = pmin(earliest, as.Date(paste0(year(earliest), "-01-01"))))

power_cache <- here("Datasets", "nasapower_daily_cache.rds")
if (file.exists(power_cache)) {
  power_data <- readRDS(power_cache)
} else {
  power_data <- list()
}

n_locs <- nrow(loc_date_ranges)
for (i in seq_len(n_locs)) {
  row <- loc_date_ranges[i, ]
  loc_key <- paste0(round(row$lat, 4), "_", round(row$lon, 4))
  if (loc_key %in% names(power_data)) next

  result <- tryCatch({
    years_needed <- seq(year(row$earliest), year(row$latest))
    yearly_dfs <- list()
    for (yr in years_needed) {
      start_d <- max(row$earliest, as.Date(paste0(yr, "-01-01")))
      end_d   <- min(row$latest, as.Date(paste0(yr, "-12-31")))
      d <- get_power(community = "ag", lonlat = c(row$lon, row$lat),
                     pars = c("T2M", "T2M_MAX", "T2M_MIN"),
                     dates = c(as.character(start_d), as.character(end_d)),
                     temporal_api = "daily")
      yearly_dfs[[as.character(yr)]] <- d %>%
        select(LON, LAT, YYYYMMDD, T2M, T2M_MAX, T2M_MIN) %>%
        rename(date = YYYYMMDD, t_mean = T2M, t_max = T2M_MAX, t_min = T2M_MIN)
    }
    bind_rows(yearly_dfs)
  }, error = function(e) { warning(e$message); NULL })

  if (!is.null(result) && nrow(result) > 0) {
    power_data[[loc_key]] <- result
  }

  if (i %% 25 == 0) { saveRDS(power_data, power_cache) }
  Sys.sleep(0.3)
}
saveRDS(power_data, power_cache)

daily_all <- bind_rows(power_data, .id = "loc_key") %>%
  mutate(
    lat = as.numeric(str_extract(loc_key, "^[^_]+")),
    lon = as.numeric(str_extract(loc_key, "[^_]+$")),
    year = year(date), month = month(date)
  ) %>%
  select(-loc_key, -LON, -LAT)

annual_metrics <- daily_all %>%
  group_by(lat, lon, year) %>%
  summarise(
    T_annual_mean = mean(t_mean, na.rm = TRUE),
    T_max = max(t_max, na.rm = TRUE),
    T_days30 = sum(t_max > 30, na.rm = TRUE),
    T_95pct = quantile(t_mean, 0.95, na.rm = TRUE),
    T_variability = sd(t_mean, na.rm = TRUE),
    .groups = "drop"
  )

amplitude <- daily_all %>%
  group_by(lat, lon, year, month) %>%
  summarise(monthly_mean = mean(t_mean, na.rm = TRUE), .groups = "drop") %>%
  group_by(lat, lon, year) %>%
  summarise(T_amplitude = max(monthly_mean) - min(monthly_mean), .groups = "drop")

annual_metrics <- left_join(annual_metrics, amplitude, by = c("lat", "lon", "year"))

daily_split <- daily_all %>%
  mutate(loc_key = paste0(round(lat, 4), "_", round(lon, 4))) %>%
  group_by(loc_key) %>% group_split() %>%
  set_names(map_chr(., ~ .x$loc_key[1]))

sample_dates <- samples %>%
  mutate(loc_key = paste0(round(lat, 4), "_", round(lon, 4))) %>%
  distinct(loc_key, collection_date, .keep_all = TRUE) %>%
  filter(!is.na(collection_date))

windowed <- tibble(loc_key = sample_dates$loc_key,
                   collection_date = sample_dates$collection_date,
                   T_sampling = NA_real_, T_30d = NA_real_,
                   T_90d = NA_real_, T_365d = NA_real_)

for (i in seq_len(nrow(windowed))) {
  key <- windowed$loc_key[i]
  sdate <- windowed$collection_date[i]
  if (!key %in% names(daily_split)) next
  city <- daily_split[[key]]
  if (is.null(city) || nrow(city) == 0) next

  idx <- which(city$date == sdate)
  if (length(idx) > 0) windowed$T_sampling[i] <- city$t_mean[idx[1]]

  w30  <- city$t_mean[city$date >= (sdate - 30) & city$date <= sdate]
  w90  <- city$t_mean[city$date >= (sdate - 90) & city$date <= sdate]
  w365 <- city$t_mean[city$date >= (sdate - 365) & city$date <= sdate]

  if (length(w30) > 0)  windowed$T_30d[i]  <- mean(w30, na.rm = TRUE)
  if (length(w90) > 0)  windowed$T_90d[i]  <- mean(w90, na.rm = TRUE)
  if (length(w365) > 0) windowed$T_365d[i] <- mean(w365, na.rm = TRUE)
}

final <- samples %>%
  mutate(loc_key = paste0(round(lat, 4), "_", round(lon, 4)))

final <- final %>%
  left_join(windowed, by = c("loc_key", "collection_date"))

final <- final %>%
  mutate(lat_r = round(lat, 4), lon_r = round(lon, 4)) %>%
  left_join(annual_metrics %>% mutate(lat_r = round(lat, 4), lon_r = round(lon, 4)) %>%
              select(-lat, -lon),
            by = c("lat_r", "lon_r", "year")) %>%
  select(-lat_r, -lon_r)

bio_join <- locations_bio %>%
  mutate(lat_r = round(lat, 4), lon_r = round(lon, 4)) %>%
  select(-city, -country, -lat, -lon) %>%
  distinct(lat_r, lon_r, .keep_all = TRUE)

final <- final %>%
  mutate(lat_r = round(lat, 4), lon_r = round(lon, 4)) %>%
  left_join(bio_join, by = c("lat_r", "lon_r")) %>%
  select(-lat_r, -lon_r, -loc_key)

stopifnot(n_distinct(final$genepid) == nrow(final))

saveRDS(final, here("Datasets", "climate_metrics.rds"))
