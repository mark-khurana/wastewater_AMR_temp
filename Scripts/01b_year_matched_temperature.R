# Extract year-matched annual mean temperature from NASA POWER and add WHO regions.
# Outputs: updated Datasets/analysis_ready.rds

library(tidyverse)
library(httr)
library(jsonlite)
library(here)

analysis <- readRDS(here("Datasets", "analysis_ready.rds")) %>%
  select(-any_of(c("T2M_annual", "WHO_region", "lat_r", "lon_r")))

loc_year <- analysis %>%
  distinct(lat, lon, year) %>%
  mutate(loc_id = paste(lat, lon, year, sep = "_"))

cache_path <- here("Datasets", "nasa_power_annual_temp.rds")
if (file.exists(cache_path)) {
  temp_data <- readRDS(cache_path)
  existing <- temp_data %>% mutate(loc_id = paste(lat, lon, year, sep = "_"))
  needed <- loc_year %>% filter(!loc_id %in% existing$loc_id)
} else {
  temp_data <- tibble()
  needed <- loc_year
}

get_nasa_annual_temp <- function(lat, lon, year) {
  url <- sprintf(
    "https://power.larc.nasa.gov/api/temporal/monthly/point?parameters=T2M&community=AG&longitude=%f&latitude=%f&start=%d&end=%d&format=JSON&header=false",
    lon, lat, year, year
  )

  tryCatch({
    resp <- GET(url, timeout(30))
    if (status_code(resp) == 200) {
      data <- fromJSON(content(resp, "text", encoding = "UTF-8"))
      annual_key <- paste0(year, "13")
      t2m <- data$properties$parameter$T2M[[annual_key]]
      if (is.null(t2m) || t2m == -999) return(NA_real_)
      return(as.numeric(t2m))
    } else {
      return(NA_real_)
    }
  }, error = function(e) {
    return(NA_real_)
  })
}

if (nrow(needed) > 0) {
  results <- tibble(
    lat = needed$lat,
    lon = needed$lon,
    year = needed$year,
    T2M_annual = NA_real_
  )

  for (i in seq_len(nrow(needed))) {
    results$T2M_annual[i] <- get_nasa_annual_temp(
      needed$lat[i], needed$lon[i], needed$year[i]
    )
    Sys.sleep(0.3)
  }

  temp_data <- bind_rows(temp_data, results)
  saveRDS(temp_data, cache_path)
}

temp_data_round <- temp_data %>%
  mutate(lat_r = round(lat, 4), lon_r = round(lon, 4)) %>%
  select(lat_r, lon_r, year, T2M_annual)

temp_data_round <- temp_data_round %>%
  distinct(lat_r, lon_r, year, .keep_all = TRUE)

analysis <- analysis %>%
  mutate(lat_r = round(lat, 4), lon_r = round(lon, 4)) %>%
  left_join(temp_data_round, by = c("lat_r", "lon_r", "year")) %>%
  select(-lat_r, -lon_r)

who_region_map <- tribble(
  ~country, ~WHO_region,
  "Algeria", "AFRO", "Angola", "AFRO", "Benin", "AFRO", "Botswana", "AFRO",
  "Burkina Faso", "AFRO", "Burundi", "AFRO", "Cabo Verde", "AFRO",
  "Cameroon", "AFRO", "Central African Republic", "AFRO", "Chad", "AFRO",
  "Comoros", "AFRO", "Congo", "AFRO",
  "Democratic Republic of the Congo", "AFRO", "Cote d'Ivoire", "AFRO",
  "Equatorial Guinea", "AFRO", "Eritrea", "AFRO", "Eswatini", "AFRO",
  "Ethiopia", "AFRO", "Gabon", "AFRO", "Gambia", "AFRO", "Ghana", "AFRO",
  "Guinea", "AFRO", "Guinea-Bissau", "AFRO", "Kenya", "AFRO",
  "Lesotho", "AFRO", "Liberia", "AFRO", "Madagascar", "AFRO",
  "Malawi", "AFRO", "Mali", "AFRO", "Mauritania", "AFRO",
  "Mauritius", "AFRO", "Mozambique", "AFRO", "Namibia", "AFRO",
  "Niger", "AFRO", "Nigeria", "AFRO", "Rwanda", "AFRO",
  "Sao Tome and Principe", "AFRO", "Senegal", "AFRO", "Seychelles", "AFRO",
  "Sierra Leone", "AFRO", "South Africa", "AFRO", "South Sudan", "AFRO",
  "Tanzania", "AFRO", "Togo", "AFRO", "Uganda", "AFRO",
  "Zambia", "AFRO", "Zimbabwe", "AFRO",
  "United Republic of Tanzania", "AFRO",
  "Antigua and Barbuda", "AMRO", "Argentina", "AMRO", "Bahamas", "AMRO",
  "Barbados", "AMRO", "Belize", "AMRO", "Bolivia", "AMRO", "Brazil", "AMRO",
  "Canada", "AMRO", "Chile", "AMRO", "Colombia", "AMRO", "Costa Rica", "AMRO",
  "Cuba", "AMRO", "Dominica", "AMRO", "Dominican Republic", "AMRO",
  "Ecuador", "AMRO", "El Salvador", "AMRO", "Grenada", "AMRO",
  "Guatemala", "AMRO", "Guyana", "AMRO", "Haiti", "AMRO",
  "Honduras", "AMRO", "Jamaica", "AMRO", "Mexico", "AMRO",
  "Nicaragua", "AMRO", "Panama", "AMRO", "Paraguay", "AMRO",
  "Peru", "AMRO", "Saint Kitts and Nevis", "AMRO",
  "Saint Lucia", "AMRO", "Saint Vincent and the Grenadines", "AMRO",
  "Suriname", "AMRO", "Trinidad and Tobago", "AMRO",
  "United States of America", "AMRO", "Uruguay", "AMRO", "Venezuela", "AMRO",
  "Bolivia (Plurinational State of)", "AMRO",
  "Venezuela (Bolivarian Republic of)", "AMRO",
  "Bangladesh", "SEARO", "Bhutan", "SEARO", "India", "SEARO",
  "Indonesia", "SEARO", "Maldives", "SEARO", "Myanmar", "SEARO",
  "Nepal", "SEARO", "Sri Lanka", "SEARO", "Thailand", "SEARO",
  "Timor-Leste", "SEARO", "North Korea", "SEARO",
  "Democratic People's Republic of Korea", "SEARO",
  "Albania", "EURO", "Andorra", "EURO", "Armenia", "EURO",
  "Austria", "EURO", "Azerbaijan", "EURO", "Belarus", "EURO",
  "Belgium", "EURO", "Bosnia and Herzegovina", "EURO", "Bulgaria", "EURO",
  "Croatia", "EURO", "Cyprus", "EURO", "Czechia", "EURO",
  "Denmark", "EURO", "Estonia", "EURO", "Finland", "EURO",
  "France", "EURO", "Georgia", "EURO", "Germany", "EURO",
  "Greece", "EURO", "Hungary", "EURO", "Iceland", "EURO",
  "Ireland", "EURO", "Israel", "EURO", "Italy", "EURO",
  "Kazakhstan", "EURO", "Kyrgyzstan", "EURO", "Latvia", "EURO",
  "Lithuania", "EURO", "Luxembourg", "EURO", "Malta", "EURO",
  "Monaco", "EURO", "Montenegro", "EURO", "Netherlands", "EURO",
  "North Macedonia", "EURO", "Norway", "EURO", "Poland", "EURO",
  "Portugal", "EURO", "Moldova", "EURO", "Romania", "EURO",
  "Russia", "EURO", "Russian Federation", "EURO", "San Marino", "EURO",
  "Serbia", "EURO", "Slovakia", "EURO", "Slovenia", "EURO",
  "Spain", "EURO", "Sweden", "EURO", "Switzerland", "EURO",
  "Tajikistan", "EURO", "Turkey", "EURO", "Turkmenistan", "EURO",
  "Ukraine", "EURO", "United Kingdom", "EURO", "Uzbekistan", "EURO",
  "Republic of Moldova", "EURO",
  "United Kingdom of Great Britain and Northern Ireland", "EURO",
  "T\u00fcrkiye", "EURO",
  "Afghanistan", "EMRO", "Bahrain", "EMRO", "Djibouti", "EMRO",
  "Egypt", "EMRO", "Iran", "EMRO", "Iran (Islamic Republic of)", "EMRO",
  "Iraq", "EMRO", "Jordan", "EMRO", "Kuwait", "EMRO",
  "Lebanon", "EMRO", "Libya", "EMRO", "Morocco", "EMRO",
  "Oman", "EMRO", "Pakistan", "EMRO", "Palestine", "EMRO",
  "Qatar", "EMRO", "Saudi Arabia", "EMRO", "Somalia", "EMRO",
  "Sudan", "EMRO", "Syria", "EMRO", "Syrian Arab Republic", "EMRO",
  "Tunisia", "EMRO", "United Arab Emirates", "EMRO", "Yemen", "EMRO",
  "Australia", "WPRO", "Brunei", "WPRO", "Brunei Darussalam", "WPRO",
  "Cambodia", "WPRO", "China", "WPRO", "Cook Islands", "WPRO",
  "Fiji", "WPRO", "Japan", "WPRO", "Kiribati", "WPRO",
  "Lao People's Democratic Republic", "WPRO", "Laos", "WPRO",
  "Malaysia", "WPRO", "Marshall Islands", "WPRO",
  "Micronesia", "WPRO", "Mongolia", "WPRO", "Nauru", "WPRO",
  "New Zealand", "WPRO", "Niue", "WPRO", "Palau", "WPRO",
  "Papua New Guinea", "WPRO", "Philippines", "WPRO",
  "South Korea", "WPRO", "Republic of Korea", "WPRO",
  "Samoa", "WPRO", "Singapore", "WPRO", "Solomon Islands", "WPRO",
  "Tonga", "WPRO", "Tuvalu", "WPRO", "Vanuatu", "WPRO",
  "Viet Nam", "WPRO", "Vietnam", "WPRO",
  "Taiwan", "WPRO", "Hong Kong", "WPRO",
  "Greenland", "EURO", "Bosnia and Herz.", "EURO",
  "C\u00f4te d'Ivoire", "AFRO", "Cote d'Ivoire", "AFRO",
  "Dem. Rep. Congo", "AFRO", "Kosovo", "EURO",
  "C\u221a\u00a5te d'Ivoire", "AFRO"
)

who_region_map <- who_region_map %>% distinct(country, .keep_all = TRUE)

analysis <- analysis %>%
  left_join(who_region_map, by = "country")

analysis$country[analysis$country == "C\u221a\u00a5te d'Ivoire"] <- "Cote d'Ivoire"

analysis$Region <- analysis$WHO_region

analysis <- analysis %>% distinct(genepid, .keep_all = TRUE)

saveRDS(analysis, here("Datasets", "analysis_ready.rds"))
