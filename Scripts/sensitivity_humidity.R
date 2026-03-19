# Sensitivity analysis: temperature x relative humidity interaction on resistome composition

library(tidyverse)
library(vegan)
library(nasapower)
library(showtext)
library(patchwork)
library(here)

showtext_auto()
font_add_google("Noto Sans", "notosans")

theme_lph <- theme_bw(base_size = 8, base_family = "notosans") +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.25, colour = "grey92"),
    legend.position  = "bottom",
    legend.key.size  = unit(0.3, "cm"),
    legend.text      = element_text(size = 6),
    legend.title     = element_text(size = 7, face = "bold"),
    plot.title       = element_text(face = "bold.italic", size = 8),
    axis.text        = element_text(size = 6.5),
    axis.title       = element_text(size = 7.5),
    plot.tag         = element_text(face = "bold", size = 10)
  )
theme_set(theme_lph)

pal_resistome <- c("Functional" = "#CC6677", "Acquired" = "#6699CC")
pal_rh <- c("Dry"      = "#E69F00",
            "Moderate"  = "#56B4E9",
            "Humid"     = "#009E73")

analysis  <- readRDS(here("Datasets", "analysis_ready.rds"))
resistome <- readRDS(here("Datasets", "resistome_matrices.rds"))

rh_cache_path <- here("Datasets", "nasapower_rh_cache.rds")

if (file.exists(rh_cache_path)) {
  rh_data <- readRDS(rh_cache_path)
} else {
  rh_data <- list()
}

samples <- analysis %>% filter(!is.na(collection_date), !is.na(lat), !is.na(lon))
loc_date_ranges <- samples %>%
  group_by(lat, lon) %>%
  summarise(
    earliest = min(collection_date, na.rm = TRUE) - 30,
    latest   = max(collection_date, na.rm = TRUE),
    .groups  = "drop"
  )

n_locs <- nrow(loc_date_ranges)

n_new <- 0
for (i in seq_len(n_locs)) {
  row <- loc_date_ranges[i, ]
  loc_key <- paste0(round(row$lat, 4), "_", round(row$lon, 4))
  if (loc_key %in% names(rh_data)) next

  result <- tryCatch({
    years_needed <- seq(year(row$earliest), year(row$latest))
    yearly_dfs <- list()
    for (yr in years_needed) {
      start_d <- max(row$earliest, as.Date(paste0(yr, "-01-01")))
      end_d   <- min(row$latest,   as.Date(paste0(yr, "-12-31")))
      d <- get_power(community = "ag",
                     lonlat    = c(row$lon, row$lat),
                     pars      = "RH2M",
                     dates     = c(as.character(start_d), as.character(end_d)),
                     temporal_api = "daily")
      yearly_dfs[[as.character(yr)]] <- d %>%
        select(LON, LAT, YYYYMMDD, RH2M) %>%
        rename(date = YYYYMMDD, rh = RH2M)
    }
    bind_rows(yearly_dfs)
  }, error = function(e) { warning(e$message); NULL })

  if (!is.null(result) && nrow(result) > 0) {
    rh_data[[loc_key]] <- result
    n_new <- n_new + 1
  }

  if (n_new %% 25 == 0 && n_new > 0) saveRDS(rh_data, rh_cache_path)
  Sys.sleep(0.3)
}
saveRDS(rh_data, rh_cache_path)

daily_rh <- bind_rows(rh_data, .id = "loc_key") %>%
  mutate(
    lat = as.numeric(str_extract(loc_key, "^[^_]+")),
    lon = as.numeric(str_extract(loc_key, "[^_]+$"))
  ) %>%
  select(-loc_key, -LON, -LAT)

rh_windowed <- samples %>%
  mutate(loc_key = paste0(round(lat, 4), "_", round(lon, 4))) %>%
  distinct(genepid, loc_key, collection_date) %>%
  mutate(RH2M_30d = NA_real_)

daily_rh_split <- daily_rh %>%
  mutate(loc_key = paste0(round(lat, 4), "_", round(lon, 4))) %>%
  group_by(loc_key) %>% group_split() %>%
  set_names(map_chr(., ~ .x$loc_key[1]))

for (i in seq_len(nrow(rh_windowed))) {
  key   <- rh_windowed$loc_key[i]
  sdate <- rh_windowed$collection_date[i]
  if (!key %in% names(daily_rh_split)) next
  city <- daily_rh_split[[key]]
  w30  <- city$rh[city$date >= (sdate - 30) & city$date <= sdate]
  if (length(w30) > 0) rh_windowed$RH2M_30d[i] <- mean(w30, na.rm = TRUE)
}

analysis <- analysis %>%
  left_join(rh_windowed %>% select(genepid, RH2M_30d), by = "genepid")

run_humidity_permanova <- function(clr_mat, climate_df, label) {
  clr_mat    <- clr_mat    %>% mutate(genepid = as.character(genepid))
  climate_df <- climate_df %>% mutate(genepid = as.character(genepid))
  shared     <- intersect(clr_mat$genepid, climate_df$genepid)

  clr_sub  <- clr_mat    %>% filter(genepid %in% shared) %>% arrange(genepid)
  clim_sub <- climate_df %>% filter(genepid %in% shared) %>% arrange(genepid)

  keep <- complete.cases(clim_sub[, c("T_30d", "RH2M_30d", "Region",
                                       "gdp_pcap_ppp", "sanitation",
                                       "health_exp_gdp", "oop_health_exp",
                                       "immunization_dpt", "animal_amc_mgkg",
                                       "pop_density", "year")])
  clr_sub  <- clr_sub[keep, ]
  clim_sub <- clim_sub[keep, ]

  dist_mat <- dist(clr_sub %>% select(-genepid))

  clim_sub <- clim_sub %>%
    mutate(log_gdp         = log(gdp_pcap_ppp),
           log_pop_density = log(pop_density + 1),
           log_animal_amc  = log(animal_amc_mgkg + 1),
           year_f          = factor(year),
           RH_tertile      = cut(RH2M_30d,
                                 breaks = quantile(RH2M_30d, c(0, 1/3, 2/3, 1)),
                                 labels = c("Dry", "Moderate", "Humid"),
                                 include.lowest = TRUE))

  perm_main <- adonis2(dist_mat ~ T_30d + RH2M_30d +
                         log_gdp + sanitation + health_exp_gdp +
                         oop_health_exp + immunization_dpt +
                         log_pop_density + log_animal_amc + Region + year_f,
                       data = clim_sub, permutations = 999, by = "margin")

  perm_int <- adonis2(dist_mat ~ T_30d * RH2M_30d +
                        log_gdp + sanitation + health_exp_gdp +
                        oop_health_exp + immunization_dpt +
                        log_pop_density + log_animal_amc + Region + year_f,
                      data = clim_sub, permutations = 999, by = "margin")

  perm_ref <- adonis2(dist_mat ~ T_30d +
                        log_gdp + sanitation + health_exp_gdp +
                        oop_health_exp + immunization_dpt +
                        log_pop_density + log_animal_amc + Region + year_f,
                      data = clim_sub, permutations = 999, by = "margin")

  strat_results <- list()
  for (tert in c("Dry", "Moderate", "Humid")) {
    idx <- which(clim_sub$RH_tertile == tert)
    if (length(idx) < 20) next
    dist_sub <- as.dist(as.matrix(dist_mat)[idx, idx])
    perm_s <- adonis2(dist_sub ~ T_30d + log_gdp + sanitation +
                        health_exp_gdp + oop_health_exp + immunization_dpt +
                        log_pop_density + log_animal_amc + Region + year_f,
                      data = clim_sub[idx, ], permutations = 999, by = "margin")
    strat_results[[tert]] <- list(
      perm = perm_s,
      n    = length(idx),
      r2   = perm_s["T_30d", "R2"],
      p    = perm_s["T_30d", "Pr(>F)"]
    )
  }

  pcoa <- cmdscale(dist_mat, k = 2, eig = TRUE)
  pcoa_df <- tibble(
    genepid    = clr_sub$genepid,
    PCoA1      = pcoa$points[, 1],
    PCoA2      = pcoa$points[, 2],
    T_30d      = clim_sub$T_30d,
    RH2M_30d   = clim_sub$RH2M_30d,
    RH_tertile = clim_sub$RH_tertile
  )
  eig_pct <- 100 * pcoa$eig / sum(pcoa$eig[pcoa$eig > 0])

  list(
    perm_ref     = perm_ref,
    perm_main    = perm_main,
    perm_int     = perm_int,
    stratified   = strat_results,
    pcoa_df      = pcoa_df,
    eig_pct      = eig_pct,
    n            = nrow(clr_sub),
    label        = label
  )
}

res_fg  <- run_humidity_permanova(resistome$fg_cluster_clr,  analysis, "Functional (FG)")
res_acq <- run_humidity_permanova(resistome$acq_cluster_clr, analysis, "Acquired")

saveRDS(list(fg = res_fg, acq = res_acq),
        here("Results", "sensitivity_humidity.rds"))

make_panel <- function(res, panel_tag, resistome_col) {
  df <- res$pcoa_df
  eig1 <- round(res$eig_pct[1], 1)

  strat <- res$stratified
  ann_text <- map_dfr(names(strat), function(tert) {
    tibble(
      RH_tertile = tert,
      label = sprintf("R\u00B2 = %.3f\n%s",
                      strat[[tert]]$r2,
                      ifelse(strat[[tert]]$p <= 0.001, "p \u2264 0\u00B7001",
                             sprintf("p = %.3f", strat[[tert]]$p)))
    )
  })

  int_p <- res$perm_int["T_30d:RH2M_30d", "Pr(>F)"]
  int_r2 <- res$perm_int["T_30d:RH2M_30d", "R2"]
  int_label <- sprintf("Interaction: R\u00B2 = %.4f, %s",
                       int_r2,
                       ifelse(int_p <= 0.001, "p \u2264 0\u00B7001",
                              sprintf("p = %0.3f", int_p)))

  df$RH_tertile <- factor(df$RH_tertile, levels = c("Dry", "Moderate", "Humid"))

  ggplot(df, aes(x = T_30d, y = PCoA1, colour = RH_tertile)) +
    geom_point(alpha = 0.5, size = 0.8) +
    geom_smooth(method = "lm", se = TRUE, linewidth = 0.8, alpha = 0.15) +
    facet_wrap(~ RH_tertile, nrow = 1) +
    scale_colour_manual(values = pal_rh, guide = "none") +
    labs(
      tag   = panel_tag,
      title = res$label,
      x     = "30-day mean temperature (\u00B0C)",
      y     = sprintf("Resistome PCoA axis 1 (%.1f%%)", eig1),
      subtitle = int_label
    ) +
    theme(
      strip.background = element_rect(fill = "grey95", colour = NA),
      strip.text       = element_text(face = "bold", size = 7),
      plot.subtitle    = element_text(size = 6.5, colour = "grey40")
    )
}

p_fg  <- make_panel(res_fg,  "A", pal_resistome["Functional"])
p_acq <- make_panel(res_acq, "B", pal_resistome["Acquired"])

strat_summary <- bind_rows(
  map_dfr(names(res_fg$stratified), function(t) {
    tibble(Resistome = "Functional", Tertile = t,
           R2 = res_fg$stratified[[t]]$r2,
           p  = res_fg$stratified[[t]]$p,
           n  = res_fg$stratified[[t]]$n)
  }),
  map_dfr(names(res_acq$stratified), function(t) {
    tibble(Resistome = "Acquired", Tertile = t,
           R2 = res_acq$stratified[[t]]$r2,
           p  = res_acq$stratified[[t]]$p,
           n  = res_acq$stratified[[t]]$n)
  })
) %>%
  mutate(Tertile   = factor(Tertile, levels = c("Dry", "Moderate", "Humid")),
         Resistome = factor(Resistome, levels = c("Functional", "Acquired")),
         sig       = ifelse(p <= 0.05, "bold", "plain"))

p_summary <- ggplot(strat_summary, aes(x = Tertile, y = R2 * 100,
                                        fill = Resistome)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6, alpha = 0.85) +
  geom_text(aes(label = sprintf("%.2f%%", R2 * 100)),
            position = position_dodge(width = 0.7),
            vjust = -0.5, size = 2.2) +
  geom_text(aes(label = sprintf("n=%d", n), y = 0),
            position = position_dodge(width = 0.7),
            vjust = 1.5, size = 1.8, colour = "grey50") +
  scale_fill_manual(values = pal_resistome) +
  labs(tag = "C",
       x = "Relative humidity tertile",
       y = "Temperature R\u00B2 (%)",
       fill = NULL) +
  theme(legend.position = "bottom")

fig <- (p_fg / p_acq) | p_summary
fig <- fig + plot_layout(widths = c(3, 1))

ggsave(here("Figures", "eFig_humidity_interaction.pdf"),
       fig, width = 10, height = 6)
ggsave(here("Figures", "eFig_humidity_interaction.png"),
       fig, width = 10, height = 6, dpi = 300)
