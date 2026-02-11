rm(list = ls())

# ── Libraries ──────────────────────────────────────────────
library(tidyverse)
library(plm)
library(lmtest)
library(sandwich)
library(psych)
library(car)
library(broom)
library(xtable)
library(stargazer)
library(ggplot2)
library(ggpubr)
library(readxl)
library(corrplot)
library(knitr)
library(tseries)
library(dynlm)
library(vars)
library(rstan)
library(brms)
library(maptools)
library(spdep)
library(INLA)
library(RColorBrewer)

# ── Load data ──────────────────────────────────────────────
file_dataImport <- "./regression_table_growth_rate_based_onset_weeks_drivers_vitD_w29-34_2016-2019.csv"
modData <- read_csv(file_dataImport)

mydata_1 <- pdata.frame(modData, index = c("GEO_ID"))

# ── Derived variables ─────────────────────────────────────
mydata_1$temp_diff2 <- abs(as.numeric(mydata_1$temp_diff))

# ── Scale variables ────────────────────────────────────────
scale_col <- function(x) scale(as.numeric(x), center = TRUE, scale = TRUE)

# Variables scaled as-is
vars_to_scale <- c(
  "temp_diff", "temp_diff2", "temp_diff_v2",
  "infected_air_passanger_log",
  "fluxes_norm_log",
  "precip_mm_NOAA", "windspeed",
  "humidity_diff", "unusual_humidity",
  "prop_fan", "prop_fan_v2",
  "conf_flu_adjusted_norm", "conf_flu_adjusted_norm_v2",
  "conf_flu_adjusted_norm_log", "conf_flu_adjusted_norm_log_v2",
  "inprop_outprop_centered",
  "feelslike_C", "vitD_def_prop",
  "density", "pop2020",
  "IECC",
  "est_prop_houses_built_1950_2009_interpol",
  "MEAN_DWELL_TIME_MIN_HOME",
  "est_prop_multi_units_interpol",
  "est_occ_101_over_interpol",
  "est_mean_household_size_interpol",
  "AR",
  "RPL_THEMES", "RPL_THEME4",
  "Poverty_Index",
  "est_median_year_built_interpol"
)

for (v in vars_to_scale) {
  mydata_1[[v]] <- scale_col(mydata_1[[v]])
}

# Variables scaled with sign reversal (×-1 then scaled)
vars_to_reverse <- c(
  "start_week", "start_week_v2",
  "IECC21_cat", "IECC21",
  "yearly_cdd_cat", "yearly_cdd_cat_2"
)

for (v in vars_to_reverse) {
  mydata_1[[v]] <- scale_col(-1 * as.numeric(mydata_1[[v]]))
}

# ── Panel regression ──────────────────────────────────────
pdata <- pdata.frame(mydata_1, index = c("GEO_ID", "year1"))

plm_model <- plm(
  onset_idx ~
    start_week_v2 +
    conf_flu_adjusted_norm_log +
    fluxes_norm_log +
    inprop_outprop_centered +
    yearly_cdd_cat +
    est_prop_houses_built_1950_2009_interpol +
    est_prop_multi_units_interpol +
    est_mean_household_size_interpol +
    prop_fan +
    temp_diff2 +
    density +
    infected_air_passanger +
    RPL_THEME4,
  data = pdata,
  index = c("GEO_ID", "year1"),
  model = "pooling"
)

summary(plm_model)
vif(plm_model)

# ── Driscoll-Kraay spatial HAC standard errors ────────────
vcov_dk <- vcovSCC(plm_model, type = "HC1")
coefs_spatial_hac <- coeftest(plm_model, vcov. = vcov_dk)

coefs_df <- data.frame(
  Estimate  = coefs_spatial_hac[, "Estimate"],
  Std.Error = coefs_spatial_hac[, "Std. Error"],
  t.value   = coefs_spatial_hac[, "t value"],
  Pr.t      = coefs_spatial_hac[, "Pr(>|t|)"]
)
rownames(coefs_df) <- rownames(coefs_spatial_hac)

write.csv(coefs_df, file = "glm_results_onset_spatial_plm.csv", row.names = TRUE)
