################################################################################
## This script combines the model based and empirical standard errors into    ##
## one nice table. This table is then saved and turned into a latex table.    ##
################################################################################

## =============================================================================
## Run option a) when analysing with recaptures, b) for only the first capture,
## and c) for only the last capture.
## =============================================================================

## option a) Old version of loading data
## =====================================

load("source/result_summaries/with_recaptures/empirical_sd_with_recaptures.RData")
emp_sd_df <- sd_df
load("source/result_summaries/with_recaptures/estimated_sd_with_recaptures.RData")
est_sd_df <- sd_df
rm(sd_df)
load("data/simulation_study/simple/simulation_1000_schemes_all_scenarios_fit_results_sim=all_no_growth.RData")
simple_fits <- scenario_fits
rm(scenario_fits)

## option b) New version of loading data without recaptures (only first capture)
## =============================================================================
load("source/result_summaries/only_first_capture/empirical_sd_only_first_capture.RData")
emp_sd_df <- sd_df
load("source/result_summaries/only_first_capture/estimated_sd_only_first_capture.RData")
est_sd_df <- sd_df
rm(sd_df)
load("data/simulation_study/only_first_capture/fit_results_simple_only_first_capture.RData")
simple_fits <- simple_fits_only_first_capture
rm(simple_fits_only_first_capture)

## option c) New version of loading data without recaptures (only last capture)
## =============================================================================
load("source/result_summaries/only_last_capture/empirical_sd_only_last_capture.RData")
emp_sd_df <- sd_df
load("source/result_summaries/only_last_capture/estimated_sd_only_last_capture.RData")
est_sd_df <- sd_df
rm(sd_df)
load("data/simulation_study/only_last_capture/fit_results_simple_only_last_capture.RData")
simple_fits <- simple_fits_only_last_capture
rm(simple_fits_only_last_capture)

## Extract failed fit scenarios
conv <- sapply(simple_fits, function(scen) {
  all(sapply(scen, function(fit) fit$message) %in% c("relative convergence (4)", 
                                                     "both X-convergence and relative convergence (5)"))
})

scenarios_to_keep <- (1:25)[conv] # convert TRUE/FALSE to indices

scen_names <- paste0(rep(paste0(rep("ME", 5), c("-67", "-33", "+0", "+33", "+67")),
                         each = 5), ":",
                     rep(paste0(rep("GC", 5), c("-10", "-5", "+0", "+5", "+10")),
                         times = 5))



## Tables for LaTeX
library("kableExtra")

caption <- ""

## TWO TABLE VERSION WITH DIFFERENCES
## ==================================

## Create correct tables for situation without recaptures
sd_df <- est_sd_df
sd_simple_true <- emp_sd_df[, 2:3]
sd_complex_true <- emp_sd_df[, 4:5]

## Simple species
df_simple <- data.frame(sd_df[, 1], # scenario
                        sd_df[, 2], # estimated standard error (male)
                        sd_simple_true[, 1], # empirical standard error (male)
                        sd_df[, 2] - sd_simple_true[, 1], # real difference (male)
                        (sd_df[, 2] - sd_simple_true[, 1]) / sd_simple_true[, 1] * 100, # relative difference (male)
                        sd_df[, 3], # estimated standard error (female)
                        sd_simple_true[, 2], # empirical standard error (female)
                        sd_df[, 3] - sd_simple_true[, 2], # real difference (female)
                        (sd_df[, 3] - sd_simple_true[, 2]) / sd_simple_true[, 2] * 100) # relative difference (female)
kable(df_simple,
      booktabs = T, 
      format = "latex", 
      digits = 2,
      col.names = c("Scenario", 
                    rep(c("Est SE", "Emp SE", "Est SE - Emp SE", 
                          "Delta SE / Emp SE"), times = 2)),
      row.names = F, linesep = "", caption = caption) %>%
  add_header_above( c(" " = 1, "Males" = 4, "Females" = 4)) 

## Complex species
df_complex <- data.frame(sd_df[, 1], # scenario
                         sd_df[, 4], # estimated standard error (male)
                         sd_complex_true[, 1], # empirical standard error (male)
                         sd_df[, 4] - sd_complex_true[, 1], # real difference (male)
                         (sd_df[, 4] - sd_complex_true[, 1]) / sd_complex_true[, 1] * 100, # relative difference (male)
                         sd_df[, 5], # estimated standard error (female)
                         sd_complex_true[, 2], # empirical standard error (female)
                         sd_df[, 5] - sd_complex_true[, 2], # real difference (female)
                         (sd_df[, 5] - sd_complex_true[, 2]) / sd_complex_true[, 2] * 100) # relative difference (female)
kable(df_complex,
      booktabs = T, 
      format = "latex", 
      digits = 2,
      col.names = c("Scenario", 
                    rep(c("Est SE", "Emp SE", "Est SE - Emp SE", 
                          "Delta SE / Emp SE"), times = 2)),
      row.names = F, linesep = "", caption = caption) %>%
  add_header_above( c(" " = 1, "Males" = 4, "Females" = 4)) 

## SINGLE TABLE VERSIONS WITH DIFFERENCES
## ======================================

# option 1 [preferred option]
kable(data.frame(sd_df[, 1], 
                 sd_df[ ,2], sd_simple_true[ ,1], 
                 sd_df[ ,3], sd_simple_true[ ,2],
                 sd_df[ ,4], sd_complex_true[ ,1], 
                 sd_df[ ,5], sd_complex_true[ ,2]), 
      booktabs = T, 
      format = "latex", 
      digits = 2,
      col.names = c("Scenario", 
                    rep(c("Adult males", "Adult females"), 4)),
      row.names = F, linesep = "", caption = caption) %>%
  add_header_above( c(" " = 1, "Males" = 2, "Females" = 2,
                      "Males" = 2, "Females" = 2)) %>% 
  add_header_above( c(" " = 1, "Simple species" = 4, "Complex species" = 4))

# option 2
kable(data.frame(scen_names[scenarios_to_keep], 
            sd_simple_est, sd_simple_true,
            sd_complex_est, sd_complex_true), 
      booktabs = T, 
      format = "latex", 
      digits = 2,
      col.names = c("Scenario", 
                    rep(c("Adult males", "Adult females"), 4)),
                    row.names = F, linesep = "", caption = caption) %>%
        add_header_above( c(" " = 1, "Estimated" = 2, "Empirical" = 2,
                            "Estimated" = 2, "Empirical" = 2)) %>% 
        add_header_above( c(" " = 1, "Simple species" = 4, "Complex species" = 4))



