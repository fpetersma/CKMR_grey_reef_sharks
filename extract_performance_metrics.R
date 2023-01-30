## =============================================================================
##  Name: extract_performance_metrics.R
##
##  Description: This script extracts a variety of performance metrics from a
##    set of fitted models. The metrics include at least: 
##      + true variance
##      + monte carlo variance
##      + mean error
##      + mean absolute error
##      + mean square error
##      + coefficient of variation
##    These metrics will be derived for male and female abundance in y0 = 2014,
##    the growth rates, and any other parameters of interest. This information 
##    will be stored in a big data.frame, which can be used to create tables 
##    for the manuscript.
## =============================================================================

## =============================================================================
## 1. LOAD THE DATA AND LIBRARIES
## =============================================================================
library(tidyverse)
library(Rfast)

load("data/1_population_multiple_sampling_schemes/vanillus/simulation_100_schemes_scenario_1-49_fit_results_sim=144.RData")
load("data/1_population_multiple_sampling_schemes/vanillus/single_combined_data_i=144.RData")

## =============================================================================
## 2. CREATE THE MASTER DATA FRAME
## =============================================================================

## Extract failed fit scenarios
conv <- sapply(scenario_fits, function(scen) {
  all(sapply(scen, function(fit) fit$message) == "relative convergence (4)")
})
names(conv) <- 1:49

## Extract the truths
true_N_m <- single_combined_data$N_hist[, "N_m"] # mature male abundance
true_N_f <- single_combined_data$N_hist[, "N_f"] # mature female abundance

true_N_m_y0 <- true_N_m[100]
true_N_f_y0 <- true_N_f[100]

true_r_m <- (true_N_m[100] - true_N_m[1]) ^ (1 / 100) / 100 + 1
true_r_f <- (true_N_f[100] - true_N_f[1]) ^ (1 / 100) / 100 + 1

## Extract the estimates
est_N_m_y0 <- sapply(scenario_fits, function(scen) { # mature male abundance
  sapply(scen, function(fit) return(exp(fit$par["N_t0_m"])))
}) 
est_N_m_y0[, !conv] <- NA # set failed convergence estimates to NA

est_N_f_y0 <- sapply(scenario_fits, function(scen) { # mature female abundance
  sapply(scen, function(fit) return(exp(fit$par["N_t0_f"])))
})  
est_N_f_y0[, !conv] <- NA # set failed convergence estimates to NA

est_r_m <- sapply(scenario_fits, function(scen) { # mature male growth rate
  sapply(scen, function(fit) return(exp(fit$par["r_m"])))
})
est_r_m[, !conv] <- NA # set failed convergence estimates to NA

est_r_f <- sapply(scenario_fits, function(scen) { # mature female growth rate
  sapply(scen, function(fit) return(exp(fit$par["r_f"])))
}) 
est_r_f[, !conv] <- NA # set failed convergence estimates to NA

## Mean error
error_N_m_y0 <- est_N_m_y0 - true_N_m_y0
error_N_f_y0 <- est_N_f_y0 - true_N_f_y0
error_r_m <- est_r_m - true_r_m
error_r_f <- est_r_f - true_r_f

median_error <- data.frame(N_f_y0 = colMedians(error_N_f_y0), 
                           N_m_y0 = colMedians(error_N_m_y0),
                           r_f = colMedians(error_r_f),
                           r_m = colMedians(error_r_m))

## Median absolute error
median_abs_error <- data.frame(N_f_y0 = colMedians(abs(error_N_f_y0)), 
                           N_m_y0 = colMedians(abs(error_N_m_y0)),
                           r_f = colMedians(abs(error_r_f)),
                           r_m = colMedians(abs(error_r_m)))

## Median squared error
median_sqrd_error <- data.frame(N_f_y0 = colMedians((error_N_f_y0) ^ 2), 
                               N_m_y0 = colMedians((error_N_m_y0) ^ 2),
                               r_f = colMedians((error_r_f) ^ 2),
                               r_m = colMedians((error_r_m) ^ 2))

## standard deviation and variance
sd <- data.frame(N_f_y0 = apply(est_N_f_y0, 2, sd),
                 N_m_y0 = apply(est_N_m_y0, 2, sd), 
                 r_f = apply(est_r_f, 2, sd), 
                 r_m = apply(est_r_m, 2, sd))
variance <- sd ^ 2

## Coefficient of variation
cv <- data.frame(N_f_y0 = apply(est_N_f_y0, 2, sd) / colmeans(est_N_f_y0),
                 N_m_y0 = apply(est_N_m_y0, 2, sd) / colmeans(est_N_m_y0), 
                 r_f = apply(est_r_f, 2, sd) / colmeans(est_r_f), 
                 r_m = apply(est_r_m, 2, sd) / colmeans(est_r_m))


## =============================================================================
## 3. CREATE FIGURES AND TABLES FOR MANUSCRIPT
## =============================================================================
