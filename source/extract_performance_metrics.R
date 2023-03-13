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

load("data/1_population_multiple_sampling_schemes/gestatii/simulation_100_schemes_scenario_1-49_fit_results_sim=555.RData")
load("data/1_population_multiple_sampling_schemes/gestatii/single_combined_data_i=555.RData")

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

true_r_m <- (true_N_m[100] - true_N_m[1]) ^ (1 / 99) / 100 + 1
true_r_f <- (true_N_f[100] - true_N_f[1]) ^ (1 / 99) / 100 + 1

## Extract the estimated abundances in y0
est_N_m_y0 <- sapply(scenario_fits, function(scen) { # mature male abundance
  sapply(scen, function(fit) return(exp(fit$par["N_t0_m"])))
}) 
est_N_m_y0[, !conv] <- NA # set failed convergence estimates to NA

est_N_f_y0 <- sapply(scenario_fits, function(scen) { # mature female abundance
  sapply(scen, function(fit) return(exp(fit$par["N_t0_f"])))
})  
est_N_f_y0[, !conv] <- NA # set failed convergence estimates to NA

## Extract the estimated growth rates
est_r_m <- sapply(scenario_fits, function(scen) { # mature male growth rate
  sapply(scen, function(fit) return(exp(fit$par["r_m"])))
})
est_r_m[, !conv] <- NA # set failed convergence estimates to NA

est_r_f <- sapply(scenario_fits, function(scen) { # mature female growth rate
  sapply(scen, function(fit) return(exp(fit$par["r_f"])))
}) 
est_r_f[, !conv] <- NA # set failed convergence estimates to NA

## ------------------------------------------------------
## Calculate the abundance averages for the past 10 years 
## ------------------------------------------------------
## The male abundance
est_N_m_10 <- sapply(scenario_fits, function(scen) { # mature female abundance
  vec <- sapply(scen, function(fit) {
    N_y0 <- exp(fit$par["N_t0_m"])
    r_10 <- exp(fit$par["r_m"]) ^ (-9:0)
    N_10 <- N_y0 * r_10
    return(mean(N_10))
  })
  ## set to NA if there was failed convergence
  if (!all(sapply(scen, function(fit) fit$message) == "relative convergence (4)")) {
    vec <- matrix(NA, nrow = 1, ncol = 100)
  }
  return(vec)
}, simplify = "array")  

error_N_m_10 <- sapply(scenario_fits, function(scen) { # mature female abundance
  ## set to NA if there was failed convergence
  if (!all(sapply(scen, function(fit) fit$message) == "relative convergence (4)")) {
    mat <- rep(NA, 100)
  } else {
    mat <- sapply(scen, function(fit) {
      N_y0 <- exp(fit$par["N_t0_m"])
      r_10 <- exp(fit$par["r_m"]) ^ (-9:0)
      N_10 <- N_y0 * r_10
      
      N_10_diff <- median(N_10 - true_N_m[91:100])
      return(N_10_diff)
    })
  }
  return(mat)
})  

## The female abundance
est_N_f_10 <- sapply(scenario_fits, function(scen) { # mature female abundance
  vec <- sapply(scen, function(fit) {
    N_y0 <- exp(fit$par["N_t0_f"])
    r_10 <- exp(fit$par["r_f"]) ^ (-9:0)
    N_10 <- N_y0 * r_10
    return(mean(N_10))
  })
  ## set to NA if there was failed convergence
  if (!all(sapply(scen, function(fit) fit$message) == "relative convergence (4)")) {
    vec <- matrix(NA, nrow = 1, ncol = 100)
  }
  return(vec)
}, simplify = "array") 

error_N_f_10 <- sapply(scenario_fits, function(scen) { # mature female abundance
  ## set to NA if there was failed convergence
  if (!all(sapply(scen, function(fit) fit$message) == "relative convergence (4)")) {
    mat <- rep(NA, 100)
  } else {
    mat <- sapply(scen, function(fit) {
      N_y0 <- exp(fit$par["N_t0_f"])
      r_10 <- exp(fit$par["r_f"]) ^ (-9:0)
      N_10 <- N_y0 * r_10
      # cat("N_10:", N_10, "\n")
      
      N_10_diff <- median(N_10 - true_N_f[91:100])
      return(N_10_diff)
    })
  }
  return(mat)
}) 

## Mean error
error_N_m_y0 <- est_N_m_y0 - true_N_m_y0
error_N_f_y0 <- est_N_f_y0 - true_N_f_y0
error_r_m <- est_r_m - true_r_m
error_r_f <- est_r_f - true_r_f

## Mean estimates
median_est <- data.frame(N_f_y0 = colMedians(est_N_f_y0), 
                             N_m_y0 = colMedians(est_N_m_y0),
                             N_f_10 = colMedians(est_N_f_10),
                             N_m_10 = colMedians(est_N_m_10),
                             r_f = colMedians(est_r_f),
                             r_m = colMedians(est_r_m))
mean_est <- data.frame(N_f_y0 = colmeans(est_N_f_y0), 
                         N_m_y0 = colmeans(est_N_m_y0),
                         N_f_10 = colmeans(est_N_f_10),
                         N_m_10 = colmeans(est_N_m_10),
                         r_f = colmeans(est_r_f),
                         r_m = colmeans(est_r_m))


## Average error
median_error <- data.frame(N_f_y0 = colMedians(error_N_f_y0), 
                           N_m_y0 = colMedians(error_N_m_y0),
                           N_f_10 = colMedians(error_N_f_10),
                           N_m_10 = colMedians(error_N_m_10),
                           r_f = colMedians(error_r_f),
                           r_m = colMedians(error_r_m))

mean_error <- data.frame(N_f_y0 = colmeans(error_N_f_y0), 
                         N_m_y0 = colmeans(error_N_m_y0),
                         N_f_10 = colmeans(error_N_f_10),
                         N_m_10 = colmeans(error_N_m_10),
                         r_f = colmeans(error_r_f),
                         r_m = colmeans(error_r_m))


## Median absolute error
median_abs_error <- data.frame(N_f_y0 = colMedians(abs(error_N_f_y0)), 
                               N_m_y0 = colMedians(abs(error_N_m_y0)),
                               N_f_10 = colMedians(abs(error_N_f_10)), 
                               N_m_10 = colMedians(abs(error_N_m_10)),
                               r_f = colMedians(abs(error_r_f)),
                               r_m = colMedians(abs(error_r_m)))

## Median squared error
median_sqrd_error <- data.frame(N_f_y0 = colMedians((error_N_f_y0) ^ 2), 
                                N_m_y0 = colMedians((error_N_m_y0) ^ 2),
                                N_f_10 = colMedians((error_N_f_10) ^ 2), 
                                N_m_10 = colMedians((error_N_m_10) ^ 2),
                                r_f = colMedians((error_r_f) ^ 2),
                                r_m = colMedians((error_r_m) ^ 2))

## standard deviation and variance
SD <- data.frame(N_f_y0 = apply(est_N_f_y0, 2, sd),
                 N_m_y0 = apply(est_N_m_y0, 2, sd), 
                 N_f_10 = apply(est_N_f_10, 2, sd),
                 N_m_10 = apply(est_N_m_10, 2, sd),
                 r_f = apply(est_r_f, 2, sd), 
                 r_m = apply(est_r_m, 2, sd))
variance <- sd ^ 2

## Coefficient of variation (works as the support for these estimates is 
## strictly positive)
CV <- data.frame(N_f_y0 = apply(est_N_f_y0, 2, raster::cv),
                 N_m_y0 = apply(est_N_m_y0, 2, raster::cv), 
                 N_f_10 = apply(est_N_f_10, 2, raster::cv),
                 N_m_10 = apply(est_N_m_10, 2, raster::cv), 
                 r_f = apply(est_r_f, 2, raster::cv), 
                 r_m = apply(est_r_m, 2, raster::cv))

## =============================================================================
## 3. CREATE FIGURES AND TABLES FOR MANUSCRIPT
## =============================================================================

## -------------------------
## Estimates and uncertainty
## -------------------------
