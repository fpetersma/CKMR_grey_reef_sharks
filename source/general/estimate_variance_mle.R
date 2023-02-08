## =============================================================================
##  Name: estimate_variance_mle.R
## 
##  Description: Estimate variance through the variance-covariance matrix.
##    We will use the numDeriv package for this. 
## =============================================================================

## =============================================================================
## 1. LOAD PACKAGES AND DATA
## =============================================================================

## Packages
library(numDeriv)
library(CKMRcpp)

## Data
load("data/1_population_multiple_sampling_schemes/gestatii/simulation_100_schemes_scenario_1-49_fit_results_sim=555.RData")
load("data/1_population_multiple_sampling_schemes/gestatii/100_schemes_dfs_suff_unique_combos_sim=555.RData")

## Add data to fit results
true_error <- 2.89          # what is the true measurement error

a_0 <- -8.27                # theoretical age at length zero              
k <- 0.0554                 # growth rate
l_inf <- 163                # asymptotic length

errors <- c(1e-8,                         # 1e-8 instead of 0 to avoid overflow
            true_error*1/3,               # 1/3 of true error
            true_error*2/3,               # 2/3 of true error
            true_error,                   # true error
            true_error*4/3,               # 4/3 of true error
            true_error*5/3,               # 5/3 of true error
            true_error*6/3)               # twice the true error
## vbgf comes from l.100 onwards in 'explore_vbgf_bias_uncertainty.R'
step_l_inf <- 0.05 * l_inf  # 8.15 is 5% of 163
a0_options <- seq(from = l_inf - step_l_inf * 3, 
                  to = l_inf + step_l_inf * 3, 
                  by = step_l_inf)
vbgfs <- matrix(c(CKMRcpp::a_0_j(l_inf, a0_options, k, a_0), 
                  a0_options), ncol = 2)



## Combine every error with every vbgf parameter combo
pars <- do.call(rbind, lapply(seq_along(errors), function(i) {
  error <- errors[i]
  out <- cbind(vbgfs, rep(error, nrow(vbgfs)))
  colnames(out) <- c("a_0", "l_inf", "sigma_l")
  return(out)
}))

## Number of fitted models for every scenario
n <- 100

## :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
## Fit all models
## :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
for (i in 1:nrow(pars)) {
  a0 <- pars[i, "a_0"]
  l_inf <- pars[i, "l_inf"]
  sigma_l <- pars[i, "sigma_l"]
  
  for (j in 1:n) {
    df_select <- dfs_suff[[j]]
    ## Create the data object for nlminb()
    dat <- list(alpha_m = 17,                         # maturity age for males (van=10, ges=17)
                alpha_f = 19,                         # maturity age for females (van=10, ges=19)
                
                # r = log(1.0000),                      # the population growth rate
                sigma_l = log(sigma_l),               # the measurement error on length
                phi = boot::logit(1 - 0.1113),        # phi is the survival rate (van=0.1535, ges=0.1113)
                
                ESTIMATE_R = 2,                       # is r fixed (0), estimated (1), 
                # or sex specific (2)?
                
                max_age = 63,                         # the maximum age to be considered (van=19, ges=63)
                max_length = 200,                     # at least the maximum length in the data
                t0 = 2014,                            # a reference year for the abundance estimation
                vbgf_l_inf = l_inf,                   # asymptotic length for VBGF
                vbgf_k = 0.0554,                      # growth coefficient for VBGF
                vbgf_a0 = a0,                         # theoretical age at length 0 for vbgf
                s_i = df_select$indiv_1_sex,
                s_j = df_select$indiv_2_sex,
                c_i = df_select$indiv_1_capture_year,
                c_j = df_select$indiv_2_capture_year,
                a_i = df_select$indiv_1_capture_age,
                a_j = df_select$indiv_2_capture_age,
                l_i = df_select$indiv_1_length,
                l_j = df_select$indiv_2_length,
                kinship = df_select$kinship,
                cov_combo_freq = df_select$covariate_combo_freq,
                n = nrow(df_select))
    
    scenario_fits[[i]][[j]]$dat <- dat
    
    ## Derive the hessian and add the fit object
    hess <- numDeriv::hessian(func=CKMRcpp::nllPOPCKMRcppAgeUnknownGestation,
                              x=scenario_fits[[i]][[j]]$par,
                              dat=dat)
    scenario_fits[[i]][[j]]$hessian <- hess
    
  }
}

