################################################################################
##                                                                            ##
##  Script that fits models with incorrect measurement error and              ##
##  models with misspecified VBGFs, and all combinations.                     ##
##                                                                            ##
##  Make sure the naming is correct (e.g., with_recaptures or               ##
##  only_first_capture)                                                       ##
##                                                                            ##
##  12/03/24                                                                  ##
################################################################################

load("data/simulation_study/with_recaptures/sufficient_dfs_with_recaptures.RData")

## :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
## Load data and libraries
## :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

library(pbapply)
library(parallel)

## :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
## Set parameter values
## :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
true_error <- 2.89          # what is the true measurement error

a_0 <- -8.27                # theoretical age at length zero              
k <- 0.0554                 # growth rate
l_inf <- 163                # asymptotic length

errors <- c(
  # 1e-8,                         # 1e-8 instead of 0 to avoid overflow
  true_error*1/3,               # 1/3 of true error
  true_error*2/3,               # 2/3 of true error
  true_error,                   # true error
  true_error*4/3,               # 4/3 of true error
  true_error*5/3 #,               # 5/3 of true error
  # true_error*6/3
)               # twice the true error
## vbgf comes from l.100 onwards in 'explore_vbgf_bias_uncertainty.R'
step_l_inf <- 0.05 * l_inf  # 8.15 is 5% of 163
a0_options <- seq(from = l_inf - step_l_inf * 2, 
                  to = l_inf + step_l_inf * 2, 
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
n <- 1000

## :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
## Fit all models
## :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

## Simple fits
## ==============
simple_fits_with_recaptures <- lapply(1:nrow(pars), function(i) {
  cat("Fitting the CKMR model in scenario:" , i, "\n")
  
  a0 <- pars[i, "a_0"]
  l_inf <- pars[i, "l_inf"]
  sigma_l <- pars[i, "sigma_l"]
  
  n_cores <- 16
  cl <- makeCluster(n_cores)
  clusterExport(cl = cl, list("a0", "l_inf", "sigma_l"), envir = environment())
  results <- pblapply(simple_suff[1:n], function(df_select) {
    ## Create parameter object for nlminb()
    par <- list(
      N_t0_m = log(500),                            # number of reproductive males
      N_t0_f = log(500))                            # number of reproductive females
    
    ## Create the data object for nlminb()
    dat <- list(alpha_m = 10,                 # maturity age for males (sim=10, com=17)
                alpha_f = 10,                 # maturity age for females (sim=10, com=19)
                
                r = log(1.0000),                      # the population growth rate
                sigma_l = log(sigma_l),               # the measurement error on length
                phi = boot::logit(1 - 0.1535),        # phi is the survival rate (sim=0.1535, com=0.1113)
                
                ESTIMATE_R = 0,                       # is r fixed (0), estimated (1), 
                # or sex specific (2)?
                
                max_age = 19,                         # the maximum age to be considered (sim=19, com=63)
                max_length = 200,                     # at least the maximum length in the data
                t0 = 2014,                            # a reference year for the abundance estimation, needs to match the data
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
    
    ## Start the optimisation (make sure the trace and relative tolerance values are correct)
    # system.time({
    res <- nlminb(start = par, 
                  objective = CKMRcpp::nllPOPCKMRcppAgeUnknown, 
                  dat = dat, 
                  control = list(trace = 1, rel.tol = 1e-7))
    # })
    
    res$dat <- dat
    
    hess <- numDeriv::hessian(func = CKMRcpp::nllPOPCKMRcppAgeUnknown,
                              x = res$par,
                              dat = dat)
    res$hess <- hess
    return(res)
  }, cl = cl); stopCluster(cl) # end of pblapply() in parallel
  return(results)
})

## Save the scenario_fits object
save(list = "simple_fits_with_recaptures",
     file = "data/simulation_study/fit_results_simple_with_recaptures.RData")


## Complex fits
## ==============
complex_fits_with_recaptures <- lapply(1:nrow(pars), function(i) {
  cat("Fitting the CKMR model in scenario:" , i, "\n")
  
  a0 <- pars[i, "a_0"]
  l_inf <- pars[i, "l_inf"]
  sigma_l <- pars[i, "sigma_l"]
  
  n_cores <- 40
  cl <- makeCluster(n_cores)
  clusterExport(cl = cl, list("a0", "l_inf", "sigma_l"), envir = environment())
  results <- pblapply(complex_suff[1:n], function(df_select) {
    ## Create parameter object for nlminb()
    par <- list(
      N_t0_m = log(500),                            # number of reproductive males
      N_t0_f = log(500))                            # number of reproductive females
    
    ## Create the data object for nlminb()
    dat <- list(alpha_m = 17,                 # maturity age for males (sim=10, com=17)
                alpha_f = 19,                 # maturity age for females (sim=10, com=19)
                
                r = log(1.0000),                      # the population growth rate
                sigma_l = log(sigma_l),               # the measurement error on length
                phi = boot::logit(1 - 0.1113),        # phi is the survival rate (sim=0.1535, com=0.1113)
                
                ESTIMATE_R = 0,                       # is r fixed (0), estimated (1), 
                # or sex specific (2)?
                
                max_age = 63,                         # the maximum age to be considered (sim=19, com=63)
                max_length = 200,                     # at least the maximum length in the data
                t0 = 2014,                            # a reference year for the abundance estimation, needs to match the data
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
    
    ## Start the optimisation (make sure the trace and relative tolerance values
    ##                         are correct)
    # system.time({
    res <- nlminb(start = par, 
                  objective = CKMRcpp::nllPOPCKMRcppAgeUnknownGestation, 
                  dat = dat, 
                  control = list(trace = 1, rel.tol = 1e-7))
    # })
    
    res$dat <- dat
    
    hess <- numDeriv::hessian(func = CKMRcpp::nllPOPCKMRcppAgeUnknownGestation,
                              x = res$par,
                              dat = dat)
    res$hess <- hess
    return(res)
  }, cl = cl); stopCluster(cl) # end of pblapply() in parallel
  return(results)
})

## Save the scenario_fits object
save(list = "complex_fits_with_recaptures",
     file = "data/simulation_study/fit_results_complex_with_recaptures.RData")
