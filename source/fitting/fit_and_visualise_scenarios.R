################################################################################
##                                                                            ##
##  Script that fits XX models with incorrect measurement error and YY        ##
##  models with misspecified VBGFs, and all combinations.                     ##
##                                                                            ##
################################################################################

## :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
## Load data and libraries
## :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

library(pbapply)
library(parallel)
load("data/vanilla_gestation_repro=U(1,4)_sample_years_139-140/smaller_files/1000_sims_dfs_suff_length_sd=2_unique_combos.RData")

## :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
## Set parameter values
## :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

errors <- c(1e-8, 2, 4, 6, 8) # 1e-8 instead of 0 to avoid overflow
## vbgfs comes from l.100 onwards in 'explore_vbgf_bias_uncertainty.R'
step_l_inf <- 8.75  # 8.75 is 5% of 175
a0_options <- seq(from = 175 - step_l_inf * 3, 
                  to = 175 + step_l_inf * 3, 
                  by = step_l_inf)
vbgfs <- matrix(c(CKMRcpp::a_0_j(175, a0_options, 0.1, -3.5), 
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

# scenario_fits_1_7 <- lapply(1:nrow(pars), function(i) {
scenario_fits_1_7 <- lapply(1:7, function(i) {
  a0 <- pars[i, "a_0"]
  l_inf <- pars[i, "l_inf"]
  sigma_l <- pars[i, "sigma_l"]
  
  n_cores <- 20 # 50 fits per core
  cl <- makeCluster(n_cores)
  clusterExport(cl = cl, list("a0", "l_inf", "sigma_l"), envir = environment())
  results <- pblapply(dfs_suff[1:n], function(df) {
    ## Create parameter object for nlminb()
    par <- list(
      # phi = boot::logit(0.87), # same as plogis(0.9) -- boot::inv.logit() is qlogis()
      N_t0_m = log(500),
      r = log(1.0002),
      # sigma_l = log(0.01),
      # phi = boot::logit(1 - 0.153),
      N_t0_f = log(500))
    
    ## Take a subset of the data if required
    df_select <- df[, ]
    
    ## Create the data object for nlminb()
    dat <- list(alpha_m = 10,                         # maturity age for males
                alpha_f = 12,                         # maturity age for females
                
                # r = log(1.0002),
                sigma_l = log(sigma_l),
                phi = boot::logit(1 - 0.097),         # phi is the surival rate
                
                fixed_r = 0,                          # is r fixed or estimated?
                
                max_age = 19,
                max_length = 250,
                t0 = 140,
                vbgf_l_inf = l_inf,
                vbgf_k = 0.1,
                vbgf_a0 = a0,
                s1 = df_select$indiv_1_sex,
                s2 = df_select$indiv_2_sex,
                c1 = df_select$indiv_1_capture_year,
                c2 = df_select$indiv_2_capture_year,
                a1 = df_select$indiv_1_capture_age,
                a2 = df_select$indiv_2_capture_age,
                l1 = df_select$indiv_1_length,
                l2 = df_select$indiv_2_length,
                kinship = df_select$kinship,
                cov_combo_freq = df_select$covariate_combo_freq,
                n = nrow(df_select))
    
    # system.time(test_new <- nllPOPCKMRcppAgeUnknown(dat, par))
    # system.time(test_old <- CKMRcpp::nllPOPCKMRcppAgeUnknown(dat, par))
    
    ## Start the optimisation (make sure the trace and relative tolerance values
    ##                         are correct)
    # system.time({
    res <- nlminb(start = par, 
                  objective = CKMRcpp::nllPOPCKMRcppAgeUnknownGestation, 
                  dat = dat, 
                  control = list(trace = 0, rel.tol = 1e-10))
    # })
    
    return(res)
  }, cl = cl); stopCluster(cl) # end of pblapply() in parallel
  return(results)
})

## Save the environment
save.image("data/simulation_scenarios_results.RData")

## :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
## Visualise the results.
## 
## What should I visualise? The abundance, the 
## :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
scenario_fits[[1]][[1]]$

## Set plot dimensions
mfrow()

test <- scenario_fits[[8]]
est <- t(sapply(test, function(x) x$par))
its <- sapply(scenario_fits, function(x) {
  mean(sapply(x, function(y) y$iterations))
})
conv <- sapply(scenario_fits, function(x) {
  mean(sapply(x, function(y) y$convergence))
})

scenarios <- as.data.frame(pars)
scenarios$mean_iterations <- its
scenarios$mean_convergence <- conv

## UPDATE: when sigma_l = 0, it fails. This makes sense, as it results in a
## probability density of positive infinity (due to division by 0). 
## Change this to 1e-8.