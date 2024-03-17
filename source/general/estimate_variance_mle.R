## =============================================================================
##  Name: estimate_variance_mle.R
## 
##  Description: Estimate variance through the variance-covariance matrix.
##    We will use the numDeriv package for this.
##  When X~N(mu, sd) and Y = logX, then var(Y) = (exp(sd^2)- 1)*exp(2mu+sd^2) 
##  See: https://blogs.sas.com/content/iml/2014/06/04/simulate-lognormal-data-with-specified-mean-and-variance.html
##       or the documentation of in R of stats::lnorm()
## =============================================================================

## =============================================================================
## 1. LOAD PACKAGES AND DATA
## =============================================================================

## Packages
# library(numDeriv)
# library(CKMRcpp)
library(pbapply)
library(parallel)
library(dplyr)
library(kableExtra)

NO_GROWTH <- TRUE

scen_names <- paste0(rep(paste0(rep("ME", 5), c("-67", "-33", "+0", "+33", "+67")), 
                         each = 5), ":",
                     rep(paste0(rep("GC", 5), c("-10", "-5", "+0", "+5", "+10")), 
                         times = 5))

## Load simulation data and fit results for simple species
load("data/simulation_study/simple/simulation_1000_schemes_all_scenarios_fit_results_sim=all_no_growth.RData")
load("data/simulation_study/simple/1000_schemes_combined_data_with_N_hist_sim=all.RData")
simple_fits <- scenario_fits
simple_sims <- combined_data

## Load simulation data and fit results for complex species
load("data/simulation_study/complex/simulation_1000_schemes_all_scenarios_fit_results_sim=all_no_growth.RData")
load("data/simulation_study/complex/1000_schemes_combined_data_with_N_hist_sim=all.RData")
complex_fits <- scenario_fits
complex_sims <- combined_data

## Clean environment
rm(scenario_fits, combined_data)

## ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
## Load data WITHOUT recaptures below
## ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
## combined data
load("data/simulation_study/combined_data_only_last_capture.RData")
complex_sims <- complex_combined
simple_sims <- simple_combined

rm(complex_combined, simple_combined)

## simple fits
load("data/simulation_study/fit_results_simple_only_last_capture.RData")
load("data/simulation_study/fit_results_complex_only_last_captures.RData")

simple_fits <- simple_fits_without_recaptures
complex_fits <- complex_fits_without_recaptures

rm(simple_fits_without_recaptures, 
   complex_fits_without_recaptures)

## -----------------------------------------------------------------------------
## Extract the variances
## -----------------------------------------------------------------------------

## Extract failed fit scenarios
conv <- sapply(simple_fits, function(scen) {
  all(sapply(scen, function(fit) fit$message) %in% c("relative convergence (4)", 
                                                     "both X-convergence and relative convergence (5)"))
})

scenarios_to_keep <- (1:25)[conv] # convert TRUE/FALSE to indices

## -----------------------------------------
## Create data frame for the simple species
## -----------------------------------------
if (NO_GROWTH) {
  sd_simple <- matrix(NA, nrow = length(scenarios_to_keep), ncol = 2)
  cv_simple <- matrix(NA, nrow = length(scenarios_to_keep), ncol = 2)
} else {
  sd_simple <- matrix(NA, nrow = length(scenarios_to_keep), ncol = 4)
  cv_simple <- matrix(NA, nrow = length(scenarios_to_keep), ncol = 4)
}

for (i in seq_along(scenarios_to_keep)) {
  sd_simple[i, ] <- rowMeans(sapply(1:1000, function(j) {
    
    hess <- simple_fits[[scenarios_to_keep[i]]][[j]]$hess
    par <- exp(simple_fits[[scenarios_to_keep[i]]][[j]]$par) # N_male, N_female
    var <- diag(solve(hess))
    
    ## NEW AND CORRECT var_real derivation
    var_real = (exp(var) - 1) * exp(2 * log(par) + var) 
    ## OLD var_real (incorrect)
    # var_real <- par ^ 2 * (exp(var) - 1) # if X=logY, then var(Y)=E(Y)^2(e^Var(X)-1)
    
    ## There is one fit that resulted that did not have positive definite hessina (and thus negative variance). this one was ignored (i=19 j = 1)
    if (any(simple_fits[[scenarios_to_keep[i]]][[j]]$var < 0)) print(paste(" STOP", i, j))
    
    return(sqrt(var_real))
  }), na.rm = TRUE)
  
  cv_simple[i, ] <- rowMeans(sapply(1:1000, function(j) {
    hess <- simple_fits[[scenarios_to_keep[i]]][[j]]$hess
    par <- exp(simple_fits[[scenarios_to_keep[i]]][[j]]$par) # N_male, N_female
    var <- diag(solve(hess))
    
    ## NEW AND CORRECT var_real derivation
    var_real = (exp(var) - 1) * exp(2 * log(par) + var) 
    ## OLD var_real (incorrect)
    # var_real <- par ^ 2 * (exp(var) - 1) # if X=logY, then var(Y)=E(Y)^2(e^Var(X)-1)
    
    ## There is one fit that resulted that did not have positive definite hessina (and thus negative variance). this one was ignored (i=19 j = 1)
    if (any(simple_fits[[scenarios_to_keep[i]]][[j]]$var < 0)) print(paste(" STOP", i, j))
    
    std_dev <- sqrt(var_real)
    
    return(std_dev / par * 100)
  }), na.rm = TRUE)
}

## -----------------------------------------
## Create data frame for the complex species
## -----------------------------------------
if (NO_GROWTH) {
  sd_complex <- matrix(NA, nrow = length(scenarios_to_keep), ncol = 2)
  cv_complex <- matrix(NA, nrow = length(scenarios_to_keep), ncol = 2)
} else {
  sd_complex <- matrix(NA, nrow = length(scenarios_to_keep), ncol = 4)
  cv_complex <- matrix(NA, nrow = length(scenarios_to_keep), ncol = 4)
}
for (i in seq_along(scenarios_to_keep)) {
  sd_complex[i, ] <- rowMeans(sapply(1:1000, function(j) {
    
    hess <- complex_fits[[scenarios_to_keep[i]]][[j]]$hess
    par <- exp(complex_fits[[scenarios_to_keep[i]]][[j]]$par) # N_male, N_female
    var <- diag(solve(hess))
    
    ## NEW AND CORRECT var_real derivation
    var_real = (exp(var) - 1) * exp(2 * log(par) + var) 
    ## OLD var_real (incorrect)
    # var_real <- par ^ 2 * (exp(var) - 1) # if X=logY, then var(Y)=E(Y)^2(e^Var(X)-1)
    
    if (any(complex_fits[[scenarios_to_keep[i]]][[j]]$var < 0)) print(paste(" STOP", i, j))
    
    return(sqrt(var_real))
  }), na.rm = TRUE)
  
  cv_complex[i, ] <- rowMeans(sapply(1:1000, function(j) {
    hess <- complex_fits[[scenarios_to_keep[i]]][[j]]$hess
    par <- exp(complex_fits[[scenarios_to_keep[i]]][[j]]$par) # N_male, N_female
    var <- diag(solve(hess))
    
    ## NEW AND CORRECT var_real derivation
    var_real = (exp(var) - 1) * exp(2 * log(par) + var) 
    ## OLD var_real (incorrect)
    # var_real <- par ^ 2 * (exp(var) - 1) # if X=logY, then var(Y)=E(Y)^2(e^Var(X)-1)
    
    if (any(complex_fits[[scenarios_to_keep[i]]][[j]]$var < 0)) print(paste(" STOP", i, j))
    
    std_dev <- sqrt(var_real)
    
    return(std_dev / par * 100)
  }), na.rm = TRUE)
}
## -----------------------------
## Combine into data frames
## -----------------------------
if (NO_GROWTH) {
  sd_df <- cbind(data.frame(scen_names[scenarios_to_keep]), 
                 sd_simple[, c(1,2)], sd_complex[, c(1,2)])
  cv_df <- cbind(data.frame(scen_names[scenarios_to_keep]),  
                 cv_simple[, c(1,2)], cv_complex[, c(1,2)])
} else {
  sd_df <- cbind(data.frame(scen_names[scenarios_to_keep]), 
                 sd_simple[, c(1,2,3,4)], sd_complex[, c(1,2,3,4)])
  cv_df <- cbind(data.frame(scen_names[scenarios_to_keep]),  
                 cv_simple[, c(1,2,3,4)], cv_complex[, c(1,2,3,4)])
}
colnames(sd_df) <- c("scenario", 
                     # "simple_r_m", "simple_r_f" , 
                     "simple_N_y0_m" , "simple_N_y0_f",
                     # "complex_r_m", "complex_r_f" , 
                     "complex_N_y0_m" , "complex_N_y0_f")
colnames(cv_df) <- c("scenario", 
                     # "simple_r_m", "simple_r_f" , 
                     "simple_N_y0_m" , "simple_N_y0_f",
                     # "complex_r_m", "complex_r_f" , 
                     "complex_N_y0_m" , "complex_N_y0_f") 

## Save if needed
save(list = "sd_df", file = "source/result_summaries/estimated_sd_only_last_capture.RData")

## Create tables for latex

labels <- c("Scenario", 
            # "$r_{\male}$", 
            # "$r_{\female}$", 
            "$N^A_{male, 100}$", 
            "$N^A_{female, 100}$", 
            # "$r_{\male}$", 
            # "$r_{\female}$",
            "$N^A_{male, 100}$", 
            "$N^A_{female, 100}$")

caption <- "" 
kable(sd_df, booktabs = T, 
      format = "latex", 
      digits = 2,
      col.names = labels, 
      row.names = F, linesep = "", caption = caption) %>%
  add_header_above( c(" " = 1, "Simple species" = 4, "Complex species" = 4))

caption <- "" 
kable(cv_df, booktabs = T, 
      format = "latex", 
      digits = 2,
      col.names = labels, 
      row.names = F, linesep = "", caption = caption) %>% 
  add_header_above( c(" " = 1, "Simple species" = 4, "Complex species" = 4))

## =============================================================================
## DERIVE THE CI COVERAGE USING A NORMAL DISTR ON log(N) (SUGGESTED BY REVIEWER)
## =============================================================================

## Create the data frames to store the estimated standard deviations
if (NO_GROWTH) {
  uncertainty_complete_simple <- array(NA, dim = c(length(scenarios_to_keep), 12, 1000))
  uncertainty_complete_complex <- array(NA, dim = c(length(scenarios_to_keep), 12, 1000))
} 

## Extract the estimated standard deviations from the Hessian matrices on the link scale
for (i in seq_along(scenarios_to_keep)) {
  uncertainty_complete_simple[i, , ] <- sapply(1:1000, function(j) {
    hess <- simple_fits[[scenarios_to_keep[i]]][[j]]$hess
    N_hat <- simple_fits[[scenarios_to_keep[i]]][[j]]$par
    var <- diag(solve(hess))
    ## There is one fit that resulted that did not have positive definite hessian (and thus negative variance). this one was ignored (i=19 j = 1)
    if (any(simple_fits[[scenarios_to_keep[i]]][[j]]$var < 0)) print(paste(" STOP", i, j))
    ci_m <- exp(qnorm(c(0.025, 0.975), N_hat[1], sqrt(var[1]))) # CI for male abundance
    ci_f <- exp(qnorm(c(0.025, 0.975), N_hat[2], sqrt(var[2]))) # CI for female abundance
    
    ## Extract true N and check if in CI
    true_N <- simple_sims[[j]]$N_hist[100, ] # extract true N in year 100
    N_m_in_ci <- true_N[1] > ci_m[1] & true_N[1] < ci_m[2] # is male abundance in CI?
    N_f_in_ci <- true_N[2] > ci_f[1] & true_N[2] < ci_f[2] # is female abundance in CI?
    return(c(N_m_link = N_hat[1], sd_m_link = sqrt(var[1]), 
             ci_m_lower = ci_m[1], ci_m_upper = ci_m[2],
             N_m_true = true_N[1], N_m_in_ci = as.numeric(N_m_in_ci),
             N_f_link = N_hat[2], sd_f_link = sqrt(var[2]),
             ci_f_lower = ci_f[1], ci_f_upper = ci_f[2],
             N_f_true = true_N[2], N_m_in_ci = as.numeric(N_f_in_ci))
    )
  })
  uncertainty_complete_complex[i, , ] <- sapply(1:1000, function(j) {
    hess <- complex_fits[[scenarios_to_keep[i]]][[j]]$hess
    N_hat <- complex_fits[[scenarios_to_keep[i]]][[j]]$par
    var <- diag(solve(hess))
    ## There is one fit that resulted that did not have positive definite hessian (and thus negative variance). this one was ignored (i=19 j = 1)
    if (any(complex_fits[[scenarios_to_keep[i]]][[j]]$var < 0)) print(paste(" STOP", i, j))
    ci_m <- exp(qnorm(c(0.025, 0.975), N_hat[1], sqrt(var[1]))) # CI for male abundance
    ci_f <- exp(qnorm(c(0.025, 0.975), N_hat[2], sqrt(var[2]))) # CI for female abundance
    
    ## Extract true N and check if in CI
    true_N <- complex_sims[[j]]$N_hist[100, ] # extract true N in year 100
    N_m_in_ci <- true_N[1] > ci_m[1] & true_N[1] < ci_m[2] # is male abundance in CI?
    N_f_in_ci <- true_N[2] > ci_f[1] & true_N[2] < ci_f[2] # is female abundance in CI?
    return(c(N_m_link = N_hat[1], sd_m_link = sqrt(var[1]), 
             ci_m_lower = ci_m[1], ci_m_upper = ci_m[2],
             N_m_true = true_N[1], N_m_in_ci = as.numeric(N_m_in_ci),
             N_f_link = N_hat[2], sd_f_link = sqrt(var[2]),
             ci_f_lower = ci_f[1], ci_f_upper = ci_f[2],
             N_f_true = true_N[2], N_m_in_ci = as.numeric(N_f_in_ci))
    )
  })
}

## Check how many times the true abundance was in this CI.
# Simple species first
coverage <- as.matrix(apply(uncertainty_complete_simple, c(1), function(df) {
  return(c(male = sum(df[6, ]) / 1000 * 100, 
           female = sum(df[12, ]) / 1000 * 100))
}))
ci_coverage_simple <- data.frame(scen = scen_names[scenarios_to_keep], 
                                 ci_coverage_male = coverage[1, ],
                                 ci_coverage_female = coverage[2, ])

# Complex species second
coverage <- as.matrix(apply(uncertainty_complete_complex, c(1), function(df) {
  return(c(male = sum(df[6, ]) / 1000 * 100, 
           female = sum(df[12, ]) / 1000 * 100))
}))
ci_coverage_complex <- data.frame(scen = scen_names[scenarios_to_keep], 
                                  ci_coverage_male = coverage[1, ],
                                  ci_coverage_female = coverage[2, ])

## Save results
save(file = "source/result_summaries/95_CI_coverage_only_last_capture.RData", 
     list = c("ci_coverage_complex", "ci_coverage_simple"))

## Tables for LaTeX
caption <- ""

kable(cbind(ci_coverage_simple, ci_coverage_complex[, -1]), 
      booktabs = T, 
      format = "latex", 
      digits = 2,
      col.names = c("Scenario", 
                    "Adult males", "Adult females", 
                    "Adult males", "Adult females"), 
      row.names = F, linesep = "", caption = caption) %>%
  add_header_above( c(" " = 1, "Simple species" = 2, "Complex species" = 2))

################################################################################
## A version based on Distance and Burnham 1987 is presented below            ##
## This one uses a different way to compute log normal intervals              ##
##
## First tries seem to suggest that the intervals are the same as previous
################################################################################

## Create the data frames to store the estimated standard deviations
if (NO_GROWTH) {
  uncertainty_burnham_simple <- array(NA, dim = c(length(scenarios_to_keep), 12, 1000))
  uncertainty_burnham_complex <- array(NA, dim = c(length(scenarios_to_keep), 12, 1000))
} 

## Extract the estimated standard deviations from the Hessian matrices on the real scale
for (i in seq_along(scenarios_to_keep)) {
  uncertainty_burnham_simple[i, , ] <- sapply(1:1000, function(j) {
    hess <- simple_fits[[scenarios_to_keep[i]]][[j]]$hess
    N_hat <- exp(simple_fits[[scenarios_to_keep[i]]][[j]]$par)
    var <- diag(solve(hess))
    var_real <- N_hat ^ 2 * (exp(var) - 1) # if X=logY, then var(Y)=E(Y)^2(e^Var(X)-1)
    ## There is one fit that resulted that did not have positive definite hessian (and thus negative variance). this one was ignored (i=19 j = 1)
    if (any(simple_fits[[scenarios_to_keep[i]]][[j]]$var < 0)) print(paste(" STOP", i, j))
    
    # Derive C from Distance book (Buckland et al., 1993; page 118), for males and females
    C <- exp(1.96 * sqrt(log(1 + (sqrt(var_real) / N_hat) ^ 2)))
    
    ci_m <- c(lower = N_hat[1] / C[1], upper = N_hat[1] * C[1]) # CI for male abundance
    ci_f <- c(lower = N_hat[2] / C[2], upper = N_hat[2] * C[2]) # CI for female abundance
    
    ## Extract true N and check if in CI
    true_N <- simple_sims[[j]]$N_hist[100, ] # extract true N in year 100
    N_m_in_ci <- true_N[1] > ci_m[1] & true_N[1] < ci_m[2] # is male abundance in CI?
    N_f_in_ci <- true_N[2] > ci_f[1] & true_N[2] < ci_f[2] # is female abundance in CI?
    return(c(N_m_real = N_hat[1], sd_m_real = sqrt(var[1]), 
             ci_m_lower = ci_m[1], ci_m_upper = ci_m[2],
             N_m_true = true_N[1], N_m_in_ci = as.numeric(N_m_in_ci),
             N_f_real = N_hat[2], sd_f_real = sqrt(var[2]),
             ci_f_lower = ci_f[1], ci_f_upper = ci_f[2],
             N_f_true = true_N[2], N_m_in_ci = as.numeric(N_f_in_ci))
    )
  })
  uncertainty_burnham_complex[i, , ] <- sapply(1:1000, function(j) {
    hess <- complex_fits[[scenarios_to_keep[i]]][[j]]$hess
    N_hat <- exp(complex_fits[[scenarios_to_keep[i]]][[j]]$par)
    var <- diag(solve(hess))
    var_real <- N_hat ^ 2 * (exp(var) - 1) # if X=logY, then var(Y)=E(Y)^2(e^Var(X)-1)
    ## There is one fit that resulted that did not have positive definite hessian (and thus negative variance). this one was ignored (i=19 j = 1)
    if (any(complex_fits[[scenarios_to_keep[i]]][[j]]$var < 0)) print(paste(" STOP", i, j))
    
    # Derive C from Distance book (Buckland et al., 1993; page 118), for males and females
    C <- exp(1.96 * sqrt(log(1 + (sqrt(var_real) / N_hat) ^ 2)))
    
    ci_m <- c(lower = N_hat[1] / C[1], upper = N_hat[1] * C[1]) # CI for male abundance
    ci_f <- c(lower = N_hat[2] / C[2], upper = N_hat[2] * C[2]) # CI for female abundance
    
    ## Extract true N and check if in CI
    true_N <- complex_sims[[j]]$N_hist[100, ] # extract true N in year 100
    N_m_in_ci <- true_N[1] > ci_m[1] & true_N[1] < ci_m[2] # is male abundance in CI?
    N_f_in_ci <- true_N[2] > ci_f[1] & true_N[2] < ci_f[2] # is female abundance in CI?
    return(c(N_m_real = N_hat[1], sd_m_real = sqrt(var[1]), 
             ci_m_lower = ci_m[1], ci_m_upper = ci_m[2],
             N_m_true = true_N[1], N_m_in_ci = as.numeric(N_m_in_ci),
             N_f_real = N_hat[2], sd_f_real = sqrt(var[2]),
             ci_f_lower = ci_f[1], ci_f_upper = ci_f[2],
             N_f_true = true_N[2], N_m_in_ci = as.numeric(N_f_in_ci))
    )
  })
}

## Check how many times the true abundance was in this CI.
# Simple species first
coverage <- as.matrix(apply(uncertainty_burnham_simple, c(1), function(df) {
  return(c(male = sum(df[6, ]) / 1000, 
           female = sum(df[12, ]) / 1000))
}))
ci_coverage_burnham_simple <- data.frame(scen = scen_names[scenarios_to_keep], 
                                 ci_coverage_male = coverage[1, ],
                                 ci_coverage_female = coverage[2, ])

# Complex species second
coverage <- as.matrix(apply(uncertainty_burnham_complex, c(1), function(df) {
  return(c(male = sum(df[6, ]) / 1000, 
           female = sum(df[12, ]) / 1000))
}))
ci_coverage_burnham_complex <- data.frame(scen = scen_names[scenarios_to_keep], 
                                  ci_coverage_male = coverage[1, ],
                                  ci_coverage_female = coverage[2, ])

## Save results
save(file = "data/simulation_study/95_CI_coverage_burnham.RData", 
     list = c("ci_coverage_burnham_complex", "ci_coverage_burnham_simple"))

## Tables for LaTeX
caption <- ""

kable(cbind(ci_coverage_burnham_simple, ci_coverage_burnham_complex[, -1]), 
      booktabs = T, 
      format = "latex", 
      digits = 2,
      col.names = c("Scenario", 
                    "Adult males", "Adult females", 
                    "Adult males", "Adult females"), 
      row.names = F, linesep = "", caption = caption) %>%
  add_header_above( c(" " = 1, "Simple species" = 2, "Complex species" = 2))
# ## Add data to fit results
# true_error <- 2.89          # what is the true measurement error
# 
# a_0 <- -8.27                # theoretical age at length zero              
# k <- 0.0554                 # growth rate
# l_inf <- 163                # asymptotic length
# 
# errors <- c(1e-8,                         # 1e-8 instead of 0 to avoid overflow
#             true_error*1/3,               # 1/3 of true error
#             true_error*2/3,               # 2/3 of true error
#             true_error,                   # true error
#             true_error*4/3,               # 4/3 of true error
#             true_error*5/3,               # 5/3 of true error
#             true_error*6/3)               # twice the true error
# ## vbgf comes from l.100 onwards in 'explore_vbgf_bias_uncertainty.R'
# step_l_inf <- 0.05 * l_inf  # 8.15 is 5% of 163
# a0_options <- seq(from = l_inf - step_l_inf * 3, 
#                   to = l_inf + step_l_inf * 3, 
#                   by = step_l_inf)
# vbgfs <- matrix(c(CKMRcpp::a_0_j(l_inf, a0_options, k, a_0), 
#                   a0_options), ncol = 2)
# 
# 
# 
# ## Combine every error with every vbgf parameter combo
# pars <- do.call(rbind, lapply(seq_along(errors), function(i) {
#   error <- errors[i]
#   out <- cbind(vbgfs, rep(error, nrow(vbgfs)))
#   colnames(out) <- c("a_0", "l_inf", "sigma_l")
#   return(out)
# }))
# 
# ## Number of fitted models for every scenario
# n <- 100
# 
# ## :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# ## Loop through the converged models and add the data and the hessian
# ## :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# ## Create list to store the extra info in
# extra_info_list <- list()
# 
# scenarios_to_keep <- (1:49)[conv]
# 
# ## Start the loopin'
# for (i in scenarios_to_keep) {
#   a0 <- pars[i, "a_0"]
#   l_inf <- pars[i, "l_inf"]
#   sigma_l <- pars[i, "sigma_l"]
#   
#   n_cores <- 25 # 2 fits per core
#   cl <- makeCluster(n_cores)
#   clusterExport(cl = cl, list("a0", "l_inf", "sigma_l", "dfs_suff", 
#                               "scenario_fits", "i"), 
#                 envir = environment())
#   # for (j in 1:n) {
#   extra_info <- pblapply(1:n, function(j)  {
#     df_select <- dfs_suff[[j]]
#     ## Create the data object for nlminb()
#     dat <- list(alpha_m = 10,                         # maturity age for males (van=10, ges=17)
#                 alpha_f = 10,                         # maturity age for females (van=10, ges=19)
#                 
#                 # r = log(1.0000),                      # the population growth rate
#                 sigma_l = log(sigma_l),               # the measurement error on length
#                 phi = boot::logit(1 - 0.1535),        # phi is the survival rate (van=0.1535, ges=0.1113)
#                 
#                 ESTIMATE_R = 2,                       # is r fixed (0), estimated (1), both(2)
#                 # or sex specific (2)?
#                 
#                 max_age = 19,                         # the maximum age to be considered (van=19, ges=63)
#                 max_length = 200,                     # at least the maximum length in the data
#                 t0 = 2014,                            # a reference year for the abundance estimation
#                 vbgf_l_inf = l_inf,                   # asymptotic length for VBGF
#                 vbgf_k = 0.0554,                      # growth coefficient for VBGF
#                 vbgf_a0 = a0,                         # theoretical age at length 0 for vbgf
#                 s_i = df_select$indiv_1_sex,
#                 s_j = df_select$indiv_2_sex,
#                 c_i = df_select$indiv_1_capture_year,
#                 c_j = df_select$indiv_2_capture_year,
#                 a_i = df_select$indiv_1_capture_age,
#                 a_j = df_select$indiv_2_capture_age,
#                 l_i = df_select$indiv_1_length,
#                 l_j = df_select$indiv_2_length,
#                 kinship = df_select$kinship,
#                 cov_combo_freq = df_select$covariate_combo_freq,
#                 n = nrow(df_select))
#     
#     scenario_fits[[i]][[j]]$dat <- dat
#     
#     ## Derive the hessian and add the fit object
#     hess <- numDeriv::hessian(func=CKMRcpp::nllPOPCKMRcppAgeUnknownGestation,
#                               x=scenario_fits[[i]][[j]]$par,
#                               dat=dat)
#     scenario_fits[[i]][[j]]$hessian <- hess
#     
#     return(list(dat = dat, hessian = hess))
#     
#   }, cl = cl); stopCluster(cl);
#   extra_info_list[[i]] <- extra_info
# }
# save(list="extra_info_list", file="extra_info_for_fits_i=144.RData")
# 
# ## Repeating the data file is pretty redundant, but do add the variances
# for(i in 1:25) {
#   if (i %in% scenarios_to_keep) {
#     for (j in 1:1000) {
#       hess <- extra_info_list[[i]][[j]]$hess
#       var <- diag(solve(hess))
#       par <- exp(scenario_fits[[i]][[j]]$par)
#       var_real <- par^2 * (exp(var) - 1) # if X=logY, then var(Y)=E(Y)^2(e^Var(X)-1)
# 
#       scenario_fits[[i]][[j]]$hess <- hess
#       scenario_fits[[i]][[j]]$var <- var
#       scenario_fits[[i]][[j]]$var_real <- var_real
#     }
#   }
#   else {
#     for (j in 1:100) {
#       scenario_fits[[i]][[j]]$hess <- NA
#       scenario_fits[[i]][[j]]$var <- NA
#     }
# 
#   }
# }
# 
# save(list="scenario_fits", file="data/1_population_multiple_sampling_schemes/vanillus/simulation_100_schemes_scenario_1-49_fit_results_sim=555_with_estimated_variance.RData")

