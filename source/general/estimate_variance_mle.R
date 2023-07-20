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
# library(numDeriv)
# library(CKMRcpp)
library(pbapply)
library(parallel)

NO_GROWTH <- TRUE

## Data
load("data/simulation_study/simple/simulation_1000_schemes_all_scenarios_fit_results_sim=all_no_growth.RData")
simple_fits <- scenario_fits
load("data/simulation_study/complex/simulation_1000_schemes_all_scenarios_fit_results_sim=all_no_growth.RData")
complex_fits <- scenario_fits
rm(scenario_fits)

## -----------------------------------------------------------------------------
## Extract the variances
## -----------------------------------------------------------------------------

## Extract failed fit scenarios
conv <- sapply(simple_fits, function(scen) {
  all(sapply(scen, function(fit) fit$message) == "relative convergence (4)")
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
    par <- exp(simple_fits[[scenarios_to_keep[i]]][[j]]$par)
    var <- diag(solve(hess))
    var_real <- par ^ 2 * (exp(var) - 1) # if X=logY, then var(Y)=E(Y)^2(e^Var(X)-1)
    ## There is one fit that resulted that did not have positive definite hessina (and thus negative variance). this one was ignored (i=19 j = 1)
    if (any(simple_fits[[scenarios_to_keep[i]]][[j]]$var < 0)) print(paste(" STOP", i, j))
    
    return(sqrt(var_real))
  }), na.rm = TRUE)
  
  cv_simple[i, ] <- rowMeans(sapply(1:1000, function(j) {
    hess <- simple_fits[[scenarios_to_keep[i]]][[j]]$hess
    par <- exp(simple_fits[[scenarios_to_keep[i]]][[j]]$par)
    var <- diag(solve(hess))
    var_real <- par^2 * (exp(var) - 1) # if X=logY, then var(Y)=E(Y)^2(e^Var(X)-1)
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
    par <- exp(complex_fits[[scenarios_to_keep[i]]][[j]]$par)
    var <- diag(solve(hess))
    var_real <- par^2 * (exp(var) - 1) # if X=logY, then var(Y)=E(Y)^2(e^Var(X)-1)
    
    if (any(complex_fits[[scenarios_to_keep[i]]][[j]]$var < 0)) print(paste(" STOP", i, j))
    
    return(sqrt(var_real))
  }), na.rm = TRUE)
  
  cv_complex[i, ] <- rowMeans(sapply(1:1000, function(j) {
    hess <- complex_fits[[scenarios_to_keep[i]]][[j]]$hess
    par <- exp(complex_fits[[scenarios_to_keep[i]]][[j]]$par)
    var <- diag(solve(hess))
    var_real <- par^2 * (exp(var) - 1) # if X=logY, then var(Y)=E(Y)^2(e^Var(X)-1)
    
    if (any(complex_fits[[scenarios_to_keep[i]]][[j]]$var < 0)) print(paste(" STOP", i, j))
    
    std_dev <- sqrt(var_real)
    
    return(std_dev / par * 100)
  }), na.rm = TRUE)
}
## -----------------------------
## Combine into data frames
## -----------------------------
if (NO_GROWTH) {
  sd_df <- cbind(data.frame(paste0(rep(1:5, each=5), "-", rep(1:5, rep=5))[scenarios_to_keep]), 
                 sd_simple[, c(2,1)], sd_complex[, c(2,1)])
  cv_df <- cbind(data.frame(paste0(rep(1:5, each=5), "-", rep(1:5, rep=5))[scenarios_to_keep]),  
                 cv_simple[, c(2,1)], cv_complex[, c(2,1)])
} else {
  sd_df <- cbind(data.frame(paste0(rep(1:5, each=5), "-", rep(1:5, rep=5))[scenarios_to_keep]), 
                 sd_simple[, c(2,1,3,4)], sd_complex[, c(2,1,3,4)])
  cv_df <- cbind(data.frame(paste0(rep(1:5, each=5), "-", rep(1:5, rep=5))[scenarios_to_keep]),  
                 cv_simple[, c(2,1,3,4)], cv_complex[, c(2,1,3,4)])
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

## Create tables for latex
library(kableExtra)

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

