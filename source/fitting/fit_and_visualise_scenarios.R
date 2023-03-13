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
data_folder <- "data/1_population_multiple_sampling_schemes/"
load(paste0(data_folder, 
            "gestatii/100_schemes_dfs_suff_unique_combos_sim=555.RData"))

## :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
## Set parameter values
## :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
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

scenario_fits <- lapply(1:nrow(pars), function(i) {
# scenario_fits <- lapply(25, function(i) {
  cat("Fitting the CKMR model in scenario:" , i, "\n")
  
  a0 <- pars[i, "a_0"]
  l_inf <- pars[i, "l_inf"]
  sigma_l <- pars[i, "sigma_l"]
  
  n_cores <- 50 # 2 fits per core
  cl <- makeCluster(n_cores)
  clusterExport(cl = cl, list("a0", "l_inf", "sigma_l"), envir = environment())
  results <- pblapply(dfs_suff[1:n], function(df) {
    ## Create parameter object for nlminb()
    par <- list(
      # r = log(1.000),                             # same r for both sexes
      r_f = log(1.000),                             # female growth rate
      r_m = log(1.000),                             # male growth rate
      # phi = boot::logit(0.87), # same as plogis(0.9) -- boot::inv.logit() is qlogis()
      N_t0_m = log(500),                            # number of reproductive males
      # sigma_l = log(0.01),
      # phi = boot::logit(1 - 0.153),
      N_t0_f = log(500))                            # number of reproductive females
    
    ## Take a subset of the data if required
    df_select <- df[, ]
    
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
    
    ## Start the optimisation (make sure the trace and relative tolerance values
    ##                         are correct)
    # system.time({
    res <- nlminb(start = par, 
                  objective = CKMRcpp::nllPOPCKMRcppAgeUnknownGestation, 
                  dat = dat, 
                  control = list(trace = 1, rel.tol = 1e-7))
    # })
    res$dat <- dat
    
    return(res)
  }, cl = cl); stopCluster(cl) # end of pblapply() in parallel
  return(results)
})


## Save the scenario_fits object
save(list = c("scenario_fits"),
     file = paste0(data_folder, "simulation_100_schemes_scenario_25_fit_results_sim=.RData"))

## :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
## Visualise the results.
## 
## What should I visualise? The abundance, the 
## :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# scenario_fits[[1]][[1]]$
# 
# ## Set plot dimensions
# mfrow()
# 
# test <- scenario_fits[[8]]
# est <- t(sapply(test, function(x) x$par))
# its <- sapply(scenario_fits, function(x) {
#   mean(sapply(x, function(y) y$iterations))
# })
# conv <- sapply(scenario_fits, function(x) {
#   mean(sapply(x, function(y) y$convergence))
# })
# 
# scenarios <- as.data.frame(pars)
# scenarios$mean_iterations <- its
# scenarios$mean_convergence <- conv

## UPDATE: when sigma_l = 0, it fails. This makes sense, as it results in a
## probability density of positive infinity (due to division by 0). 
## Change this to 1e-8.
result_list <- scenario_fits[[25]]
# result_list <- scenario_fits_25[[1]]
## Look at the estimates
N_est <- t(sapply(result_list, function(res){
  return(c(N_male = exp(res$par["N_t0_m"]), N_female = exp(res$par["N_t0_f"])))
}))

r_est <- sapply(result_list, function(res){
  exp(res$par["r"])
})

summary(cbind(N_est, r_est))
CKMRcpp::plotCKMRabundance(result_list, c(-30, 0), y0=2014, med=T)
lines(1984:2014, combined_data[[1]]$N_hist[70:100, 1])
lines(1984:2014, combined_data[[1]]$N_hist[70:100, 2])


save(list = "result_list", 
     file = paste0(data_folder, "simulation_scenarios_results_100_schemes_i=",
                   sim_i, ".RData"))

MCE_data <- t(sapply(seq(from = 1, to = 901, by = 100), function(x) {
  mean_male <- mean(N_est[x:(x+99), 1])
  mean_female <- mean(N_est[x:(x+99), 2])
  mean_r <- mean(r_est[x:(x+99)])
  
  return(c(mean_male = mean_male, mean_female = mean_female, mean_r = mean_r))
}))

MCE_100 <- apply(MCE_data, 2, sd)
MCE_1000 <- MCE_100 / sqrt(10)
MCE_100
MCE_1000
## :::::::::::::::::::::::::::::::::::::::::::::::::::
## Plot multiple scenarios
## :::::::::::::::::::::::::::::::::::::::::::::::::::
data_folder <- "data/1_population_multiple_sampling_schemes/"
# load(paste0(data_folder, "vanillus/simulation_100_schemes_scenario_1-49_fit_results_sim=144.RData"))
# load(paste0(data_folder, "vanillus/single_combined_data_i=144.RData"))

load(paste0(data_folder, "gestatii/simulation_100_schemes_scenario_1-49_fit_results_sim=555.RData"))
load(paste0(data_folder, "gestatii/single_combined_data_i=555.RData"))

scenarios_to_drop <- c(1:7,
                       seq(from=8, to=36, by=7), 
                       seq(from=14, to=42, by=7),
                       43:49)
scenarios_to_keep <- c(1:49)[-scenarios_to_drop]

## Create male plot
p_m <- CKMRcpp::plotCKMRabundance(scenario_fits[scenarios_to_keep],
                               c(-20, 0),
                               y0 = 2014, 
                               sex = "male", 
                                y_axis = "Number of mature males",
                               truth = single_combined_data$N_hist)
## Create female plot
p_f <- CKMRcpp::plotCKMRabundance(scenario_fits[scenarios_to_keep],
                                  year_lim = c(-20, 0),
                                  y0 = 2014, 
                                  sex = "female", 
                                  y_axis = "Number of mature females",
                                  truth = single_combined_data$N_hist)

## Print the plots
p_m
p_f

## Save the plots
save(list = c("p_m", "p_f", "single_combined_data", "scenario_fits"), 
     file = "data/vanilla_abundance_plots_5x5.RData")

## Save the plots as svg with dimensions 1100x700