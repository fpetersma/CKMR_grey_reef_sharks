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
data_folder <- "data/1_population_multiple_sampling_schemes/vanillus/"
load(paste0(data_folder, "100_schemes_dfs_suff_unique_combos_sim=144.RData"))
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
vbgfs <- matrix(c(CKMRcpp::a_0_j(l_inf, a0_options, 0.0554, -8.27), 
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
# scenario_fits_25 <- lapply(25, function(i) {
  cat("Fitting the CKMR model in scenario:" , i, "\n")
  
  a0 <- pars[i, "a_0"]
  l_inf <- pars[i, "l_inf"]
  sigma_l <- pars[i, "sigma_l"]
  
  n_cores <- 25 # 50 fits per core
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
    dat <- list(alpha_m = 10,                         # maturity age for males
                alpha_f = 10,                         # maturity age for females
                
                # r = log(1.0000),                      # the population growth rate
                sigma_l = log(sigma_l),               # the measurement error on length
                phi = boot::logit(1 - 0.1535),        # phi is the survival rate
                
                ESTIMATE_R = 2,                       # is r fixed (0), estimated (1), 
                                                      # or sex specific (2)?
    
                max_age = 19,                         # the maximum age to be considered
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
                  objective = CKMRcpp::nllPOPCKMRcppAgeUnknown, 
                  dat = dat, 
                  control = list(trace = 1, rel.tol = 1e-7))
    # })
    
    return(res)
  }, cl = cl); stopCluster(cl) # end of pblapply() in parallel
  return(results)
})


## Save the environment
save(list = c("scenario_fits"),
     file = "data/simulation_100_schemes_scenario_1-49_fit_results_sim=144.RData")

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
load(paste0(data_folder, "simulation_100_schemes_scenario_1-49_fit_results_sim=144.RData"))
load(paste0(data_folder, "100_schemes_combined_data_with_N_hist_sim=144.RData"))

p_m <- CKMRcpp::plotCKMRabundancePretty(scenario_fits[c(9:13, 16:20, 23:27, 30:34, 37:41)],
                               c(-20, 0),
                               y0 = 2014, 
                               sex = "male", 
                               truth = combined_data[[1]]$N_hist)

p_f <- CKMRcpp::plotCKMRabundancePretty(scenario_fits[c(9:13, 16:20, 23:27, 30:34, 37:41)],
                               c(-20, 0),
                               y0 = 2014, 
                               sex = "female", 
                               truth = combined_data[[1]]$N_hist)
p_m
p_f


plotCKMRabundancePretty(fits_list = scenario_fits[c(30)],
                                 year_lim = c(-20, 0),
                                 y0 = 2014, 
                                 sex = "female", 
                                 truth = combined_data[[1]]$N_hist)


p <- egg::ggarrange(p_f +  theme(axis.text.x = element_blank(),
                                 axis.ticks.x = element_blank(),
                                 axis.title.x = element_blank()),
                    p_m , nrow = 2)
p
## Done. Probably best to keep the two plots separate, not combined.
## Next things: plot with the different vbgf curves, and summary statistics of 
## bias, error, things like that. 

# BELOW IS OLD STUFF
# plots <- list()
# for (i in 1:49) {
#   # readline(prompt="Press [Enter] for next plot, or [Esc] to exit.")
# 
#   plots[[i]] <- plotCKMRabundancePretty(fits = scenario_fits[[i]],
#                                    year_lim = c(-20,0),
#                                    max_y_axis = 3000,
#                                    truth = combined_data[[1]]$N_hist,
#                                    y0 = 2014)
#   # p1 <- CKMRcpp::plotCKMRabundance(fits = scenario_fits[[i]],
#   #                                  year_lim = c(-20,0),
#   #                                  max_y_axis = 3000,
#   #                                  fixed_r = NULL,
#   #                                  med = T,
#   #                                  truth = combined_data[[1]]$N_hist,
#   #                                  y0 = 2014)
# 
# 
# }
# gridExtra::grid.arrange(grobs = plots, ncol = 7, nrow = 7)
# 
# ## What if I only look at a 5x5 grid
# fits_5by5 <- scenario_fits[c(9:13, 16:20, 23:27, 30:34, 37:41)]
# plots <- list()
# for (i in 1:25) {
#   # readline(prompt="Press [Enter] for next plot, or [Esc] to exit.")
#   
#   plots[[i]] <- plotCKMRabundancePretty(fits = fits_5by5[[i]],
#                                         year_lim = c(-20,0),
#                                         max_y_axis = 3000,
#                                         truth = combined_data[[1]]$N_hist,
#                                         y0 = 2014)
#   
#   
# }
# gridExtra::grid.arrange(grobs = plots, ncol = 5)
# egg::ggarrange(plots = plots, ncol = 5)
# 
# # And now for male
# plots <- list()
# for (i in 1:25) {
#   # readline(prompt="Press [Enter] for next plot, or [Esc] to exit.")
#   
#   plots[[i]] <- plotCKMRabundancePretty(fits = fits_5by5[[i]],
#                                         year_lim = c(-20,0),
#                                         max_y_axis = 3000,
#                                         female = FALSE,
#                                         truth = combined_data[[1]]$N_hist,
#                                         y0 = 2014)
#   
# }
# gridExtra::grid.arrange(grobs = plots, ncol = 5)
