################################################################################
##
##                      Test nllCKMRcppVanilla()
##
################################################################################

##================================================================
##                      One thousand models                     =
##================================================================

library(Rcpp)
library(parallel)
library(pbapply)
# library(CKMRcpp)
## BE CAREFULL WHICH DATA TO LOAD
## ==============================

file_folder <- "data/vanilla_variable_reproduction_sample_years_136-140_sample_size_150/"

load(paste0(file_folder, "smaller_files/1000_sims_dfs_suff_length_sd=2_unique_combos.RData"))

## MAKE SURE THE LATEST VERSION OF THE LIKELIHOOD IS COMPILED
## ==========================================================
# sourceCpp("source/fitting/nllCKMRcppVanilla_unknown_age.cpp")
# Rcpp::sourceCpp("source/fitting/nllCKMRcppVanilla_gestation.cpp")

# ## Convert ages to length without noise
# dfs_suff <- pblapply(dfs_suff, function(df) {
#   df$indiv_1_length <- round(vbgf(df$indiv_1_capture_age))
#   df$indiv_2_length <- round(vbgf(df$indiv_2_capture_age))
#   return(df)
# })

## The unknown age version is a lot slower, so run in parallel
n_cores <- 20 # 50 fits per core
cl <- makeCluster(n_cores)
result_list <- pblapply(dfs_suff[1:1000], function(df) {
  ## Source cpp file if not using the CKMRcpp package
  # Rcpp::sourceCpp("source/fitting/nllCKMRcppVanilla_gestation.cpp")
  
  ## Create parameter object for nlminb()
  par <- list(
    # phi = boot::logit(0.87), # same as plogis(0.9) -- boot::inv.logit() is qlogis()
    N_t0_m = log(500), 
    # r = log(1.2),
    # sigma_l = log(0.01),
    # phi = boot::logit(1 - 0.153),
    N_t0_f = log(500))
  
  ## Take a subset of the data if required
  df_select <- df[, ]
  
  ## Create the data object for nlminb()
  dat <- list(alpha_m = 10,
              alpha_f = 10,
              
              ## MAKE SURE sigma_l AND phi ARE CORRECT!
              r = log(1.0002),
              sigma_l = log(2),
              phi = boot::logit(1 - 0.1675),
              
              max_age = 19,
              max_length = 250,
              t0 = 140,
              vbgf_l_inf = 175,
              vbgf_k = 0.1,
              vbgf_a0 = -3.5,
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
                  objective = CKMRcpp::nllPOPCKMRcppAgeUnknown, 
                  dat = dat, 
                  control = list(trace = 0, rel.tol = 1e-10))
  # })
  
  return(res)
}, cl = cl); stopCluster(cl); # end of pblapply() in parallel


# save(list = c("result_list"), file = paste0(file_folder, "result_list_1-1000_unknown_age_sd=10.RData"))

## Look at the estimates
N_est <- t(sapply(result_list, function(res){
  return(c(N_male = exp(res$par["N_t0_m"]), N_female = exp(res$par["N_t0_f"])))
}))

r_est <- sapply(result_list, function(res){
  exp(res$par["r"])
})

# phi_est <- sapply(result_list, function(res){
#   boot::inv.logit(res$par["phi"])
# })

# sigma_l_est <- sapply(result_list, function(res){
#   exp(res$par["sigma_l"])
# })

summary(cbind(N_est, r_est))

par(mfrow = c(3, 1))

hist(N_est[, 1], main = "", xlab = "Number of mature males"); abline(v = c(mean(N_est[, 1]), median(N_est[, 1])), col = "red", lty = c(1, 2))
hist(N_est[, 2], main = "",  xlab = "Number of mature females"); abline(v = c(mean(N_est[, 2]), median(N_est[, 2])), col = "red", lty = c(1, 2))
hist(r_est, main = "", xlab = "Yearly growth rate"); abline(v = c(mean(r_est), median(r_est)), col = "red", lty = c(1, 2))
# hist(phi_est, main = "", xlab = "Yearly surival rate"); abline(v = c(mean(phi_est), median(phi_est)), col = "red", lty = c(1, 2))

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

