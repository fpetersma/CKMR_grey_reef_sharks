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

## BE CAREFULL WHICH DATA TO LOAD
## ==============================

load("data/vanilla_sample_years_139-140_sample_size_375/1000_sims_dfs_length_unique_combos.RData")

## MAKE SURE THE LATEST VERSION OF THE LIKELIHOOD IS COMPILED
## ==========================================================
sourceCpp("source/fitting/nllCKMRcppVanilla.cpp")

result_list <- pblapply(dfs[1:10], function(df) {
  par <- list(
    # phi = boot::logit(0.87), # same as plogis(0.9) -- boot::inv.logit() is qlogis()
    N_t0_m = log(500), 
    r = log(1.00),
    N_t0_f = log(500))
  # sigma_vbgf = log(2))
  
  df_select <- df[, ]
  
  dat <- list(alpha_m = 10,
              alpha_f = 10,
              
              # r = log(1.00),
              # sigma_vbgf = log(0.000001),
              phi = boot::logit(1 - 0.153),
              
              max_age = 19,
              t0 = 140,
              # vbgf_l_inf = 175,
              # vbgf_k = 0.1,
              # vbgf_t0 = -3.5,
              s1 = df_select$indiv_1_sex,
              s2 = df_select$indiv_2_sex,
              c1 = df_select$indiv_1_capture_year,
              c2 = df_select$indiv_2_capture_year,
              a1 = df_select$indiv_1_capture_age,
              a2 = df_select$indiv_2_capture_age,
              pair_found = df_select$pop_found,
              cov_combo_freq = df_select$covariate_combo_freq,
              n = nrow(df_select))
  
  res <- nlminb(start = par, 
                objective = nllPOPCKMRcppAgeKnown, 
                dat = dat, 
                control = list(trace = 0))
  
  return(res)
  
})

## Look at the estimates
N_est_suff<- t(sapply(result_list, function(res){
  return(c(N_male = exp(res$par["N_t0_m"]), N_female = exp(res$par["N_t0_f"])))
}))

summary(N_est)

par(mfrow = c(3, 1))

hist(N_est[, 1], main = "", xlab = "Number of mature males"); abline(v = c(mean(N_est[, 1]), median(N_est[, 1])), col = "red", lty = c(1, 2))
hist(N_est[, 2], main = "",  xlab = "Number of mature females"); abline(v = c(mean(N_est[, 2]), median(N_est[, 2])), col = "red", lty = c(1, 2))

r_est <- sapply(result_list, function(res){
  exp(res$par["r"])
})

hist(r_est, main = "", xlab = "Yearly growth rate"); abline(v = c(mean(r_est), median(r_est)), col = "red", lty = c(1, 2))
