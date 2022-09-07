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

load("data/vanilla_sample_years_136-140_sample_size_150/1000_sims_dfs_suff_age_unique_combos.RData")

## MAKE SURE THE LATEST VERSION OF THE LIKELIHOOD IS COMPILED
## ==========================================================
sourceCpp("source/fitting/nllCKMRcppVanilla.cpp")

result_list <- pblapply(dfs_suff[1:1000], function(df) {
  par <- list(
    # phi = boot::logit(0.87), # same as plogis(0.9) -- boot::inv.logit() is qlogis()
    N_t0_m = log(500), 
    r = log(1.2),
    # phi = boot::logit(1 - 0.153),
    N_t0_f = log(500))
  
  # sigma_vbgf = log(2))
  
  df_select <- df[, ]
  
  dat <- list(alpha_m = 10,
              alpha_f = 10,
              
              # r = log(1.00),
              # sigma_vbgf = log(0.000001),
              phi = boot::logit(1 - 0.153),
              
              max_age = 19,
              t0 = 130,
              # vbgf_l_inf = 175,
              # vbgf_k = 0.1,
              # vbgf_t0 = -3.5,
              s1 = df_select$indiv_1_sex,
              s2 = df_select$indiv_2_sex,
              c1 = df_select$indiv_1_capture_year,
              c2 = df_select$indiv_2_capture_year,
              a1 = df_select$indiv_1_capture_age,
              a2 = df_select$indiv_2_capture_age,
              kinship = df_select$kinship,
              cov_combo_freq = df_select$covariate_combo_freq,
              n = nrow(df_select))
  
  res <- nlminb(start = par, 
                objective = nllPOPCKMRcppAgeKnown, 
                dat = dat, 
                control = list(trace = 0))
  
  return(res)
  
})

## Look at the estimates
N_est<- t(sapply(result_list, function(res){
  return(c(N_male = exp(res$par["N_t0_m"]), N_female = exp(res$par["N_t0_f"])))
}))

r_est <- sapply(result_list, function(res){
  exp(res$par["r"])
})

phi_est <- sapply(result_list, function(res){
  boot::inv.logit(res$par["phi"])
})

summary(cbind(N_est, r_est, phi_est))

par(mfrow = c(3, 1))

hist(N_est[, 1], main = "", xlab = "Number of mature males"); abline(v = c(mean(N_est[, 1]), median(N_est[, 1])), col = "red", lty = c(1, 2))
hist(N_est[, 2], main = "",  xlab = "Number of mature females"); abline(v = c(mean(N_est[, 2]), median(N_est[, 2])), col = "red", lty = c(1, 2))
# hist(r_est, main = "", xlab = "Yearly growth rate"); abline(v = c(mean(r_est), median(r_est)), col = "red", lty = c(1, 2))
hist(phi_est, main = "", xlab = "Yearly surival rate"); abline(v = c(mean(phi_est), median(phi_est)), col = "red", lty = c(1, 2))

MCE_data <- t(sapply(seq(from = 1, to = 901, by = 100), function(x) {
  mean_male <- mean(N_est[x:(x+99), 1])
  mean_female <- mean(N_est[x:(x+99), 2])
  mean_r <- mean(r_est[x:(x+99)])
  
  return(c(mean_male = mean_male, mean_female = mean_female, mean_r = mean_r))
}))

MCE_100 <- apply(MCE_data, 2, sd)
MCE_1000 <- MCE_100 / sqrt(10)
