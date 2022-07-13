#' ========================================================================== ##
#' Test nllCKMRcppAgeKnown()
#' ========================================================================== ##

library(Rcpp)

load("data/df_sufficient_age_known.RData")

sourceCpp("source/fitting/nllCKMRcppAgeKnown.cpp")

par <- list(
  # phi = boot::logit(0.5), # same as plogis(0.9) -- boot::inv.logit() is qlogis()
  N_t0_m = log(500), 
  N_t0_f = log(500))
# r = log(1.03)) 
# sigma_vbgf = log(2))

df_select <- df[, ]
df_select <- df_sufficient[, ]

dat <- list(alpha_m = 10, 
            alpha_f = 12,
            
            r = log(1.02),
            # sigma_vbgf = log(0.000001),
            # phi = boot::logit(0.87),
            
            max_age = 19,
            t0 = 40,
            # vbgf_l_inf = 175,
            # vbgf_k = 0.1,
            # vbgf_t0 = -3.5,
            s1 = df_select$indiv_1_sex,
            s2 = df_select$indiv_2_sex,
            c1 = df_select$indiv_1_capture_year,
            c2 = df_select$indiv_2_capture_year,
            a1 = df_select$indiv_1_age,
            a2 = df_select$indiv_2_age,
            pair_found = df_select$pop_found,
            cov_combo_freq = df_select$covariate_combo_freq,
            n = nrow(df_select))

sourceCpp("source/fitting/nllCKMRcppAgeKnown.cpp")

system.time({
  result_Rcpp <- nllPOPCKMRcppAgeKnown(dat = dat, par = par)
})

res <- nlminb(start = par, 
              objective = nllPOPCKMRcppAgeKnown, 
              dat = dat, 
              control = list(trace = 1))

sum(indiv$Sex == "F" & is.na(indiv$DeathY) & indiv$AgeLast >= 12)
sum(indiv$Sex == "M" & is.na(indiv$DeathY) & indiv$AgeLast >= 10)

# res <- optim(par = par, 
#              fn = nllPOPCKMRcpp, 
#              method = "BFGS",
#              dat = dat, 
#              control = list(trace = 0))

