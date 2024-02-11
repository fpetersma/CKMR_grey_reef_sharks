#' ========================================================================== ##
#' Test nllCKMRcpp()
#' ========================================================================== ##

library(Rcpp)

load("data/df_sufficient_age_known.RData")

sourceCpp("source/fitting/nllCKMRcpp.cpp")

par <- list(
  #phi = boot::logit(0.87), # same as plogis(0.9) -- boot::inv.logit() is qlogis()
            N_t0_m = log(600), 
            N_t0_f = log(900))
            # r = log(1.03)) 
            # sigma_vbgf = log(2))

df_select <- df[, ]
df_select <- df_sufficient[, ]

dat <- list(alpha_m = 10, 
            alpha_f = 12,
            
            r = log(1.00),
            sigma_vbgf = log(0.00001),
            phi = boot::logit(0.87),
            
            max_age = 19,
            t0 = 140,
            vbgf_l_inf = 175,
            vbgf_k = 0.1,
            vbgf_t0 = -3.5,
            s1 = df_select$indiv_1_sex,
            s2 = df_select$indiv_2_sex,
            c1 = df_select$indiv_1_capture_year,
            c2 = df_select$indiv_2_capture_year,
            l1 = df_select$indiv_1_length,
            l2 = df_select$indiv_2_length,
            pair_found = df_select$pop_found,
            cov_combo_freq = df_select$covariate_combo_freq,
            n = nrow(df_select))

sourceCpp("source/fitting/nllCKMRcpp.cpp")

system.time({
  result_Rcpp <- nllPOPCKMRcpp(dat = dat, par = par)
})

res <- nlminb(start = par, 
              objective = nllPOPCKMRcpp, 
              dat = dat, 
              control = list(trace = 1))



# res <- optim(par = par, 
#              fn = nllPOPCKMRcpp, 
#              method = "BFGS",
#              dat = dat, 
#              control = list(trace = 0))

