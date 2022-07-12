#' ========================================================================== ##
#' Test nllCKMRcpp()
#' ========================================================================== ##

library(Rcpp)

load("C:/Users/felix/OneDrive - University of St Andrews/Documents/University of St Andrews/PhD/Sharks/Close-kin mark-recapture/simulation_software/df_and_df_sufficient.RData")

sourceCpp("C:/Users/felix/OneDrive - University of St Andrews/Documents/University of St Andrews/PhD/Sharks/Close-kin mark-recapture/simulation_software/nllCKMRcpp.cpp")

par <- list(phi = boot::logit(0.855), # same as plogis(0.9) -- boot::inv.logit() is qlogis()
            N_t0_m = log(900), 
            N_t0_f = log(600)) #,
            # r = log(0.002), 
            # sigma_vbgf = log(2))

df_select <- df[, ]
df_select <- df_sufficient[, ]

dat <- list(alpha_m = 10, 
            alpha_f = 12,
            
            r = log(0.002),
            sigma_vbgf = log(2),
            
            max_age = 19,
            t0 = 80,
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

