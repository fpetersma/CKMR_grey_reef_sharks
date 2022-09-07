################################################################################
##
##                      Test nllCKMRcppParentAgeKnown()
##
################################################################################

##================================================================
##                          Single model                         =
##================================================================

library(Rcpp)

load("data/100_sims_dfs_suff_correct.RData")

sourceCpp("source/fitting/nllCKMRcppParentAgeKnown.cpp")

par <- list(
  # phi = boot::logit(0.87), # same as plogis(0.9) -- boot::inv.logit() is qlogis()
  N_t0_m = log(500), 
  r = log(1.01),
  N_t0_f = log(500))
# sigma_vbgf = log(2))

df_select <- dfs_suff[[1]]
POPs <- df_select[df_select$pop_found == TRUE, c("comb_id", "indiv_1_id", "indiv_2_id")]
POPs$comb_id_rev <- paste(POPs$indiv_2_id, POPs$indiv_1_id, sep = "_")
df_select$comb_id_rev <- paste(df_select$indiv_2_id, df_select$indiv_1_id, sep = "_")
df_select$pop_found <- df_select$comb_id %in% POPs$comb_id_rev |
  df_select$comb_id %in% POPs$comb_id 

dat <- list(alpha_m = 10, 
            alpha_f = 12,
            
            # r = log(1.01),
            sigma_vbgf = log(1),
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
            a1 = df_select$indiv_1_age,
            a2 = df_select$indiv_2_age,
            pair_found = df_select$pop_found,
            cov_combo_freq = df_select$covariate_combo_freq,
            n = nrow(df_select))

sourceCpp("source/fitting/nllCKMRcppParentAgeKnown.cpp")

system.time({
  result_Rcpp <- nllPOPCKMRcppParentAgeKnown(dat = dat, par = par)
})

res <- nlminb(start = par, 
              objective = nllPOPCKMRcppParentAgeKnown, 
              dat = dat, 
              control = list(trace = 1))

sum(indiv$Sex == "F" & is.na(indiv$DeathY) & indiv$AgeLast >= 12)
sum(indiv$Sex == "M" & is.na(indiv$DeathY) & indiv$AgeLast >= 10)

# res <- optim(par = par, 
#              fn = nllPOPCKMRcpp, 
#              method = "BFGS",
#              dat = dat, 
#              control = list(trace = 0))



##================================================================
##                      One hundred models                       =
##================================================================

library(Rcpp)
library(parallel)
library(pbapply)

load("data/100_sims_dfs_suff.RData")

sourceCpp("source/fitting/nllCKMRcppParentAgeKnown.cpp")

result_list <- pblapply(dfs_suff, function(df) {
  par <- list(
    # phi = boot::logit(0.87), # same as plogis(0.9) -- boot::inv.logit() is qlogis()
    N_t0_m = log(500), 
    r = log(1.03),
    N_t0_f = log(500))
  # sigma_vbgf = log(2))
  
  df_select <- df[, ]
  df_select <- df_sufficient[, ]
  
  dat <- list(alpha_m = 10, 
              alpha_f = 12,
              
              # r = log(1.00),
              sigma_vbgf = log(1),
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
              a1 = df_select$indiv_1_age,
              a2 = df_select$indiv_2_age,
              pair_found = df_select$pop_found,
              cov_combo_freq = df_select$covariate_combo_freq,
              n = nrow(df_select))
  
  res <- nlminb(start = par, 
                objective = nllPOPCKMRcppAgeKnown, 
                dat = dat, 
                control = list(trace = 1))
  
  return(res)
  
})

N_male <- sapply(result_list, function(res){
  res$par["N_t0_m"]
})

N_female <- sapply(result_list, function(res){
  res$par["N_t0_f"]
})


system.time({
  result_Rcpp <- nllPOPCKMRcppAgeKnown(dat = dat, par = par)
})



sum(indiv$Sex == "F" & is.na(indiv$DeathY) & indiv$AgeLast >= 12)
sum(indiv$Sex == "M" & is.na(indiv$DeathY) & indiv$AgeLast >= 10)

# res <- optim(par = par, 
#              fn = nllPOPCKMRcpp, 
#              method = "BFGS",
#              dat = dat, 
#              control = list(trace = 0))
