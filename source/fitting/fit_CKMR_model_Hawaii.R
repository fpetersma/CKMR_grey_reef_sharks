## ========================================================================== ##
## This script fits a CKMR model to the POP data simulated in 
## 'GR_sims_Hawaii.R'. The complete environment was saved in 
## 'pop_Hawaii.RRdata'.  
## 
## The intention is to have most functions in separate scripts, but initially
## everything will be there. 
## ========================================================================== ##

## LOAD LIBRARIES, DATA AND PREP -----------------------------------------------
## =============================================================================
library(dplyr) # for clever df operations
library(pbapply) # add progress bar to apply, and allows for parallel 
library(parallel) # parellel processing

source("~/University of St Andrews/PhD/Sharks/Close-kin mark-recapture/simulation_software/CKMR_functions.R")

load("~/University of St Andrews/PhD/Sharks/Close-kin mark-recapture/simulation_software/pop_Hawaii.RData")

## View the POP data
# View(POPs)

## Extract sampled individuals from the population
sampled_indiv <- indiv[!is.na(indiv$SampY),] # a total of 1252 sampled 
                                             # individuals = 1566252 comparisons

## Keep relevant information and rename columns to match theory
sampled_indiv <- subset(sampled_indiv, select = c(Me, Sex, AgeLast, SampY))
colnames(sampled_indiv) <- c("id", "sex", "age", "capture_year")

## In practice age is latent, but length is not. Convert age to length through 
## the VBGF, and add a normal error. Roughly based on other studies.
## Derive the length corresponding to the age, and add error from N(0, 2)
sampled_indiv$length <- vbgf(sampled_indiv$age)
sampled_indiv$length <- sampled_indiv$length + rnorm(nrow(sampled_indiv), 0, 2)

## Create matrix of all combinations with information
df <- subset(expand.grid(sampled_indiv$id, sampled_indiv$id), 
            Var1 != Var2)
colnames(df) <- c("indiv_1_id", "indiv_2_id")

## Add columns with individual capture info, to be filled in in the loop below
info_cols <- c("sex", "capture_year", "length")
df[, paste0("indiv_1_" , info_cols)] <- NA
df[, paste0("indiv_2_" , info_cols)] <- NA

## Loop through the unique ids and add capture information to main df
for (id in unique(sampled_indiv$id)) {
  ## Extract info for individual with id
  info <- sampled_indiv[sampled_indiv$id == id, info_cols]
  
  ## Add info when individual id is first in the comparison
  df[df$indiv_1_id == id, paste0("indiv_1_", info_cols)] <- info

  ## Add info when individual id is second in the comparison
  df[df$indiv_2_id == id, paste0("indiv_2_", info_cols)] <- info
}

## Round lengths to nearest integer, as this is how it is recorded.
df$indiv_1_length <- round(df$indiv_1_length)
df$indiv_2_length <- round(df$indiv_2_length)

## Make sure the right columsn are numeric
df$indiv_1_capture_year <- as.numeric(df$indiv_1_capture_year)
df$indiv_1_length <- as.numeric(df$indiv_1_length)
df$indiv_2_capture_year <- as.numeric(df$indiv_2_capture_year)
df$indiv_2_length <- as.numeric(df$indiv_2_length)

## Add unique id for the combination
df$comb_id <- paste0(df$indiv_1_id, "_", df$indiv_2_id)
POPs$comb_id <- paste0(POPs$Var1, "_", POPs$Var2)

## Add whether a POP was found, or not
df$pop_found <- df$comb_id %in% POPs$comb_id

## Find unique cov
df$covariate_combo_id <- apply(df, 1, function(row) {
  return(paste0(row[3], row[4], row[5], row[6], row[7], row[8], row[10]))
})

df$covariate_combo_id <- as.numeric(as.factor(df$covariate_combo_id))

## Good news: only 126,782 unique probabilities to derive, instead of 1,566,252

## Add frequency of every identical covariate combination
df <- transform(df, covariate_combo_freq = ave(seq(nrow(df)),
                                               covariate_combo_id,
                                               FUN = length))

## Only keep the first of every identical covariate combo
df_sufficient <- df %>% 
  group_by(covariate_combo_id) %>% 
  slice(1) %>% 
  ungroup()

## Create a look-up table for length vs age probabilities
ages <- 0:20
lengths <- 50:200
age_length_prob_long <- expand.grid(ages, lengths)

## A long table format version
age_length_prob_long$prob <- dnorm(x = age_length_prob_long$Var2, 
                                   mean = vbgf(age_length_prob_long$Var1), 
                                   sd = 2)
## A matrix version with length on the rows and age on the columns
age_length_prob_matrix <- reshape2::acast(data = age_length_prob_long, 
                                          formula = Var2 ~ Var1, 
                                          value.var = "prob")

## DERIVE THE LOG-LIKELIHOOD ---------------------------------------------------
## =============================================================================

## Running on a single core: 01m 12s on 1 boosted core
result_single <- sum(pbapply(df_sufficient[26009, ], 1, function(obs) {
  pair_prob <- pairProb(pair = "PO", 
                        s1 = obs["indiv_1_sex"], 
                        c1 = as.numeric(obs["indiv_1_capture_year"]),
                        l1 = as.numeric(obs["indiv_1_length"]),
                        s2 = obs["indiv_2_sex"],
                        c2 = as.numeric(obs["indiv_2_capture_year"]),
                        l2 = as.numeric(obs["indiv_2_length"]),
                        alpha_m = 10, 
                        alpha_f = 12, 
                        max_age = 20,
                        phi = 0.9,
                        N_t0_m = 3000,
                        N_t0_f = 3000, 
                        t0 = 80, 
                        r = 0.02, 
                        sigma_vbgf = 2) # , 
                        # al_prob_matrix = age_length_prob_matrix)
  if (!as.logical(obs["pop_found"])) {
    pair_prob <- 1 - pair_prob
  }
  
  return(log(pair_prob) * as.numeric(obs["covariate_combo_freq"]))
}))

## Run in parallel
cl <- makeCluster(6)
clusterExport(cl, c("vbgf", "pairProb"))

## Cluster version: 14s on 6 boosted cores. 
result_cluster <- sum(pbapply(df_sufficient[26009, ], 1, function(obs) {
  pair_prob <- pairProb(pair = "PO", 
                        s1 = obs["indiv_1_sex"], 
                        c1 = as.numeric(obs["indiv_1_capture_year"]),
                        l1 = as.numeric(obs["indiv_1_length"]),
                        s2 = obs["indiv_2_sex"],
                        c2 = as.numeric(obs["indiv_2_capture_year"]),
                        l2 = as.numeric(obs["indiv_2_length"]),
                        alpha_m = 10, 
                        alpha_f = 12, 
                        max_age = 20,
                        phi = 0.9,
                        N_t0_m = 3000,
                        N_t0_f = 3000, 
                        t0 = 80, 
                        r = 0.02, 
                        sigma_vbgf = 2) # , 
                        # al_prob_matrix = age_length_prob_matrix)
  if (!as.logical(obs["pop_found"])) {
    pair_prob <- 1 - pair_prob
  }
  
  # return(log(pair_prob))
  return(log(pair_prob) * as.numeric(obs["covariate_combo_freq"]))
}, cl = cl))

stopCluster(cl)

#' CONCLUSION 08/07/2022
#' =====================
#' The clustered version using all data returns a llk of -3295.087.
#' The clustered version using sufficient data returns a llk of -3295.087.
#' So it seems to work? Rcpp version is way faster, but returns the wrong 
#' likelihood... Why is this?
#' 
#' Observation 26009 in df_sufficient yields a likelihood of:
#' R:     -30.6924334643983
#' Rcpp:   16.9353381021967
#' That is roughly half... 
#' 
#' Turns out there was a typo in line 171 of R function:
#' l.171 >> prob <- sum(output_offspring_level)
#' This should have been, and has been changed to:
#' l.171 >> prob <- sum(output_parent_level)
#' 
#' The Rcpp version was actually correct (or at least, more correct) than R, lol
#'