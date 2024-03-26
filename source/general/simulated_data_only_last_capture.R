################################################################################
## Reduce simulated data sets and save
## 12/03/2024
################################################################################

data_folder <- "data/simulation_study/simple/" 

## Load the correct simple or complex data file
load(file = paste0(data_folder, "1000_simulated_data_sets.RData"))

## Below turns 12.8 GB to 61.2 MB
simple_sims_reduced <- lapply(simulated_data_sets, function(x) {
  return(x[x$no_samples > 0, ]) # only keep individuals that were sampled
})
rm(simulated_data_sets)

data_folder <- "data/simulation_study/complex/"

## Load the correct simple or complex data file
load(file = paste0(data_folder, "1000_simulated_data_sets.RData"))

## Below turns 9.8 GB to 67.2 MB
complex_sims_reduced <- lapply(simulated_data_sets, function(x) {
  return(x[x$no_samples > 0, ]) # only keep individuals that were sampled
})
rm(simulated_data_sets)

save(list = c("simple_sims_reduced", "complex_sims_reduced"),
     file = "data/simulation_study/simulated_data_reduced_both_species.RData")

## Phase 2
## ==============

## Only keep last sampling year
simple_sims_only_last_capture <- lapply(simple_sims_reduced, function(x) {
  freq <- sum(x$no_samples == 2)
  x[x$no_samples == 2, 10:12] <- matrix(c(2014, NA, 1), nrow = freq, 
                                        ncol = 3, byrow = TRUE)
  return(x)
})

complex_sims_only_last_capture <- lapply(complex_sims_reduced, function(x) {
  freq <- sum(x$no_samples == 2)
  x[x$no_samples == 2, 11:13] <- matrix(c(2014, NA, 1), nrow = freq, 
                                        ncol = 3, byrow = TRUE)
  return(x)
})

save(list = c("simple_sims_only_last_capture", "complex_sims_only_last_capture"),
     file = "data/simulation_study/simulated_data_only_last_capture_both_species.RData")

## Find pairs
## ===============

## If the lines before have been run before, just load the RData file
# load("data/simulation_study/simulated_data_only_last_capture_both_species.RData")
library(pbapply)
library(parallel)

## Find the pairs in parallel. 
n_cores <- 22
cl <- makeCluster(n_cores)
simple_pairs <- pblapply(simple_sims_only_last_capture, function(indiv) {
  return(CKMRcpp::findRelativesCustom(indiv = indiv, verbose = FALSE,
                                      sampled = TRUE))
}, cl = cl); stopCluster(cl)

simple_POPs <- pblapply(simple_pairs, function(pairs) {
  return(pairs[pairs$OneTwo == 1, ]) ## Parent-Offspring pairs)
})
rm(simple_pairs)

# Find the complex pairs in parallel. 
n_cores <- 22
cl <- makeCluster(n_cores)
complex_pairs <- pblapply(complex_sims_only_last_capture, function(indiv) {
  return(CKMRcpp::findRelativesCustom(indiv = indiv, verbose = FALSE,
                                      sampled = TRUE))
}, cl = cl); stopCluster(cl)

complex_POPs <- pblapply(complex_pairs, function(pairs) {
  return(pairs[pairs$OneTwo == 1, ]) ## Parent-Offspring pairs)
})

rm(complex_pairs)

save(list = c("simple_POPs", "complex_POPs"),
     file = "data/simulation_study/POPs_only_last_capture.RData")

## Add population histories
## If population size is the same, just do this once and repeat for every sampling scheme

## Load the correct simple or complex data file
load(file = paste0("data/simulation_study/simple/1000_simulated_data_sets.RData"))
simple_sims_complete <- simulated_data_sets

n_cores <- 22
cl <- makeCluster(n_cores)
first_breed_male <- 10
first_breed_female <- 10
clusterExport(cl, varlist = c("first_breed_male", "first_breed_female"))

simple_N_hist <- pblapply(simulated_data_sets, function(indiv) {
  N_hist <- t(sapply(min(indiv$DeathY, na.rm = T):max(indiv$DeathY, na.rm = T), 
                     function(year) {
                       N_male <- sum(CKMRcpp::extractTheLiving(indiv, year, TRUE, 
                                                               first_breed_male, "male", TRUE))
                       N_female <- sum(CKMRcpp::extractTheLiving(indiv, year, TRUE, 
                                                                 first_breed_female, "female", TRUE))
                       return(c(N_m = N_male, N_f = N_female))
                     }))
  row.names(N_hist) <- min(indiv$DeathY, na.rm = T):max(indiv$DeathY, na.rm = T)
  return(N_hist)
}, cl = cl); stopCluster(cl)

## Load the correct simple or complex data file
load(file = paste0("data/simulation_study/complex/1000_simulated_data_sets.RData"))
complex_sims_complete <- simulated_data_sets

n_cores <- 22
cl <- makeCluster(n_cores)

first_breed_male <- 17          # applies only to males
first_breed_female <- 19        # applies only to females
clusterExport(cl, varlist = c("first_breed_male", "first_breed_female"))

complex_N_hist <- pblapply(simulated_data_sets, function(indiv) {
  N_hist <- t(sapply(min(indiv$DeathY, na.rm = T):max(indiv$DeathY, na.rm = T), 
                     function(year) {
                       N_male <- sum(CKMRcpp::extractTheLiving(indiv, year, TRUE, 
                                                               first_breed_male, "male", TRUE))
                       N_female <- sum(CKMRcpp::extractTheLiving(indiv, year, TRUE, 
                                                                 first_breed_female, "female", TRUE))
                       return(c(N_m = N_male, N_f = N_female))
                     }))
  row.names(N_hist) <- min(indiv$DeathY, na.rm = T):max(indiv$DeathY, na.rm = T)
  return(N_hist)
}, cl = cl); stopCluster(cl)

rm(simulated_data_sets)

## Create combined data objects in reduced format
simple_combined <- pblapply(1:length(simple_sims_only_last_capture), function(i) {
  out <- list(POPs = simple_POPs[[i]],
              N_hist = simple_N_hist[[i]],
              indiv = simple_sims_only_last_capture[[i]])
  return(out)
})
complex_combined <- pblapply(1:length(complex_sims_only_last_capture), function(i) {
  out <- list(POPs = complex_POPs[[i]],
              N_hist = complex_N_hist[[i]],
              indiv = complex_sims_only_last_capture[[i]])
  return(out)
})

save(list = c("simple_combined", "complex_combined"),
     file = "data/simulation_study/combined_data_only_last_capture.RData")

## Prepare final data frames for analysis
## Simple species first
## ==================================
n_cores <- 16
cl <- makeCluster(n_cores)
simple_dfs <- pblapply(simple_combined, function(x) {
  POPs <- x$POPs
  indiv <- x$indiv
  
  sampled_indiv <- indiv
  sampled_indiv$SampY[sampled_indiv$SampY == "2013_2014"] <- 2014
  sampled_indiv$SampY <- as.integer(sampled_indiv$SampY)
  
  ## Keep relevant information and rename columns to match theory
  sampled_indiv$SampAge <- as.integer(sampled_indiv$SampY - sampled_indiv$BirthY)
  sampled_indiv <- subset(sampled_indiv, select = c(Me, Sex, SampAge, SampY))
  
  sampled_indiv$Sex <- as.integer(sampled_indiv$Sex)
  
  colnames(sampled_indiv) <- c("id", "sex", "capture_age", "capture_year")
  
  ## In practice age is latent, but length is not. Convert age to length through 
  ## the VBGF, and add a normal error. Roughly based on other studies.
  sampled_indiv$length <- as.integer(CKMRcpp::vbgf(a = sampled_indiv$capture_age,
                                                   a_0 = -8.27, 
                                                   k = 0.0554, 
                                                   l_inf = 163))
  
  df_ids <- t(combn(1:nrow(sampled_indiv), 2)) ## only unique comparisons
  df <- as.data.frame(matrix(sampled_indiv$id[df_ids], nrow = nrow(df_ids), 
                             ncol = 2, byrow = FALSE))
  colnames(df) <- c("indiv_1_id", "indiv_2_id")
  
  ## Add columns with individual capture info, to be filled in in the loop below
  info_cols <- c("sex", "capture_year", "length", "capture_age")
  df[, paste0("indiv_1_" , info_cols)] <- NA
  df[, paste0("indiv_2_" , info_cols)] <- NA
  
  ## Loop through the unique ids and add capture information to main df
  ## [This is the slow part]
  for (id in unique(sampled_indiv$id)) {
    ## Extract info for individual with id
    info <- sampled_indiv[sampled_indiv$id == id, info_cols][1, ] # if multiple rows (in case of recaptures) only keep the first row
    
    ## Add info when individual id is first in the comparison
    df[df$indiv_1_id == id, paste0("indiv_1_", info_cols)] <- info
    
    ## Add info when individual id is second in the comparison
    df[df$indiv_2_id == id, paste0("indiv_2_", info_cols)] <- info
  }
  
  ## Round lengths to nearest integer, as this is how it is recorded.
  df$indiv_1_length <- as.integer(round(df$indiv_1_length))
  df$indiv_2_length <- as.integer(round(df$indiv_2_length))
  
  ## Make sure the right columns are numeric
  df$indiv_1_capture_year <- as.integer(df$indiv_1_capture_year)
  df$indiv_1_length <- as.integer(df$indiv_1_length)
  df$indiv_1_capture_age <- as.integer(df$indiv_1_capture_age)
  df$indiv_2_capture_year <- as.integer(df$indiv_2_capture_year)
  df$indiv_2_length <- as.integer(df$indiv_2_length)
  df$indiv_2_capture_age <- as.integer(df$indiv_2_capture_age)
  
  ## Add unique id for the combination
  df$comb_id <- as.numeric(paste0(df$indiv_1_id, df$indiv_2_id))
  POPs$comb_id <- as.numeric(paste0(POPs$Var1, POPs$Var2))
  POPs$comb_id_rev <- as.numeric(paste0(POPs$Var2, POPs$Var1))
  
  ## Add whether a POP was found, or not. as direction is unknown, check for PO and OP
  df$pop_found <- as.integer(df$comb_id %in% POPs$comb_id |
                               df$comb_id %in% POPs$comb_id_rev)
  df$kinship <- ifelse(df$comb_id %in% POPs$comb_id | df$comb_id %in% POPs$comb_id_rev,
                       yes = 1L, no = 0L)   # 0=U, 1=PO/OP, S=2
  df$kinship[df$indiv_1_id == df$indiv_2_id] <- 2L
  
  return(df)
  
}, cl = cl); stopCluster(cl)

## :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
## Run below to add noise the length with the preferred uncertainty
## :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
set.seed(2023) # seed for reproducibility

simple_dfs <- pblapply(simple_dfs, function(x) {
  x$indiv_1_length <- 
    as.integer(
      round(CKMRcpp::vbgf(a = x$indiv_1_capture_age,
                          a_0 = -8.27, 
                          k = 0.0554, 
                          l_inf = 163) + rnorm(nrow(x), 0, 2.89)))
  x$indiv_2_length <- 
    as.integer(
      round(CKMRcpp::vbgf(a = x$indiv_2_capture_age,
                          a_0 = -8.27, 
                          k = 0.0554, 
                          l_inf = 163) + rnorm(nrow(x), 0, 2.89)))
  return(x)
})

## Save dfs
save(list = "simple_dfs", 
     file = "data/simulation_study/simple_dfs_only_last_capture.RData")

## Complex species now
## ===============================
n_cores <- 16
cl <- makeCluster(n_cores)
complex_dfs <- pblapply(complex_combined, function(x) {
  POPs <- x$POPs
  indiv <- x$indiv
  
  sampled_indiv <- indiv
  sampled_indiv$SampY[sampled_indiv$SampY == "2013_2014"] <- 2013
  sampled_indiv$SampY <- as.integer(sampled_indiv$SampY)
  
  ## Keep relevant information and rename columns to match theory
  sampled_indiv$SampAge <- as.integer(sampled_indiv$SampY - sampled_indiv$BirthY)
  sampled_indiv <- subset(sampled_indiv, select = c(Me, Sex, SampAge, SampY))
  
  sampled_indiv$Sex <- as.integer(sampled_indiv$Sex)
  
  colnames(sampled_indiv) <- c("id", "sex", "capture_age", "capture_year")
  
  ## In practice age is latent, but length is not. Convert age to length through 
  ## the VBGF, and add a normal error. Roughly based on other studies.
  ## Derive the length corresponding to the age, and add error from N(0, 2)
  sampled_indiv$length <- as.integer(CKMRcpp::vbgf(a = sampled_indiv$capture_age,
                                                   a_0 = -8.27, 
                                                   k = 0.0554, 
                                                   l_inf = 163))
  
  df_ids <- t(combn(1:nrow(sampled_indiv), 2)) ## only unique comparisons
  df <- as.data.frame(matrix(sampled_indiv$id[df_ids], nrow = nrow(df_ids), 
                             ncol = 2, byrow = FALSE))
  colnames(df) <- c("indiv_1_id", "indiv_2_id")
  
  ## Add columns with individual capture info, to be filled in in the loop below
  info_cols <- c("sex", "capture_year", "length", "capture_age")
  df[, paste0("indiv_1_" , info_cols)] <- NA
  df[, paste0("indiv_2_" , info_cols)] <- NA
  
  ## Loop through the unique ids and add capture information to main df
  ## [This is the slow part]
  for (id in unique(sampled_indiv$id)) {
    ## Extract info for individual with id
    info <- sampled_indiv[sampled_indiv$id == id, info_cols][1, ] # if multiple rows (in case of recaptures) only keep the first row
    
    ## Add info when individual id is first in the comparison
    df[df$indiv_1_id == id, paste0("indiv_1_", info_cols)] <- info
    
    ## Add info when individual id is second in the comparison
    df[df$indiv_2_id == id, paste0("indiv_2_", info_cols)] <- info
  }
  
  ## Round lengths to nearest integer, as this is how it is recorded.
  df$indiv_1_length <- as.integer(round(df$indiv_1_length))
  df$indiv_2_length <- as.integer(round(df$indiv_2_length))
  
  ## Make sure the right columns are numeric
  df$indiv_1_capture_year <- as.integer(df$indiv_1_capture_year)
  df$indiv_1_length <- as.integer(df$indiv_1_length)
  df$indiv_1_capture_age <- as.integer(df$indiv_1_capture_age)
  df$indiv_2_capture_year <- as.integer(df$indiv_2_capture_year)
  df$indiv_2_length <- as.integer(df$indiv_2_length)
  df$indiv_2_capture_age <- as.integer(df$indiv_2_capture_age)
  
  ## Add unique id for the combination
  df$comb_id <- as.numeric(paste0(df$indiv_1_id, df$indiv_2_id))
  POPs$comb_id <- as.numeric(paste0(POPs$Var1, POPs$Var2))
  POPs$comb_id_rev <- as.numeric(paste0(POPs$Var2, POPs$Var1))
  
  ## Add whether a POP was found, or not. as direction is unknown, check for PO and OP
  df$pop_found <- as.integer(df$comb_id %in% POPs$comb_id |
                               df$comb_id %in% POPs$comb_id_rev)
  df$kinship <- ifelse(df$comb_id %in% POPs$comb_id | df$comb_id %in% POPs$comb_id_rev,
                       yes = 1L, no = 0L)   # 0=U, 1=PO/OP, S=2
  df$kinship[df$indiv_1_id == df$indiv_2_id] <- 2L
  
  return(df)
  
}, cl = cl); stopCluster(cl)

## :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
## Run below to add noise the length with the preferred uncertainty
## :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
set.seed(2024) # seed for reproducibility

complex_dfs <- pblapply(complex_dfs, function(x) {
  x$indiv_1_length <- 
    as.integer(
      round(CKMRcpp::vbgf(a = x$indiv_1_capture_age,
                          a_0 = -8.27, 
                          k = 0.0554, 
                          l_inf = 163) + rnorm(nrow(x), 0, 2.89)))
  x$indiv_2_length <- 
    as.integer(
      round(CKMRcpp::vbgf(a = x$indiv_2_capture_age,
                          a_0 = -8.27, 
                          k = 0.0554, 
                          l_inf = 163) + rnorm(nrow(x), 0, 2.89)))
  return(x)
})

## Save dfs
save(list = "complex_dfs", 
     file = "data/simulation_study/complex_dfs_only_last_capture.RData")

## :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
## Run below to add unique covariate combo ids and return sufficient dfs
## :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

## Simple
n_cores <- 32
cl <- makeCluster(n_cores)
# clusterExport(cl, c("vbgf"))
simple_suff <- pblapply(simple_dfs, function(x) {
  library(dplyr)
  ## Find unique covariate combination using length or age (make sure this is correct!)
  x$covariate_combo_id <- apply(x, 1, function(row) {
    return(paste0(row["indiv_1_sex"], row["indiv_1_capture_year"], 
                  # row["indiv_1_capture_age"],
                  row["indiv_1_length"],
                  row["indiv_2_sex"], row["indiv_2_capture_year"], 
                  # row["indiv_2_capture_age"],
                  row["indiv_2_length"],
                  row["kinship"]))
  })
  
  x$covariate_combo_id <- as.numeric(as.factor(x$covariate_combo_id))
  
  ## Add frequency of every identical covariate combination
  x <- transform(x, covariate_combo_freq = ave(seq(nrow(x)),
                                               covariate_combo_id,
                                               FUN = length))
  
  df_sufficient <- x %>% 
    group_by(covariate_combo_id) %>% 
    slice(1) %>% 
    ungroup() %>% 
    select(-c(comb_id, indiv_1_id, indiv_2_id, covariate_combo_id, pop_found))
  return(as.data.frame(df_sufficient))
}, cl = cl); stopCluster(cl)

## Complex
n_cores <- 32
cl <- makeCluster(n_cores)
# clusterExport(cl, c("vbgf"))
complex_suff <- pblapply(complex_dfs, function(x) {
  library(dplyr)
  ## Find unique covariate combination using length or age (make sure this is correct!)
  x$covariate_combo_id <- apply(x, 1, function(row) {
    return(paste0(row["indiv_1_sex"], row["indiv_1_capture_year"], 
                  # row["indiv_1_capture_age"],
                  row["indiv_1_length"],
                  row["indiv_2_sex"], row["indiv_2_capture_year"], 
                  # row["indiv_2_capture_age"],
                  row["indiv_2_length"],
                  row["kinship"]))
  })
  
  x$covariate_combo_id <- as.numeric(as.factor(x$covariate_combo_id))
  
  ## Add frequency of every identical covariate combination
  x <- transform(x, covariate_combo_freq = ave(seq(nrow(x)),
                                               covariate_combo_id,
                                               FUN = length))
  
  df_sufficient <- x %>% 
    group_by(covariate_combo_id) %>% 
    slice(1) %>% 
    ungroup() %>% 
    select(-c(comb_id, indiv_1_id, indiv_2_id, covariate_combo_id, pop_found))
  return(as.data.frame(df_sufficient))
}, cl = cl); stopCluster(cl)

save(list = c("simple_suff", "complex_suff"),
     file = "data/simulation_study/sufficient_dfs_only_last_capture.RData")
