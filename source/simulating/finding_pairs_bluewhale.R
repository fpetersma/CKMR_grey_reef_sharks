################################################################################
##  Extract POP from a indiv object. This script is mainly to run stuff       ##
##  on the many cores of the bluewhale server.                                ##
##
##  Note: when rerunning this, it might suffice to load some preprocessed
##  Rdata files. This will safe a lot of time. I repeat, a *lot*.
################################################################################

## Load packages and source custom functions ===================================
library(fishSim)
library(ids)
library(pbapply)
library(doParallel)
library(CKMRcpp)

# source("source/fitting/CKMR_functions.R")
# source("source/simulating/custom_functions_fishSim.R")

data_folder <- "data/vanilla_variable_reproduction_sample_years_139-140_sample_size_375/" 

## Load the correct 100_sims_vanilla file
load(file = paste0(data_folder, "1000_sims_vanilla_repro=U(1,4)_no_gestation.RData"))

## Find the pairs in parallel. 
n_cores <- 30
cl <- makeCluster(n_cores)
clusterExport(cl, c("findRelativesCustom"))
pairs_list <- pblapply(simulated_data_sets, function(indiv) {
  indiv_samp <- indiv[!is.na(indiv$SampY), ]
  return(findRelativesCustom(indiv = indiv, 
                             sampled = TRUE))
}, cl = cl); stopCluster(cl);

## Save data
# save(file = paste0(data_folder, "1000_sims_vanilla_pairs.RData"), list = "pairs_list")

POPs_list <- pblapply(pairs_list, function(pairs) {
  return(pairs[pairs$OneTwo == 1, ]) ## Parent-Offspring pairs)
})

## Save data
# save(list = c("POPs_list"), file = paste0(data_folder, "100_sims_vanilla_POPs.RData"))
# load(paste0(data_folder, "100_sims_vanilla_POPs.RData"))

## Self-captures
selfie_list <- pblapply(simulated_data_sets, function(indiv) {
  self_captures <- indiv[indiv$no_samples > 1, ]
})

## Prep data ready for analysis
combined_data <- lapply(1:length(simulated_data_sets), function(i) {
  out <- list(POPs = POPs_list[[i]], 
              self = selfie_list[[i]],
              indiv = simulated_data_sets[[i]])
  return(out)
})

# save(list = c("combined_data"), file = paste0(data_folder, "1000_sims_vanilla_combined_data.RData"))

n_cores <- 24
cl <- makeCluster(n_cores)
clusterExport(cl, c("vbgf"))
dfs <- pblapply(combined_data, function(x) {
  POPs <- x$POPs
  self <- x$self
  indiv <- x$indiv
  
  ## Extract sampled individuals from the population
  sampled_indiv <- indiv[!is.na(indiv$SampY),] # a total of 1258 sampled 
  # individuals = 1581306 comparisons
  
  ## Add the recaptures [This is new can can be removed if it does not work]
  self <- self[rep(seq_len(nrow(self)), self$no_samples), ]
  for (id in unique(self$Me)) {
    ## Get all samplings of the same individual with id
    selfs <- self[self$Me == id, ]
    ## Extract all the sampling years of this individual
    samp_years <- as.numeric(unlist(strsplit(selfs[1, "SampY"], split = "_")))
    ## Assign one sampling year to every sampling occassion. 
    for (i in 1:nrow(selfs)) {
      selfs[i, "SampY"] <- samp_years[i]
    }
    
    ## Updated the rows in 'self'
    self[self$Me == id, ] <- selfs
  }
  sampled_indiv <- rbind(sampled_indiv[sampled_indiv$no_samples == 1, ], 
                         self)
  sampled_indiv$SampY <- as.numeric(sampled_indiv$SampY)
  
  
  ## Keep relevant information and rename columns to match theory
  sampled_indiv$SampAge <- sampled_indiv$SampY - sampled_indiv$BirthY
  
  sampled_indiv <- subset(sampled_indiv, select = c(Me, Sex, SampAge, SampY))
  colnames(sampled_indiv) <- c("id", "sex", "capture_age", "capture_year")
  
  ## In practice age is latent, but length is not. Convert age to length through 
  ## the VBGF, and add a normal error. Roughly based on other studies.
  ## Derive the length corresponding to the age, and add error from N(0, 2)
  sampled_indiv$length <- vbgf(sampled_indiv$capture_age)
  # sampled_indiv$length <- sampled_indiv$length + rnorm(nrow(sampled_indiv), 0, 2)
  
  ## Create matrix of all combinations with information 
  df_ids <- as.matrix(subset(expand.grid(1:nrow(sampled_indiv), 
                                         1:nrow(sampled_indiv)), 
                             Var1 != Var2))
  df <- as.data.frame(matrix(sampled_indiv$id[df_ids], nrow = nrow(df_ids), 
                             ncol = 2, byrow = FALSE))
  ## OLD expand grid below, no longer works with self captures ------------>>>>>
  # df_ids <- subset(expand.grid(sampled_indiv$id, sampled_indiv$id), ## all comparisons
  # Var1 != Var2)
  df_ids <- t(combn(1:nrow(sampled_indiv), 2)) ## only unique comparisons
  df <- as.data.frame(matrix(sampled_indiv$id[df_ids], nrow = nrow(df_ids), 
                             ncol = 2, byrow = FALSE))
  ## <<<<<< --------------------------------------------------------------------
  colnames(df) <- c("indiv_1_id", "indiv_2_id")
  
  ## Add columns with individual capture info, to be filled in in the loop below
  info_cols <- c("sex", "capture_year", "length", "capture_age")
  df[, paste0("indiv_1_" , info_cols)] <- NA
  df[, paste0("indiv_2_" , info_cols)] <- NA
  
  ## Loop through the unique ids and add capture information to main df
  for (id in unique(sampled_indiv$id)) {
    ## Extract info for individual with id
    info <- sampled_indiv[sampled_indiv$id == id, info_cols][1, ] # if multiple rows (in case of recaptures) only keep the first row
    
    ## Add info when individual id is first in the comparison
    df[df$indiv_1_id == id, paste0("indiv_1_", info_cols)] <- info
    
    ## Add info when individual id is second in the comparison
    df[df$indiv_2_id == id, paste0("indiv_2_", info_cols)] <- info
  }

  ## Round lengths to nearest integer, as this is how it is recorded.
  df$indiv_1_length <- round(df$indiv_1_length)
  df$indiv_2_length <- round(df$indiv_2_length)
  
  ## Make sure the right columns are numeric
  df$indiv_1_capture_year <- as.numeric(df$indiv_1_capture_year)
  df$indiv_1_length <- as.numeric(df$indiv_1_length)
  df$indiv_1_capture_age <- as.numeric(df$indiv_1_capture_age)
  df$indiv_2_capture_year <- as.numeric(df$indiv_2_capture_year)
  df$indiv_2_length <- as.numeric(df$indiv_2_length)
  df$indiv_2_capture_age <- as.numeric(df$indiv_2_capture_age)
  
  ## Add unique id for the combination
  df$comb_id <- paste0(df$indiv_1_id, "_", df$indiv_2_id)
  POPs$comb_id <- paste0(POPs$Var1, "_", POPs$Var2)
  POPs$comb_id_rev <- paste0(POPs$Var2, "_", POPs$Var1)
  
  ## Add whether a POP was found, or not. as direction is unknown, check for PO and OP
  df$pop_found <- df$comb_id %in% POPs$comb_id |
    df$comb_id %in% POPs$comb_id_rev
  df$kinship <- ifelse(df$comb_id %in% POPs$comb_id | df$comb_id %in% POPs$comb_id_rev,
                       yes = "PO/OP", no = "U")
  df$kinship[df$indiv_1_id == df$indiv_2_id] <- "S"
  
  # ## Find unique covariate combination using length or age (make sure this is correct!)
  # df$covariate_combo_id <- apply(df, 1, function(row) {
  #   return(paste0(row["indiv_1_sex"], row["indiv_1_capture_year"], 
  #                 row["indiv_1_capture_age"],
  #                 # row["indiv_1_length"],
  #                 row["indiv_2_sex"], row["indiv_2_capture_year"], 
  #                 row["indiv_2_capture_age"],
  #                 # row["indiv_2_length"],
  #                 row["kinship"]))
  # })
  # 
  # df$covariate_combo_id <- as.numeric(as.factor(df$covariate_combo_id))
  # 
  # ## Good news: only 138,425 unique probabilities to derive, instead of 1,573,770
  # 
  # ## Add frequency of every identical covariate combination
  # df <- transform(df, covariate_combo_freq = ave(seq(nrow(df)),
  #                                                covariate_combo_id,
  #                                                FUN = length))
  
  return(df)
  
}, cl = cl); stopCluster(cl);

## :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
## Code below can be run to remove covariate_combo_id and covariate_combo_freq
## :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
library(dplyr)
dfs <- pblapply(dfs, function(x) {
  # x <- x[, -c(ncol(x) - 1, ncol(x))]
  x <- dplyr::select(x, -c(covariate_combo_id, covariate_combo_freq))
  return(x)
})

## :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
## Run below to add noise the length with the preferred uncertainty
## :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
dfs <- pblapply(dfs, function(x) {
  x$indiv_1_length <- round(vbgf(x$indiv_1_capture_age) + rnorm(nrow(x), 0, 10))
  x$indiv_2_length <- round(vbgf(x$indiv_2_capture_age) + rnorm(nrow(x), 0, 10))
  return(x)
})

## :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
## Run below to add unique covariate combo ids and return sufficient dfs
## :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
n_cores <- 30
cl <- makeCluster(n_cores)
# clusterExport(cl, c("vbgf"))
dfs_suff <- pblapply(dfs, function(x) {
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
    ungroup()
  return(df_sufficient)
}, cl = cl); stopCluster(cl);


# save(list = c("dfs"), file = "data/test_data_dfs.RData")
# save(list = c("dfs_suff"), file = "data/test_data_dfs_suff.RData")

save(list = c("dfs"), file = paste0(data_folder, "1000_sims_dfs.RData"))
save(list = c("dfs_suff"), file = paste0(data_folder, "1000_sims_dfs_suff.RData"))

# ## Looking up relationship between captured pairs
# pairs <- findRelativesPar(indiv = indiv, 
#                           sampled = TRUE, 
#                           nCores = 20)
# POPs <- pairs[pairs$OneTwo == 1, ] ## Parent-Offspring pairs
# 
# ## Save data
# save.image("data/data_for_testing_estimator.RData")
# 
# ## Extracting simulated abundances
# N_true <- t(sapply(simulated_data_sets, function(x) {
#   mature_females <- sum(is.na(x$DeathY) & x$Sex == "F" & x$AgeLast >= 12)
#   mature_males <- sum(is.na(x$DeathY) & x$Sex == "M" & x$AgeLast >= 10)
#   return(c(N_m = mature_males, N_f = mature_females))
# }))
# 
# hist(N_true[, 1], main = "", xlab = "Number of mature males"); abline(v = c(mean(N_true[, 1]), median(N_true[, 1])), col = "red", lty = c(1, 2))
# 
# hist(N_true[, 2], main = "",  xlab = "Number of mature females"); abline(v = c(mean(N_true[, 2]), median(N_true[, 2])), col = "red", lty = c(1, 2))
