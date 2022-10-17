## Load the required file

## Remove unique ids and turn numeric into integers
dfs_min <- lapply(dfs_suff, function(x) {
  x <- as.data.frame(dplyr::select(x, -c(indiv_1_id,
                                         indiv_2_id,
                                         comb_id,
                                         covariate_combo_id,
                                         pop_found)))
  
  x$indiv_1_capture_year <- as.integer(x$indiv_1_capture_year)
  x$indiv_1_capture_age <- as.integer(x$indiv_1_capture_age)
  x$indiv_1_length <- as.integer(x$indiv_1_length)
  
  x$indiv_2_capture_year <- as.integer(x$indiv_2_capture_year)
  x$indiv_2_capture_age <- as.integer(x$indiv_2_capture_age)
  x$indiv_2_length <- as.integer(x$indiv_2_length)
  
  x$kinship <- as.integer(gsub("U", 0, 
                               gsub("PO/OP", 1, 
                                    gsub("S", 2, x$kinship))))
  
  x$indiv_1_sex <- as.integer(x$indiv_1_sex == "F")
  x$indiv_2_sex <- as.integer(x$indiv_2_sex == "F")
  
  x$covariate_combo_freq <- as.integer(x$covariate_combo_freq)
  
  # x$covariate_combo_id <- as.integer(x$covariate_combo_id)

  # x$pop_found <- as.integer(x$pop_found)
  
  return(x)
})

dfs_suff <- dfs_min

save(list = "dfs_suff", file = "data/vanilla_variable_reproduction_sample_years_139-140_sample_size_375/smaller_files/1000_sims_dfs_suff_length_sd=2_unique_combos.RData")
