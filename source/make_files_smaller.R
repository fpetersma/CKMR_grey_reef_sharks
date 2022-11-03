## Load the required file
library(pbapply)


## Remove unique ids and turn numeric into integers
dfs_min <- pblapply(dfs_suff, function(x) {
  # x <- as.data.frame(dplyr::select(x, -c(indiv_1_id,
  #                                        indiv_2_id,
  #                                        comb_id,
  #                                        pop_found)))
  # x <- as.data.frame(dplyr::select(x, -c(covariate_combo_id)))
  # 
  # x$indiv_1_capture_year <- as.integer(x$indiv_1_capture_year)
  # x$indiv_1_capture_age <- as.integer(x$indiv_1_capture_age)
  # x$indiv_1_length <- as.integer(x$indiv_1_length)
  # 
  # x$indiv_2_capture_year <- as.integer(x$indiv_2_capture_year)
  # x$indiv_2_capture_age <- as.integer(x$indiv_2_capture_age)
  # x$indiv_2_length <- as.integer(x$indiv_2_length)
  
  x$kinship <- as.integer(gsub("U", 0,
                               gsub("PO/OP", 1,
                                    gsub("S", 2, x$kinship))))

  # x$indiv_1_sex <- as.integer(x$indiv_1_sex == "F")
  # x$indiv_2_sex <- as.integer(x$indiv_2_sex == "F")
  # 
  # x$covariate_combo_freq <- as.integer(x$covariate_combo_freq)
  
  return(x)
})

dfs_suff <- dfs_min

save(list = "dfs_suff", file = "data/vanilla_gestation_repro=U(1,4)_sample_years_136-140/smaller_files/1000_sims_dfs_suff_length_sd=2_unique_combos.RData")
