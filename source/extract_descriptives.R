## First load the correct dfs and simulated_data_sets files. 

load("data/vanilla_sample_years_136-140_sample_size_150/100_sims_mix.RData")
load("data/vanilla_sample_years_136-140_sample_size_150/100_sims_dfs_unique_combos.RData")

## Then run the rest:

descriptives <- matrix(NA, nrow = 100, ncol = 6)
colnames(descriptives) <- c("ID", "n_alive", "n_moms", "n_dads", "n_sampled", "n_pop")

for (i in 1:100) {
  sim_dat <- simulated_data_sets[[i]]
  n_alive <- sum(is.na(sim_dat$DeathY))
  n_moms <- sum(sim_dat$Sex == "F" & is.na(sim_dat$DeathY) & sim_dat$AgeLast >= 10)
  n_dads <- sum(sim_dat$Sex == "M" & is.na(sim_dat$DeathY) & sim_dat$AgeLast >= 10)
  
  df <- dfs[[i]]
  n_sampled <- length(unique(df$indiv_1_id))
  n_pop <- sum(df$pop_found)
  
  descriptives[i, ] <- c(i, n_alive, n_moms, n_dads, n_sampled, n_pop)
}


