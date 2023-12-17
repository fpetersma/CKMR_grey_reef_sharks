## Extract descriptive statistics from simulation 

## Vanilla
load("D:/felix/CKMR_grey_reef_sharks/data/1_population_multiple_sampling_schemes/vanilla/100_schemes_combined_data_with_N_hist_sim=144.RData")

n_POPs <- sapply(combined_data, function(x) nrow(x$POPs))
n_selfs <- sapply(combined_data, function(x) nrow(x$self))

## Complex
load("D:/felix/CKMR_grey_reef_sharks/data/1_population_multiple_sampling_schemes/complex/100_schemes_combined_data_with_N_hist_sim=555.RData")

n_POPs <- c(n_POPs, sapply(combined_data, function(x) nrow(x$POPs)))
n_selfs <- c(n_selfs, sapply(combined_data, function(x) nrow(x$self)))
summary(n_POPs)
summary(n_selfs)
