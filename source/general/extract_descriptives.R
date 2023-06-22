## =============================================================================
##  Name: extract_descriptives.R
##
##  Description: This script extracts a variety of descriptive statistics from 
##    the simulations of the complex and vanilla species. Descriptives of 
##    interest are:
##      - total number of POPs
##      - total number of recaptures
##      - total adult (male and female) abundance in year 2014
##      - realised growth rates
## =============================================================================

## =============================================================================
## 1. LOAD THE DATA AND LIBRARIES
## =============================================================================
library(tidyverse)
library(Rfast)
library(kableExtra)
library(CKMRcpp)

## Vanilla
load("D:/felix/CKMR_grey_reef_sharks/data/simulation_study/vanilla/1000_schemes_combined_data_with_N_hist_sim=all.RData")
combined_data_vanilla <- combined_data

n_POPs <- sapply(combined_data_vanilla, function(x) nrow(x$POPs))

n_selfs <- sapply(combined_data_vanilla, function(x) nrow(x$self))

true_N_m <- t(sapply(combined_data_vanilla, function(x) return(x$N_hist[, "N_m"]))) # mature male abundance
true_N_f <- t(sapply(combined_data_vanilla, function(x) return(x$N_hist[, "N_f"]))) # mature female abundance

N_male_2014 <- true_N_m[, 100]
N_female_2014 <- true_N_f[, 100]

r_male <- (true_N_m[, 100] / true_N_m[, 1]) ^ (1 / 99) + 1
r_female <- (true_N_f[, 100] / true_N_m[, 1]) ^ (1 / 99) + 1

summary(n_POPs) # 25-76
summary(n_selfs) # 4-25
summary(N_male_2014) # 600-1019
summary(N_female_2014) # 630-992
summary(r_male) - 1 # -0.3% to +0.3%
summary(r_female) - 1 #-0.3% to +0.3%


## Complex
load("D:/felix/CKMR_grey_reef_sharks/data/simulation_study/complex/1000_schemes_combined_data_with_N_hist_sim=all.RData")
combined_data_complex <- combined_data

n_POPs <- sapply(combined_data_complex, function(x) nrow(x$POPs))

n_selfs <- sapply(combined_data_complex, function(x) nrow(x$self))

true_N_m <- t(sapply(combined_data_complex, function(x) return(x$N_hist[, "N_m"]))) # mature male abundance
true_N_f <- t(sapply(combined_data_complex, function(x) return(x$N_hist[, "N_f"]))) # mature female abundance

N_male_2014 <- true_N_m[, 100]
N_female_2014 <- true_N_f[, 100]

r_male <- (true_N_m[, 100] / true_N_m[, 1]) ^ (1 / 99) + 1
r_female <- (true_N_f[, 100] / true_N_m[, 1]) ^ (1 / 99) + 1

summary(n_POPs) # 25-76
summary(n_selfs) # 4-25
summary(N_male_2014) # 600-1019
summary(N_female_2014) # 630-992
summary(r_male) - 1 # -0.3% to +0.3%
summary(r_female) - 1 #-0.3% to +0.3%
  