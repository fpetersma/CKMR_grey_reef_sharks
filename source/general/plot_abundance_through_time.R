## =============================================================================
## 1. LOAD THE DATA AND LIBRARIES
## =============================================================================
library(tidyverse)
library(Rfast)
library(ggplot2)

## Run 2 lines below for simple species
# load("data/simulation_study/simple/simulation_1000_schemes_all_scenarios_fit_results_sim=all.RData")
# load("data/simulation_study/simple/1000_schemes_combined_data_with_N_hist_sim=all.RData")

## Run 2 lines below for complex species
load("data/simulation_study/complex/simulation_1000_schemes_all_scenarios_fit_results_sim=all.RData")
load("data/simulation_study/complex/1000_schemes_combined_data_with_N_hist_sim=all.RData")

## =============================================================================
## 2. CREATE THE MASTER DATA FRAME
## =============================================================================

## Set constants
year_lim <- c(-19, 0)
n_years <- year_lim[2] - year_lim[1] + 1
y0 <- 100
max_y_axis <- 1000

scen_names <- paste0(rep(paste0(rep("ME", 5), c("-67", "-33", "+0", "+33", "+67")), 
                         each = 5), ":",
                     rep(paste0(rep("GC", 5), c("-10", "-5", "+0", "+5", "+10")), 
                         times = 5))


## -----------------------
## Create long format data
## -----------------------

## Extract the true abundances for male and female respectively
true_N_m <- t(sapply(combined_data, function(dat) return(dat$N_hist[, "N_m"])))
true_N_f <- t(sapply(combined_data, function(dat) return(dat$N_hist[, "N_f"]))) 

## Loop over the fits
for (i in 1:length(scenario_fits)) {
  fits <- scenario_fits[[i]]
  
  ## Extract estimates for selected parameter
  est <- t(sapply(fits, function(x) x$par))
  conv <- sapply(fits, function(x) x$message)
  objective <- sapply(fits, function(x) x$objective)
  
  ## Derive years and create matrices for abundance estimates
  years <- year_lim[1]:year_lim[2]
  est_N_f <- matrix(NA, nrow = length(fits), ncol = length(years))
  est_N_m <- matrix(NA, nrow = length(fits), ncol = length(years))
  
  ## For every fit, derive the abundance for all years
  for (j in 1:length(fits)) {
    # est_N_f[j, ] <- exp(est[j, "N_t0_f"])
    # est_N_m[j, ] <- exp(est[j, "N_t0_m"]) 
    est_N_f[j, ] <- exp(est[j, "N_t0_f"]) * exp(est[j, "r_f"]) ^ years
    est_N_m[j, ] <- exp(est[j, "N_t0_m"]) * exp(est[j, "r_m"]) ^ years
  }
  
  ## Subtract the true abundances to get the error
  error_N_f <- est_N_f - true_N_f[, (100-ncol(est_N_f)+1):100]
  error_N_m <- est_N_m - true_N_m[, (100-ncol(est_N_m)+1):100]
  
  ## Turn the data in long format so ggplot2 knows what to do
  error_f_df <- cbind(data.frame(year = y0 + years), 
                      t(error_N_f),
                      "1001" = Rfast::colMedians(error_N_f),# median
                      "1002" = 0)                        # correct = no error
  error_f_long <- reshape2::melt(error_f_df, value.name = "N_error", id.vars = "year")
  error_f_long$sim_id <- i
  error_f_long$sex <- "F"
  
  ## Turn the data in long format so ggplot2 knows what to do
  error_m_df <- cbind(data.frame(year = y0 + years), 
                      t(error_N_m),
                      "1001" = Rfast::colMedians(error_N_m),# median
                      "1002" = 0)                        # correct = no error
  error_m_long <- reshape2::melt(error_m_df, value.name = "N_error", id.vars = "year")
  error_m_long$sim_id <- i
  error_m_long$sex <- "M"
  
  # if (any(conv != "relative convergence (4)")) {
  if (any(is.infinite(objective))) {
    error_m_long$conv <- "failed"    
    error_f_long$conv <- "failed" 
  } else {
    error_m_long$conv <- "successful"    
    error_f_long$conv <- "successful" 
  }
  
  # Combine the male and female data sets
  # If it's the first simulation scenario, create error_long; else, append error_long
  if(i == 1) {
    error_long <- rbind(error_f_long, error_m_long)
  } else {
    error_long <- rbind(error_long, error_f_long, error_m_long)
  }
}

## Set failed abundance estimates to NA
error_long$N_error[error_long$conv == "failed"] <- NA

dimension <- sqrt(max(error_long$sim_id))

error_long$sim_id <- as.factor(error_long$sim_id)
levels(error_long$sim_id) <- scen_names

## Add the scenario descriptions
error_long$Measurment_scenario <- rep(c("-67% of true length measurement error (2.89)",
                                          "-33% of true length measurement error (2.89)",
                                          "The true length measurement error (2.89)",
                                          "+33% of true length measurement error (2.89)",
                                          "+67% of true length measurement error (2.89)"),
                                        each = n_years * 1002 * 2 * dimension)
error_long$VBGF_scenario <- rep(rep(c("Growth curve shifted vertically by -10%",
                                        "Growth curve shifted vertically by -5%",
                                        "The true growth curve",
                                        "Growth curve shifted vertically by +5%",
                                        "Growth curve shifted vertically by +10%"),
                                      each = n_years * 1002 * 2), 
                                  times = dimension)

## =============================================================================
## 3 START THE PLOTTING
## =============================================================================

p_m <- ggplot(subset(error_long, sex == "M")) +
  geom_line(
    mapping = aes(x = year, y = N_error, colour = variable),
    show.legend = F, linewidth = 1) +
  scale_color_manual(values = c(rep(alpha("darkgrey", 0.1), 1000),
                                alpha("black", 0.5),
                                alpha("red", 0.4))) +
  theme_bw() +
  ylab("Error in adult abundance estimate") +
  xlab("Year") +
  facet_wrap(~ sim_id , nrow = dimension, labeller = label_parsed) + 
  coord_cartesian(ylim=c(-1000,3000)) +
  scale_x_continuous(breaks = seq(from = min(years) + y0, to = max(years) + y0, by = 5),
                     labels = seq(from = min(years) + y0, to = max(years) + y0, by = 5))

p_f <- ggplot(subset(error_long, sex == "F")) +
  geom_line(
    mapping = aes(x = year, y = N_error, colour = variable),
    show.legend = F, linewidth = 1) +
  scale_color_manual(values = c(rep(alpha("darkgrey", 0.1), 1000),
                                alpha("black", 0.5),
                                alpha("red", 0.4))) +
  theme_bw() +
  ylab("Error in adult abundance estimate") +
  xlab("Year") +
  facet_wrap(~ sim_id , nrow = dimension, labeller = label_parsed) + 
  coord_cartesian(ylim=c(-1000,3000)) +
  scale_x_continuous(breaks = seq(from = min(years) + y0, to = max(years) + y0, by = 5),
                     labels = seq(from = min(years) + y0, to = max(years) + y0, by = 5))
 
p_m
p_f

## Save the plots as svg with dimensions 1100x700 [manually]
## name: "[sex]_error_trend_[species].svg" (eg "female_error_trend_simple.svg")

hist(sapply(scenario_fits[[13]], function(fit) fit$par[,3]))













