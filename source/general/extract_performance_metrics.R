## =============================================================================
##  Name: extract_performance_metrics.R
##
##  Description: This script extracts a variety of performance metrics from a
##    set of fitted models. The metrics include at least: 
##      + true variance
##      + monte carlo variance
##      + mean error
##      + mean absolute error
##      + mean square error
##      + coefficient of variation
##    These metrics will be derived for male and female abundance in y0 = 2014,
##    the growth rates, and any other parameters of interest. This information 
##    will be stored in a big data.frame, which can be used to create tables 
##    for the manuscript.
## =============================================================================

## =============================================================================
## 1. LOAD THE DATA AND LIBRARIES
## =============================================================================
library(tidyverse)
library(Rfast)
# library(viridis)

load("data/1_population_multiple_sampling_schemes/vanillus/simulation_100_schemes_scenario_1-49_fit_results_sim=144.RData")
load("data/1_population_multiple_sampling_schemes/vanillus/single_combined_data_i=144.RData")

load("data/1_population_multiple_sampling_schemes/gestatii/simulation_100_schemes_scenario_1-49_fit_results_sim=555.RData")
load("data/1_population_multiple_sampling_schemes/gestatii/single_combined_data_i=555.RData")

## =============================================================================
## 2. CREATE THE MASTER DATA FRAME
## =============================================================================

## Extract failed fit scenarios
conv <- sapply(scenario_fits, function(scen) {
  all(sapply(scen, function(fit) fit$message) == "relative convergence (4)")
})
names(conv) <- 1:49

## Extract the truths
true_N_m <- single_combined_data$N_hist[, "N_m"] # mature male abundance
true_N_f <- single_combined_data$N_hist[, "N_f"] # mature female abundance

true_N_m_y0 <- true_N_m[100]
true_N_f_y0 <- true_N_f[100]

true_N_m_10 <- mean(true_N_m[91:100])
true_N_f_10 <- mean(true_N_f[91:100])

true_r_m <- (true_N_m[100] - true_N_m[1]) ^ (1 / 99) / 100 + 1
true_r_f <- (true_N_f[100] - true_N_f[1]) ^ (1 / 99) / 100 + 1

## Extract the estimated abundances in y0
est_N_m_y0 <- sapply(scenario_fits, function(scen) { # mature male abundance
  sapply(scen, function(fit) return(exp(fit$par["N_t0_m"])))
}) 
est_N_m_y0[, !conv] <- NA # set failed convergence estimates to NA

est_N_f_y0 <- sapply(scenario_fits, function(scen) { # mature female abundance
  sapply(scen, function(fit) return(exp(fit$par["N_t0_f"])))
})  
est_N_f_y0[, !conv] <- NA # set failed convergence estimates to NA

## Extract the estimated growth rates
est_r_m <- sapply(scenario_fits, function(scen) { # mature male growth rate
  sapply(scen, function(fit) return(exp(fit$par["r_m"])))
})
est_r_m[, !conv] <- NA # set failed convergence estimates to NA

est_r_f <- sapply(scenario_fits, function(scen) { # mature female growth rate
  sapply(scen, function(fit) return(exp(fit$par["r_f"])))
}) 
est_r_f[, !conv] <- NA # set failed convergence estimates to NA

## ------------------------------------------------------
## Calculate the abundance averages for the past 10 years 
## ------------------------------------------------------
## The male abundance
est_N_m_10 <- sapply(scenario_fits, function(scen) { # mature female abundance
  vec <- sapply(scen, function(fit) {
    N_y0 <- exp(fit$par["N_t0_m"])
    r_10 <- exp(fit$par["r_m"]) ^ (-9:0)
    N_10 <- N_y0 * r_10
    return(mean(N_10))
  })
  ## set to NA if there was failed convergence
  if (!all(sapply(scen, function(fit) fit$message) == "relative convergence (4)")) {
    vec <- matrix(NA, nrow = 1, ncol = 100)
  }
  return(vec)
}, simplify = "array")  

error_N_m_10 <- sapply(scenario_fits, function(scen) { # mature female abundance
  ## set to NA if there was failed convergence
  if (!all(sapply(scen, function(fit) fit$message) == "relative convergence (4)")) {
    mat <- rep(NA, 100)
  } else {
    mat <- sapply(scen, function(fit) {
      N_y0 <- exp(fit$par["N_t0_m"])
      r_10 <- exp(fit$par["r_m"]) ^ (-9:0)
      N_10 <- N_y0 * r_10
      
      N_10_diff <- median(N_10 - true_N_m[91:100])
      return(N_10_diff)
    })
  }
  return(mat)
})  

## The female abundance
est_N_f_10 <- sapply(scenario_fits, function(scen) { # mature female abundance
  vec <- sapply(scen, function(fit) {
    N_y0 <- exp(fit$par["N_t0_f"])
    r_10 <- exp(fit$par["r_f"]) ^ (-9:0)
    N_10 <- N_y0 * r_10
    return(mean(N_10))
  })
  ## set to NA if there was failed convergence
  if (!all(sapply(scen, function(fit) fit$message) == "relative convergence (4)")) {
    vec <- matrix(NA, nrow = 1, ncol = 100)
  }
  return(vec)
}, simplify = "array") 

error_N_f_10 <- sapply(scenario_fits, function(scen) { # mature female abundance
  ## set to NA if there was failed convergence
  if (!all(sapply(scen, function(fit) fit$message) == "relative convergence (4)")) {
    mat <- rep(NA, 100)
  } else {
    mat <- sapply(scen, function(fit) {
      N_y0 <- exp(fit$par["N_t0_f"])
      r_10 <- exp(fit$par["r_f"]) ^ (-9:0)
      N_10 <- N_y0 * r_10
      # cat("N_10:", N_10, "\n")
      
      N_10_diff <- median(N_10 - true_N_f[91:100])
      return(N_10_diff)
    })
  }
  return(mat)
}) 

## Mean error
error_N_m_y0 <- est_N_m_y0 - true_N_m_y0
error_N_f_y0 <- est_N_f_y0 - true_N_f_y0
error_r_m <- est_r_m - true_r_m
error_r_f <- est_r_f - true_r_f

## Mean estimates
median_est <- data.frame(N_f_y0 = colMedians(est_N_f_y0), 
                             N_m_y0 = colMedians(est_N_m_y0),
                             N_f_10 = colMedians(est_N_f_10),
                             N_m_10 = colMedians(est_N_m_10),
                             r_f = colMedians(est_r_f),
                             r_m = colMedians(est_r_m))
mean_est <- data.frame(N_f_y0 = colmeans(est_N_f_y0), 
                         N_m_y0 = colmeans(est_N_m_y0),
                         N_f_10 = colmeans(est_N_f_10),
                         N_m_10 = colmeans(est_N_m_10),
                         r_f = colmeans(est_r_f),
                         r_m = colmeans(est_r_m))


## Average error
median_error <- data.frame(N_f_y0 = colMedians(error_N_f_y0), 
                           N_m_y0 = colMedians(error_N_m_y0),
                           N_f_10 = colMedians(error_N_f_10),
                           N_m_10 = colMedians(error_N_m_10),
                           r_f = colMedians(error_r_f),
                           r_m = colMedians(error_r_m))

mean_error <- data.frame(N_f_y0 = colmeans(error_N_f_y0), 
                         N_m_y0 = colmeans(error_N_m_y0),
                         N_f_10 = colmeans(error_N_f_10),
                         N_m_10 = colmeans(error_N_m_10),
                         r_f = colmeans(error_r_f),
                         r_m = colmeans(error_r_m))


## Mean absolute error
mean_abs_error <- data.frame(N_f_y0 = colmeans(abs(error_N_f_y0)), 
                               N_m_y0 = colmeans(abs(error_N_m_y0)),
                               N_f_10 = colmeans(abs(error_N_f_10)), 
                               N_m_10 = colmeans(abs(error_N_m_10)),
                               r_f = colmeans(abs(error_r_f)),
                               r_m = colmeans(abs(error_r_m)))

## Mean squared error
mean_sqrd_error <- data.frame(N_f_y0 = colmeans((error_N_f_y0) ^ 2), 
                                N_m_y0 = colmeans((error_N_m_y0) ^ 2),
                                N_f_10 = colmeans((error_N_f_10) ^ 2), 
                                N_m_10 = colmeans((error_N_m_10) ^ 2),
                                r_f = colmeans((error_r_f) ^ 2),
                                r_m = colmeans((error_r_m) ^ 2))

## standard deviation and variance
std_dev <- data.frame(N_f_y0 = apply(est_N_f_y0, 2, sd),
                 N_m_y0 = apply(est_N_m_y0, 2, sd), 
                 N_f_10 = apply(est_N_f_10, 2, sd),
                 N_m_10 = apply(est_N_m_10, 2, sd),
                 r_f = apply(est_r_f, 2, sd), 
                 r_m = apply(est_r_m, 2, sd))
variance <- std_dev ^ 2

## Coefficient of variation (works as the support for these estimates is 
## strictly positive)
coeff_var <- data.frame(N_f_y0 = apply(est_N_f_y0, 2, raster::cv),
                 N_m_y0 = apply(est_N_m_y0, 2, raster::cv), 
                 N_f_10 = apply(est_N_f_10, 2, raster::cv),
                 N_m_10 = apply(est_N_m_10, 2, raster::cv), 
                 r_f = apply(est_r_f, 2, raster::cv), 
                 r_m = apply(est_r_m, 2, raster::cv))

## =============================================================================
## 3. CREATE FIGURES AND TABLES FOR MANUSCRIPT
## =============================================================================

scenarios_to_drop <- c(1:7,
                       seq(from=8, to=36, by=7), 
                       seq(from=14, to=42, by=7),
                       43:49)
scenarios_to_keep <- c(1:49)[-scenarios_to_drop]

## :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
## Create a big table with all the results for the supplemenatry materials
## :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

summarise_performance_metrics <- function(column, truth) {
  df <- data.frame(
    true = truth, 
    median_est = median_est[, column], 
    mean_est = mean_est[, column], 
    std_dev = std_dev[, column], 
    coeff_var = coeff_var[, column], 
    median_error = median_error[, column], 
    mean_error = mean_error[, column], 
    mean_abs_error = mean_abs_error[, column], 
    mean_sqrd_error = mean_sqrd_error[, column])
  df <- df[scenarios_to_keep, ]
  df$scenario <- paste0(rep(1:5, each = 5), "-", 1:5)
  df <- na.omit(df[, c(10, 1:9)])
  return(df)
}

summary_N_f_y0 <- summarise_performance_metrics("N_f_y0", true_N_f_y0)

summary_N_m_y0 <- summarise_performance_metrics("N_m_y0", true_N_m_y0)

summary_N_f_10 <- summarise_performance_metrics("N_f_10", true_N_f_10)

summary_N_m_10 <- summarise_performance_metrics("N_m_10", true_N_m_10)

summary_r_f <- summarise_performance_metrics("r_f", true_r_f)

summary_r_m <- summarise_performance_metrics("r_m", true_r_m)

labels <- c("Scen.", 
            "Mdn. Est.", 
            "Mean Est.", 
            "Std. Dev.", 
            "CV", 
            "Mdn. Err.", 
            "Mean Err.", 
            "MAE", 
            "MSE")

caption <- "Various performance metrics for the parameter [quantity], extracted from 100 simulations of the gestating shark population. The columns are, from left to right: scenario, median estimate, mean estimate, standard deviation, coefficient of variation (\\%), median error, mean error, mean absolute error, and mean squared error. Scenarios 1-1, 1-2, 1-3, 2-1, 2-2, and 3-1 were not included as (most of) the simulations did not fit successfully. Scenario 3-3 uses the correct measurement error (2.89) and growth curve specification."
kable(summary_N_f_y0[, -2], booktabs = T, format = "latex", digits = 2,
      col.names = labels, row.names = F, linesep = "", caption = caption)

kable(summary_N_m_y0[, -2], booktabs = T, format = "latex", digits = 2,
      col.names = labels, row.names = F, linesep = "", caption = caption)

kable(summary_N_f_10[, -2], booktabs = T, format = "latex", digits = 2,
      col.names = labels, row.names = F, linesep = "", caption = caption)

kable(summary_N_m_10[, -2], booktabs = T, format = "latex", digits = 2,
      col.names = labels, row.names = F, linesep = "", caption = caption)

kable(summary_r_f[, -2], booktabs = T, format = "latex", digits = 3,
      col.names = labels, row.names = F, linesep = "", caption = caption)

kable(summary_r_m[, -2], booktabs = T, format = "latex", digits = 3,
      col.names = labels, row.names = F, linesep = "", caption = caption)


## :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
## Boxplots for the main manuscript
## :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

## -----------------------------------------
## The abundance in y0 estimates violin plot
## -----------------------------------------
est_N_f_y0_long <- est_N_f_y0[, scenarios_to_keep]
colnames(est_N_f_y0_long) <- paste0(rep(1:5, each = 5), "-", 1:5)
est_N_f_y0_long <- reshape2::melt(est_N_f_y0_long,
                                    value.name = "Abundance", 
                                    na.rm = TRUE)[, -1] # Removes NA scenarios
est_N_f_y0_long$Sex = "Female"

est_N_m_y0_long <- est_N_m_y0[, scenarios_to_keep]
colnames(est_N_m_y0_long) <- paste0(rep(1:5, each = 5), "-", 1:5)
est_N_m_y0_long <- reshape2::melt(est_N_m_y0_long,
                                    value.name = "Abundance", 
                                    na.rm = TRUE)[, -1]
est_N_m_y0_long$Sex = "Male"

est_N_y0_long <- rbind(est_N_m_y0_long, est_N_f_y0_long)
colnames(est_N_y0_long) <- c("Scenario", "Abundance", "Sex")
est_N_y0_long$Scenario <- as.factor(as.character(est_N_y0_long$Scenario))

y0_plot <- ggplot(data=est_N_y0_long, aes(fill=Sex, x=Scenario, y=Abundance)) +
  geom_violin(position="dodge", alpha=0.5, draw_quantiles = c(0.5)) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  geom_hline(yintercept = true_N_f_y0, 
             colour = scales::hue_pal()(2)[1], 
             alpha = 0.5, linewidth = 1)+
  geom_hline(yintercept = true_N_m_y0, 
             colour = scales::hue_pal()(2)[2], 
             alpha = 0.5, linewidth = 1)+
  # coord_flip() +  # comment out for normal plot
  # scale_x_discrete(limits = rev(levels(est_N_y0_long$Scenario))) +
  scale_fill_discrete(name = "") 
  # scale_fill_viridis(discrete=T, name="")

## ----------------------------------------------------------
## The mean abundance over year 1 to 10 estimates violin plot
## ----------------------------------------------------------
est_N_f_10_long <- est_N_f_10[, scenarios_to_keep]
colnames(est_N_f_10_long) <- paste0(rep(1:5, each = 5), "-", 1:5)
est_N_f_10_long <- reshape2::melt(est_N_f_10_long,
                                  value.name = "Abundance", 
                                  na.rm = TRUE)[, -1] # Removes NA scenarios
est_N_f_10_long$Sex = "Female"

est_N_m_10_long <- est_N_m_10[, scenarios_to_keep]
colnames(est_N_m_10_long) <- paste0(rep(1:5, each = 5), "-", 1:5)
est_N_m_10_long <- reshape2::melt(est_N_m_10_long,
                                  value.name = "Abundance", 
                                  na.rm = TRUE)[, -1]
est_N_m_10_long$Sex = "Male"

est_N_10_long <- rbind(est_N_m_10_long, est_N_f_10_long)
colnames(est_N_10_long) <- c("Scenario", "Abundance", "Sex")
est_N_10_long$Scenario <- as.factor(as.character(est_N_10_long$Scenario))

mean_y_10_plot <- ggplot(data=est_N_10_long, aes(fill=Sex, x=Scenario, y=Abundance)) +
  geom_violin(position="dodge", alpha=0.5, draw_quantiles = c(0.5)) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  geom_hline(yintercept = true_N_f_10, 
             colour = scales::hue_pal()(2)[1], 
             alpha = 0.5, linewidth = 1)+
  geom_hline(yintercept = true_N_m_10, 
             colour = scales::hue_pal()(2)[2], 
             alpha = 0.5, linewidth = 1)+
  # coord_flip() +  # comment out for normal plot
  # scale_x_discrete(limits = rev(levels(est_N_10_long$Scenario))) +
  ylim(c(0, 3000)) +
  scale_fill_discrete(name = "") 

## --------------------------------
## The growth estimates violin plot
## --------------------------------
est_r_f_long <- est_r_f[, scenarios_to_keep]
colnames(est_r_f_long) <- paste0(rep(1:5, each = 5), "-", 1:5)
est_r_f_long <- reshape2::melt(est_r_f_long,
                                  value.name = "Growth", 
                                  na.rm = TRUE)[, -1] # Removes NA scenarios
est_r_f_long$Sex = "Female"

est_r_m_long <- est_r_m[, scenarios_to_keep]
colnames(est_r_m_long) <- paste0(rep(1:5, each = 5), "-", 1:5)
est_r_m_long <- reshape2::melt(est_r_m_long,
                                  value.name = "Growth", 
                                  na.rm = TRUE)[, -1]
est_r_m_long$Sex = "Male"

est_r_long <- rbind(est_r_m_long, est_r_f_long)
colnames(est_r_long) <- c("Scenario", "Growth", "Sex")
est_r_long$Scenario <- as.factor(as.character(est_r_long$Scenario))

r_plot <- ggplot(data=est_r_long, aes(fill=Sex, x=Scenario, y=Growth)) +
  geom_violin(position="dodge", alpha=0.5, draw_quantiles = c(0.5)) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  geom_hline(yintercept = true_r_f, 
             colour = scales::hue_pal()(2)[1], 
             alpha = 0.5, linewidth = 1)+
  geom_hline(yintercept = true_r_m, 
             colour = scales::hue_pal()(2)[2], 
             alpha = 0.5, linewidth = 1)+
  # coord_flip() +  # comment out for normal plot, together with line below
  # scale_x_discrete(limits = rev(levels(est_r_long$Scenario))) +
  scale_fill_discrete(name = "") 

## -------------------------------
## Combine all plots into one plot
## -------------------------------

library(gridExtra)
# library(svglite)
p <- grid.arrange(y0_plot + 
                    theme(legend.position = "none") +
                    ylab("Abundance in 2014") +
                    xlab(element_blank()) +
                    coord_cartesian(ylim=c(0, 3000)), # truncate at 3000 without removing values outside
                  mean_y_10_plot + 
                    theme(legend.position = "none") +
                    ylab("Mean abundance (2005-2014)") +
                    xlab(element_blank()) +
                    coord_cartesian(ylim=c(0, 3000)), # truncate at 3000 without removing values outside
                  r_plot + 
                    ylab("Yearly growth rate"),
                  nrow=3, ncol=1); p;
# now export it with width 1000 and height 750


## :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
## Now  the errors [currently commented out]
## ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# error_N_f_y0_long <- error_N_f_y0[, scenarios_to_keep]
# colnames(error_N_f_y0_long) <- paste0(rep(1:5, each = 5), "-", 1:5)
# error_N_f_y0_long <- reshape2::melt(error_N_f_y0_long,
#                                     value.name = "Abundance", 
#                                     na.rm = TRUE)[, -1] # Removes NA scenarios
# error_N_f_y0_long$Sex = "Female"
# 
# error_N_m_y0_long <- error_N_m_y0[, scenarios_to_keep]
# colnames(error_N_m_y0_long) <- paste0(rep(1:5, each = 5), "-", 1:5)
# error_N_m_y0_long <- reshape2::melt(error_N_m_y0_long,
#                                     value.name = "Abundance", 
#                                     na.rm = TRUE)[, -1]
# error_N_m_y0_long$Sex = "Male"
# 
# error_N_y0_long <- rbind(error_N_m_y0_long, error_N_f_y0_long)
# colnames(error_N_y0_long) <- c("Scenario", "Abundance", "Sex")
# error_N_y0_long$Scenario <- as.factor(as.character(error_N_y0_long$Scenario))
# 
# ggplot(data=error_N_y0_long, aes(fill=Sex, x=Scenario, y=Abundance)) +
#   geom_violin(position="dodge", alpha=0.5, draw_quantiles = c(0.5)) +
#   theme_minimal() +
#   geom_hline(yintercept = 0, colour = "red", alpha = 0.5)+
#   scale_x_discrete(limits = rev(levels(error_N_y0_long$Scenario))) +
#   scale_fill_viridis(discrete=T, name="") +
#   coord_flip()


