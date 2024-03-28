## =============================================================================
##  Name: extract_performance_metrics_LINK.R
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
library(kableExtra)

## Load two lines below for the simple species
# load("data/simulation_study/simple/simulation_1000_schemes_all_scenarios_fit_results_sim=all_no_growth.RData")
# load("data/simulation_study/simple/1000_schemes_combined_data_with_N_hist_sim=all.RData")

## Load the two lines below for the complex species
load("data/simulation_study/complex/simulation_1000_schemes_all_scenarios_fit_results_sim=all_no_growth.RData")
load("data/simulation_study/complex/1000_schemes_combined_data_with_N_hist_sim=all.RData")

## =============================================================================
## 2. CREATE THE MASTER DATA FRAME
## =============================================================================

## Extract failed fit scenarios
conv <- sapply(scenario_fits, function(scen) {
  all(sapply(scen, function(fit) fit$message) == "relative convergence (4)")
})
names(conv) <- 1:25

## Extract the truths
true_N_m <- t(sapply(combined_data, function(dat) return(dat$N_hist[, "N_m"]))) # mature male abundance
true_N_f <- t(sapply(combined_data, function(dat) return(dat$N_hist[, "N_f"]))) # mature male abundance

true_N_m_y0 <- log(true_N_m[, 100])
true_N_f_y0 <- log(true_N_f[, 100])

true_N_m_10 <- log(rowMeans(true_N_m[, 91:100]))
true_N_f_10 <- log(rowMeans(true_N_f[, 91:100]))

true_r_m <- log((true_N_m[, 100] / true_N_m[, 1]) ^ (1 / 99))
true_r_f <- log((true_N_f[, 100] / true_N_f[, 1]) ^ (1 / 99))

## Extract the estimated abundances in y0
est_N_m_y0 <- sapply(scenario_fits, function(scen) { # mature male abundance
  sapply(scen, function(fit) return(fit$par["N_t0_m"]))
}) 
est_N_m_y0[, !conv] <- NA # set failed convergence estimates to NA

est_N_f_y0 <- sapply(scenario_fits, function(scen) { # mature female abundance
  sapply(scen, function(fit) return(fit$par["N_t0_f"]))
})  
est_N_f_y0[, !conv] <- NA # set failed convergence estimates to NA

## Extract the estimated growth rates
est_r_m <- sapply(scenario_fits, function(scen) { # mature male growth rate
  sapply(scen, function(fit) return(fit$par["r_m"]))
})
est_r_m[, !conv] <- NA # set failed convergence estimates to NA

est_r_f <- sapply(scenario_fits, function(scen) { # mature female growth rate
  sapply(scen, function(fit) return(fit$par["r_f"]))
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
    return(log(mean(N_10)))
  })
  ## set to NA if there was failed convergence
  if (!all(sapply(scen, function(fit) fit$message) == "relative convergence (4)")) {
    vec <- matrix(NA, nrow = 1, ncol = 1000)
  }
  return(vec)
}, simplify = "array")  

error_N_m_10 <- est_N_m_10 - matrix(true_N_m_10, ncol = 25, nrow = 1000)

bias_N_m_10 <- error_N_m_10 / matrix(true_N_m_10, ncol = 25, nrow = 1000) * 100

## The female abundance
est_N_f_10 <- sapply(scenario_fits, function(scen) { # mature female abundance
  vec <- sapply(scen, function(fit) {
    N_y0 <- exp(fit$par["N_t0_f"])
    r_10 <- exp(fit$par["r_f"]) ^ (-9:0)
    N_10 <- N_y0 * r_10
    return(log(mean(N_10)))
  })
  ## set to NA if there was failed convergence
  if (!all(sapply(scen, function(fit) fit$message) == "relative convergence (4)")) {
    vec <- matrix(NA, nrow = 1, ncol = 1000)
  }
  return(vec)
}, simplify = "array")  

error_N_f_10 <- est_N_f_10 - matrix(true_N_f_10, ncol = 25, nrow = 1000)

bias_N_f_10 <- error_N_f_10 / matrix(true_N_f_10, ncol = 25, nrow = 1000) * 100

## Mean error
error_N_m_y0 <- est_N_m_y0 - matrix(true_N_m_y0, nrow = 1000, ncol = 25, byrow = F) 
error_N_f_y0 <- est_N_f_y0 - matrix(true_N_f_y0, nrow = 1000, ncol = 25, byrow = F) 
error_r_m <- est_r_m - matrix(true_r_m, nrow = 1000, ncol = 25, byrow = F) 
error_r_f <- est_r_f - matrix(true_r_f, nrow = 1000, ncol = 25, byrow = F) 
# new: bias
bias_N_m_y0 <- error_N_m_y0 / matrix(true_N_m_y0, nrow = 1000, ncol = 25, byrow = F) * 100
bias_N_f_y0 <- error_N_f_y0 / matrix(true_N_f_y0, nrow = 1000, ncol = 25, byrow = F) * 100
bias_r_m <- error_r_m / matrix(true_r_m, nrow = 1000, ncol = 25, byrow = F) * 100
bias_r_f <- error_r_f / matrix(true_r_f, nrow = 1000, ncol = 25, byrow = F) * 100

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


## Average bias
median_bias <- data.frame(N_f_y0 = colMedians(bias_N_f_y0), 
                          N_m_y0 = colMedians(bias_N_m_y0),
                          N_f_10 = colMedians(bias_N_f_10),
                          N_m_10 = colMedians(bias_N_m_10),
                          r_f = colMedians(bias_r_f),
                          r_m = colMedians(bias_r_m))

mean_bias <- data.frame(N_f_y0 = colmeans(bias_N_f_y0), 
                        N_m_y0 = colmeans(bias_N_m_y0),
                        N_f_10 = colmeans(bias_N_f_10),
                        N_m_10 = colmeans(bias_N_m_10),
                        r_f = colmeans(bias_r_f),
                        r_m = colmeans(bias_r_m))

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

## standard deviation and variance of estimates
std_dev <- data.frame(N_f_y0 = apply(est_N_f_y0, 2, sd),
                      N_m_y0 = apply(est_N_m_y0, 2, sd), 
                      N_f_10 = apply(est_N_f_10, 2, sd),
                      N_m_10 = apply(est_N_m_10, 2, sd),
                      r_f = apply(est_r_f, 2, sd), 
                      r_m = apply(est_r_m, 2, sd))
variance <- std_dev ^ 2

## standard deviation and variance of error
std_dev_error <- data.frame(N_f_y0 = apply(error_N_f_y0, 2, sd),
                            N_m_y0 = apply(error_N_m_y0, 2, sd), 
                            N_f_10 = apply(error_N_f_10, 2, sd),
                            N_m_10 = apply(error_N_m_10, 2, sd),
                            r_f = apply(error_r_f, 2, sd), 
                            r_m = apply(error_r_m, 2, sd))
variance_error <- std_dev_error ^ 2

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

## :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
## Boxplots for the main manuscript 
##
## Initially we had estimated abundance, but these no longer work with the new
## simulation strategy. Therefore, we now make similar plots, but based the 
## error. 
## :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

## ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: 
## Below are the new bias plots
## :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
## -----------------------------------------
## The error in y0 estimates violin plot
## -----------------------------------------
bias_N_f_y0_long <- bias_N_f_y0
colnames(bias_N_f_y0_long) <- paste0(rep(1:5, each = 5), "-", 1:5)
bias_N_f_y0_long <- reshape2::melt(bias_N_f_y0_long,
                                   # na.rm = TRUE
                                   value.name = "Abundance")[, -1]

bias_N_f_y0_long$Measurment_scenario <- rep(c("-67% of true length measurement error (2.89)",
                                              "-33% of true length measurement error (2.89)",
                                              "The true length measurement error (2.89)",
                                              "+33% of true length measurement error (2.89)",
                                              "+67% of true length measurement error (2.89)"),
                                            each = 5000)
bias_N_f_y0_long$VBGF_scenario <- rep(rep(c("Growth curve shifted vertically by -10%",
                                            "Growth curve shifted vertically by -5%",
                                            "The true growth curve",
                                            "Growth curve shifted vertically by +5%",
                                            "Growth curve shifted vertically by +10%"),
                                          each = 1000), 
                                      times = 5)
bias_N_f_y0_long$Sex = "Female"
bias_N_f_y0_long <- na.omit(bias_N_f_y0_long)

bias_N_m_y0_long <- bias_N_m_y0
colnames(bias_N_m_y0_long) <- paste0(rep(1:5, each = 5), "-", 1:5)
bias_N_m_y0_long <- reshape2::melt(bias_N_m_y0_long,
                                   value.name = "Abundance")[, -1]

bias_N_m_y0_long$Measurment_scenario <- rep(c("-67% of true length measurement error (2.89)",
                                              "-33% of true length measurement error (2.89)",
                                              "The true length measurement error (2.89)",
                                              "+33% of true length measurement error (2.89)",
                                              "+67% of true length measurement error (2.89)"),
                                            each = 5000)
bias_N_m_y0_long$VBGF_scenario <- rep(rep(c("Growth curve shifted vertically by -10%",
                                            "Growth curve shifted vertically by -5%",
                                            "The true growth curve",
                                            "Growth curve shifted vertically by +5%",
                                            "Growth curve shifted vertically by +10%"),
                                          each = 1000), 
                                      times = 5)
bias_N_m_y0_long$Sex = "Male"
bias_N_m_y0_long <- na.omit(bias_N_m_y0_long)

bias_N_y0_long <- rbind(bias_N_m_y0_long[, c(1, 3, 4, 2, 5)], 
                        bias_N_f_y0_long[, c(1, 3, 4, 2, 5)])
colnames(bias_N_y0_long)[1] <- c("Scenario")
bias_N_y0_long$Scenario <- as.factor(as.character(bias_N_y0_long$Scenario))

y0_plot <- ggplot(data=bias_N_y0_long, aes(fill=Sex, x=Scenario, y=Abundance)) +
  geom_boxplot(position="dodge", 
               # draw_quantiles = c(0.5), # for violin plots
               coef = 1.5,   # for boxplots, how many times the IQR values are included (defaults to 1.5)
               alpha=0.5) +
  stat_summary(mapping = aes(fill=Sex), position=position_dodge(0.75),
               fun = mean, geom="point", shape=21, color="black") +
  theme_minimal() +
  theme(legend.position = "bottom") +
  geom_hline(yintercept = 0, 
             colour = "black", 
             alpha = 0.5)+
  ylab("Relative error in abundance (%)") +
  # coord_flip() +  # comment out for normal plot
  # scale_x_discrete(limits = rev(levels(est_N_10_long$Scenario))) +
  # coord_cartesian(ylim = c(-100, 500)) +
  scale_fill_discrete(name = "")

## ----------------------------------------------------------
## The mean abundance over year 1 to 10 estimates violin plot
## ----------------------------------------------------------
bias_N_f_10_long <- bias_N_f_10
colnames(bias_N_f_10_long) <- paste0(rep(1:5, each = 5), "-", 1:5)
bias_N_f_10_long <- reshape2::melt(bias_N_f_10_long,
                                   # na.rm = TRUE
                                   value.name = "Abundance")[, -1]

bias_N_f_10_long$Measurment_scenario <- rep(c("-67% of true length measurement error (2.89)",
                                              "-33% of true length measurement error (2.89)",
                                              "The true length measurement error (2.89)",
                                              "+33% of true length measurement error (2.89)",
                                              "+67% of true length measurement error (2.89)"),
                                            each = 5000)
bias_N_f_10_long$VBGF_scenario <- rep(rep(c("Growth curve shifted vertically by -10%",
                                            "Growth curve shifted vertically by -5%",
                                            "The true growth curve",
                                            "Growth curve shifted vertically by +5%",
                                            "Growth curve shifted vertically by +10%"),
                                          each = 1000), 
                                      times = 5)
bias_N_f_10_long$Sex = "Female"
bias_N_f_10_long <- na.omit(bias_N_f_10_long)

bias_N_m_10_long <- bias_N_m_10
colnames(bias_N_m_10_long) <- paste0(rep(1:5, each = 5), "-", 1:5)
bias_N_m_10_long <- reshape2::melt(bias_N_m_10_long,
                                   value.name = "Abundance")[, -1]

bias_N_m_10_long$Measurment_scenario <- rep(c("-67% of true length measurement error (2.89)",
                                              "-33% of true length measurement error (2.89)",
                                              "The true length measurement error (2.89)",
                                              "+33% of true length measurement error (2.89)",
                                              "+67% of true length measurement error (2.89)"),
                                            each = 5000)
bias_N_m_10_long$VBGF_scenario <- rep(rep(c("Growth curve shifted vertically by -10%",
                                            "Growth curve shifted vertically by -5%",
                                            "The true growth curve",
                                            "Growth curve shifted vertically by +5%",
                                            "Growth curve shifted vertically by +10%"),
                                          each = 1000), 
                                      times = 5)
bias_N_m_10_long$Sex = "Male"
bias_N_m_10_long <- na.omit(bias_N_m_10_long)

bias_N_10_long <- rbind(bias_N_m_10_long[, c(1, 3, 4, 2, 5)], 
                        bias_N_f_10_long[, c(1, 3, 4, 2, 5)])
colnames(bias_N_10_long)[1] <- c("Scenario")
bias_N_10_long$Scenario <- as.factor(as.character(bias_N_10_long$Scenario))

mean_y_10_plot <- ggplot(data=bias_N_10_long, aes(fill=Sex, x=Scenario, y=Abundance)) +
  geom_boxplot(position="dodge", 
               # draw_quantiles = c(0.5), # for violin plots
               coef = 1.5,   # for boxplots, how many times the IQR values are included (defaults to 1.5)
               alpha=0.5) +
  stat_summary(mapping = aes(fill=Sex), position=position_dodge(0.75),
               fun = mean, geom="point", shape=21, color="black") +
  theme_minimal() +
  theme(legend.position = "bottom") +
  geom_hline(yintercept = 0, 
             colour = "black", 
             alpha = 0.5)+
  # coord_flip() +  # comment out for normal plot
  # scale_x_discrete(limits = rev(levels(est_N_10_long$Scenario))) +
  # coord_cartesian(ylim = c(-100, 500)) +
  scale_fill_discrete(name = "") 


## --------------------------------
## The growth estimates violin plot
## --------------------------------
bias_r_f_long <- bias_r_f
colnames(bias_r_f_long) <- paste0(rep(1:5, each = 5), "-", 1:5)
bias_r_f_long <- reshape2::melt(bias_r_f_long,
                                # na.rm = TRUE
                                value.name = "Growth")[, -1]

bias_r_f_long$Measurment_scenario <- rep(c("-67% of true length measurement error (2.89)",
                                           "-33% of true length measurement error (2.89)",
                                           "The true length measurement error (2.89)",
                                           "+33% of true length measurement error (2.89)",
                                           "+67% of true length measurement error (2.89)"),
                                         each = 5000)
bias_r_f_long$VBGF_scenario <- rep(rep(c("Growth curve shifted vertically by -10%",
                                         "Growth curve shifted vertically by -5%",
                                         "The true growth curve",
                                         "Growth curve shifted vertically by +5%",
                                         "Growth curve shifted vertically by +10%"),
                                       each = 1000), 
                                   times = 5)
bias_r_f_long$Sex = "Female"
bias_r_f_long <- na.omit(bias_r_f_long)

bias_r_m_long <- bias_r_m
colnames(bias_r_m_long) <- paste0(rep(1:5, each = 5), "-", 1:5)
bias_r_m_long <- reshape2::melt(bias_r_m_long,
                                value.name = "Growth")[, -1]

bias_r_m_long$Measurment_scenario <- rep(c("-67% of true length measurement error (2.89)",
                                           "-33% of true length measurement error (2.89)",
                                           "The true length measurement error (2.89)",
                                           "+33% of true length measurement error (2.89)",
                                           "+67% of true length measurement error (2.89)"),
                                         each = 5000)
bias_r_m_long$VBGF_scenario <- rep(rep(c("Growth curve shifted vertically by -10%",
                                         "Growth curve shifted vertically by -5%",
                                         "The true growth curve",
                                         "Growth curve shifted vertically by +5%",
                                         "Growth curve shifted vertically by +10%"),
                                       each = 1000), 
                                   times = 5)
bias_r_m_long$Sex = "Male"
bias_r_m_long <- na.omit(bias_r_m_long)

bias_r_long <- rbind(bias_r_m_long[, c(1, 3, 4, 2, 5)], 
                     bias_r_f_long[, c(1, 3, 4, 2, 5)])
colnames(bias_r_long)[1] <- c("Scenario")
bias_r_long$Scenario <- as.factor(as.character(bias_r_long$Scenario))

r_plot <- ggplot(data=bias_r_long, aes(fill=Sex, x=Scenario, y=Growth)) +
  geom_boxplot(position="dodge", 
               # draw_quantiles = c(0.5), # for violin plots
               coef = 1.5,   # for boxplots, how many times the IQR values are included (defaults to 1.5)
               alpha=0.5) +
  stat_summary(mapping = aes(fill=Sex), position=position_dodge(0.75),
               fun = mean, geom="point", shape=21, color="black") +
  theme_minimal() +
  theme(legend.position = "bottom") +
  geom_hline(yintercept = 0, 
             colour = "black", 
             alpha = 0.5)+
  ylab("Relative error in growth rate") +
  # coord_flip() +  # comment out for normal plot
  # scale_x_discrete(limits = rev(levels(est_N_10_long$Scenario))) +
  scale_fill_discrete(name = "") 
## -------------------------------
## Combine all plots into one plot
## -------------------------------

library(gridExtra)
# library(svglite)
lay <- rbind(c(1),
             c(1),
             c(2))

p <- grid.arrange(y0_plot + 
                    theme(legend.position = "none") +
                    ylab("Rel. error abundance 2014") +
                    xlab(element_blank()) +
                    coord_cartesian(ylim=c(-100, 1000)), # truncate without removing values outside
                  # mean_y_10_plot +
                  #   theme(legend.position = "none") +
                  #   ylab("Rel. error abundance (2005-2014)") +
                  #   xlab(element_blank()) +
                  #   coord_cartesian(ylim=c(-100, 1000)), # truncate  without removing values outside
                  r_plot + 
                    coord_cartesian(ylim=c(-50, 50))+
                    ylab("Rel. error yearly growth rate"),
                  layout_matrix = lay); p

## now export it with width 1000 and height 750


## :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
## Create a big table with all the results for the supplementary materials
## :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

summarise_performance_metrics <- function(column, truth) {
  df <- data.frame(
    # true = truth, 
    # median_est = median_est[, column],
    # mean_est = mean_est[, column],
    # std_dev = std_dev[, column],
    # coeff_var = coeff_var[, column],
    bias = median_bias[, column],
    median_error = median_error[, column], 
    mean_error = mean_error[, column], 
    mean_abs_error = mean_abs_error[, column], 
    mean_sqrd_error = mean_sqrd_error[, column],
    error_std_dev = std_dev_error[, column])
  
  df$scenario <- paste0(rep(1:5, each = 5), "-", 1:5)
  df <- na.omit(df[, c(ncol(df), 1:(ncol(df)-1))])
  return(df)
}

summary_N_f_y0 <- summarise_performance_metrics("N_f_y0", true_N_f_y0)

summary_N_m_y0 <- summarise_performance_metrics("N_m_y0", true_N_m_y0)

summary_N_f_10 <- summarise_performance_metrics("N_f_10", true_N_f_10)

summary_N_m_10 <- summarise_performance_metrics("N_m_10", true_N_m_10)

summary_r_f <- summarise_performance_metrics("r_f", true_r_f)

summary_r_m <- summarise_performance_metrics("r_m", true_r_m)

labels <- c("Scen.", 
            # "Mdn. Est.", 
            # "Mean Est.", 
            # "Std. Dev.", 
            # "CV", 
            "Bias",
            "Mdn. Err.", 
            "Mean Err.", 
            "MAE", 
            "MSE",
            "Err. Std. Dev.")

caption <- "Various performance metrics for the parameter [quantity], extracted from 100 simulations of the gestating shark population. The columns are, from left to right: scenario, median estimate, mean estimate, standard deviation, coefficient of variation (\\%), median error, mean error, mean absolute error, and mean squared error. Scenarios 1-1, 1-2, 1-3, 2-1, 2-2, and 3-1 were not included as (most of) the simulations did not fit successfully. Scenario 3-3 uses the correct measurement error (2.89) and growth curve specification."
kable(summary_N_f_y0, booktabs = T, format = "latex", digits = 2,
      col.names = labels, row.names = F, linesep = "", caption = caption)

kable(summary_N_m_y0, booktabs = T, format = "latex", digits = 2,
      col.names = labels, row.names = F, linesep = "", caption = caption)

kable(summary_N_f_10, booktabs = T, format = "latex", digits = 2,
      col.names = labels, row.names = F, linesep = "", caption = caption)

kable(summary_N_m_10, booktabs = T, format = "latex", digits = 2,
      col.names = labels, row.names = F, linesep = "", caption = caption)

kable(summary_r_f, booktabs = T, format = "latex", digits = 3,
      col.names = labels, row.names = F, linesep = "", caption = caption)

kable(summary_r_m, booktabs = T, format = "latex", digits = 3,
      col.names = labels, row.names = F, linesep = "", caption = caption)
