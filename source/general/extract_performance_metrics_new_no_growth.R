## =============================================================================
##  Name: extract_performance_metrics_new.R
##
##  Description: This script extracts a variety of performance metrics from a
##    set of fitted models. The metrics include at least: 
##      + true variance
##      + monte carlo variance
##      + mean error
##      + mean absolute error
##      + mean square error
##      + coefficient of variation
##    These metrics will be derived for male and female abundance in y0 = 100,
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

NO_GROWTH <- TRUE

scen_names <- paste0(rep(paste0(rep("ME", 5), c("-67", "-33", "+0", "+33", "+67")), 
                         each = 5), ":",
                     rep(paste0(rep("GC", 5), c("-10", "-5", "+0", "+5", "+10")), 
                         times = 5))

## Load simulation data and fit results for simple species
load("data/simulation_study/simple/simulation_1000_schemes_all_scenarios_fit_results_sim=all_no_growth.RData")
load("data/simulation_study/simple/1000_schemes_combined_data_with_N_hist_sim=all.RData")
simple_fits <- scenario_fits
simple_sims <- combined_data

## Load simulation data and fit results for complex species
load("data/simulation_study/complex/simulation_1000_schemes_all_scenarios_fit_results_sim=all_no_growth.RData")
load("data/simulation_study/complex/1000_schemes_combined_data_with_N_hist_sim=all.RData")
complex_fits <- scenario_fits
complex_sims <- combined_data

## Clean environment
rm(scenario_fits, combined_data)

## =============================================================================
## 2. CREATE THE MASTER DATA FRAME
## =============================================================================

############################
 ##### SIMPLE SPECIES #####
############################

## Extract failed fit scenarios
conv <- sapply(simple_fits, function(scen) {
  all(sapply(scen, function(fit) fit$message) == "relative convergence (4)")
})

scenarios_to_keep <- (1:25)[conv] # convert TRUE/FALSE to indices

## Extract the truths
true_N_m <- t(sapply(simple_sims, function(dat) return(dat$N_hist[, "N_m"]))) # mature male abundance
true_N_f <- t(sapply(simple_sims, function(dat) return(dat$N_hist[, "N_f"]))) # mature female abundance

true_N_m_y0 <- true_N_m[, 100]
true_N_f_y0 <- true_N_f[, 100]

## Extract the estimated abundances in y0
est_N_m_y0 <- sapply(simple_fits, function(scen) { # mature male abundance
  sapply(scen, function(fit) return(exp(fit$par["N_t0_m"])))
}) 
est_N_m_y0[, !conv] <- NA # set failed convergence estimates to NA

est_N_f_y0 <- sapply(simple_fits, function(scen) { # mature female abundance
  sapply(scen, function(fit) return(exp(fit$par["N_t0_f"])))
})  
est_N_f_y0[, !conv] <- NA # set failed convergence estimates to NA

## Mean error
error_N_m_y0 <- est_N_m_y0 - matrix(true_N_m_y0, nrow = 1000, 
                                    ncol = ncol(est_N_m_y0), byrow = F) 
error_N_f_y0 <- est_N_f_y0 - matrix(true_N_f_y0, nrow = 1000, 
                                    ncol = ncol(est_N_f_y0), byrow = F) 

# new: bias
bias_N_m_y0 <- error_N_m_y0 / matrix(true_N_m_y0, nrow = 1000, ncol = ncol(est_N_m_y0), byrow = F) * 100
bias_N_f_y0 <- error_N_f_y0 / matrix(true_N_f_y0, nrow = 1000, ncol = ncol(est_N_f_y0), byrow = F) * 100

## Mean estimates
median_est <- data.frame(N_f_y0 = colMedians(est_N_f_y0), 
                         N_m_y0 = colMedians(est_N_m_y0))
mean_est <- data.frame(N_f_y0 = colmeans(est_N_f_y0), 
                       N_m_y0 = colmeans(est_N_m_y0))

## Average error
median_error <- data.frame(N_f_y0 = colMedians(error_N_f_y0), 
                           N_m_y0 = colMedians(error_N_m_y0))

mean_error <- data.frame(N_f_y0 = colmeans(error_N_f_y0), 
                         N_m_y0 = colmeans(error_N_m_y0))


## Average bias
median_bias <- data.frame(N_f_y0 = colMedians(bias_N_f_y0), 
                          N_m_y0 = colMedians(bias_N_m_y0))

mean_bias <- data.frame(N_f_y0 = colmeans(bias_N_f_y0), 
                        N_m_y0 = colmeans(bias_N_m_y0))

## Mean absolute error
mean_abs_error <- data.frame(N_f_y0 = colmeans(abs(error_N_f_y0)), 
                             N_m_y0 = colmeans(abs(error_N_m_y0)))

## Mean squared error
mean_sqrd_error <- data.frame(N_f_y0 = colmeans((error_N_f_y0) ^ 2), 
                              N_m_y0 = colmeans((error_N_m_y0) ^ 2))

## standard deviation and variance of error
std_dev_error <- data.frame(N_f_y0 = apply(error_N_f_y0, 2, sd),
                            N_m_y0 = apply(error_N_m_y0, 2, sd))
variance_error <- std_dev_error ^ 2

## standard deviation and variance of estimates
std_dev <- data.frame(N_f_y0 = apply(est_N_f_y0, 2, sd),
                      N_m_y0 = apply(est_N_m_y0, 2, sd))
variance <- std_dev ^ 2

## Coefficient of variation (works as the support for these estimates is 
## strictly positive)
## This one is the correct one: std_dev_error comes from the error, mean_est 
## will standardise it
coeff_var <- std_dev_error / mean_est * 100

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
colnames(bias_N_f_y0_long) <- scen_names
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
colnames(bias_N_m_y0_long) <- scen_names
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
# bias_N_m_y0_long <- na.omit(bias_N_m_y0_long)

bias_N_y0_long <- rbind(bias_N_m_y0_long[, c(1, 3, 4, 2, 5)], 
                        bias_N_f_y0_long[, c(1, 3, 4, 2, 5)])
colnames(bias_N_y0_long)[1] <- c("Scenario")
bias_N_y0_long$Scenario <- factor(bias_N_y0_long$Scenario, levels = scen_names)

## A SINGLE PLOT
# y0_plot <- ggplot(data=bias_N_y0_long, aes(fill=Sex, x=Scenario, y=Abundance)) +
#   geom_boxplot(position="dodge", 
#                outlier.alpha = 0.1,
#                # outlier.shape = NA,
#                # draw_quantiles = c(0.5), # for violin plots
#                coef = 5,   # for boxplots, how many times the IQR values are included (defaults to 1.5)
#                alpha=0.5) +
#   stat_summary(mapping = aes(fill=Sex), position=position_dodge(0.75),
#                fun = mean, geom="point", shape=21, color="black") +
#   theme_minimal() +
#   theme(legend.position = "bottom") +
#   geom_hline(yintercept = 0, 
#              colour = "black", 
#              alpha = 0.5)+
#   ylab("Relative error in abundance") +
#   # coord_flip() +  # comment out for normal plot
#   # scale_x_discrete(limits = rev(levels(est_N_5_long$Scenario))) +
#   # coord_cartesian(ylim = c(-100, 500)) +
#   scale_fill_discrete(name = ""); y0_plot

## A 5x5 PLOT
y0_plot <- ggplot(data=bias_N_y0_long, aes(fill=Sex, x=1, y=Abundance)) +
  geom_boxplot(position="dodge", 
               outlier.alpha = 0.1,
               # outlier.shape = NA,
               # draw_quantiles = c(0.5), # for violin plots
               coef = 5,   # for boxplots, how many times the IQR values are included (defaults to 1.5)
               alpha=0.5) +
  stat_summary(mapping = aes(fill=Sex), position=position_dodge(0.75),
               fun = mean, geom="point", shape=21, color="black") +
  theme_minimal() +
  theme(legend.position = "bottom", 
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank()) +
  geom_hline(yintercept = 0, 
             colour = "black", 
             alpha = 0.5)+
  facet_wrap(~ Scenario , nrow = 5, labeller = label_parsed) +
  ylab("Relative error in abundance (%)") +
  # coord_flip() +  # comment out for normal plot
  # scale_x_discrete(limits = rev(levels(est_N_5_long$Scenario))) +
  # coord_cartesian(ylim = c(-100, 500)) +
  scale_fill_discrete(name = ""); y0_plot
# Export this plot in 800x1000

## -------------------------------
## Combine all plots into one plot
## -------------------------------

library(gridExtra)
# library(svglite)
lay <- rbind(c(1),
             # c(1),
             c(2),
             # c(2),
             c(3))

p <- grid.arrange(y0_plot + 
                    theme(legend.position = "none",
                          axis.text.x = element_blank()) +
                    ylab("Relative error\nabundance (100)") +
                    xlab(element_blank()) +
                    coord_cartesian(ylim=c(-100, 500)), # truncate without removing values outside
                  mean_y_5_plot +
                    theme(legend.position = "none",
                          axis.text.x = element_blank()) +
                    ylab("Relative error\nmean abundance (96-100)") +
                    xlab(element_blank()) +
                    coord_cartesian(ylim=c(-100, 500)), # truncate  without removing values outside
                  r_plot + 
                    # coord_cartesian(ylim=c(-100, 100))+ # for simple
                    coord_cartesian(ylim=c(-50, 50))+ # for complex
                    ylab("Relative error\nyearly growth rate"),
                  layout_matrix = lay); p

## now export it with width 1000 and height 750

## Alternative for the no growth version:
lay <- rbind(c(1))

p <- grid.arrange(y0_plot + 
                    # theme(legend.position = "none",
                    #       axis.text.x = element_blank()) +
                    ylab("Relative error\nabundance (100)") +
                    # xlab(element_blank()) +
                    coord_cartesian(ylim=c(-100, 500)),
                  layout_matrix = lay); p

## now export it with width 1000 and height 500

## Some custom code to combine the two plots, as there is no need to keep them
## separate when only N is estimated
lay <- rbind(c(1), c(2))

p <- grid.arrange(y0_plot_simple + 
                    theme(legend.position = "none",
                          axis.text.x = element_blank()) +
                    ylab("Simple species") +
                    xlab(element_blank()) +
                    coord_cartesian(ylim=c(-100, 500)),
                  y0_plot_complex + 
                    # theme(legend.position = "none",
                    #       axis.text.x = element_blank()) +
                    ylab("Complex species") +
                    coord_cartesian(ylim=c(-100, 500)),
                  layout_matrix = lay, 
                  left = "Relative error in adult abundance estimate"); p


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
    median_rel_error = median_bias[, column],
    mean_rel_error = mean_bias[, column],
    median_error = median_error[, column], 
    mean_error = mean_error[, column], 
    mean_abs_error = mean_abs_error[, column], 
    root_mean_sqrd_error = sqrt(mean_sqrd_error[, column]))
  # error_std_dev = std_dev_error[, column])
  
  df$scenario <- scen_names
  df <- na.omit(df[, c(ncol(df), 1:(ncol(df)-1))])
  return(df)
}

summary_N_f_y0 <- summarise_performance_metrics("N_f_y0", true_N_f_y0)

summary_N_m_y0 <- summarise_performance_metrics("N_m_y0", true_N_m_y0)

labels <- c("Scen.", 
            # "Mdn. Est.", 
            # "Mean Est.", 
            # "Std. Dev.", 
            # "CV", 
            "Mnd. Rel. Err.",
            "Mean Rel. Err.",
            "Mdn. Err.", 
            "Mean Err.", 
            "MAE", 
            "RMSE")

caption <- "Various performance metrics for the parameter [quantity], extracted from 1000 simulations of the [species] shark population. The columns are, from left to right: scenario, median estimate, mean estimate, standard deviation, coefficient of variation (\\%), median error, mean error, mean absolute error, and mean squared error. Scenarios 1-1, 1-2, 1-3, 2-1, 2-2, and 3-1 were not included as (most of) the simulations did not fit successfully. Scenario 3-3 uses the correct measurement error (2.89) and growth curve specification."
kable(summary_N_f_y0, booktabs = T, format = "latex", digits = 2,
      col.names = labels, row.names = F, linesep = "", caption = caption)

kable(summary_N_m_y0, booktabs = T, format = "latex", digits = 2,
      col.names = labels, row.names = F, linesep = "", caption = caption)

## =============================================================================
## 4. CREATE STANDARD DEVIATION AND CV TABLES FOR APPENDIX D
## =============================================================================

## Create tables for CV and SD similar to "estimate_variance_mle.R" 

if (NO_GROWTH) {
  ## Run section 1 and 2 above first for simple, then run lines below
  sd_simple <- std_dev_error[conv, c(2, 1)]
  cv_simple <- coeff_var[conv, c(2, 1)]
  
  ## Now run section 1 and 2 for the complex species, then run lines below
  sd_complex <- std_dev_error[conv, c(2, 1)]
  cv_complex <- coeff_var[conv, c(2, 1)]
  
  ## Now create tables using lines below
  sd_df <- cbind(data.frame(scen_names[conv]), sd_simple, sd_complex)
  cv_df <- cbind(data.frame(scen_names[conv]), cv_simple, cv_complex)
  
  colnames(sd_df) <- c("scenario", 
                       "simple_N_y0_m" , "simple_N_y0_f",
                       "complex_N_y0_m" , "complex_N_y0_f")
  colnames(cv_df) <- c("scenario", 
                       "simple_N_y0_m" , "simple_N_y0_f",
                       "complex_N_y0_m" , "complex_N_y0_f") 
  
  ## Create tables for latex
  library(kableExtra)
  
  labels <- c("Scenario", 
              "$N^A_{male, 100}$", 
              "$N^A_{female, 100}$", 
              "$N^A_{male, 100}$", 
              "$N^A_{female, 100}$")
  
  caption <- "" 
  kable(sd_df, booktabs = T, 
        format = "latex", 
        digits = 2,
        col.names = labels, 
        row.names = F, linesep = "", caption = caption) %>%
    add_header_above( c(" " = 1, "simple species" = 4, "Complex species" = 4))
  
  caption <- "" 
  kable(cv_df, booktabs = T, 
        format = "latex", 
        digits = 2,
        col.names = labels, 
        row.names = F, linesep = "", caption = caption) %>% 
    add_header_above( c(" " = 1, "simple species" = 4, "Complex species" = 4))
} else {
  ## Run section 1 and 2 above first for simple, then run lines below
  sd_simple <- std_dev_error[conv, c(6, 5, 2, 1)]
  cv_simple <- coeff_var[conv, c(6, 5, 2, 1)]
  
  ## Now run section 1 and 2 for the complex species, then run lines below
  sd_complex <- std_dev_error[conv, c(6, 5, 2, 1)]
  cv_complex <- coeff_var[conv, c(6, 5, 2, 1)]
  
  ## Now create tables using lines below
  sd_df <- cbind(data.frame(scen_names[conv]), sd_simple, sd_complex)
  cv_df <- cbind(data.frame(scen_names[conv]), cv_simple, cv_complex)
  
  colnames(sd_df) <- c("scenario", 
                       "simple_r_m", "simple_r_f" , "simple_N_y0_m" , "simple_N_y0_f",
                       "complex_r_m", "complex_r_f" , "complex_N_y0_m" , "complex_N_y0_f")
  colnames(cv_df) <- c("scenario", 
                       "simple_r_m", "simple_r_f" , "simple_N_y0_m" , "simple_N_y0_f",
                       "complex_r_m", "complex_r_f" , "complex_N_y0_m" , "complex_N_y0_f") 
  
  ## Create tables for latex
  library(kableExtra)
  
  labels <- c("Scenario", 
              "$r_{male}$", 
              "$r_{female}$", 
              "$N^A_{male, 2014}$", 
              "$N^A_{female, 2014}$", 
              "$r_{male}$", 
              "$r_{female}$",
              "$N^A_{male, 2014}$", 
              "$N^A_{female, 2014}$")
  
  caption <- "" 
  kable(sd_df, booktabs = T, 
        format = "latex", 
        digits = 2,
        col.names = labels, 
        row.names = F, linesep = "", caption = caption) %>%
    add_header_above( c(" " = 1, "simple species" = 4, "Complex species" = 4))
  
  caption <- "" 
  kable(cv_df, booktabs = T, 
        format = "latex", 
        digits = 2,
        col.names = labels, 
        row.names = F, linesep = "", caption = caption) %>% 
    add_header_above( c(" " = 1, "simple species" = 4, "Complex species" = 4))
}


