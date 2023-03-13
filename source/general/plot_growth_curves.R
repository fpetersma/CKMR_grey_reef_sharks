## ========================================================================== ##
##  Name: plot_growth_curves.R                                                ##
##                                                                            ##
##  Description: In this script we plot the growth curves for the vanillus    ##
##    and gestatii population. It is a combination of the scripts             ##
##    fit_and_visualise_scenarios.R and explore_vbgf_bias_uncertainty.R.      ##
## ========================================================================== ##

## ------------------------------------
## Set the values for the growth curves
## ------------------------------------

a_0 <- -8.27                # theoretical age at length zero              
k <- 0.0554                 # growth rate
l_inf <- 163                # asymptotic length

## vbgf comes from l.100 onwards in 'explore_vbgf_bias_uncertainty.R'
step_l_inf <- 0.05 * l_inf  # 8.15 is 5% of 163
## only do this for 5 scenarios, as we are only using those in the manuscript
l_inf_options <- seq(from = l_inf - step_l_inf * 2, 
                  to = l_inf + step_l_inf * 2, 
                  by = step_l_inf)
vbgfs <- matrix(c(CKMRcpp::a_0_j(l_inf, l_inf_options, k, a_0), 
                  l_inf_options), ncol = 2)

## ----------------------------------
## Create a matrix of values and plot
ages <- 0:63
a_0_l_inf_matrix <- matrix(NA, nrow = length(ages), ncol = nrow(vbgfs))
for (i in 1:nrow(vbgfs)) {
  a_0 <- vbgfs[i, 1]
  l_inf <- vbgfs[i, 2]
  
  a_0_l_inf_matrix[, i] <- CKMRcpp::vbgf(a = ages, l_inf = l_inf, a_0 = a_0, k = 0.1)
}

matplot(ages, a_0_l_inf_matrix, type = "l", 
        col = c(rep("black", 2), "red", rep("black", 2)),
        lty = 2, lwd = c(1,1,2,1,1),
        ylim = c(0, 200))

## Can I plot this a little nicer?
library(ggplot2)
a_0_l_inf_long <- reshape2::melt(a_0_l_inf_matrix)
colnames(a_0_l_inf_long) <- c("Age", "Case", "Length")
a_0_l_inf_long$Age <- a_0_l_inf_long$Age - 1
a_0_l_inf_long$Case <- as.factor(a_0_l_inf_long$Case)

## Add labels
labels <- c("-10%", "-5%", "0%", "+5%", "+10%")
a_0_l_inf_long$label <- ifelse(a_0_l_inf_long$Age == max(a_0_l_inf_long$Age), 
                               labels[as.integer(a_0_l_inf_long$Case)], 
                               NA_character_)

growth_curve_plot <- ggplot(a_0_l_inf_long, aes(x = Age, y = Length, colour = Case, linetype = Case)) +
  geom_line(#linewidth = 1, 
            alpha = 1, 
            show.legend = F) +                        
  scale_color_manual(values=c("black", "black", "red", "black", "black")) +     # make truth red and wrong cases black
  scale_linetype_manual(values=c(2,2,1,2,2)) +
  theme_bw() +
  xlim(0, 70) +
  ylab("Length (cm)") +
  xlab("Age (years)") +
  ggrepel::geom_label_repel(aes(label = label), ## Add labels     
                   nudge_x = 3, #c(-30, -20, -10, 0, 10), 
                   na.rm = TRUE, show.legend = FALSE,
                   label.padding = 0.25,
                   box.padding = 0.0,
                   alpha = 1,
                   label.size = 0.5
                   ) 

## --------------------------------------------------
## Now create a plot of the various measurement error
## --------------------------------------------------
## Enter the tidyverse
library(tidyverse)

## Set true error and deviations
true_error <- 2.89                        # what is the true measurement error

errors <- c(1e-8,                         # 1e-8 instead of 0 to avoid overflow
            true_error*1/3,               # 1/3 of true error
            true_error*2/3,               # 2/3 of true error
            true_error,                   # true error
            true_error*4/3,               # 4/3 of true error
            true_error*5/3,               # 5/3 of true error
            true_error*6/3)               # twice the true error

errors <- errors[2:6]                     # only keep middle 5 scenarios

## Start on the plot
x_values <- seq(from=-20, to=20, length.out = 10000)

labels <- c("-67%", "-33%", "0%", "+33%", "+67%")

error_df <- data.frame(x = rep(x_values, times = length(errors)), 
                       e = rep(errors, each = length(x_values)),
                       case = factor(rep(1:5, each=length(x_values)),  labels=labels)) %>% 
  mutate(y = dnorm(x, 0, e))

error_plot <- ggplot(error_df) +
  geom_line(aes(x = x, y = y, group = case, 
                colour = case,
                # linetype = case
                ),  
            linewidth = 1, 
            alpha = 0.5) +
   theme_bw() +
  # scale_color_manual(values=c("black", "black", "red", "black", "black")) +     # make truth red and wrong cases black
  # scale_linetype_manual(values=c(2,2,1,2,2)) +
  labs(colour = "Deviation\ntrue error", 
       x = "Measurement error on length",
       y = "Probability density")

## :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
## Combine both plots into one plot
## :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
library(gridExtra)

p <- grid.arrange(growth_curve_plot,
                  error_plot,
                  nrow=1, ncol=2); p

## Save as svg with width 1000 and height 400