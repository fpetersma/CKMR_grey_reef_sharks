################################################################################
##                                                                            ##
##  Explore the different Von Bertalanffy growth functions.                   ##
##                                                                            ##
##  Misspecified uncertainty in the length measurements is relatively easy,   ##
##  but the a 'shift' in the curve is a bit harder. Let's see if I can        ##
##  achieve what I want to achieve.                                           ##
##                                                                            ##
################################################################################

## :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
## Load packages
## :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
library(tidyverse)
library(CKMRcpp)

## :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
## Create a default vbgf function with plotting capabilities
## :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
ages <- 0:19
lengths <- vbgf(ages) # default for vbgf() are t_0=-3.5, k = 0.1, l_inf=175
plot(ages, lengths, type = "l", ylim = c(0, 200))

## Try to do shift it
t_0_options <- 0:-5
k_options <- seq(from=0.02, to=0.18, by=0.02)
l_inf_options <- seq(from=160, to=200, by=5)

## t_0
t_0_matrix <- matrix(NA, nrow = 20, ncol = length(t_0_options))
for (i in seq_along(t_0_options)) {
  t_0 <- t_0_options[i]
  t_0_matrix[, i] <- vbgf(a = ages, t_0 = t_0)
}

matplot(ages, t_0_matrix, type = "l")

## k
k_matrix <- matrix(NA, nrow = 20, ncol = length(k_options))
for (i in seq_along(k_options)) {
  k <- k_options[i]
  k_matrix[, i] <- vbgf(a = ages, k = k)
}

matplot(ages, k_matrix, type = "l")

## l_inf
l_inf_matrix <- matrix(NA, nrow = 20, ncol = length(l_inf_options))
for (i in seq_along(l_inf_options)) {
  l_inf <- l_inf_options[i]
  l_inf_matrix[, i] <- vbgf(a = ages, l_inf = l_inf)
}

matplot(ages, l_inf_matrix, type = "l")

## It seems that mainly t_0 and to a lesser extent l_inf are important
## equation 0.5 step to 5 works quite well, but 0.5 to 6.5 is better
t_0_l_inf_options <- matrix(c(c(-2, -2.5, -3, -3.5, -4, -4.5, -5), 
                              seq(from=155.5, to=194.5, by=6.5)), ncol = 2)
t_0_l_inf_matrix <- matrix(NA, nrow = 20, ncol = nrow(t_0_l_inf_options))
for (i in 1:nrow(t_0_l_inf_options)) {
  t_0 <- t_0_l_inf_options[i, 1]
  l_inf <- t_0_l_inf_options[i, 2]
  
  t_0_l_inf_matrix[, i] <- vbgf(a = ages, l_inf = l_inf, t_0 = t_0, k = 0.15)
}

matplot(ages, t_0_l_inf_matrix, type = "l", 
        col = c(rep("black", 3), "red", rep("black", 3)),
        lty = 2, lwd = c(1,1,1,2,1,1,1),
        ylim = c(0, 200))

## Are the changes in the mean equal?
diff(colMeans(t_0_l_inf_matrix)) # yeah, pretty much!
sd(diff(colMeans(t_0_l_inf_matrix))) # 0.01015705
mean(diff(colMeans(t_0_l_inf_matrix))) # 7.8258 = 4.5%

## Can I achieve increments of 5% = 8.75?
t_0_l_inf_options <- matrix(c(seq(from=-2, to=-5, by=-0.5), 
                              seq(from=151.4, to=199.1, by=7.7)), ncol = 2)
t_0_l_inf_matrix <- matrix(NA, nrow = 20, ncol = nrow(t_0_l_inf_options))
for (i in 1:nrow(t_0_l_inf_options)) {
  t_0 <- t_0_l_inf_options[i, 1]
  l_inf <- t_0_l_inf_options[i, 2]
  
  t_0_l_inf_matrix[, i] <- vbgf(a = ages, l_inf = l_inf, t_0 = t_0, k = 0.15)
}

matplot(ages, t_0_l_inf_matrix, type = "l", 
        col = c(rep("black", 3), "red", rep("black", 3)),
        lty = 2, lwd = c(1,1,1,2,1,1,1),
        ylim = c(0, 200))

## Are the changes in the mean equal?
diff(colMeans(t_0_l_inf_matrix)) # yeah, pretty much, althought not as good as before.
sd(diff(colMeans(t_0_l_inf_matrix))) # 0.06975593
mean(diff(colMeans(t_0_l_inf_matrix))) # 8.769889
diff(colMeans(t_0_l_inf_matrix)) / 175 # Pretty damn close to 5%.
# closer to 6.3% if looked at the mean of the colMeans, which is probably better
diff(colMeans(t_0_l_inf_matrix)) / mean(colMeans(t_0_l_inf_matrix)) 

## :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
## Simulations are at k = 0.1, so can I find the right balance there?
## :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

## Write a fucntion to calculate a_0_j (also added to the package CKMRcpp)
a_0_j <- function(l_inf_i, l_inf_j, k, a_0_i) {
  d <- l_inf_i - l_inf_j
  a_0 <- log(1 - (CKMRcpp::vbgf(0, a_0_i, k, l_inf_i) - d) / l_inf_j) / k
  return(a_0)
}
a0_options <- seq(from=175-8.75*3, to=175+8.75*3, by=8.75)
t_0_l_inf_options <- matrix(c(a_0(175, a0_options, 0.1, -3.5), 
                              a0_options), 
                              ncol = 2)
t_0_l_inf_matrix <- matrix(NA, nrow = 20, ncol = nrow(t_0_l_inf_options))
for (i in 1:nrow(t_0_l_inf_options)) {
  t_0 <- t_0_l_inf_options[i, 1]
  l_inf <- t_0_l_inf_options[i, 2]
  
  t_0_l_inf_matrix[, i] <- vbgf(a = ages, l_inf = l_inf, t_0 = t_0, k = 0.1)
}

matplot(ages, t_0_l_inf_matrix, type = "l", 
        col = c(rep("black", 3), "red", rep("black", 3)),
        lty = 2, lwd = c(1,1,1,2,1,1,1),
        ylim = c(0, 200))

## Are the changes in the mean equal?
diff(colMeans(t_0_l_inf_matrix)) # yeah, pretty much, although not as good as before.
sd(diff(colMeans(t_0_l_inf_matrix))) # 0.06975593
mean(diff(colMeans(t_0_l_inf_matrix))) # 8.769889
diff(colMeans(t_0_l_inf_matrix)) / mean(colMeans(t_0_l_inf_matrix)) # Pretty damn close to 5%. 
