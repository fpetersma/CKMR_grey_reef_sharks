################################################################################
##  Extract POP from a indiv object. This script is mainly to run stuff       ##
##  on the many cores of the bluewhale server.                                ##
##
##
################################################################################

## Load packages and source custom functions ===================================
library(fishSim)
library(ids)

## Load data
load(file = "data/data_for_testing_estimator.RData")

## Looking up relationship between captured pairs
pairs <- findRelativesPar(indiv = indiv, 
                          sampled = TRUE, 
                          nCores = 30)
POPs <- pairs[pairs$OneTwo == 1, ] ## Parent-Offspring pairs

## Save data
save.image("data/data_for_testing_estimator.RData")

