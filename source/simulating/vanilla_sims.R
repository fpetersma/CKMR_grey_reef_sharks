## ========================================================================== ##
## vanilla_sims.R
## 
## In this script I will simulate the vanilla population data using a mix of 
## information on life  history trait estimates from the GBR, Hawaii and Palmyra 
## study. This is a very basic population and details can be found in the 
## Overleaf document.
## ========================================================================== ##

## Load packages and source custom functions ===================================
library(fishSim)
library(parallel)
library(pbapply)
library(ids)
source("source/simulating/custom_functions_fishSim.R")

## Initialise parameters =======================================================
## Set life history parameters
first_breed_male <- 10          # applies only to males
first_breed_female <- 10        # applies only to females
first_litter <- 2               # first litter 
batch_size <- 2                 # all possible values for random draw
max_clutch <- Inf               # maximum litter sizes
max_age <- 19                   # maximum age
mort_rate <- 0.153              # flat mortality rate (0.131 gives stable population)
surv_rate <- 1 - mort_rate      # survival is 1 minus mortality
force_Y1 <- NA                  # first-year mortality (not always used)
female_curve <- c(rep(0, 10), 
                  rep(1, 90))   # female maturity curve for categories 0:maxAge
male_curve <- c(rep(0, 10), 
                rep(1, 90))     # male maturity curve for categories 0:maxAge
mate_type <- "ageSex"           # type of maturity/mating structure
mort_type <- "flat"             # type of mortality structure
fecundity_dist <- "uniform"     # uniform distribution (custom)
sex_ratio <- c(0.5, 0.5)        # the offspring male/female sex ratio
single_paternity <- FALSE        # is there a single father to a litter?

## Set simulation and sampling parameters
years <- 101:140             # number of years to run simulation
sampling_years <- c((max(years) - 1):max(years))     # years in which sampling occurs
sample_size <- c(rep(0, sampling_years[1] - 1), 
                 rep(375, 2))   # number of sampled individuals in each year
lethal_sampling <- FALSE        # is sampling lethal?
retrospective_sampling <- TRUE

## Store simdata sets
n_cores <- 30
cl <- makeCluster(n_cores)
clusterExport(cl, c(ls(), "makeFounders", "mort", "birthdays", "uuid"))
simulated_data_sets <- pblapply(1:1000, function(i) {
  
  ## Set a seed for reproducibility
  set.seed(170593 + 121 * i)
  
  ## Create initial population =================================================
  ## Set initial population in year 0
  indiv <-  makeFounders(
    pop = 8500,                   # starting pop. size 8500 from literature
    osr = sex_ratio,
    stocks = c(1), # a single stock
    maxAge = max_age,             #
    survCurv = surv_rate ^ (1:(max_age )) / sum(surv_rate ^ (1:(max_age)))
  )
  
  # ## Add marker to pregnant females
  # indiv <- addPregnancy(
  #   indiv = indssssiv,
  #   matingAges = c(12, 14, 16, 18))
  # 
  
  ## Start the simulation ========================================================
  for (y in years) {
    # cat("Starting simulation of year:", y, "...")
    
    ## 1. Mating/birth
    indiv <- mateOrBirth(
      indiv = indiv,
      batchSize = batch_size,
      fecundityDist = fecundity_dist, # uniform is custom
      osr = sex_ratio,
      year = y,
      firstBreedFemale = first_breed_female,
      firstBreedMale = first_breed_male,
      firstLitter = first_litter, # firstLitter is custom
      type = mate_type, 
      maxClutch = max_clutch,
      singlePaternity = single_paternity,  
      maleCurve = male_curve,   
      femaleCurve = female_curve,
      no_gestation = TRUE,
    )
    
    # ## 2. Sampling
    # if (y %in% sampling_years && !retrospective_sampling) {
    #   indiv <- captureOnlyFirst(
    #     indiv = indiv, 
    #     n = sample_size[y], 
    #     year = y, 
    #     fatal = lethal_sampling
    #   )
    # }
    
    ## 3. Survival
    indiv <- mort(
      indiv = indiv, 
      year = y, 
      type = mort_type, 
      mortRate = mort_rate,
      maxAge = max_age - 1
    )
    
    ## 4. Age incrementation
    indiv <- birthdays(indiv)
  
    # cat(" Cycle completed!\n")
  }
  # cat(paste0("Simulation ", i, " completed!\n"))
  # simulated_data_sets[[i]] <- indiv
  return(indiv)
}, cl = cl)

## How many individuals are still alive in 'indiv'?
hist(sapply(simulated_data_sets, function(x) {nrow(x[is.na(x$DeathY), ])}), xlab = "indivs")

## Set sampling parameters
n_sample_year <- 2
sampling_years <- c((max(years) - n_sample_year + 1):max(years))     # years in which sampling occurs
sample_size <- c(rep(0, sampling_years[1] - 1), 
                 rep(375, n_sample_year))   # number of sampled individuals in each year
lethal_sampling <- FALSE        # is sampling lethal?
retrospective_sampling <- TRUE

## Remove previous sampling if required
simulated_data_sets <- lapply(simulated_data_sets, function(indiv) {
  indiv$SampY <- NA
  return(indiv)
})
all(is.na(simulated_data_sets[[1]]$SampY)) # should be TRUE before continuing

## Add sampling
if (retrospective_sampling) {
  simulated_data_sets <- lapply(simulated_data_sets, function(indiv) {
    for (year in sampling_years) {
      n <- sample_size[year]
      indiv <- retroCapture(indiv, n = n, year = year, fatal = lethal_sampling)
    }
    return(indiv)
  })
}
all(is.na(simulated_data_sets[[1]]$SampY))  # should be FALSE
unique(simulated_data_sets[[1]]$SampY)      # check if this seems correct
sum(!is.na(simulated_data_sets[[1]]$SampY)) # seem correct as well?

save.image(file = "data/vanilla_sample_years_139-140_sample_size_375/1000_sims_mix.RData")


## Create summary stats
par(mfrow = c(3, 1))
## Extracting simulated abundances
N_true <- t(sapply(simulated_data_sets, function(x) {
  mature_females <- sum(is.na(x$DeathY) & x$Sex == "F" & x$AgeLast >= 10)
  mature_males <- sum(is.na(x$DeathY) & x$Sex == "M" & x$AgeLast >= 10)
  return(c(N_m = mature_males, N_f = mature_females))
}))

r_true <- sapply(simulated_data_sets, function(x) {
  alive <- sum(is.na(x$DeathY))
  return((alive / 8500) ^ (1 / 40))
})
summary(r_true)

hist(N_true[, 1], main = "", xlab = "Number of mature males"); abline(v = c(mean(N_true[, 1]), median(N_true[, 1])), col = "red", lty = c(1, 2))

hist(N_true[, 2], main = "",  xlab = "Number of mature females"); abline(v = c(mean(N_true[, 2]), median(N_true[, 2])), col = "red", lty = c(1, 2))

hist(r_true, main = "", xlab = "Yearly growth rate"); abline(v = c(mean(r_true), median(r_true)), col = "red", lty = c(1, 2))
