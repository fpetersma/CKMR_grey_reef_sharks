## ========================================================================== ##
## vanilla_with_gestation_and_variable_reproduction_sims.R
## 
## In this script I will simulate the vanilla population data using a mix of 
## information on life  history trait estimates from the GBR, Hawaii and Palmyra 
## study. This is a very basic population and details can be found in the 
## Overleaf document.
##
## Three things have been added to "vanilla_sims.R":
##  - Reproduction for sharks is rarely as predictable as for (big) mammals. In
##    our scenario, we are going to assume that litter size is random from the 
##    discretised U(1, 4), thus ERRO = 2.5. 
##  - Females have a gestation of one year, which means that the male population
##    of the year BEFORE the birth of the offspring should be considered, and 
##    potential fathers only have to survive until y_offspring - 1. 
##    Females, on the other hand, have to have been alive in y_offspring AND
##    y_offspring - 1. Not sure yet how to implement this in the likelihood.
##  - Females mature later than males, at 12 years of age instead of 10. 
##
## I will likely first implement the variable reproduction in the likelihood, as 
## that is relatively straightforward; the gestation is a bit trickier. Last, I 
## will add the difference in maturity age.
## ========================================================================== ##

## Load packages and source custom functions ===================================
library(fishSim)
library(parallel)
library(pbapply)
library(ids)
library(CKMRcpp)

## Initialise parameters =======================================================
## Set life history parameters
first_breed_male <- 10          # applies only to males
first_breed_female <- 10        # applies only to females
first_litter <- 1:4             # first litter 
batch_size <- 1:4               # all possible values for random draw
max_clutch <- Inf               # maximum litter sizes
max_age <- 19                   # maximum age
mort_rate <- 0.1675             # flat mortality rate (0.1675 gives stable population)
surv_rate <- 1 - mort_rate      # survival is 1 minus mortality
force_Y1 <- NA                  # first-year mortality (not always used)
female_curve <- c(rep(0, first_breed_female), 
                  rep(1, 90))   # female maturity curve for categories 0:maxAge
male_curve <- c(rep(0, first_breed_male), 
                rep(1, 90))     # male maturity curve for categories 0:maxAge
mate_type <- "ageSex"           # type of maturity/mating structure
mort_type <- "age"             # type of mortality structure
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
rm(simulated_data_sets)
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
  #   indiv = indiv,
  #   matingAges = c(12, 14, 16, 18))
  # 
  
  ## Start the simulation ========================================================
  for (y in years) {
    # cat("Starting simulation of year:", y, "...")
    
    ## 1. Mating/birth
    indiv <- CKMRcpp::mateOrBirth(
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
      no_gestation = TRUE
    )
    
    ## 2. Survival
    indiv <- mort(
      indiv = indiv, 
      year = y, 
      type = mort_type, 
      # mortRate = mort_rate,
      ageMort = rep(mort_rate, max_age + 1),
      maxAge = max_age - 1 
    )
    
    ## 3s. Age incrementation
    indiv <- birthdays(indiv)
    
    # cat(" Cycle completed!\n")
  }
  # cat(paste0("Simulation ", i, " completed!\n"))
  # simulated_data_sets[[i]] <- indiv
  return(indiv)
}, cl = cl); stopCluster(cl);

# load("data/vanilla_sample_years_139-140_sample_size_375/1000_sims_mix.RData")

## How many individuals are still alive in 'indiv'?
hist(sapply(simulated_data_sets, function(x) {nrow(x[is.na(x$DeathY), ])}), xlab = "indivs")
summary(sapply(simulated_data_sets, function(x) {nrow(x[is.na(x$DeathY), ])}), xlab = "indivs")

## Set sampling parameters
n_samples <- 750
n_sample_year <- 5
sampling_years <- c((max(years) - n_sample_year + 1):max(years))     # years in which sampling occurs
sample_size <- c(rep(0, sampling_years[1] - 1), 
                 rep(n_samples / n_sample_year, n_sample_year))   # number of sampled individuals in each year
lethal_sampling <- FALSE        # is sampling lethal?
retrospective_sampling <- TRUE

## Remove previous sampling if required
simulated_data_sets <- lapply(simulated_data_sets, function(indiv) {
  indiv$SampY <- NA
  return(indiv[, 1:9]) # only return first 9 columns (before sampling columns)
})
all(is.na(simulated_data_sets[[1]]$SampY)) # should be TRUE before continuing
ncol(simulated_data_sets[[1]]) == 9 # should also be TRUE

## Add sampling
if (retrospective_sampling) {
  simulated_data_sets <- pblapply(simulated_data_sets, function(indiv) {
    for (year in sampling_years) {
      n <- sample_size[year]
      indiv <- retroCapture2(indiv, n = n, year = year, fatal = lethal_sampling)
    }
    indiv$no_samples <- rowSums(!is.na(indiv[, 10:ncol(indiv)]))
    return(indiv)
  })
}

all(is.na(simulated_data_sets[[1]]$SampY))  # should be FALSE
unique(simulated_data_sets[[1]]$SampY)      # check if this seems correct
sum(!is.na(simulated_data_sets[[1]]$SampY)) # seem correct as well?

save.image(file = "data/vanilla_variable_reproduction_sample_years_136-140_sample_size_150/1000_sims_mix.RData")

## Create summary stats
par(mfrow = c(3, 1))
## Extracting simulated abundances
N_true <- t(sapply(simulated_data_sets, function(x) {
  mature_females <- sum(is.na(x$DeathY) & x$Sex == "F" & x$AgeLast >= 10)
  mature_males <- sum(is.na(x$DeathY) & x$Sex == "M" & x$AgeLast >= 10)
  return(c(N_m = mature_males, N_f = mature_females))
}))

# r_true <- sapply(simulated_data_sets, function(x) {
#   alive_120 <- sum(!is.na(x$DeathY) & x$DeathY >= 120 & x$BirthY <= 111)
#   alive_140 <- sum(is.na(x$DeathY) & x$BirthY <= 131)
#   return((alive_140 / alive_120) ^ (1 / 21))
# })
## Derive growht rate based on animals alive AT THE END OF YEAR
r_true <- sapply(simulated_data_sets, function(x) {
  alive_101 <- sum(!is.na(x$DeathY) & x$DeathY > 101 & x$BirthY <= 101)
  alive_140 <- sum(is.na(x$DeathY) | x$DeathY > 140 & x$BirthY <= 140)
  return((alive_140 / alive_101) ^ (1 / 40))
})
summary(data.frame(N_true, r_true))

hist(N_true[, 1], main = "", xlab = "Number of mature males"); abline(v = c(mean(N_true[, 1]), median(N_true[, 1])), col = "red", lty = c(1, 2))

hist(N_true[, 2], main = "",  xlab = "Number of mature females"); abline(v = c(mean(N_true[, 2]), median(N_true[, 2])), col = "red", lty = c(1, 2))

hist(r_true, main = "", xlab = "Yearly growth rate"); abline(v = c(mean(r_true), median(r_true)), col = "red", lty = c(1, 2))
