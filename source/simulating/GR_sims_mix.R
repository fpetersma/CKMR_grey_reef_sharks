## ========================================================================== ##
## GR_sims_mix.R
## 
## In this script I will simulate GR population data using a mix of information 
## on life  history trait estimates from the GBR, Hawaii and Palmyra study. 
## 
## Here, I will simulate accurately, where sharks mate in year t and give birth
## in year t+1; in previous scripts (GR_sims_Hawaii.R and GR_sims_Palmyra.R)
## this was incorrect, where sharks mated and birthed in the same  year. 
## ========================================================================== ##

## Load packages and source custom functions ===================================
library(fishSim)
library(ids)
source("source/simulating/custom_functions_fishSim.R")

## Initialise parameters =======================================================
## Set life history parameters
gestation <- 1                  # gestation years of a mother
first_breed_male <- 10          # applies only to males
first_breed_female <- 12        # applies only to females
first_litter <- 3:6             # first litter 
batch_size <- 3:6               # all possible values for random draw
max_clutch <- Inf               # maximum litter sizes
max_age <- 19                   # maximum age
mort_rate <- 0.131              # flat mortality rate (0.131 gives stable population)
surv_rate <- 1 - mort_rate      # survival is 1 minus mortality
force_Y1 <- NA                  # first-year mortality (not always used)
female_curve <- c(rep(0, 12), 
                  rep(1, 88))   # female maturity curve for categories 0:maxAge
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
                 rep(350, 2))   # number of sampled individuals in each year
lethal_sampling <- FALSE        # is sampling lethal?
retrospective_sampling <- TRUE

## Store simdata sets
simulated_data_sets <- list()

for (i in 1:100) {
  
  ## Create initial population =================================================
  ## Set initial population in year 0
  indiv <-  makeFounders(
    pop = 8500,                   # starting pop. size 8500 from literature
    osr = sex_ratio,
    stocks = c(1),                # a single stock
    maxAge = max_age + 1,             #
    survCurv = surv_rate ^ (1:(max_age + 1)) / sum(surv_rate ^ (1:(max_age + 1)))
  )
  
  ## Remove 1 year from the ages in indiv to match my system 
  indiv$AgeLast <- indiv$AgeLast - 1
  
  ## Add marker to pregnant females
  indiv <- addPregnancy(
    indiv = indiv,
    matingAges = c(12, 14, 16, 18))
  
  
  ## Start the simulation ========================================================
  for (y in years) {
    # cat("Starting simulation of year:", y, "...")
    
    ## 1. Sampling
    if (y %in% sampling_years && !retrospective_sampling) {
      indiv <- captureOnlyFirst(
        indiv = indiv, 
        n = sample_size[y], 
        year = y, 
        fatal = lethal_sampling
      )
    }
    
    ## 2. Survival
    indiv <- mort(
      indiv = indiv, 
      year = y, 
      type = mort_type, 
      mortRate = mort_rate,
      maxAge = max_age - 1
    )
    
    ## 3. Age incrementation
    indiv <- birthdays(indiv)
    
    ## 4. Mating/birth
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
      femaleCurve = female_curve
    )
  
    # cat(" Cycle completed!\n")
  }
  cat(paste0("Simulation ", i, " completed!\n"))
  simulated_data_sets[[i]] <- indiv
}

## Set sampling parameters
n_sample_year <- 2
sampling_years <- c((max(years) - n_sample_year + 1):max(years))     # years in which sampling occurs
sample_size <- c(rep(0, sampling_years[1] - 1), 
                 rep(750, n_sample_year))   # number of sampled individuals in each year
lethal_sampling <- FALSE        # is sampling lethal?
retrospective_sampling <- TRUE

## Remove previous sampling if required
simulated_data_sets <- lapply(simulated_data_sets, function(indiv) {
  indiv$SampY <- NA
  return(indiv)
})
all(is.na(simulated_data_sets[[1]]$SampY)) # should be TRUE

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

save.image(file = "data/vanilla_sample_years_139-140_sample_size_375/100_sims_mix.RData")


## =============================================================================
## Create initial population (single population, outdated code)
## =============================================================================
## Set initial population in year 0
indiv <-  makeFounders(
  pop = 8500,                   # starting pop. size 8500 from literature
  osr = sex_ratio,
  stocks = c(1),                # a single stock
  maxAge = max_age + 1,             #
  survCurv = surv_rate ^ (1:(max_age + 1)) / sum(surv_rate ^ (1:(max_age + 1)))
)

## Remove 1 year from the ages in indiv to match my system
indiv$AgeLast <- indiv$AgeLast - 1

## Add marker to pregnant females
indiv <- addPregnancy(
  indiv = indiv,
  matingAges = c(12, 14, 16, 18))


## Start the simulation ========================================================
for (y in years) {
  cat("Starting simulation of year:", y, "...")
  
  ## 1. Survival
  indiv <- mort(
    indiv = indiv, 
    year = y, 
    type = mort_type, 
    mortRate = mort_rate,
    maxAge = max_age - 1
  )
  
  ## 2. Age incrementation
  indiv <- birthdays(indiv)
  
  ## 3. Mating/birth
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
    femaleCurve = female_curve
  )

  ## 4. Sampling
  if (y %in% sampling_years) {
    indiv <- capture(
      indiv = indiv, 
      n = sample_size[y], 
      year = y, 
      fatal = lethal_sampling
    )
  }
  cat(" Cycle completed!\n")
}

## How many individuals are still alive in 'indiv'?
nrow(indiv[is.na(indiv$DeathY), ])
# View(indiv[is.na(indiv$DeathY), ])

## Save data for fitting
save.image(file = "data/data_for_testing_estimator.RData")

## Looking up relationship between captured pairs, or run on Bluewhale
pairs <- findRelativesPar(indiv = indiv, 
                          sampled = TRUE, 
                          nCores = 6)
POPs <- pairs[pairs$OneTwo == 1,] ## Parent-Offspring pairs



## Start 100  simulations ======================================================
indiv_alive <- rep(NA, 10)
for (sim_i in 1:10) {
  ## Set initial population in year 0
  indiv <-  makeFounders(
    pop = 8500,                   # starting pop. size 8500 from literature
    osr = sex_ratio,
    stocks = c(1),                # a single stock
    maxAge = max_age,             #
    survCurv = surv_rate ^ (1:max_age) / sum(surv_rate ^ (1:max_age))
  )
  
  ## Add marker to pregnant females
  indiv <- addPregnancy(
    indiv = indiv,
    matingAges = c(12, 14, 16, 18)
  )
  
  
  for (y in years) {
    cat("Starting simulation of year:", y, "...")
    # 1. Mating
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
      femaleCurve = female_curve
    )
    
    ## 2. Mortality
    indiv <- mort(
      indiv = indiv, 
      year = y, 
      type = mort_type, 
      mortRate = mort_rate,
      maxAge = max_age
    )
    
    ## 3. Birthdays
    indiv <- birthdays(indiv)
    
    ## 4. Sampling
    if (y %in% sampling_years) {
      indiv <- capture(
        indiv = indiv, 
        n = sample_size[y], 
        year = y, 
        fatal = lethal_sampling
      )
    }
    cat(" Cycle completed!\n")
  }
  ## How many individuals are still alive in 'indiv'?
  indiv_alive[sim_i] <- nrow(indiv[is.na(indiv$DeathY), ])
}

summary(indiv_alive)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 7781    8439    8630    8614    8768    9122 
