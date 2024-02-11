## ========================================================================== ##
## GR_sims_Palmyra.R
## 
## In this script I will simulate GR population data using information life 
## history trait estimates from a Hawaii study. The growth curve is very 
## different, and so is the maximum age.
## ========================================================================== ##

## Load packages and source custom functions ===================================
library(fishSim)
library(ids)
source("C:/Users/felix/OneDrive - University of St Andrews/Documents/University of St Andrews/PhD/Sharks/Close-kin mark-recapture/simulation_software/custom_functions_fishSim.R")

## Initialise parameters =======================================================
## Set life history parameters
recovery_time <- 2              # years until next breeding/litter
first_breed <- 18               # applies only to females
first_litter <- 3:6             # first litter 
batch_size <- 3:6               # all possible values for random draw
max_clutch <- Inf               # maximum litter sizes
max_age <- 60                   # maximum age
mort_rate <- 0.115               # flat mortality rate (0.115 gives stable)
surv_rate <- 1 - mort_rate      # survival is 1 minus mortality
force_Y1 <- NA                # first-year mortality (not always used)
female_curve <- c(rep(0, 19), 
                  rep(c(0, 1), 44))   # female maturity curve for categories 0:maxAge
male_curve <- c(rep(0, 17), 
                rep(1, 90))     # male maturity curve for categories 0:maxAge
mate_type <- "ageSex"           # type of maturity/mating structure
mort_type <- "flat"             # type of mortality structure
fecundity_dist <- "uniform"     # uniform distribution (custom)
sex_ratio <- c(0.5, 0.5)        # the offspring male/female sex ratio
single_paternity <- TRUE        # is there a single father to a litter?

## Set simulation and sampling parameters
years <- 1:10                  # number of years to run simulation, makes sure all starting ind are dead by the end
sampling_years <- c(max(years) - 1, max(years))     # years in which sampling occurs
sample_size <- c(rep(0, sampling_years[1] - 1), 
                 600, 700)   # number of sampled individuals in each year
lethal_sampling <- FALSE        # is sampling lethal?

## Create initial population ===================================================
## Set initial population in year 0
indiv <-  makeFounders(
  pop = 8500,                   # starting pop. size 8500 from literature
  osr = sex_ratio,
  stocks = c(1),                # a single stock
  maxAge = max_age,             #
  survCurv = surv_rate ^ (1:max_age) / sum(surv_rate ^ (1:max_age))
)

## Add birth recovery remaining
indiv <- addBirthRecovery(
  indiv = indiv,
  recoveryTime = recovery_time,
  maturityAge = first_breed,
  random = FALSE,
  breedingYears <- c(seq(from = 18, to = 60, by = 2))
)

## Start the simulation ========================================================

for (y in years) {
  cat("Starting simulation of year:", y, "...")
  # 1. Mating
  indiv <- mateWithRecovery(
    indiv = indiv,
    recovery = recovery_time,
    batchSize = batch_size,
    fecundityDist = fecundity_dist, # uniform is custom
    osr = sex_ratio,
    year = y,
    firstBreed = first_breed,
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
  
  ## 4. Recovery
  indiv <- recover(indiv)
  
  # 5. Sampling
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

## Start 100  simulations ======================================================
indiv_alive_Palmyra <- rep(NA, 100)
for (sim_i in 1:100) {
  ## Set initial population in year 0
  indiv <-  makeFounders(
    pop = 8500,                   # starting pop. size 8500 from literature
    osr = sex_ratio,
    stocks = c(1),                # a single stock
    maxAge = max_age,             #
    survCurv = surv_rate ^ (1:max_age) / sum(surv_rate ^ (1:max_age))
  )
  
  ## Add birth recovery remaining
  indiv <- addBirthRecovery(
    indiv = indiv,
    recoveryTime = recovery_time,
    maturityAge = first_breed,
    random = FALSE,
    breedingYears <- c(seq(from = 18, to = 60, by = 2))
  )
  
  
  for (y in years) {
    cat("Starting simulation of year:", y, "...")
    # 1. Mating
    indiv <- mateWithRecovery(
      indiv = indiv,
      recovery = recovery_time,
      batchSize = batch_size,
      fecundityDist = fecundity_dist, # uniform is custom
      osr = sex_ratio,
      year = y,
      firstBreed = first_breed,
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
    
    ## 4. Recovery
    indiv <- recover(indiv)
    
    # 5. Sampling
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
  indiv_alive_Palmyra[sim_i] <- nrow(indiv[is.na(indiv$DeathY), ])
}

summary(indiv_alive_Palmyra)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 8070    8434    8608    8619    8777    9153 

## Rerunning the above with repeated c(1,0) maturity curve and recoverytime = 1,
## instead of recoverytime = 2, gave the following
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 7475    8105    8256    8252    8412    8768 
