## ========================================================================== ##
## GR_sims_Palmyra.R
## 
## In this script I will simulate GR population data using information life 
## history trait estimates from the Palmyra study (Bradley et al, 2017).
## Also see the script "checking_vbgf_Bradley_2017.R".
## The one thing I do not follow from the paper is the survival, as this is
## what I adjust to achieve a constant population size through time. 
## ========================================================================== ##

## Load packages and source custom functions ===================================
# library(fishSim)
library(parallel)
library(pbapply)
# library(ids)

## Initialise parameters =======================================================
## Set life history parameters
first_breed_male <- 17          # applies only to males
first_breed_female <- 19        # applies only to females
first_litter <- 1:4            # first litter 
batch_size <- 1:4               # all possible values for random draw
max_clutch <- Inf               # maximum litter sizes
max_age <- 63                   # maximum age
mort_rate <- 0.092               # flat mortality rate (0.XXX gives stable)
surv_rate <- 1 - mort_rate      # survival is 1 minus mortality
force_Y1 <- NA                # first-year mortality (not always used)
female_curve <- c(rep(0, first_breed_female),  # female maturity curve 
                  rep(1, 100 - first_breed_female))  
male_curve <- c(rep(0, first_breed_male), # male maturity curve 
                rep(1, 100 - first_breed_male))
mate_type <- "ageSex"           # type of maturity/mating structure
mort_type <- "age"             # type of mortality structure
fecundity_dist <- "uniform"     # uniform distribution (custom)
sex_ratio <- c(0.5, 0.5)        # the offspring male/female sex ratio
single_paternity <- FALSE        # is there a single father to a litter?
no_gestation <- FALSE

## Set simulation and sampling parameters
years <- 1915:2014             # number of years to run simulation

## Store simdata sets
n_cores <- 25
cl <- makeCluster(n_cores)
rm(simulated_data_sets)
clusterExport(cl, c(ls()))
simulated_data_sets <- pblapply(1:1000, function(i) {
  
  ## Set a seed for reproducibility
  set.seed(170593 + 121 * i)
  
  system.time({
  ## Create initial population =================================================
  ## Set initial population in year 0
  indiv <-  CKMRcpp::createFounders(
    pop = 8500,                   # starting pop. size 8500 from literature
    osr = sex_ratio,
    stocks = c(1), # a single stock
    maxAge = max_age,             #
    survCurv = surv_rate ^ (1:(max_age )) / sum(surv_rate ^ (1:(max_age)))
  )
  
  ## Add marker to pregnant females
  indiv <- CKMRcpp::addPregnancy(
    indiv = indiv,
    matingAges = seq(from = first_breed_female, to = max_age, by = 2))
  
  ## Create a graveyard to bury the dead
  graveyard <- data.frame(
    Me = character(0),
    Sex = character(0),
    Dad = character(0),
    Mum = character(0),
    BirthY = character(0),
    DeathY = character(0),
    Stock = character(0),
    AgeLast = character(0),
    SampY = character(0),
    Pregnant = character(0)
  )
  
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
      no_gestation = no_gestation
    )
    
    ## 2. Survival
    indiv <- fishSim::mort(
      indiv = indiv, 
      year = y, 
      type = mort_type, 
      ageMort = rep(mort_rate, max_age + 1),
      maxAge = max_age - 1 
    )
    
    ## 3. Age incrementation
    indiv <- fishSim::birthdays(indiv)
    
    ## 4. Every 5 years, bury the dead, and remove them from indiv
    if (y %% 5 == 0) {
      graveyard <- rbind(graveyard, indiv[!is.na(indiv$DeathY), ])
      indiv <- indiv[is.na(indiv$DeathY), ]
    }

    # cat(" Cycle completed!\n")
  }
  
  })
  
  ## Dig up the dead and add them to indiv
  indiv <- rbind(indiv, graveyard)
  
  # nrow(indiv[is.na(indiv$DeathY), ])
  
  # cat(paste0("Simulation ", i, " completed!\n"))
  # simulated_data_sets[[i]] <- indiv
  return(indiv)
}, cl = cl); stopCluster(cl);

## How many individuals are still alive in 'indiv'?
nrow(indiv[is.na(indiv$DeathY), ])

summary(sapply(simulated_data_sets, function(x) {nrow(x[is.na(x$DeathY), ])}), xlab = "indivs")

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
