## ========================================================================== ##
## GR_sims_Palmyra.R
## 
## In this script I will simulate GR population data using information life 
##  history trait estimates from the Palmyra study (Bradley et al, 2017).
##  Also see the script "checking_vbgf_Bradley_2017.R".
##  The one thing I do not follow from the paper is the survival, as this is
##  what I adjust to achieve a constant population size through time. 

## Little note: stacking two data.frames with data.table::rbindlist() is almost
##  twice as fast as turning them into matrices and then using rbind().. 
##  MIND = BLOWN!

## Also, starting using environments as lists:
##  https://adv-r.hadley.nz/environments.html
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
first_litter <- 3:6            # first litter 
batch_size <- 3:6               # all possible values for random draw
max_clutch <- Inf               # maximum litter sizes
max_age <- 63                   # maximum age
mort_rate <- 0.1113              # flat mortality rate (0.1115 is REALLY close, slightly negative)
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

## Store simulated sets
n_cores <- 20
cl <- makeCluster(n_cores)
rm(simulated_data_sets)
clusterExport(cl, c(ls()))
simulated_data_sets <- pblapply(1:1000, function(i) {
  
  ## Set a seed for reproducibility
  set.seed(170593 + 121 * i)
  
  # system.time({
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
    Me = integer(0),
    Sex = integer(0),
    Dad = integer(0),
    Mum = integer(0),
    BirthY = integer(0),
    DeathY = integer(0),
    Stock = integer(0),
    AgeLast = integer(0),
    SampY = integer(0),
    Pregnant = integer(0)
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
      ageMort = rep(mort_rate, max_age + 1L),
      maxAge = max_age - 1L 
    )
    
    ## 3. Age incrementation
    indiv <- fishSim::birthdays(indiv)
    
    ## 4. Every 5 years, bury the dead, and remove them from indiv
    ## NOTE: always remove those who died from at least two years ago.
    ## It is important to keep the fathers who died in, as they might have died
    ## between giving mating and birth, but still need to be assigned the 
    ## fatherhood.
    if (y %% 5 == 0) {
      graveyard <- as.data.frame(
        data.table::rbindlist(list(graveyard, indiv[!is.na(indiv$DeathY) &
                                                      indiv$DeathY <= y - 2, ])))
      indiv <- indiv[is.na(indiv$DeathY) | indiv$DeathY > y - 2, ]
    }

    # cat(" Cycle completed!\n")
  }
  
  # })
  
  ## Dig up the dead and put them on a pile with the living
  indiv <- as.data.frame(data.table::rbindlist(list(indiv, graveyard)))
  
  # nrow(indiv[is.na(indiv$DeathY), ])
  
  # cat(paste0("Simulation ", i, " completed!\n"))
  # simulated_data_sets[[i]] <- indiv
  return(indiv)
}, cl = cl); stopCluster(cl);

## How many individuals are still alive in 'indiv'?
hist(sapply(simulated_data_sets, function(x) {nrow(x[is.na(x$DeathY), ])}), xlab = "indivs")
summary(sapply(simulated_data_sets, function(x) {nrow(x[is.na(x$DeathY), ])}))

## :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
## Run retrospecitve sampling
## :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

## Set sampling parameters
n_samples <- 750
n_sample_year <- 2
sampling_years <- c((max(years) - n_sample_year + 1):max(years))     # years in which sampling occurs
sample_size <- c(rep(0, sampling_years[1] - 1), 
                 rep(n_samples / n_sample_year, n_sample_year))   # number of sampled individuals in each year
lethal_sampling <- FALSE        # is sampling lethal?
retrospective_sampling <- TRUE

## Remove previous sampling if required
simulated_data_sets <- lapply(simulated_data_sets, function(indiv) {
  indiv$SampY <- NA
  return(as.data.frame(indiv[, 1:10])) # only return first 10 columns (before sampling columns)
})
all(is.na(simulated_data_sets[[1]]$SampY)) # should be TRUE before continuing
ncol(simulated_data_sets[[1]]) == 10 # should also be TRUE

## Add sampling
if (retrospective_sampling) {
  simulated_data_sets <- pblapply(simulated_data_sets, function(indiv) {
    for (year in sampling_years) {
      n <- sample_size[year]
      indiv <- CKMRcpp::retroCapture2(indiv, n = n, year = year, fatal = lethal_sampling)
    }
    indiv$no_samples <- rowSums(!is.na(indiv[, 11:ncol(indiv)]))
    return(indiv)
  })
}

all(is.na(simulated_data_sets[[1]]$SampY))  # should be FALSE
unique(simulated_data_sets[[1]]$SampY)      # check if this seems correct
sum(!is.na(simulated_data_sets[[1]]$SampY)) # seem correct as well?

## ::::::::::::::::::::::::::::::::
## Summary statistics
## ::::::::::::::::::::::::::::::::

## Extracting simulated abundances
# N_true <- t(sapply(simulated_data_sets, function(x) {
#   mature_females <- sum(is.na(x$DeathY) & x$Sex == 1 & x$AgeLast >= 19)
#   mature_males <- sum(is.na(x$DeathY) & x$Sex == 0 & x$AgeLast >= 17)
#   return(c(N_m = mature_males, N_f = mature_females))
# }))

N_true <- t(sapply(simulated_data_sets, function(x) {
  mature_females <- nrow(
    CKMRcpp::extractTheLiving(x, 2014, F, min_age = first_breed_female + 1,
                              sex = "female"))
  mature_males <- nrow(
    CKMRcpp::extractTheLiving(x, 2014, F, min_age = first_breed_male,
                              sex = "male"))
  return(c(N_m = mature_males, N_f = mature_females))
}))

## Derive growht rate based on animals alive AT THE END OF YEAR
# r_true <- sapply(simulated_data_sets, function(x) {
#   alive_1915 <- sum(!is.na(x$DeathY) & x$DeathY > 1915 & x$BirthY <= 1915)
#   alive_2014 <- sum(is.na(x$DeathY) | x$DeathY > 2014 & x$BirthY <= 2014)
#   return((alive_2014 / alive_1915) ^ (1 / 99))
# })

r_true <- sapply(simulated_data_sets, function(x) {
  alive_1915 <- nrow(CKMRcpp::extractTheLiving(x, 1915, F)) # alive at end of year
  alive_2014 <- nrow(CKMRcpp::extractTheLiving(x, 2014, F)) # alive at end of year
  return((alive_2014 / alive_1915) ^ (1 / (2014 - 1915)))
})


summary(data.frame(N_true, r_true))
