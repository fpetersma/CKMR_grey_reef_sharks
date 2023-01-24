## ========================================================================== ##
## vanilla_sims.R
## 
## In this script I will simulate the vanilla population data using a mix of 
## information on life  history trait estimates from the GBR, Hawaii and Palmyra 
## study. This is a very basic population and details can be found in the 
## Overleaf document.
## ========================================================================== ##

## Load packages and source custom functions ===================================
# library(fishSim)
library(parallel)
library(pbapply)
# library(ids)

## Initialise parameters =======================================================
## Set life history parameters
first_breed_male <- 10          # applies only to males
first_breed_female <- 10        # applies only to females
first_litter <- 2               # first litter 
batch_size <- 2                 # all possible values for random draw
max_clutch <- Inf               # maximum litter sizes
max_age <- 19                   # maximum age
mort_rate <- 0.1535              # flat mortality rate (0.1535 gives stable population, but worked with 0.153 which induced a slgiht growth)
surv_rate <- 1 - mort_rate      # survival is 1 minus mortality
force_Y1 <- NA                  # first-year mortality (not always used)
female_curve <- c(rep(0, first_breed_female), 
                  rep(1, 100 - first_breed_female))   # female maturity curve for categories 0:maxAge
male_curve <- c(rep(0, first_breed_male), 
                rep(1, 100 - first_breed_male))     # male maturity curve for categories 0:maxAge
mate_type <- "ageSex"           # type of maturity/mating structure
mort_type <- "age"             # type of mortality structure
fecundity_dist <- "uniform"     # uniform distribution (custom)
sex_ratio <- c(0.5, 0.5)        # the offspring male/female sex ratio
single_paternity <- FALSE        # is there a single father to a litter?
no_gestation <- TRUE

## Set simulation and sampling parameters
years <- 1915:2014             # number of years to run simulation

## Store simdata sets
n_cores <- 50
cl <- makeCluster(n_cores)
# rm(simulated_data_sets)
clusterExport(cl, c(ls()))
simulated_data_sets <- pblapply(1:1000, function(i) {
  
  ## Set a seed for reproducibility
  set.seed(170593 + 121 * i)
  
  ## Create initial population =================================================
  ## Set initial population in year 0
  indiv <-  CKMRcpp::createFounders(
    pop = 8500,                   # starting pop. size 8500 from literature
    osr = sex_ratio,
    stocks = c(1), # a single stock
    maxAge = max_age,         
    y0 = min(years) - 1,
    survCurv = surv_rate ^ (1:(max_age )) / sum(surv_rate ^ (1:(max_age)))
  )
  
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
    SampY = integer(0)
  )
  
  ## Start the simulation ========================================================
  for (year in years) {
    # cat("Starting simulation of year:", y, "...")
    
    ## 1. Mating/birth
    indiv <- CKMRcpp::mateOrBirth(
      indiv = indiv,
      batchSize = batch_size,
      fecundityDist = fecundity_dist, # uniform is custom
      osr = sex_ratio,
      year = year,
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
    indiv_g <- indiv
    ## 2. Survival
    indiv <- fishSim::mort(
      indiv = indiv, 
      year = year, 
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
    if (year %% 5 == 0) {
      graveyard <- as.data.frame(
        data.table::rbindlist(list(graveyard, indiv[!is.na(indiv$DeathY) &
                                                      indiv$DeathY <= year - 2, ])))
      indiv <- indiv[is.na(indiv$DeathY) | indiv$DeathY > year - 2, ]
    }
    
  }
  ## Dig up the dead and put them on a pile with the living
  indiv <- as.data.frame(data.table::rbindlist(list(indiv, graveyard)))
  # cat(paste0("Simulation ", i, " completed!\n"))
  # simulated_data_sets[[i]] <- indiv
  return(indiv)
}, cl = cl); stopCluster(cl);

## How many individuals are still alive in 'indiv'?
hist(sapply(simulated_data_sets, function(x) {nrow(x[is.na(x$DeathY), ])}), xlab = "indivs")
summary(sapply(simulated_data_sets, function(x) {nrow(x[is.na(x$DeathY), ])}))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 6542    8001    8406    8398    8798   10208 

# Create 100 repeats
simulated_populations <- simulated_data_sets
## Choose from 1, 144, 333, 800, 
simulated_data_sets <- rep(simulated_populations[144], 100)

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
  return(as.data.frame(indiv[, 1:9])) # only return first 10 columns (before sampling columns)
})
all(is.na(simulated_data_sets[[1]]$SampY)) # should be TRUE before continuing
ncol(simulated_data_sets[[1]]) == 9 # should also be TRUE

## Add sampling
if (retrospective_sampling) {
  simulated_data_sets <- pblapply(simulated_data_sets, function(indiv) {
    for (year in sampling_years) {
      n <- sample_size[year]
      indiv <- CKMRcpp::retroCapture2(indiv, n = n, year = year, fatal = lethal_sampling)
    }
    indiv$no_samples <- rowSums(!is.na(indiv[, 10:ncol(indiv)]))
    return(indiv)
  })
}

all(is.na(simulated_data_sets[[1]]$SampY))  # should be FALSE
unique(simulated_data_sets[[1]]$SampY)      # check if this seems correct
sum(!is.na(simulated_data_sets[[1]]$SampY)) # seem correct as well?


# save(simulated_data_sets, file = "data/1_population_multiple_sampling_schemes/vanillus/1000_schemes_simulated_data_set.RData")

## :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
## Create summary statistics
## :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# par(mfrow = c(3, 1))
## Extracting simulated abundances
N_true <- t(sapply(simulated_data_sets, function(x) {
  mature_females <- nrow(
    CKMRcpp::extractTheLiving(x, 2014, F, min_age = first_breed_female,
                              sex = "female"))
  mature_males <- nrow(
    CKMRcpp::extractTheLiving(x, 2014, F, min_age = first_breed_male,
                              sex = "male"))
  return(c(N_m = mature_males, N_f = mature_females))
}))


r_true <- sapply(simulated_data_sets, function(x) {
  alive_1915 <- nrow(CKMRcpp::extractTheLiving(x, 1915, F)) # alive at end of year
  alive_2014 <- nrow(CKMRcpp::extractTheLiving(x, 2014, F)) # alive at end of year
  return((alive_2014 / alive_1915) ^ (1 / (2014 - 1915)))
})

summary(data.frame(N_true, r_true))

hist(N_true[, 1], main = "", xlab = "Number of mature males"); abline(v = c(mean(N_true[, 1]), median(N_true[, 1])), col = "red", lty = c(1, 2))

hist(N_true[, 2], main = "",  xlab = "Number of mature females"); abline(v = c(mean(N_true[, 2]), median(N_true[, 2])), col = "red", lty = c(1, 2))

hist(r_true, main = "", xlab = "Yearly growth rate"); abline(v = c(mean(r_true), median(r_true)), col = "red", lty = c(1, 2))
