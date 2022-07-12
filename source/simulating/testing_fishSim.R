## ========================================================================== ##
## Script to explore the functionality of the fishSim package by SM Baylis/   ##
## https://github.com/SMBaylis/fishSim                                        ##
## ========================================================================== ##

## Install & load package, and initiate data
## ========================================================================== ##

## Install fishSim
# install.packages("fastmatch") # required by fishSim 
# devtools::install_github(repo = "SMBaylis/fishSim") 

## Load library
library(fishSim)

## 1. The relatively simple scenario
## =============================================================================
## 1.1 Setup -------------------------------------------------------------------

indiv <- makeFounders(pop = 8500, 
                      osr = c(0.5, 0.5), 
                      stocks = c(1),
                      maxAge = 18,
                      survCurv = 0.82 ^ (1:18) / sum(0.82 ^ (1:18)))

## 1.2 Mating ------------------------------------------------------------------

# mate() is based on fathers and mothers, whereas altMate() is based on only 
# mothers, and fathers are randomly selected
indiv <- altMate(
  indiv = indiv,
  batchSize = 4,
  fecundityDist = "truncPoisson",
  osr = c(0.5, 0.5),
  year = 1,
  firstBreed = 0,         # set as zero when maturity curves are specified
  type = "ageSex",        # ageSex for age-sex specific maturity curves
  maxClutch = 4,
  singlePaternity = TRUE,  # does a litter have a single father?
  # exhaustFathers = FALSE,
  # maturityCurve, # non age-specific maturity curve
  maleCurve = c(rep(0, 10), rep(1, 9)),   # maturity curves, values in [0, 1]
  femaleCurve = c(rep(0, 12), rep(1, 7)) # maturity curves, values in [0, 1]
  # important is that these curves have 0:maxAge
)

## 1.3 Mortality ---------------------------------------------------------------
indiv <- mort(
  indiv, 
  year = 1, 
  type = "flat", 
  mortRate = 0.26,
  maxAge = 18
)

## 1.4 Birthdays ---------------------------------------------------------------

indiv <- birthdays(indiv)

## 1.5 Sampling ----------------------------------------------------------------

## non-lethal sampling
indiv <- capture(
  indiv, 
  n = 10, 
  year = 1, 
  fatal = FALSE
)

## 1.6 What's my population doing? Can I make it stay the same size? -----------
check_growthrate(
  forceY1 = 0.21,
  mateType = "ageSex",
  mortType = "flat",
  batchSize = 4,
  firstBreed = 0,
  maxClutch = 4,
  osr = c(0.5, 0.5),
  # maturityCurve,
  femaleCurve = c(rep(0, 12), rep(1, 7)),
  maxAge = 18,
  mortRate = 0.16
  # ageMort,
  # stockMort,
  # ageStockMort
)

PoNG(
  mateType = "ageSex",
  mortType = "flat",
  batchSize = 4.1,
  firstBreed = 0,
  maxClutch = 6,
  osr = c(0.50, 0.50),
  femaleCurve = c(rep(0, 7), # 0% change of maturity at age 0--6
                  0.5, # 50% change of maturity at age 7
                  rep(1, 11)), # 100% change of maturity at age 8--18
  maxAge = 18,
  mortRate = 0.22
)


## 1.7 Turning processes into a simple simulation ------------------------------

## Set life history parameters
first_breed <- 0                # should be zero when female_curve is specified
batch_size <- 4.1               # expected litter size per female 
max_clutch <- 6                 # maximum litter sizes
max_age <- 18                   # maximum age
mort_rate <- 0.15               # flat mortality rate
surv_rate <- 1 - mort_rate      # survival is 1 minus mortality
forceY1 <- 0.35                 # first-year mortality
female_curve <- c(rep(0, 12), 
                  rep(1, 88))    # female maturity curve for categories 0:maxAge
male_curve <- c(rep(0, 10), 
                rep(1, 90))     # male maturity curve for categories 0:maxAge
mate_type <- "ageSex"           # type of maturity/mating structure
mort_type <- "flat"             # type of mortality structure
sex_ratio <- c(0.5, 0.5)        # the offspring male/female sex ratio
single_paternity <- TRUE        # is there a single father to a litter?

## Set simulation and sampling parameters
years <- 1:50                   # number of years to run simulation
sampling_years <- c(49, 50)     # years in which sampling occurs
sample_size <- c(rep(0, sampling_years[1] - 1), 
                 600, 700)   # number of sampled individuals in each year
lethal_sampling <- FALSE        # is sampling lethal?

## Set initial population in year 0
indiv <-  makeFounders(
  pop = 8500,                   # starting pop. size 8500 from literature
  osr = sex_ratio,
  stocks = c(1),                # a single stock
  maxAge = max_age,             #
  survCurv = surv_rate ^ (1:max_age) / sum(surv_rate ^ (1:max_age))
)

## Start the simulation
for (y in years) {
  cat("Starting simulation of year:", y, "...")
  # 1. Mating
  indiv <- altMate(
    indiv = indiv,
    batchSize = batch_size,
    fecundityDist = "truncPoisson",
    osr = sex_ratio,
    year = y,
    firstBreed = first_breed,
    type = mate_type, 
    maxClutch = max_clutch,
    singlePaternity = single_paternity,  
    maleCurve = male_curve,   
    femaleCurve = female_curve
  )
  
  # 2. Mortality
  indiv <- mort(
    indiv = indiv, 
    year = y, 
    type = mort_type, 
    mortRate = mort_rate,
    maxAge = max_age
  )
  
  # 3. Birthdays
  indiv <- birthdays(indiv)
  
  if (y %in% sampling_years) {
    # 4. Sampling
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

## 1.8 Looking up relationship between captured pairs --------------------------
pairs <- findRelativesPar(indiv = indiv, 
                          sampled = TRUE, 
                          nCores = 6)
POPs <- pairs[pairs$OneTwo == 1,] ## Parent-Offspring pairs
# HSPs <- pairs[pairs$TwoTwo == 1,] ## Half-sibling pairs

## 1.9 Creating custom functionality of recovery after breeding ----------------

## Create initial pop
indiv <- makeFounders(pop = 8500, 
                      osr = c(0.5, 0.5), 
                      stocks = c(1),
                      maxAge = 18,
                      survCurv = 0.82 ^ (1:18) / sum(0.82 ^ (1:18)))
## Add random birth recovery remaining
indiv <- addBirthRecovery(
  indiv = indiv,
  recoveryTime = 2,
  maturityAge = 12,
  random = FALSE,
  breedingYears <- c(12, 14, 16, 18)
)

## Add mating cycle
indiv <- mateWithRecovery(
  indiv = indiv,
  recovery = 2,
  batchSize = c(3, 4),
  fecundityDist = "uniform", # uniform is custom
  osr = c(0.5, 0.5),
  year = 1,
  firstBreed = 12,
  firstLitter = c(1,2), # firstLitter is custom
  type = "ageSex", 
  maxClutch = Inf,
  singlePaternity = TRUE,  
  maleCurve = c(rep(0, 10), 
                rep(1, 90)) ,   
  femaleCurve = c(rep(0, 12), 
                  rep(1, 88))
)

## Add mortality
indiv <- mort(
  indiv, 
  year = 1, 
  type = "flat", 
  mortRate = 0.18,
  maxAge = 18
)

## Add birthdays and recovery
indiv <- birthdays(indiv)
indiv <- recover(indiv)

## Add sampling
## non-lethal sampling (here we can add non-random sampling)
indiv <- capture(
  indiv, 
  n = 100, 
  year = 1, 
  fatal = FALSE
)

## CONCLUSION ==================================================================
#' This seems to be working. Functionality has been created to have biennial 
#' breeding instead of annual breeding. Moreover, breeding in the first year
#' for female produces 1-2 offspring, whereas this is 3-4 in the following years.
#' This is based on GBR data, and should probably be more similar to 
#' Hawaii data, as that is likely more similar to Palmyra, as it is closer.
#' Breeding in effect means giving birth. As the sharks mature at age 11, their
#' first mating year in reality is 11 years, so breading is 12 to reflect the
#' gestation of roughly 12 months. The sharks then rest for a year, and then 
#' mate again to give birth every two years (biennially). 
#' 
#' Order of operations: 
#' 1) give birth
#' 2) death
#' 3) birthdays and recover
#' 4) potential sampling
#' 
#' I will now write a script to mimic the Palmyra population based on Hawaii 
#' results, as I am not sure if I can use the Bradley 2017 results. This script
#' is called:
#'                       "GR_sims_Hawaii.R"
#'
#'