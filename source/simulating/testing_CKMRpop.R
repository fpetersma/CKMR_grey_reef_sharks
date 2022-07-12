## ========================================================================== ##
## Script to explore the functionality of the CKMRpop package by Eric         ##
## Anderson.                                                                  ##
## ========================================================================== ##

## Load package and initiate data
## ========================================================================== ##

library(CKMRpop)

## Life history
data(species_1_life_history) # example of life history object
gr_life_history <- list("max-age" = 18,
                        "fem-surv-probs" = rep(0.74, 18),
                        "fem-asrf" = c(rep(0, 11), # immature first 11 years
                                     0.5, # half fecundity first year
                                     rep(1, 6)), # fully mature remainder life
                        "fem-prob-repro" = c(rep(0, 11), # no repro first 11 years
                                           rep(1, 7)), # always repro, unless with baby
                        "repro-inhib" = 1, # mom will need 12-14 months for baby, so breed biennially
                                            # problem is that this only works with binary repro, which does not apply to this species
                                            # maybe I should send an email to Eric Anderson?
                        "male-surv-probs" = rep(0.74, 18),
                        "male-asrp" = c(rep(0, 9), # relative repro potential of males. identical when mature
                                      rep(1, 9)),
                        "male-prob-repro" = c(rep(0, 9), # all males will reproduce, if possible
                                            rep(1, 9)),
                        "offsp-dsn" = "binary", # I thnk this is right, but still needs some more understanding
                        "mate-fidelity" = -1, # female mate with a random male
                        "sex-ratio" = 0.5 # probability of offspring being male
                        )

SPD <- gr_life_history

# before we tell spip what the cohort sizes are, we need to 
# tell it how long we will be running the simulation
SPD$`number-of-years` <- 100

# this is our cohort size
cohort_size <- 1000

# Do some matrix algebra to compute starting values from the
# stable age distribution:
L <- leslie_from_spip(SPD, cohort_size)

# then we add those to the spip parameters
SPD$`initial-males` <- floor(L$stable_age_distro_fem)
SPD$`initial-females` <- floor(L$stable_age_distro_male)

# tell spip to use the cohort size
SPD$`cohort-size` <- paste("const", cohort_size, collapse = " ")

## Sample 5%
samp_frac <- 0.05
samp_start_year <- 80
samp_stop_year <- 100
SPD$`discard-all` <- 0
SPD$`gtyp-ppn-fem-post` <- paste(
  samp_start_year, "-", samp_stop_year, " ", 
  samp_frac, " ", samp_frac, " ", samp_frac, " ",
  paste(rep(0, SPD$`max-age` - 3), collapse = " "),
  sep = ""
)
SPD$`gtyp-ppn-male-post` <- SPD$`gtyp-ppn-fem-post`

set.seed(5)  # set a seed for reproducibility of results
spip_dir <- run_spip(pars = SPD, )  # run spip
slurped <- slurp_spip(spip_dir, 1) # read the spip output into R
