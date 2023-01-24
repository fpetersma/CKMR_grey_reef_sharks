# # ## Run code below to check the functions, and the pmf should sum to one
sum(fAgeGivenLength(a = 0:19, l = rep(148, 20), sigma_l = 5, p_geom = 0.1535,
                    max_age = 19, pmf_age = "geom"))
# # ## Below should be below one
sum(fAgeGivenLength(a = 6:8, l = rep(148, 3), sigma_l = 5, p_geom = 0.1535,
                    max_age = 19, pmf_age = "geom"))

sum(fAgeGivenLength(a = 10, l = 162, sigma_l = 3, p_geom = 0.1535,
                    max_age = 19, pmf_age = "geom"))

fAgeGivenLength <- function(a, 
                            l,
                            sigma_l,
                            p_geom,
                            max_age,
                            pmf_age = "geom", 
                            vbgf_pars = c(l_inf = 175, 
                                          k = 0.2,
                                          a_0 = -2)) {
  ## Run line below for testing
  # p_geom <- 1 - 0.8455; max_age <- 19; sigma_l <- 2; l <- c(151); a <- c(8);
  
  if (length(l) != length(a)) {
    stop("for every l there should be an a!")
  }
  
  prob_length_given_age <- fLengthGivenAge(l, a, sigma_l, vbgf_pars)
  prob_sampled_age <- fSampledAge(a, p_geom, max_age)
  ## OLD -----------------------
  # expected_length <- vbgf_pars["l_inf"] *
  #   (1 - exp(-vbgf_pars["k"] * (a - vbgf_pars["a_0"])))
  # prob_length_given_age <- 
  #  pnorm(q = l + 0.5, mean = expected_length, sd = sigma_l) - 
  #   pnorm(q = l - 0.5, mean = expected_length, sd = sigma_l)
  # prob_sampled_age <- 
  #  dgeom(a, p_geom) / pgeom(20, p_geom)
  ## ---------------------------
  
  ## requires the function fSampledLength()
  prob_sampled_length <- fSampledLength(l = l, 
                                        sigma_l = sigma_l,
                                        max_age = max_age,
                                        p_geom = p_geom,
                                        pmf_age = pmf_age,
                                        vbgf_pars = vbgf_pars)
  
  out <- prob_length_given_age * prob_sampled_age / prob_sampled_length
  
  return(out)
}

fSampledLength <- function(l, 
                           sigma_l,
                           max_age,
                           p_geom,
                           pmf_age = "geom", 
                           vbgf_pars = c(l_inf = 175, 
                                         k = 0.2,
                                         a_0 = -2)) {
  ## Run line below for testing
  # p_geom <- 1 - 0.8455; max_age <- 19; sigma_l <- 1; l <- c(135, 141, 58);
  
  if (pmf_age == "geom") {
    ## Use sapply() to allow for l to have length > 1
    prob_mass <- sapply(l, function(l_i) {
      ## Create potential ages and convert to expected lengths through the vbgf
      potential_ages <- 0:max_age
      
      return(sum(fLengthGivenAge(l_i, potential_ages, sigma_l, vbgf_pars) * 
                   fSampledAge(potential_ages, p_geom, max_age)))
    })
    
    return(prob_mass)
  } else {
    stop("Function only implemented for geometric age distribution as of yet.")
  }
}

fSampledAge <- function(a, p_geom, max_age) {
  if (any(a < 0 | a > max_age)) {
    stop("all values for a have to be non-negative and at most equal to max_age")
  }
  
  return(dgeom(a, p_geom) / pgeom(max_age + 1, p_geom))
}


fLengthGivenAge <- function(l, a, sigma_l, vbgf_pars = c(l_inf = 175, 
                                                k = 0.2,
                                                a_0 = -2)) {
  expected_length <- vbgf_pars["l_inf"] * 
    (1 - exp(-vbgf_pars["k"] * (a - vbgf_pars["a_0"])))
  
  out <- pDiscreteNorm(x = l, mu = expected_length, sigma = sigma_l)
  
  return(out)
}
  
pDiscreteNorm <- function(x, mu, sigma) {
  return(pnorm(x + 0.5, mean = mu, sd = sigma) - 
           pnorm(x - 0.5, mean = mu, sd = sigma))
}


## OLD fSampledLength()
# fSampledLength <- function(l, 
#                            sigma_l,
#                            max_age,
#                            p_geom,
#                            pmf_age = "geom", 
#                            vbgf_pars = c(l_inf = 175, 
#                                          k = 0.2,
#                                          a_0 = -2)) {
#   ## Run line below for testing
#   # p_geom <- 1 - 0.8455; max_age <- 19; sigma_l <- 1; l <- c(135, 141, 58);
#   
#   if (pmf_age == "geom") {
#     ## Use sapply() to allow for l to have length > 1
#     prob_mass <- sapply(l, function(l_i) {
#       ## Create potential ages and convert to expected lengths through the vbgf
#       potential_ages <- 0:max_age
#       expected_lengths <- vbgf_pars["l_inf"] * 
#         (1 - exp(-vbgf_pars["k"] * (potential_ages - vbgf_pars["a_0"])))
#       
#       return(sum(
#         (pnorm(q = l_i + 0.5, mean = expected_lengths, sd = sigma_l) - 
#            pnorm(q = l_i - 0.5, mean = expected_lengths, sd = sigma_l)) * # discrete normal pmf
#           dgeom(potential_ages, p_geom) / # pmf geom part
#           pgeom(20, p_geom))) # cdf geom part
#     })
#     
#     return(prob_mass)
#   } else {
#     stop("Function only implemented for geometric age distribution as of yet.")
#   }
# }

