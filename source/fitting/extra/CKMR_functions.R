## ========================================================================== ##
## A collections of functions required for fitting a CKMR model.
## The following functions are included:
##  - 
##  -
##  -
## ========================================================================== ##
nllR <- function(dat, par) {
  ## ===========================================================================
  ## Derive look up table for age-at-length probabilities, if required 
  ## ---------------------------------------------------------------------------
  ages <- 0:dat$max_age
  lengths <- 1:dat$max_length
  age_length_prob_long <- expand.grid(ages, lengths)
  
  ## A long table format version
  age_length_prob_long$prob <- dnorm(x = age_length_prob_long$Var2, 
                                     mean = vbgf(age_length_prob_long$Var1), 
                                     sd = 2)
  ## A matrix version with length on the rows and age on the columns
  age_length_prob_matrix <- reshape2::acast(data = age_length_prob_long, 
                                            formula = Var2 ~ Var1, 
                                            value.var = "prob")
  
  ## ===========================================================================
  ## Estimate the negative log-likelihood in Rcpp 
  ## ---------------------------------------------------------------------------
  sourceCpp("C:/Users/felix/OneDrive - University of St Andrews/Documents/University of St Andrews/PhD/Sharks/Close-kin mark-recapture/simulation_software/nllCKMRcpp.cpp")
  
  nll <- nll_CKMR_Rcpp(dat, par)
  
  return(nll)
  
}

#' vbgf()
#' A function that derives the length based on age and the vbgf parameters
#'
#' @param a Age (numeric)
#' @param t_0 Theoretical age when length is zero (numeric)
#' @param k Growth parameter (numeric)
#' @param l_inf Asymptotic length (numeric)
#'
#' @return
#' @export
#'
#' @examples
vbgf <- function(a, t_0 = -3.5, k = 0.1, l_inf = 175) {
  return(l_inf * (1 - exp(-k * (a - t_0))))
}

#' invvbgf()
#' A function that derives the length based on age and the vbgf parameters
#'
#' @param a Age (numeric)
#' @param t_0 Theoretical age when length is zero (numeric)
#' @param k Growth parameter (numeric)
#' @param l_inf Asymptotic length (numeric)
#'
#' @return
#' @export
#'
#' @examples
invvbgf <- function(l, t_0 = -3.5, k = 0.1, l_inf = 175) {
  return(t_0 - log(1 - l/l_inf) / k)
}

#' Title
#'
#' @param pair A string indicating the type of pair/kinship (character)
#' @param s1 Sex of individual 1 (character)
#' @param c1 Capture year of individual 1 (numeric)
#' @param l1 Length when captured of individual 1 (numeric)
#' @param s2 Sex of individual 2 (character)
#' @param c2 Capture year of individual 2 (numeric)
#' @param l2 Length when captured of individual 2 (numeric)
#' @param alpha_m Male age of maturity (numeric)
#' @param alpha_f Female age of maturity (numeric)
#' @param max_age Maximum age (numeric)
#' @param phi Survival probability from year t to t+1 (numeric)
#' @param N_t0_f Female population size in year t0 (numeric)
#' @param N_t0_m Male population size in year t0 (numeric)
#' @param t0 Reference year (numeric)
#' @param r Population growth rate (numeric)
#' @param sigma_vbgf Standard error on the age-at-length relationship (numeric)
#'
#' @return
#' @export
#'
#' @examples
#' observation_prob <- 1 - pairProb(pair = "PO", 
#'                                  s1 = df$indiv_1_sex[1], 
#'                                  c1 = df$indiv_1_capture_year[1],
#'                                  l1 = df$indiv_1_length[1],
#'                                  s2 = df$indiv_2_sex[1],
#'                                  c2 = df$indiv_2_capture_year[1],
#'                                  l2 = df$indiv_2_length[1],
#'                                  alpha_m = 10, 
#'                                  alpha_f = 12, 
#'                                  max_age = 20,
#'                                  phi = 0.9,
#'                                  N_t0 = c(female = 3000, male = 3000), # named vector
#'                                  t0 = 80, 
#'                                  r = 0.02, 
#'                                  sigma_vbgf = 1.5)
#' 
pairProb <- function(pair, s1, c1, l1, s2, c2, l2, alpha_m, alpha_f, max_age,
                     phi, N_t0_f, N_t0_m, t0, r, sigma_vbgf, 
                     al_prob_matrix = NULL) {
  
  ## Derive potential ages for parent
  ages_parent <- 0:max_age
  if (pair == "PO") {
    # if (s1 == "F") {
    
    output_parent_level <- rep(NA, length(ages_parent))
    ## Loop over ages for parent
    for (index_age_parent in seq_along(ages_parent)) {
      a1 <- ages_parent[index_age_parent]
      y1 <- c1 - a1
      
      if (s1 == "F") {
        ## Set maximum age for the offspring
        max_age_offspring <- c2 - c1 - alpha_f + a1
      }
      if (s1 == "M") {
        ## Set maximum age for the offspring
        max_age_offspring <- c2 - c1 - alpha_m + a1
      }
      print(max_age_offspring)
      
      ## Probability equals zero if max_age_offspring is negative
      if (max_age_offspring < 0) {
        output_parent_level[index_age_parent] <- 0
      } else {
        ## Derive potential ages for offspring
        ages_offspring <- 0:max_age_offspring
        
        ## Create output for prob for every potential age of the offspring
        output_offspring_level <- rep(NA, length(ages_offspring))
        
        ## Loop over ages for offspring
        for (index_age_offspring in seq_along(ages_offspring)) {
          a2 <- ages_offspring[index_age_offspring]
          y2 <- c2 - a2
          
          if (s1 == "F") {
            out <- 1 / (N_t0_f * (1 + r) ^ (y2 - t0))
          }
          if (s1 == "M") {
            out <- 1 / (N_t0_m * (1 + r) ^ (y2 - t0))
          }
          ## Account for survival of parent i if j was born after c1 was captured
          if (c1 < y2) {
            out <- out * phi ^ (c1 - y2)
          }
          
          ## Multiply above by probability density of offspring age
          out <- out * dnorm(x = l2, mean = vbgf(a2), sd = sigma_vbgf)
          # out <- out * al_prob_matrix[as.character(l2), as.character(a2)]
          
          
          ## Add to output_offspring_level
          output_offspring_level[index_age_offspring] <- out
          
          print(out)
        }
        ## Sum the probabilities for every potential offspring age
        output_offspring_level <- sum(output_offspring_level)
        
        ## Multiply above by probability density of parent age
        output_offspring_level <- output_offspring_level * dnorm(x = l1,
                                                                 mean = vbgf(a1),
                                                                 sd = sigma_vbgf)
        ## Look up table version is slower, so don't use the line below
        # output_offspring_level <- output_offspring_level * 
        #   al_prob_matrix[as.character(l1), as.character(a1)]
        
        ## Add to output_parent_level
        output_parent_level[index_age_parent] <- output_offspring_level
        
        print(output_offspring_level)
      }
    }
    ## Sum the probabilities for every potential offspring age
    prob <- sum(output_parent_level)
    # } else if (s1 == "M") {
    #   output_parent_level <- rep(NA, length(ages_parent))
    #   ## Loop over ages for parent
    #   for (index_age_parent in seq_along(ages_parent)) {
    #     a1 <- ages_parent[index_age_parent]
    #     y1 <- c1 - a1
    #     
    #     ## Set maximum age for the offspring
    #     max_age_offspring <- c2 - c1 - alpha_m + a1
    #     
    #     ## Probability equals zero if max_age_offspring is negative
    #     if (max_age_offspring < 0) {
    #       output_parent_level[index_age_parent] <- 0
    #     } else {
    #       ## Derive potential ages for offspring
    #       ages_offspring <- 0:max_age_offspring
    #       
    #       ## Create output for prob for every potential age of the offspring
    #       output_offspring_level <- rep(NA, length(ages_offspring))
    #       
    #       ## Loop over ages for offspring
    #       for (index_age_offspring in seq_along(ages_offspring)) {
    #         a2 <- ages_offspring[index_age_offspring]
    #         y2 <- c2 - a2
    #         
    #         out <- 1 / (N_t0["male"] * (1 + r) ^ (c2 - a2 - t0))
    #         
    #         ## Account for survival of parent i if j was born after c1 was captured
    #         if (c1 < y2) {
    #           out <- out * phi ^ (c1 - y2)
    #         }
    #         
    #         ## Multiply above by probability density of offspring age
    #         out <- out * dnorm(x = l2, mean = vbgf(a2), sd = sigma_vbgf)
    #         # out <- out * al_prob_matrix[as.character(l2), as.character(a2)]
    #         
    #         ## Add to output_offspring_level
    #         output_offspring_level[index_age_offspring] <- out
    #       }
    #       ## Sum the probabilities for every potential offspring age
    #       output_offspring_level <- sum(output_offspring_level)
    #       
    #       ## Multiply above by probability density of parent age
    #       output_offspring_level <- output_offspring_level * dnorm(x = l1,
    #                                                                mean = vbgf(a1),
    #                                                                sd = sigma_vbgf)
    #       ## Look up table version is slower, so don't use the line below
    #       # output_offspring_level <- output_offspring_level * 
    #       #   al_prob_matrix[as.character(l1), as.character(a1)]
    #       
    #       ## Add to output_parent_level
    #       output_parent_level[index_age_parent] <- output_offspring_level
    #     }
    #   }
    #   ## Sum the probabilities for every potential offspring age
    #   prob <- sum(output_offspring_level)
    # }
  } else if (pair == "HS") {
    ## Not implemented yet
  }
  print(paste0("prob = ", prob))
  
  return(prob)
}