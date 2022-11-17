## ========================================================================== ##
## A new file with custom functions to add to fishSim to make it more similar ##
## to the grey reef shark case study from Palmyra. Recovered after            ##
## accidentally replacing the file. Very stupid.                              ##  
## ========================================================================== ##

#' Imagine two VBGFs with equal k, function i and curve j.
#' If we specife all parameters for function i, then this function gives you
#' the a_0 for j such that the difference between l_inf_i and l_inf_j is 
#' identical to the difference between the lengths at a = 0. 
#' 
#' Uses CKMRcpp::vbgf()
#' 
#' 
a_0_j <- function(l_inf_i, l_inf_j, k, a_0_i) {
  d <- l_inf_i - l_inf_j
  a_0 <- log(1 - (CKMRcpp::vbgf(0, a_0_i, k, l_inf_i) - d) / l_inf_j) / k
  return(a_0)
}

addBirthRecovery <- function(indiv, 
                             recoveryTime, 
                             maturityAge,
                             random = TRUE,
                             breedingYears = NA) {
  if (random) {
    potential_recovery <- 0:(recoveryTime - 1) # -1 to match the order
    
    indiv$Recovery <- 0
    # Add random years since offspring for mature females
    indiv$Recovery[indiv$Sex == "F" & indiv$AgeLast >= maturityAge] <- 
      sample(potential_recovery, 
             sum(indiv$Sex == "F" & indiv$AgeLast >= maturityAge), 
             TRUE)
  } else {
    indiv$Recovery <- 0
    # Add random years since offspring for mature females
    indiv$Recovery[indiv$Sex == "F" & 
                     indiv$AgeLast >= maturityAge &
                     indiv$AgeLast %in% (breedingYears + 1)] <- 1
  }
  return(indiv)
}

addPregnancy <- function(indiv, matingAges) {
  
  is_pregnant <- indiv$Sex == 1 & indiv$AgeLast %in% (matingAges + 1)
  
  indiv$Pregnant <- as.integer(is_pregnant)
  
  return(indiv)
}

captureOnlyFirst <- function(indiv, 
                             n = 1, 
                             year = "-1", 
                             fatal = FALSE, 
                             sex = NULL, 
                             age = NULL) {
  if (!is.null(sex)) {
    is.sex <- indiv[, 2] == sex
  } else is.sex <- TRUE
  
  if (!is.null(age)) {
    is.age <- indiv[, 8] %in% age
  } else is.age <- TRUE
  
  is.alive <- is.na(indiv[, 6]) & is.sex & is.age 
  is.dead <- !is.alive
  n.alive <- sum(is.alive)
  n <- min(n, n.alive)
  if (n > 0) {
    sample.loc <- sample.int(n.alive, size = n)
    ## The line below is different from fishSim::capture(). 
    ## Ensure that we only consider individuals that have not been
    ## captured before. Basically, from the sampled individuals, only updated
    ## the sample year of the ones that don't have a sample year yet
    sample_loc_new <- sample.loc[is.na(indiv[is.alive, ][sample.loc, 9])]
    indiv[is.alive, ][sample_loc_new, 9] <- year
    if (fatal) {
      indiv[is.alive, ][sample_loc_new, 6] <- year
    }
  }
  return(indiv)
}

#' createFounders() is an almost literal copy of fishSim::makeFounders, but
#' instead calls the function uuid() explicitly through ids::uuid().
createFounders <- function (pop = 1000, 
                            osr = c(0.5, 0.5), 
                            stocks = c(0.3, 0.3, 0.4), 
                            maxAge = 20, 
                            survCurv = 0.7^(1:maxAge)/sum(0.7^(1:maxAge))) {
  if (sum(osr) != 1) 
    warning("osr does not sum to 1")
  if (sum(stocks) != 1) 
    warning("stocks do not sum to 1")
  if (sum(survCurv) != 1) 
    warning("survCurv does not sum to 1")
  if (length(survCurv) != maxAge) 
    warning("survCurv and maxAge imply different maximum ages")
  indiv <- data.frame(Me = integer(pop), Sex = integer(pop), 
                      Dad = integer(pop), Mum = integer(pop), BirthY = integer(pop), 
                      DeathY = integer(pop), Stock = integer(pop), AgeLast = integer(pop), 
                      SampY = integer(pop))
  if (pop > 0) {
    indiv[, 1] <- 1:nrow(indiv) # much less memory than ids::uuid()
    indiv[, 2] <- sample(c(0, 1), pop, TRUE, # 0 is male, 1 is female
                         prob = osr)
    indiv[, 3] <- c(rep(-1, pop)) # -1 indicates founder
    indiv[, 4] <- c(rep(-1, pop)) # -1 indicates founder
    indiv[, 6] <- as.integer(c(rep(NA, pop)))
    indiv[, 7] <- as.integer(sample(1:length(stocks), pop, TRUE, 
                                    prob = stocks))
    indiv[, 8] <- sample.int(maxAge, pop, TRUE, prob = survCurv)
    indiv[, 5] <- 1L - indiv[, 8]
    indiv[, 9] <- as.integer(c(rep(NA, pop)))
  }
  
  return(indiv)
}

extractSampledIndiv <- function(indiv) {
  ## Extract the self captures
  self <- indiv[indiv$no_samples > 1, ]
  
  ## Extract sampled individuals from the population
  sampled_indiv <- indiv[!is.na(indiv$SampY),]
  
  ## Add the recaptures 
  self <- self[rep(seq_len(nrow(self)), self$no_samples), ]
  
  ## Seperate the sample years and allocate one to every occasion that 
  ## individual was sampled. 
  for (id in unique(self$Me)) {
    ## Get all samplings of the same individual with id
    selfs <- self[self$Me == id, ]
    
    ## Extract all the sampling years of this individual
    samp_years <- as.numeric(unlist(strsplit(selfs[1, "SampY"], split = "_")))
    
    ## Assign one sampling year to every sampling occassion. 
    for (i in 1:nrow(selfs)) {
      selfs[i, "SampY"] <- samp_years[i]
    }
    
    ## Updated the rows in 'self'
    self[self$Me == id, ] <- selfs
  }
  ## Bind the single and multiple captures together
  sampled_indiv <- rbind(sampled_indiv[sampled_indiv$no_samples == 1, ], 
                         self)
  
  ## Ensure SampY is numeric
  sampled_indiv$SampY <- as.numeric(sampled_indiv$SampY)
  
  ## Keep relevant information and rename columns to match theory
  sampled_indiv$SampAge <- sampled_indiv$SampY - sampled_indiv$BirthY
  sampled_indiv <- subset(sampled_indiv, select = c(Me, Sex, SampAge, SampY))
  colnames(sampled_indiv) <- c("id", "sex", "capture_age", "capture_year")
  
  ## Return object
  return(sampled_indiv)
}

#' extractTheLiving()
#' 
#' Returns the living at the start of the year (DEFAULT) or end of the year
#' Note that the number of animals alive at the start of year t should match
#' the living at the end of year t-1.
extractTheLiving <- function(indiv,
                             year,
                             start_of_year = TRUE,
                             min_age = 0,
                             sex = NULL) {
  if (is.null(sex)) {
    indiv <- indiv
  } else if (sex == "male") {
    indiv <- indiv[indiv$Sex == 0, ]
  } else if (sex == "female") {
    indiv <- indiv[indiv$Sex == 1, ]
  } else {
    stop("sex should be either male, female, or NULL!")
  }
  
  if (start_of_year) {
    # Born before 'year', and either still alive or dead in 'year' or later
    alive <- indiv$BirthY < year &                     # Born before 'year'?
      ((indiv$DeathY >= year & !is.na(indiv$DeathY)) | # Died in or after 'year'?
         is.na(indiv$DeathY))                          # Still alive?
    # Old enough at the start of the year
    old_enough <-  alive & indiv$BirthY < (year - min_age + 1)
  } else {
    # Born in or before 'year', and either still alive or died after 'year' 
    alive <- indiv$BirthY <= year &                    # Born in or before 'year'?
      ((indiv$DeathY > year & !is.na(indiv$DeathY)) |  # Died after 'year'?
         is.na(indiv$DeathY))                          # Still alive?
    # Old enough at the endof the year
    old_enough <-  alive & indiv$BirthY <= (year - min_age + 1)
  }
  
  living <- subset(indiv, old_enough)
  
  return(living)
}

parents <- function (ID, indiv) {
  outs <- c(rep(0, length(ID) * 2))
  for (i in 1:length(ID)) {
    if (ID[i] == -1) { # -1 is founder
      outs[(i * 2) - 1] <- -1 # -1 is founder
      outs[(i * 2)] <- -1 # -1 is founder
    }
    else {
      outs[(i * 2) - 1] <- indiv[indiv[, 1] == ID[i], 
                                 3]
      outs[(i * 2)] <- indiv[indiv[, 1] == ID[i], 4]
    }
  }
  return(outs)
}

grandparents <- function (ID, indiv) {
  parents(parents(ID, indiv), indiv)
}


#' findRelativesCustom
#' 
#' A version of fishSim::findRelatives() that only looks as far back as great-
#' grandparents.
#'
#' @param indiv 
#' @param sampled 
#' @param verbose 
#' @param delimitIndiv 
#'
#' @return
#' @export
#'
#' @examples
findRelativesCustom <- function(indiv, 
                                sampled = TRUE, 
                                verbose = FALSE,
                                delimitIndiv = TRUE) {
  if (sampled) {
    if (sum(!is.na(indiv[, 9])) == 0) 
      stop("no sampled individuals")
    if (verbose) 
      print(data.frame(table(indiv[!is.na(indiv[, 9]), 
                                   9], dnn = "Sample Year")))
    sampled <- indiv[!is.na(indiv[, 9]), 1]
  } else {
    sampled <- indiv[, 1]
  }
  if (delimitIndiv) {
    keepers <- indiv$Me %in% sampled | indiv$Me %in% indiv$Mum | 
      indiv$Me %in% indiv$Dad
    indiv <- indiv[keepers, ]
  }
  ancestors <- matrix(data = c(sampled, rep(-1, 
                                            length(sampled) * 6)), 
                      nrow = length(sampled))
  colnames(ancestors) <- c("self", "father", "mother", 
                           "FF", "FM", "MF", "MM" #, 
                           # "FFF", "FFM", "FMF", "FMM", "MFF", "MFM", 
                           # "MMF", "MMM"
  )
  for (i in 1:nrow(ancestors)) {
    ancestors[i, 2:3] <- CKMRcpp::parents(ancestors[i, 1], indiv)
    ancestors[i, 4:7] <- CKMRcpp::grandparents(ancestors[i, 1], indiv)
    # ancestors[i, 8:15] <- fishSim::great.grandparents(ancestors[i, 1], indiv)
  }
  expand.grid.unique <- function(x, y, include.equals = FALSE) {
    x <- unique(x)
    y <- unique(y)
    g <- function(i) {
      z <- setdiff(y, x[seq_len(i - include.equals)])
      if (length(z)) 
        cbind(x[i], z, deparse.level = 0)
    }
    do.call(rbind, lapply(seq_along(x), g))
  }
  
  pairs <- expand.grid.unique(ancestors[, 1], ancestors[, 1], 
                              include.equals = TRUE) # include self comparison?
  colnames(pairs) <- c("Var1", "Var2")
  related <- c(rep(NA, nrow(pairs)))
  totalRelatives <- c(rep(NA, nrow(pairs)))
  OneTwo <- c(rep(NA, nrow(pairs)))        # For POPs
  # OneThree <- c(rep(NA, nrow(pairs)))    # For GGPs
  # OneFour <- c(rep(NA, nrow(pairs)))
  TwoTwo <- c(rep(NA, nrow(pairs)))        # For HSPs and FSPs
  # TwoThree <- c(rep(NA, nrow(pairs)))
  # TwoFour <- c(rep(NA, nrow(pairs)))
  # ThreeThree <- c(rep(NA, nrow(pairs)))
  # ThreeFour <- c(rep(NA, nrow(pairs)))
  # FourFour <- c(rep(NA, nrow(pairs)))
  for (i in 1:length(related)) {
    # allAncestors <- ancestors[ancestors[, 1] == pairs[i, 1], 1:7] %in% ## gone for speed
    #   ancestors[ancestors[, 1] == pairs[i, 2], 1:7]
    indi1 <- pairs[i, 1]
    indi2 <- pairs[i, 2]
    indi1Par <- ancestors[ancestors[, 1] == pairs[i, 1], 
                          2:3]
    indi2Par <- ancestors[ancestors[, 1] == pairs[i, 2], 
                          2:3]
    indi1GP <- ancestors[ancestors[, 1] == pairs[i, 1], 4:7]
    indi2GP <- ancestors[ancestors[, 1] == pairs[i, 2], 4:7]
    # indi1GGP <- ancestors[ancestors[, 1] == pairs[i, 1], 
    # 8:15]
    # indi2GGP <- ancestors[ancestors[, 1] == pairs[i, 2], 
    # 8:15]
    # related[i] <- is.integer(any(allAncestors))## gone for speed
    # totalRelatives[i] <- sum(allAncestors) ## gone for speed
    OneTwo[i] <- sum(c(indi1 %in% indi2Par, indi2 %in% indi1Par)) # For POP
    # OneThree[i] <- sum(c(indi1 %in% indi2GP, indi2 %in% indi1GP)) # For GGP
    # OneFour[i] <- sum(c(indi1 %in% indi2GGP, indi2 %in% indi1GGP))
    
    TwoTwo[i] <- sum(c(indi1Par %in% indi2Par)) # For HSP and FSP
    # TwoThree[i] <- sum(c(indi1Par %in% indi2GP, indi2Par %in%
    #                        indi1GP))
    # TwoFour[i] <- sum(c(indi1Par %in% indi2GGP, indi2Par %in% 
    # indi1GGP))
    
    # ThreeThree[i] <- sum(c(indi1GP %in% indi2GP))
    # ThreeFour[i] <- sum(c(indi1GP %in% indi2GGP, indi2GP %in% 
    # indi1GGP))
    
    # FourFour[i] <- sum(c(indi1GGP %in% indi2GGP))
    
    ## Optimise this print function away; it just takes time.
    # if (i %% 1000 == 0) {
    #   cat("\r", i, " of ", length(related), 
    #       " comparisons", sep = "")
    #   flush.console()
    # }
  }
  pairs <- data.frame(pairs, 
                      # related, totalRelatives, 
                      OneTwo, 
                      #OneThree, #OneFour,
                      TwoTwo#, TwoThree, #TwoFour,
                      #ThreeThree#, ThreeFour, 
                      #FourFour
  )
  return(pairs)
}


findRelativesParCustom <- function(indiv, 
                                   sampled = TRUE, 
                                   verbose = TRUE,
                                   nCores = 1, 
                                   delimitIndiv = TRUE) {
  require(doParallel)
  registerDoParallel(nCores)
  if (sampled) {
    if (sum(!is.na(indiv[, 9])) == 0) 
      stop("no sampled individuals")
    if (verbose) 
      print(data.frame(table(indiv[!is.na(indiv[, 9]), 
                                   9], dnn = "Sample Year")))
    sampled <- indiv[!is.na(indiv[, 9]), 1]
  } else {
    sampled <- indiv[, 1]
  }
  if (delimitIndiv) {
    keepers <- indiv$Me %in% sampled | indiv$Me %in% indiv$Mum | 
      indiv$Me %in% indiv$Dad
    indiv <- indiv[keepers, ]
  }
  ancestors <- matrix(data = sampled, nrow = length(sampled))
  parents.o <- foreach(i = 1:nrow(ancestors), .combine = rbind) %dopar% 
    {
      fishSim::parents(ancestors[i, 1], indiv)
    }
  print(paste("parents found at ", Sys.time(), sep = ""))
  grandparents.o <- foreach(i = 1:nrow(ancestors), .combine = rbind) %dopar% 
    {
      fishSim::grandparents(ancestors[i, 1], indiv)
    }
  print(paste("grandparents found at ", Sys.time()))
  ggrandparents.o <- foreach(i = 1:nrow(ancestors), .combine = rbind) %dopar% 
    {
      fishSim::great.grandparents(ancestors[i, 1], indiv)
    }
  print(paste("great-grandparents found at ", Sys.time(), 
              sep = ""))
  
  ancestors <- cbind(ancestors, parents.o, grandparents.o, ggrandparents.o) #, 
  
  colnames(ancestors) <- c("self", "father", "mother", 
                           "FF", "FM", "MF", "MM", 
                           "FFF", "FFM", "FMF", "FMM", "MFF",  "MFM", "MMF", "MMM") #, 
  expand.grid.unique <- function(x, y, include.equals = FALSE) {
    x <- unique(x)
    y <- unique(y)
    g <- function(i) {
      z <- setdiff(y, x[seq_len(i - include.equals)])
      if (length(z)) 
        cbind(x[i], z, deparse.level = 0)
    }
    do.call(rbind, lapply(seq_along(x), g))
  }
  print("made it")
  
  pairs <- expand.grid.unique(ancestors[, 1], ancestors[, 1])
  
  
  
  colnames(pairs) <- c("Var1", "Var2")
  related <- c(rep(NA, nrow(pairs)))
  totalRelatives <- c(rep(NA, nrow(pairs)))
  OneTwo <- c(rep(NA, nrow(pairs)))
  OneThree <- c(rep(NA, nrow(pairs)))
  OneFour <- c(rep(NA, nrow(pairs)))
  OneFive <- c(rep(NA, nrow(pairs)))
  OneSix <- c(rep(NA, nrow(pairs)))
  OneSeven <- c(rep(NA, nrow(pairs)))
  TwoTwo <- c(rep(NA, nrow(pairs)))
  TwoThree <- c(rep(NA, nrow(pairs)))
  TwoFour <- c(rep(NA, nrow(pairs)))
  TwoFive <- c(rep(NA, nrow(pairs)))
  TwoSix <- c(rep(NA, nrow(pairs)))
  TwoSeven <- c(rep(NA, nrow(pairs)))
  ThreeThree <- c(rep(NA, nrow(pairs)))
  ThreeFour <- c(rep(NA, nrow(pairs)))
  ThreeFive <- c(rep(NA, nrow(pairs)))
  ThreeSix <- c(rep(NA, nrow(pairs)))
  ThreeSeven <- c(rep(NA, nrow(pairs)))
  FourFour <- c(rep(NA, nrow(pairs)))
  FourFive <- c(rep(NA, nrow(pairs)))
  FourSix <- c(rep(NA, nrow(pairs)))
  FourSeven <- c(rep(NA, nrow(pairs)))
  FiveFive <- c(rep(NA, nrow(pairs)))
  FiveSix <- c(rep(NA, nrow(pairs)))
  FiveSeven <- c(rep(NA, nrow(pairs)))
  SixSix <- c(rep(NA, nrow(pairs)))
  SixSeven <- c(rep(NA, nrow(pairs)))
  SevenSeven <- c(rep(NA, nrow(pairs)))
  for (i in 1:length(related)) {
    allAncestors <- ancestors[ancestors[, 1] == pairs[i, 1], 1:15] %in% 
      ancestors[ancestors[, 1] == pairs[i, 2], 1:15]
    
    indi1 <- pairs[i, 1]
    indi2 <- pairs[i, 2]
    indi1Par <- ancestors[ancestors[, 1] == pairs[i, 1], 
                          2:3]
    indi2Par <- ancestors[ancestors[, 1] == pairs[i, 2], 
                          2:3]
    indi1GP <- ancestors[ancestors[, 1] == pairs[i, 1], 
                         4:7]
    indi2GP <- ancestors[ancestors[, 1] == pairs[i, 2], 
                         4:7]
    indi1GGP <- ancestors[ancestors[, 1] == pairs[i, 1], 
                          8:15]
    indi2GGP <- ancestors[ancestors[, 1] == pairs[i, 2], 
                          8:15]
    
    related[i] <- any(allAncestors)
    totalRelatives[i] <- sum(allAncestors)
    OneTwo[i] <- sum(c(indi1 %in% indi2Par, indi2 %in% indi1Par))
    OneThree[i] <- sum(c(indi1 %in% indi2GP, indi2 %in% 
                           indi1GP))
    OneFour[i] <- sum(c(indi1 %in% indi2GGP, indi2 %in% 
                          indi1GGP))
    TwoTwo[i] <- sum(c(indi1Par %in% indi2Par))
    TwoThree[i] <- sum(c(indi1Par %in% indi2GP, indi2Par %in% 
                           indi1GP))
    TwoFour[i] <- sum(c(indi1Par %in% indi2GGP, indi2Par %in% 
                          indi1GGP))
    ThreeThree[i] <- sum(c(indi1GP %in% indi2GP))
    ThreeFour[i] <- sum(c(indi1GP %in% indi2GGP, indi2GP %in% 
                            indi1GGP))
    FourFour[i] <- sum(c(indi1GGP %in% indi2GGP))
    
    if (i%%1000 == 0) {
      cat("\r", i, " of ", length(related), " comparisons", 
          sep = "")
      flush.console()
    }
  }
  pairs <- data.frame(pairs, related, totalRelatives, OneTwo, 
                      OneThree, OneFour, 
                      #OneFive, OneSix, OneSeven, 
                      TwoTwo, TwoThree, TwoFour, 
                      #TwoFive, TwoSix, TwoSeven, 
                      ThreeThree, ThreeFour, 
                      # ThreeFive, ThreeSix, ThreeSeven, 
                      FourFour #, 
                      # FourFive, FourSix, FourSeven, FiveFive, FiveSix, FiveSeven, 
                      # SixSix, SixSeven, SevenSeven
  )
  stopImplicitCluster()
  return(pairs)
}

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
  # p_geom <- 1 - 0.8455; max_age <- 19; sigma_l <- 2; l <- c(151, 145); a <- c(8, 9);
  
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



fLengthGivenAge <- function(l, a, sigma_l, vbgf_pars = c(l_inf = 175, 
                                                         k = 0.2,
                                                         a_0 = -2)) {
  expected_length <- vbgf_pars["l_inf"] * 
    (1 - exp(-vbgf_pars["k"] * (a - vbgf_pars["a_0"])))
  
  out <- pDiscreteNorm(x = l, mu = expected_length, sigma = sigma_l)
}

fSampledAge <- function(a, p_geom, max_age) {
  if (any(a < 0 | a > max_age)) {
    stop("all values for a have to be non-negative and at most equal to max_age")
  }
  
  return(dgeom(a, p_geom) / pgeom(max_age + 1, p_geom))
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


mateOrBirth <- function (indiv,
                         batchSize = 0.5, 
                         fecundityDist = "uniform", 
                         osr = c(0.5, 0.5), 
                         year = "-1", 
                         firstBreedFemale = 10, # only for females
                         firstBreedMale = 10,
                         firstLitter = NA,
                         type = "flat", 
                         maxClutch = Inf, 
                         singlePaternity = TRUE, 
                         exhaustFathers = FALSE, 
                         maturityCurve, 
                         maleCurve, 
                         femaleCurve,
                         no_gestation = TRUE,
                         quiet = TRUE) {
  if (!(type %in% c("flat", "age", "ageSex"))) {
    stop("'type' must be one of 'flat', 'age', or 'ageSex'.")
  }
  if (!(fecundityDist %in% c("poisson", "truncPoisson", 
                             "binomial", "uniform"))) {
    stop("'fecundityDist' must be one of 'poisson', 'truncPoisson', 'binomial', or 'uniform.")
  }
  
  if (no_gestation) {
    ## Subset the birthing females
    mothers <- subset(indiv, indiv[, 2] == 1 &   # are they female
                        indiv[, 8] >= firstBreedFemale & # are they old enough to breed
                        is.na(indiv[, 6]))       # are they alive
    if (nrow(mothers) == 0) {
      warning("There are no carrying females in the population")
    }
    ## Subset potential fathers
    fathers <- subset(indiv, indiv[, 2] == 0 & 
                        indiv[, 8] >= firstBreedMale & # only include males that were mature a year ago
                        is.na(indiv[, 6])) # either alive 
    
    if (nrow(fathers) == 0) {
      warning("There were no mature males in the population one year ago.")
    }
  } else {
    ## Subset the birthing females
    mothers <- subset(indiv, indiv[, 2] == 1 &   # are they female
                        indiv[, 8] > firstBreedFemale & # are they old enough to breed
                        is.na(indiv[, 6]) &        # are they alive
                        indiv[, 10] == 1)          # are they pregnant?
    if (nrow(mothers) == 0) {
      warning("There are no carrying females in the population")
    }
    ## Subset the mating females
    mating_females <- subset(indiv, indiv[, 2] == 1 &   # are they female
                               indiv[, 8] >= firstBreedFemale & # are they old enough to breed
                               is.na(indiv[, 6]) &        # are they alive
                               indiv[, 10] == 0)          # are they pregnant?
    if (nrow(mating_females) == 0)  {
      warning("There are no mating females in the population")
    }
    ## Subset potential fathers = all mature males that were alive and mature ONE YEAR AGO!
    fathers <- subset(indiv, indiv[, 2] == 0 & 
                        indiv[, 8] - 1 >= firstBreedMale & # only include males that were mature a year ago
                        (is.na(indiv[, 6]) | indiv[, 6] == year)) # either alive or were still alive one year ago
    
    if (nrow(fathers) == 0) {
      warning("There were no mature males in the population one year ago.")
    }
  }
  if (type == "ageSex") {
    if (no_gestation) {
      fathers <- fathers[runif(nrow(fathers)) < maleCurve[fathers[, 8] + 1], , 
                         drop = FALSE]
      mothers <- mothers[runif(nrow(mothers)) < femaleCurve[mothers[, 8] + 1], , 
                         drop = FALSE]
    }
    else {
      fathers <- fathers[runif(nrow(fathers)) < maleCurve[fathers[, 8] + 1], , 
                         drop = FALSE]
      mating_females <- mating_females[runif(nrow(mating_females)) < 
                                         femaleCurve[mating_females[, 8] + 1], , 
                                       drop = FALSE]
    }
    
    if (fecundityDist == "uniform") {
      if (length(batchSize) == 1) {
        clutch <- rep(batchSize, nrow(mothers))
      } else {
        clutch <- sample(batchSize, size = nrow(mothers), replace = TRUE)
      }
    }
    ## Here a line is added to set first litter at 1-2 pups, not 3-4. NOT USED
    # clutch[mothers$AgeLast == firstBreedFemale + 1] <- sample(
    #   firstLitter,
    #   size = sum(mothers$AgeLast == firstBreedFemale + 1),
    #   replace = TRUE)
    
    mothers <- subset(mothers, clutch > 0)
    
    ## Add/remove pregnancy indicators
    if (!no_gestation) {
      indiv$Pregnant[indiv$Me %in% mothers$Me] <- 0 # No longer pregnant
      indiv$Pregnant[indiv$Me %in% mating_females$Me] <- 1 # Mating females now pregnant
    }
    clutch <- clutch[clutch > 0]
  }
  clutch[clutch > maxClutch] <- maxClutch
  
  sprog.m <- CKMRcpp::createFounders(pop = 0) # this creates empty pop, but without column 10
  if (!no_gestation) {
    sprog.m$Pregnant <- integer(0) # add column 10 for pregnancy
  }
  
  for (s in unique(mothers[, 7])) {
    mothersInStock <- mothers[mothers[, 7] == s, , drop = FALSE]
    clutchInStock <- clutch[mothers[, 7] == s]
    fathersInStock <- fathers[fathers[, 7] == s, , drop = FALSE]
    if (nrow(fathersInStock) == 0) {
      warning(paste("There were no mature males in stock ", 
                    s, ", so ", nrow(mothersInStock), " mature females did not produce offspring", 
                    sep = ""))
      sprog.stock <- CKMRcpp::createFounders(pop = 0)
    }
    else if (nrow(fathersInStock > 0)) {
      n.sprogs <- sum(clutchInStock)
      sprog.stock <- data.frame(Me = integer(n.sprogs), 
                                Sex = integer(n.sprogs), Dad = integer(n.sprogs), 
                                Mum = integer(n.sprogs), BirthY = integer(n.sprogs), 
                                DeathY = integer(n.sprogs), Stock = integer(n.sprogs), 
                                AgeLast = integer(n.sprogs), SampY = integer(n.sprogs),
                                Pregnant = integer(n.sprogs))
      ticker <- 1
      for (m in 1:nrow(mothersInStock)) {
        if (nrow(fathersInStock) == 0) {
          warning(paste("All fathers in stock ", 
                        s, " are exhausted.", sep = ""))
        }
        else {
          sprog.stock[ticker:(ticker + clutchInStock[m] - 
                                1), 4] <- mothersInStock[m, 1]
          if (singlePaternity == TRUE) {
            sprog.stock[ticker:(ticker + clutchInStock[m] - 
                                  1), 3] <- fathersInStock[sample(1:nrow(fathersInStock), 
                                                                  1), 1]
          }
          else if (singlePaternity == FALSE) {
            if (nrow(fathersInStock) >= clutchInStock[m]) {
              sprog.stock[ticker:(ticker + clutchInStock[m] - 
                                    1), 3] <- fathersInStock[sample(1:nrow(fathersInStock), 
                                                                    clutchInStock[m]), 1]
            }
            else {
              sprog.stock[ticker:(ticker + nrow(fathersInStock) - 
                                    1), 3] <- fathersInStock[, 1]
            }
          }
          if (exhaustFathers == TRUE) {
            fathersInStock <- fathersInStock[!fathersInStock[, 
                                                             1] %in% sprog.stock[, 3], , drop = FALSE]
          }
          ticker <- ticker + clutchInStock[m]
        }
      }
    }
    sprog.stock <- sprog.stock[!is.na(sprog.stock[, 3]), 
                               , drop = FALSE]
    sprog.stock[, 1] <- -1 ## -1 is a placeholder; assign id's later
    sprog.stock[, 2] <- sample(c(0, 1), nrow(sprog.stock), # M=0, F=1
                               TRUE, prob = osr)
    sprog.stock[, 5] <- c(rep(year, nrow(sprog.stock)))
    sprog.stock[, 6] <- c(rep(NA, nrow(sprog.stock)))
    sprog.stock[, 7] <- as.integer(rep(s), nrow(sprog.stock))
    sprog.stock[, 8] <- c(rep(0, nrow(sprog.stock)))
    sprog.stock[, 9] <- c(rep(NA, nrow(sprog.stock)))
    if (!no_gestation) {
      sprog.stock[, 10] <- c(rep(0, nrow(sprog.stock)))
    }
    sprog.m <- data.table::rbindlist(list(sprog.m, sprog.stock))
  }
  names(sprog.m) <- names(indiv)
  sprog.m[, 1] <-  max(indiv$Me) + (1:nrow(sprog.m))
  indiv <- data.table::rbindlist(list(indiv, sprog.m))
  
  if (!quiet) cat("There were", nrow(mothers), "breeding mothers this cycle.\n")
  return(as.data.frame(indiv))
}

mateWithRecovery <- function (indiv,
                              recovery = 2, 
                              batchSize = 0.5, 
                              fecundityDist = "uniform", 
                              osr = c(0.5, 0.5), 
                              year = "-1", 
                              firstBreed = 12, # only for females
                              firstLitter = NA,
                              type = "flat", 
                              maxClutch = Inf, 
                              singlePaternity = TRUE, 
                              exhaustFathers = FALSE, 
                              maturityCurve, 
                              maleCurve, 
                              femaleCurve,
                              quiet = TRUE) {
  if (!(type %in% c("flat", "age", "ageSex"))) {
    stop("'type' must be one of 'flat', 'age', or 'ageSex'.")
  }
  if (!(fecundityDist %in% c("poisson", "truncPoisson", 
                             "binomial", "uniform"))) {
    stop("'fecundityDist' must be one of 'poisson', 'truncPoisson', 'binomial', or 'uniform.")
  }
  mothers <- subset(indiv, indiv[, 2] == "F" &   # are they female
                      indiv[, 8] >= firstBreed & # are they old enough to breed
                      is.na(indiv[, 6]) &        # are they alive
                      indiv[, 10] == 0)          # recovered from previous breeding cycle?
  if (nrow(mothers) == 0) 
    warning("There are no mature females in the population")
  fathers <- subset(indiv, indiv[, 2] == "M" & 
                      # indiv[, 8] >= firstBreed & # commented out so that males only uses their maturity curve
                      is.na(indiv[, 6]))
  if (nrow(fathers) == 0) 
    warning("There are no mature males in the population")
  # if (type == "flat") {
  #   if (fecundityDist == "poisson") {
  #     clutch <- rpois(n = nrow(mothers), lambda = batchSize)
  #   }
  #   if (fecundityDist == "truncPoisson") {
  #     clutch <- rTruncPoisson(n = nrow(mothers), T = batchSize)
  #   }
  #   if (fecundityDist == "binomial") {
  #     clutch <- rbinom(nrow(mothers), 1, prob = batchSize)
  #   }
  #   mothers <- subset(mothers, clutch > 0)
  #   clutch <- clutch[clutch > 0]
  # }
  # else if (type == "age") {
  #   mothers <- mothers[runif(nrow(mothers)) < maturityCurve[mothers[, 
  #                                                                   8] + 1], , drop = FALSE]
  #   fathers <- fathers[runif(nrow(fathers)) < maturityCurve[fathers[, 
  #                                                                   8] + 1], , drop = FALSE]
  #   if (fecundityDist == "poisson") {
  #     clutch <- rpois(n = nrow(mothers), lambda = batchSize)
  #   }
  #   if (fecundityDist == "truncPoisson") {
  #     clutch <- rTruncPoisson(n = nrow(mothers), T = batchSize)
  #   }
  #   if (fecundityDist == "binomial") {
  #     clutch <- rbinom(nrow(mothers), 1, prob = batchSize)
  #   }
  #   mothers <- subset(mothers, clutch > 0)
  #   clutch <- clutch[clutch > 0]
  # }
  else if (type == "ageSex") {
    mothers <- mothers[runif(nrow(mothers)) < femaleCurve[mothers[, 8] + 1], , 
                       drop = FALSE]
    fathers <- fathers[runif(nrow(fathers)) < maleCurve[fathers[, 
                                                                8] + 1], , drop = FALSE]
    if (fecundityDist == "poisson") {
      clutch <- rpois(n = nrow(mothers), lambda = batchSize)
    }
    if (fecundityDist == "truncPoisson") {
      clutch <- rTruncPoisson(n = nrow(mothers), T = batchSize)
    }
    if (fecundityDist == "binomial") {
      clutch <- rbinom(nrow(mothers), 1, prob = batchSize)
    }
    if (fecundityDist == "uniform") {
      clutch <- sample(batchSize, size = nrow(mothers), replace = TRUE)
    }
    ## Here a line is added to set first litter at 1-2 pups, not 3-4. 
    clutch[mothers$AgeLast == firstBreed] <- sample(firstLitter,
                                                    size = sum(mothers$AgeLast == firstBreed),
                                                    replace = TRUE)
    mothers <- subset(mothers, clutch > 0)
    indiv$Recovery[indiv$Me %in% mothers$Me] <- recovery # set recovery of mothers
    clutch <- clutch[clutch > 0]
  }
  clutch[clutch > maxClutch] <- maxClutch
  
  sprog.m <- CKMRcpp::createFounders(pop = 0) # this creates empty pop, but without column 10
  sprog.m$Recovery <- character(0) # add column 10 for recovery years remaining
  
  for (s in unique(mothers[, 7])) {
    mothersInStock <- mothers[mothers[, 7] == s, , drop = FALSE]
    clutchInStock <- clutch[mothers[, 7] == s]
    fathersInStock <- fathers[fathers[, 7] == s, , drop = FALSE]
    if (nrow(fathersInStock) == 0) {
      warning(paste("There were no mature males in stock ", 
                    s, ", so ", nrow(mothersInStock), " mature females did not produce offspring", 
                    sep = ""))
      sprog.stock <- CKMRcpp::createFounders(pop = 0)
    }
    else if (nrow(fathersInStock > 0)) {
      n.sprogs <- sum(clutchInStock)
      sprog.stock <- data.frame(Me = character(n.sprogs), 
                                Sex = character(n.sprogs), Dad = character(n.sprogs), 
                                Mum = character(n.sprogs), BirthY = integer(n.sprogs), 
                                DeathY = integer(n.sprogs), Stock = integer(n.sprogs), 
                                AgeLast = integer(n.sprogs), SampY = integer(n.sprogs))
      ticker <- 1
      for (m in 1:nrow(mothersInStock)) {
        if (nrow(fathersInStock) == 0) {
          warning(paste("All fathers in stock ", 
                        s, " are exhausted.", sep = ""))
        }
        else {
          sprog.stock[ticker:(ticker + clutchInStock[m] - 
                                1), 4] <- mothersInStock[m, 1]
          if (singlePaternity == TRUE) {
            sprog.stock[ticker:(ticker + clutchInStock[m] - 
                                  1), 3] <- fathersInStock[sample(1:nrow(fathersInStock), 
                                                                  1), 1]
          }
          else if (singlePaternity == FALSE) {
            if (nrow(fathersInStock) >= clutchInStock[m]) {
              sprog.stock[ticker:(ticker + clutchInStock[m] - 
                                    1), 3] <- fathersInStock[sample(1:nrow(fathersInStock), 
                                                                    clutchInStock[m]), 1]
            }
            else {
              sprog.stock[ticker:(ticker + nrow(fathersInStock) - 
                                    1), 3] <- fathersInStock[, 1]
            }
          }
          if (exhaustFathers == TRUE) {
            fathersInStock <- fathersInStock[!fathersInStock[, 
                                                             1] %in% sprog.stock[, 3], , drop = FALSE]
          }
          ticker <- ticker + clutchInStock[m]
        }
      }
    }
    sprog.stock <- sprog.stock[!is.na(sprog.stock[, 3]), 
                               , drop = FALSE]
    sprog.stock[, 1] <- ids::uuid(n = nrow(sprog.stock), drop_hyphens = TRUE)
    sprog.stock[, 2] <- sample(c("M", "F"), nrow(sprog.stock), 
                               TRUE, prob = osr)
    sprog.stock[, 5] <- c(rep(year, nrow(sprog.stock)))
    sprog.stock[, 6] <- c(rep(NA, nrow(sprog.stock)))
    sprog.stock[, 7] <- as.integer(rep(s), nrow(sprog.stock))
    sprog.stock[, 8] <- c(rep(0, nrow(sprog.stock)))
    sprog.stock[, 9] <- c(rep(NA, nrow(sprog.stock)))
    sprog.stock[, 10] <- c(rep(0, nrow(sprog.stock)))
    sprog.m <- rbind(sprog.m, sprog.stock)
  }
  names(sprog.m) <- names(indiv)
  indiv <- rbind(indiv, sprog.m)
  
  if (!quiet) cat("There were", nrow(mothers), "breeding mothers this cycle.\n")
  return(indiv)
}

pDiscreteNorm <- function(x, mu, sigma) {
  return(pnorm(x + 0.5, mean = mu, sd = sigma) - 
           pnorm(x - 0.5, mean = mu, sd = sigma))
}

plotCKMRabundance <- function(
    fits,           # A list of fits
    # par_name,     # The name of the abundance parameter
    year_lim,       # Number of years backward and forward from reference year
    max_y_axis = 6000,
    fixed_r = NULL, # 
    y0 = 140,       # The reference year
    med = F,     # TRUE for the median, else the mean
    alpha = 0.05,    # Signifance level for confidence intervals
    truth = NULL
) {
  
  ## Extract estimates for selected parameter
  est <- t(sapply(fits, function(x) x$par))
  
  ## Derive years and create matrices for abundance estimates
  years <- year_lim[1]:year_lim[2]
  abun_f <- matrix(NA, nrow = length(fits), ncol = length(years))
  abun_m <- matrix(NA, nrow = length(fits), ncol = length(years))
  
  ## For every fit, derive the abundance for all years
  if (is.null(fixed_r)) {
    for (i in 1:length(fits)) {
      abun_f[i, ] <- exp(est[i, "N_t0_f"]) * exp(est[i, "r"]) ^ years
      abun_m[i, ] <- exp(est[i, "N_t0_m"]) * exp(est[i, "r"]) ^ years
    }
  } else {
    for (i in 1:length(fits)) {
      abun_f[i, ] <- exp(est[i, "N_t0_f"]) * fixed_r ^ years
      abun_m[i, ] <- exp(est[i, "N_t0_m"]) * fixed_r ^ years
    }
  }
  
  ## Create matrices to store the median/mean abundance, and lower/upper bounds
  stat_f <- matrix(NA, nrow = length(years), ncol = 3)
  stat_m <- matrix(NA, nrow = length(years), ncol = 3)
  
  ## Extract mean/median, and CI
  if (med) {
    for (j in 1:length(years)) {
      stat_f[j, ] <- c(median = median(abun_f[, j]), 
                       lower = quantile(abun_f[, j], alpha / 2),
                       upper = quantile(abun_f[, j], 1 - alpha / 2))
      
      stat_m[j, ] <- c(median = median(abun_m[, j]), 
                       lower = quantile(abun_m[, j], alpha / 2),
                       upper = quantile(abun_m[, j], 1 - alpha / 2))
    }
  } else {
    for (j in 1:length(years)) {
      stat_f[j, ] <- c(mean = mean(abun_f[, j]), 
                       lower = quantile(abun_f[, j], alpha / 2),
                       upper = quantile(abun_f[, j], 1 - alpha / 2))
      
      stat_m[j, ] <- c(mean = mean(abun_m[, j]), 
                       lower = quantile(abun_m[, j], alpha / 2),
                       upper = quantile(abun_m[, j], 1 - alpha / 2))
    }
  }
  
  ## Set parameters to allow for two plots in one figure
  par(mfrow = c(1, 2));
  
  ## Create plot for the female adult population
  p_f <- matplot(x = y0 + years, 
                 y = stat_f, 
                 type = "l",
                 col = c("darkgreen"),
                 lty = c(1, 3, 3),
                 ylim = c(0, max_y_axis),
                 ylab = "Female adult abudance")
  ## If the true trend is provided, add it to the plot
  if (!is.null(truth)) {
    lines(x = years + y0, 
          y = truth$N_f * truth$r ^ years, 
          col = "red", 
          lty = 2)
  }
  
  ## Create plot for the female adult population
  p_m <- matplot(x = y0 + years, 
                 y = stat_m, 
                 type = "l",
                 col = c("purple"),
                 lty = c(1, 3, 3),
                 ylim = c(0, max_y_axis),
                 ylab = "Male adult abudance")
  ## If the true trend is provided, add it to the plot
  if (!is.null(truth)) {
    lines(x = years + y0, 
          y = truth$N_m * truth$r ^ years, 
          col = "red", 
          lty = 2)
  }
}

recover <- function(indiv) 
{
  ## Reduce 'Recovery' by one year for all living individuals
  indiv[is.na(indiv[, 6]), 10] <- indiv[is.na(indiv[, 6]), 10] - 1L
  ## Set all negative 'Recovery' to 0
  indiv[indiv[, 10] < 0, 10] <- 0
  ## Set Recovery for dead animals to 0
  indiv[!is.na(indiv[, 6]), 10] <- 0
  
  return(indiv)
}

retroCapture2 <- function (indiv, 
                           n = 1, 
                           year = "-1", 
                           fatal = FALSE) {
  ## check which individuals were alive at the sampling occasion
  is_alive <- (is.na(indiv$DeathY) | indiv$DeathY >= year) & indiv$BirthY <= year
  is_dead <- !is_alive
  n_alive <- sum(is_alive)
  n <- min(n, n_alive)
  
  if (n > 0) {
    sample_loc <- sample.int(n_alive, size = n)
    ## The line below is different from fishSim::capture(). 
    ## Ensure that we add a year to the ones that haven't been sampled before,
    ## and append a year to the ones sampled before. 
    sample_loc_new <- sample_loc[is.na(indiv[is_alive, ][sample_loc, 9])]
    sample_loc_previous <- sample_loc[!is.na(indiv[is_alive, ][sample_loc, 9])]
    
    indiv[is_alive, ][sample_loc_new, 9] <- year
    indiv[is_alive, ][sample_loc_previous, 9] <- 
      paste0(indiv[is_alive, ][sample_loc_previous, 9], "_", year)
    indiv[,  paste0("Samp", year)] <- NA
    indiv[is_alive, ][sample_loc, paste0("Samp", year)] <- year
    if (fatal) {
      indiv[is_alive, ][sample_loc_new, 6] <- year
    }
  }
  return(indiv)
}

retroCapture <- function (indiv, 
                          n = 1, 
                          year = "-1", 
                          fatal = FALSE) {
  ## check which individuals were alive at the sampling occasion
  is_alive <- (is.na(indiv$DeathY) | indiv$DeathY >= year) & indiv$BirthY <= year
  is_dead <- !is_alive
  n_alive <- sum(is_alive)
  n <- min(n, n_alive)
  
  if (n > 0) {
    sample_loc <- sample.int(n_alive, size = n)
    ## The line below is different from fishSim::capture(). 
    ## Ensure that we add a year to the ones that haven't been sampled before,
    ## and append a year to the ones sampled before. 
    sample_loc_new <- sample_loc[is.na(indiv[is_alive, ][sample_loc, 9])]
    sample_loc_previous <- sample_loc[!is.na(indiv[is_alive, ][sample_loc, 9])]
    indiv[is_alive, ][sample_loc_new, 9] <- year
    indiv[is_alive, ][sample_loc_previous, 9] <- 
      paste0(indiv[is_alive, ][sample_loc_previous, 9], "_", year)
    if (fatal) {
      indiv[is_alive, ][sample_loc_new, 6] <- year
    }
  }
  return(indiv)
}

uniformCheckGrowthrate <- function (
    fecundityDist = "uniform", 
    forceY1 = NA,      # year 1 mortality?
    mateType, # mating type
    mortType = "flat", # mortality constant or variable?
    batchSize, # expected litter size
    firstBreed = 0, # first year of breeding (at least min maturity curve)
    maxClutch = Inf, # max litter size
    osr = c(0.5, 0.5), # sex ratio
    maturityCurve, # equivalent to femaleCurve, but used if mateType == "age"
    femaleCurve, # female maturity curve
    maxAge = Inf, # self explanatory
    mortRate, # rate of mortality
    ageMort, # input for the mort() call
    stockMort, # input for the mort() call
    ageStockMort) { # input for the mort() call 
  ## Check inputs
  if (!(mateType %in% c("flat", "age", "ageSex"))) {
    stop("'mateType' must be one of 'flat', 'age', or 'ageSex'.")
  }
  if (!(mortType %in% c("flat", "age", "stock", "ageStock"))) {
    stop("'mortType' must be one of 'flat', 'age', 'stock', or 'ageStock'.")
  }
  if (missing(batchSize)) 
    stop("'batchSize' must be specified.")
  ## creating "batches" from a Poisson to run a simulation to get an expected
  ## batchSize. For uniform that can be much easier, as it is just the expectation,
  ## which is the mean. 
  if (batchSize != Inf && fecundityDist == "uniform") { # this is new!
    batchSize <- mean(batchSize)
  } else if (batchSize != Inf) { # this is old
    batches <- rpois(1e+06, lambda = batchSize)
    batchSize <- mean(batches[batches <= maxClutch]) 
  }
  
  
  ## mateType == "flat", not relevant to GR
  if (mateType == "flat") {
    if (mortType == "flat") {
      mat <- matrix(data = 0, nrow = length(0:firstBreed) + 
                      1, ncol = length(0:firstBreed) + 1)
      mat[1, ((2 + firstBreed):ncol(mat))] <- batchSize * 
        osr[2]
      for (i in 1:ncol(mat)) {
        if ((i + 1) <= nrow(mat)) 
          mat[i + 1, i] <- 1 - mortRate
      }
      mat[nrow(mat), ncol(mat)] <- 1 - mortRate
    }
    else if (mortType == "age") {
      mat <- matrix(data = 0, nrow = length(ageMort) + 
                      1, ncol = length(ageMort) + 1)
      mat[1, ((2 + firstBreed):ncol(mat))] <- batchSize * 
        osr[2]
      for (i in 1:ncol(mat)) {
        if ((i + 1) <= nrow(mat)) 
          mat[i + 1, i] <- 1 - ageMort[i]
      }
      mat[nrow(mat), ncol(mat)] <- 1 - ageMort[length(ageMort)]
    }
    else if (mortType == "stock") {
      mat <- matrix(data = 0, nrow = length(0:firstBreed) + 
                      1, ncol = length(0:firstBreed) + 1)
      mat.l <- lapply(seq_len(length(stockMort)), function(X) mat)
      for (s in 1:length(stockMort)) {
        mat.l[[s]][1, ((2 + firstBreed):ncol(mat.l[[s]]))] <- batchSize * 
          osr[2]
        for (i in 1:ncol(mat.l[[s]])) {
          if ((i + 1) <= nrow(mat.l[[s]])) 
            mat.l[[s]][i + 1, i] <- 1 - stockMort[s]
        }
        mat.l[[s]][nrow(mat.l[[s]]), ncol(mat.l[[s]])] <- 1 - 
          stockMort[s]
      }
    }
    else if (mortType == "ageStock") {
      mat <- matrix(data = 0, nrow = length(ageStockMort[, 
                                                         1]) + 1, ncol = length(ageStockMort[, 1]) + 
                      1)
      mat.l <- lapply(seq_len(ncol(ageStockMort)), function(X) mat)
      for (s in 1:ncol(ageStockMort)) {
        mat.l[[s]][1, ((2 + firstBreed):ncol(mat.l[[s]]))] <- batchSize * 
          osr[2]
        for (i in 1:ncol(mat.l[[s]])) {
          if ((i + 1) <= nrow(mat.l[[s]])) 
            mat.l[[s]][i + 1, i] <- 1 - ageStockMort[i, 
                                                     s]
        }
        mat.l[[s]][nrow(mat.l[[s]]), ncol(mat.l[[s]])] <- 1 - 
          ageStockMort[nrow(ageStockMort), s]
      }
    }
  }
  ## mateType == "age" also not relevant to GR
  else if (mateType == "age") {
    if (firstBreed > 0) 
      maturityCurve[1:(firstBreed - 1)] <- 0
    if (mortType == "flat") {
      mat <- matrix(data = 0, nrow = length(maturityCurve), 
                    ncol = length(maturityCurve))
      mat[1, ] <- maturityCurve * batchSize * osr[2]
      for (i in 1:ncol(mat)) {
        if ((i + 1) <= nrow(mat)) 
          mat[i + 1, i] <- 1 - mortRate
      }
      mat[nrow(mat), ncol(mat)] <- 1 - mortRate
    }
    else if (mortType == "age") {
      mat <- matrix(data = 0, nrow = max(c(length(maturityCurve), 
                                           length(ageMort))), ncol = max(c(length(maturityCurve), 
                                                                           length(ageMort))))
      mat[1, (1:length(maturityCurve))] <- maturityCurve * 
        batchSize * osr[2]
      mat[1, (length(maturityCurve):ncol(mat))] <- maturityCurve[length(maturityCurve)] * 
        batchSize * osr[2]
      for (i in 1:length(ageMort)) {
        if ((i + 1) <= nrow(mat)) 
          mat[i + 1, i] <- 1 - ageMort[i]
      }
      for (i in length(ageMort):ncol(mat)) {
        if ((i + 1) <= nrow(mat)) 
          mat[i + 1, i] <- 1 - ageMort[length(ageMort)]
      }
      mat[nrow(mat), ncol(mat)] <- 1 - ageMort[length(ageMort)]
    }
    else if (mortType == "stock") {
      mat <- matrix(data = 0, nrow = length(maturityCurve), 
                    ncol = length(maturityCurve))
      mat.l <- lapply(seq_len(length(stockMort)), function(X) mat)
      for (s in 1:length(stockMort)) {
        mat.l[[s]][1, (1:length(maturityCurve))] <- maturityCurve * 
          batchSize * osr[2]
        mat.l[[s]][1, (length(maturityCurve):ncol(mat.l[[s]]))] <- maturityCurve[length(maturityCurve)] * 
          batchSize * osr[2]
        for (i in 1:ncol(mat.l[[s]])) {
          if ((i + 1) <= nrow(mat.l[[s]])) 
            mat.l[[s]][i + 1, i] <- 1 - stockMort[s]
        }
        mat.l[[s]][nrow(mat.l[[s]]), ncol(mat.l[[s]])] <- 1 - 
          stockMort[s]
      }
    }
    else if (mortType == "ageStock") {
      mat <- matrix(data = 0, nrow = max(c(length(maturityCurve), 
                                           nrow(ageStockMort))), ncol = max(c(length(maturityCurve), 
                                                                              nrow(ageStockMort))))
      mat.l <- lapply(seq_len(ncol(ageStockMort)), function(X) mat)
      for (s in 1:ncol(ageStockMort)) {
        mat.l[[s]][1, (1:length(maturityCurve))] <- maturityCurve * 
          batchSize * osr[2]
        mat.l[[s]][1, (length(maturityCurve):ncol(mat.l[[s]]))] <- maturityCurve[length(maturityCurve)] * 
          batchSize * osr[2]
        for (i in 1:ncol(mat.l[[s]])) {
          if ((i + 1) <= nrow(mat.l[[s]])) 
            mat.l[[s]][i + 1, i] <- 1 - ageStockMort[i, 
                                                     s]
        }
        for (i in nrow(ageStockMort):ncol(mat.l[[s]])) {
          if ((i + 1) <= nrow(mat.l[[s]])) 
            mat.l[[s]][i + 1, i] <- 1 - ageStockMort[nrow(ageMort), 
                                                     s]
        }
        mat.l[[s]][nrow(mat.l[[s]]), ncol(mat.l[[s]])] <- 1 - 
          ageStockMort[nrow(ageStockMort), s]
      }
    }
  }
  ## This is the one we are interested in!
  else if (mateType == "ageSex") {
    if (mortType == "flat") {   # mortType is assumed flat for GR
      mat <- matrix(data = 0, nrow = length(femaleCurve), 
                    ncol = length(femaleCurve))
      mat[1, ] <- femaleCurve * batchSize * osr[2]
      for (i in 1:ncol(mat)) {
        if ((i + 1) <= nrow(mat)) 
          mat[i + 1, i] <- 1 - mortRate
      }
      mat[nrow(mat), ncol(mat)] <- 1 - mortRate
    }
    else if (mortType == "age") {
      mat <- matrix(data = 0, nrow = max(c(length(femaleCurve), 
                                           length(ageMort))), ncol = max(c(length(femaleCurve), 
                                                                           length(ageMort))))
      mat[1, (1:length(femaleCurve))] <- femaleCurve * 
        batchSize * osr[2]
      mat[1, (length(femaleCurve):ncol(mat))] <- femaleCurve[length(femaleCurve)] * 
        batchSize * osr[2]
      for (i in 1:length(ageMort)) {
        if ((i + 1) <= nrow(mat)) 
          mat[i + 1, i] <- 1 - ageMort[i]
      }
      for (i in length(ageMort):ncol(mat)) {
        if ((i + 1) <= nrow(mat)) 
          mat[i + 1, i] <- 1 - ageMort[length(ageMort)]
      }
      mat[nrow(mat), ncol(mat)] <- 1 - ageMort[length(ageMort)]
    }
    else if (mortType == "stock") {
      mat <- matrix(data = 0, nrow = length(femaleCurve), 
                    ncol = length(femaleCurve))
      mat.l <- lapply(seq_len(length(stockMort)), function(X) mat)
      for (s in 1:length(stockMort)) {
        mat.l[[s]][1, (1:length(femaleCurve))] <- femaleCurve * 
          batchSize * osr[2]
        mat.l[[s]][1, (length(femaleCurve):ncol(mat.l[[s]]))] <- femaleCurve[length(femaleCurve)] * 
          batchSize * osr[2]
        for (i in 1:ncol(mat.l[[s]])) {
          if ((i + 1) <= nrow(mat.l[[s]])) 
            mat.l[[s]][i + 1, i] <- 1 - stockMort[s]
        }
        mat.l[[s]][nrow(mat.l[[s]]), ncol(mat.l[[s]])] <- 1 - 
          stockMort[s]
      }
    }
    else if (mortType == "ageStock") {
      mat <- matrix(data = 0, nrow = max(c(length(femaleCurve), 
                                           nrow(ageStockMort))), ncol = max(c(length(femaleCurve), 
                                                                              nrow(ageStockMort))))
      mat.l <- lapply(seq_len(ncol(ageStockMort)), function(X) mat)
      for (s in 1:ncol(ageStockMort)) {
        mat.l[[s]][1, (1:length(femaleCurve))] <- femaleCurve * 
          batchSize * osr[2]
        mat.l[[s]][1, (length(femaleCurve):ncol(mat.l[[s]]))] <- femaleCurve[length(femaleCurve)] * 
          batchSize * osr[2]
        for (i in 1:ncol(mat.l[[s]])) {
          if ((i + 1) <= nrow(mat.l[[s]])) 
            mat.l[[s]][i + 1, i] <- 1 - ageStockMort[i, 
                                                     s]
        }
        for (i in nrow(ageStockMort):ncol(mat.l[[s]])) {
          if ((i + 1) <= nrow(mat.l[[s]])) 
            mat.l[[s]][i + 1, i] <- 1 - ageStockMort[nrow(ageMort), 
                                                     s]
        }
        mat.l[[s]][nrow(mat.l[[s]]), ncol(mat.l[[s]])] <- 1 - 
          ageStockMort[nrow(ageStockMort), s]
      }
    }
  }
  if (!is.na(forceY1)) {
    if (mortType %in% c("flat", "age")) 
      mat[2, 1] <- 1 - forceY1
    if (mortType %in% c("stock", "ageStock")) 
      for (i in 1:length(mat.l)) mat.l[[i]][2, 1] <- 1 - 
          forceY1
  }
  if (mortType %in% c("flat", "age")) {
    return(eigen(mat)$values[1])
  }
  if (mortType %in% c("stock", "ageStock")) {
    outs <- c(rep(NA, length(mat.l)))
    for (i in 1:length(outs)) outs[i] <- eigen(mat.l[[i]])$values[1]
    return(outs)
  }
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
vbgf <- function(a, a_0 = -3.5, k = 0.1, l_inf = 175) {
  return(l_inf * (1 - exp(-k * (a - a_0))))
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
invvbgf <- function(l, a_0 = -3.5, k = 0.1, l_inf = 175) {
  return(a_0 - log(1 - l/l_inf) / k)
}
