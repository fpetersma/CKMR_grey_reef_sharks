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
                            y0 = 0,
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
    indiv[, 5] <- 1L - indiv[, 8] + y0 # add start year y0 to match the simulation
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
    ## Get all sample occassions of the same individual with id
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
                             sex = NULL,
                             only_index = FALSE) {
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
    old_enough <- alive & (indiv$BirthY < (year - min_age + 1))
  } else {
    # Born in or before 'year', and either still alive or died after 'year' 
    alive <- indiv$BirthY <= year &                    # Born in or before 'year'?
      ((indiv$DeathY > year & !is.na(indiv$DeathY)) |  # Died after 'year'?
         is.na(indiv$DeathY))                          # Still alive?
    # Old enough at the endof the year
    old_enough <-  alive & indiv$BirthY <= (year - min_age + 1)
  }

  
  ## If requested, only return the index
  if (only_index) {
    living <- old_enough
  ## Else, return the living from indiv with all information included
  } else {
    living <- subset(indiv, old_enough)
  }
  
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
  TwoTwo <- c(rep(NA, nrow(pairs)))        # For HSPs and FSPs

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
    # related[i] <- is.integer(any(allAncestors))## gone for speed
    # totalRelatives[i] <- sum(allAncestors) ## gone for speed
    OneTwo[i] <- sum(c(indi1 %in% indi2Par, indi2 %in% indi1Par)) # For POP
    TwoTwo[i] <- sum(c(indi1Par %in% indi2Par)) # For HSP and FSP

    
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
                      TwoTwo
                      )
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
    
    ## I WONDER IF WHAT I DO IS RIGHT HERE... I DONT THINK I SHOULD LOOK AT AGE LAST, BUT INSTEAD 
    ## COMPARE TO BIRTHYEAR
    fathers <- CKMRcpp::extractTheLiving(indiv, year = year - 1, start_of_year = T, 
                                         min_age = firstBreedMale, sex = "male")
    # cat("number of fathers step 1:", nrow(fathers), "\n")
    # fathers <- subset(indiv, indiv[, 2] == 0 & 
    #                     year - indiv[, 5] - 1 >= firstBreedMale &  ## Using birthyear is important, instead of AgeLast, which stops counting after death
    #                     # indiv[, 8] - 1 >= firstBreedMale & # only include males that were mature a year ago
    #                     (is.na(indiv[, 6]) | indiv[, 6] == year)) # either still alive or alive one year ago
    
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
    # cat("number of fathers step 2:", nrow(fathers), "\n")
    
    
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
      if (!no_gestation) {
        sprog.stock <- data.frame(Me = integer(n.sprogs), 
                                  Sex = integer(n.sprogs), Dad = integer(n.sprogs), 
                                  Mum = integer(n.sprogs), BirthY = integer(n.sprogs), 
                                  DeathY = integer(n.sprogs), Stock = integer(n.sprogs), 
                                  AgeLast = integer(n.sprogs), SampY = integer(n.sprogs),
                                  Pregnant = integer(n.sprogs))
      } else {
        sprog.stock <- data.frame(Me = integer(n.sprogs), 
                                  Sex = integer(n.sprogs), Dad = integer(n.sprogs), 
                                  Mum = integer(n.sprogs), BirthY = integer(n.sprogs), 
                                  DeathY = integer(n.sprogs), Stock = integer(n.sprogs), 
                                  AgeLast = integer(n.sprogs), SampY = integer(n.sprogs))
      }
      
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

pDiscreteNorm <- function(x, mu, sigma) {
  return(pnorm(x + 0.5, mean = mu, sd = sigma) - 
           pnorm(x - 0.5, mean = mu, sd = sigma))
}

plotCKMRabundance <- function(
    fits_list,           # A list of simulation objects 
    # par_name,     # The name of the abundance parameter
    year_lim,       # Number of years backward and forward from reference year
    max_y_axis = 3000,
    sex = "both",  # male, female, or both 
    y0 = 2014,       # The reference year
    truth = NULL,
    y_axis = "Abundance",
    share_r = FALSE # do males and females share the growth parameter, or not
) {
  ## load libraries
  library(tidyverse, quietly = T, warn.conflicts = F)
  
  for (i in 1:length(fits_list)) {
    fits <- fits_list[[i]]
    
    ## Extract estimates for selected parameter
    est <- t(sapply(fits, function(x) x$par))
    conv <- sapply(fits, function(x) x$message)
    
    ## Derive years and create matrices for abundance estimates
    years <- year_lim[1]:year_lim[2]
    abun_f <- matrix(NA, nrow = length(fits), ncol = length(years))
    abun_m <- matrix(NA, nrow = length(fits), ncol = length(years))
    
    ## For every fit, derive the abundance for all years
    for (j in 1:length(fits)) {
      abun_f[j, ] <- exp(est[j, "N_t0_f"]) * exp(est[j, "r_f"]) ^ years
      abun_m[j, ] <- exp(est[j, "N_t0_m"]) * exp(est[j, "r_m"]) ^ years
    }
    
    ## Turn the data in long format so ggplot2 knows what to do
    abun_f_df <- cbind(data.frame(year = y0 + years), 
                       t(abun_f),
                       "101"= Rfast::colMedians(abun_f), #median
                       "102" = truth[years + nrow(truth), 2]) #truth
    abun_f_long <- reshape2::melt(abun_f_df, value.name = "N", id.vars = "year")
    abun_f_long$sim_id <- i
    abun_f_long$sex <- "F"
    
    ## Turn the data in long format so ggplot2 knows what to do
    abun_m_df <- cbind(data.frame(year = y0 + years), 
                       t(abun_m),
                       "101"= Rfast::colMedians(abun_m), #median
                       "102" = truth[years + nrow(truth), 1]) #truth
    abun_m_long <- reshape2::melt(abun_m_df, value.name = "N", id.vars = "year")
    abun_m_long$sim_id <- i
    abun_m_long$sex <- "M"
    
    if (any(conv != "relative convergence (4)")) {
      abun_m_long$conv <- "failed"    
      abun_f_long$conv <- "failed" 
    } else {
      abun_m_long$conv <- "successful"    
      abun_f_long$conv <- "successful" 
    }
    
    # Combine the male and datasets
    # If it's the first simulation scenario, create abun_long; else, add to abun_long
    if (i == 1) {
      abun_long <- rbind(abun_f_long, abun_m_long)
    } else {
      abun_long <- rbind(abun_long, abun_f_long, abun_m_long)
    }
  }
  ## Here starts the ggplot2 magic ----------------------------------------- <<<
  
  abun_long$N[abun_long$conv == "failed"] <- NA
  
  dimension <- sqrt(max(abun_long$sim_id))
  
  abun_long$sim_id <- as.factor(abun_long$sim_id)
  levels(abun_long$sim_id) <- paste0(rep(1:5, each = 5), "-", 1:5)
  
  if (sex == "both") {
    # Create the plot with 100 fitted abundances, and a mean line
    p <- ggplot(abun_long) +
      geom_line(
        mapping = aes(x = year, y = N, colour = variable),
        show.legend = F, linewidth = 1) +
      scale_color_manual(values = c(rep(alpha("darkgrey", 0.2), 100),
                                    alpha("black", 0.5),
                                    alpha("red", 0.4))) +
      theme_bw() +
      ylab(y_axis ) +
      xlab("Year") +
      facet_wrap(~ sim_id + sex, nrow = dimension, labeller = label_parsed) + 
      coord_cartesian(ylim=c(0, max_y_axis)) +
      scale_x_continuous(breaks = seq(from = min(years) + y0, to = max(years) + y0, by = 5),
                         labels = seq(from = min(years) + y0, to = max(years) + y0, by = 5))
    
  } else if (sex == "male") {
    p <- ggplot(subset(abun_long, sex == "M")) +
      geom_line(
        mapping = aes(x = year, y = N, colour = variable),
        show.legend = F, linewidth = 1) +
      scale_color_manual(values = c(rep(alpha("darkgrey", 0.2), 100),
                                    alpha("black", 0.5),
                                    alpha("red", 0.4))) +
      theme_bw() +
      ylab(y_axis ) +
      xlab("Year") +
      facet_wrap(~ sim_id , nrow = dimension, labeller = label_parsed) + 
      coord_cartesian(ylim=c(0,max_y_axis)) +
      scale_x_continuous(breaks = seq(from = min(years) + y0, to = max(years) + y0, by = 5),
                         labels = seq(from = min(years) + y0, to = max(years) + y0, by = 5))
    
  } else {
    p <- ggplot(subset(abun_long, sex == "F")) +
      geom_line(
        mapping = aes(x = year, y = N, colour = variable),
        show.legend = F, linewidth = 1) +
      scale_color_manual(values = c(rep(alpha("darkgrey", 0.2), 100),
                                    alpha("black", 0.5),
                                    alpha("red", 0.4))) +
      theme_bw() +
      ylab(y_axis ) +
      xlab("Year") +
      facet_wrap(~ sim_id , nrow = dimension, labeller = label_parsed) + 
      coord_cartesian(ylim=c(0,max_y_axis)) +
      scale_x_continuous(breaks = seq(from = min(years) + y0, to = max(years) + y0, by = 5),
                         labels = seq(from = min(years) + y0, to = max(years) + y0, by = 5))
    
  }
  return(p)
}

retroCapture2 <- function (indiv, 
                           n = 1, 
                           year = "-1", 
                           fatal = FALSE) {
  ## check which individuals were alive at the sampling occasion
  the_living <- CKMRcpp::extractTheLiving(indiv, year, TRUE, 0, NULL)
  is_alive <- indiv$Me %in% the_living$Me
  
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
