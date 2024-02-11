## ========================================================================== ##
## A new file with custom functions to add to fishSim to make it more similar ##
## to the grey reef shark case study from Palmyra. Recovered after            ##
## accidentally replacing the file. Very stupid.                              ##  
## ========================================================================== ##

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
  
  is_pregnant <- indiv$Sex == "F" & indiv$AgeLast %in% (matingAges + 1)
  
  indiv$Pregnant <- is_pregnant
  
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
                                verbose = TRUE,
                                delimitIndiv = TRUE) {
  library(doParallel)
  
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
  ancestors <- matrix(data = c(sampled, rep("blanks", 
                                            length(sampled) * 14)), 
                      nrow = length(sampled))
  colnames(ancestors) <- c("self", "father", "mother", 
                           "FF", "FM", "MF", "MM", "FFF", 
                           "FFM", "FMF", "FMM", "MFF", "MFM", 
                           "MMF", "MMM")
  for (i in 1:nrow(ancestors)) {
    ancestors[i, 2:3] <- fishSim::parents(ancestors[i, 1], indiv)
    ancestors[i, 4:7] <- fishSim::grandparents(ancestors[i, 1], indiv)
    ancestors[i, 8:15] <- fishSim::great.grandparents(ancestors[i, 1], indiv)
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
  OneTwo <- c(rep(NA, nrow(pairs)))
  OneThree <- c(rep(NA, nrow(pairs)))
  OneFour <- c(rep(NA, nrow(pairs)))
  TwoTwo <- c(rep(NA, nrow(pairs)))
  TwoThree <- c(rep(NA, nrow(pairs)))
  TwoFour <- c(rep(NA, nrow(pairs)))
  ThreeThree <- c(rep(NA, nrow(pairs)))
  ThreeFour <- c(rep(NA, nrow(pairs)))
  FourFour <- c(rep(NA, nrow(pairs)))
  for (i in 1:length(related)) {
    allAncestors <- ancestors[ancestors[, 1] == pairs[i, 1], 1:15] %in% 
      ancestors[ancestors[, 1] == pairs[i, 2], 1:15]
    indi1 <- pairs[i, 1]
    indi2 <- pairs[i, 2]
    indi1Par <- ancestors[ancestors[, 1] == pairs[i, 1], 
                          2:3]
    indi2Par <- ancestors[ancestors[, 1] == pairs[i, 2], 
                          2:3]
    indi1GP <- ancestors[ancestors[, 1] == pairs[i, 1], 4:7]
    indi2GP <- ancestors[ancestors[, 1] == pairs[i, 2], 4:7]
    indi1GGP <- ancestors[ancestors[, 1] == pairs[i, 1], 
                          8:15]
    indi2GGP <- ancestors[ancestors[, 1] == pairs[i, 2], 
                          8:15]
    related[i] <- any(allAncestors)
    totalRelatives[i] <- sum(allAncestors)
    OneTwo[i] <- sum(c(indi1 %in% indi2Par, indi2 %in% indi1Par))
    OneThree[i] <- sum(c(indi1 %in% indi2GP, indi2 %in% indi1GP))
    OneFour[i] <- sum(c(indi1 %in% indi2GGP, indi2 %in% indi1GGP))
    
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
      cat("\r", i, " of ", length(related), 
          " comparisons", sep = "")
      flush.console()
    }
  }
  pairs <- data.frame(pairs, related, totalRelatives, OneTwo, 
                      OneThree, OneFour, 
                      TwoTwo, TwoThree, TwoFour, 
                      ThreeThree, ThreeFour, 
                      FourFour)
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
  # gggrandparents.o <- foreach(i = 1:nrow(ancestors), .combine = rbind) %dopar% 
  #   {
  #     great2.grandparents(ancestors[i, 1], indiv)
  #   }
  # print(paste("great-great-grandparents found at ", Sys.time(), 
  #             sep = ""))
  # ggggrandparents.o <- foreach(i = 1:nrow(ancestors), .combine = rbind) %dopar% 
  #   {
  #     great3.grandparents(ancestors[i, 1], indiv)
  #   }
  # print(paste("great-great-great-grandparents found at ", 
  #             Sys.time(), sep = ""))
  # gggggrandparents.o <- foreach(i = 1:nrow(ancestors), .combine = rbind) %dopar% 
  #   {
  #     great4.grandparents(ancestors[i, 1], indiv)
  #   }
  # print(paste("great-great-great-great-grandparents found at ", 
  #             Sys.time(), sep = ""))
  ancestors <- cbind(ancestors, parents.o, grandparents.o, ggrandparents.o) #, 
  # gggrandparents.o, ggggrandparents.o, 
  # gggggrandparents.o)
  colnames(ancestors) <- c("self", "father", "mother", 
                           "FF", "FM", "MF", "MM", 
                           "FFF", "FFM", "FMF", "FMM", "MFF",  "MFM", "MMF", "MMM") #, 
  # "FFFF", "FFFM", "FFMF", "FFMM", 
  # "FMFF", "FMFM", "FMMF", "FMMM", "MFFF", "MFFM", "MFMF", 
  # "MFMM", "MMFF", "MMFM", "MMMF", "MMMM", "FFFFF", "FFFFM", 
  # "FFFMF", "FFFMM", "FFMFF", "FFMFM", "FFMMF", "FFMMM", 
  # "FMFFF", "FMFFM", "FMFMF", "FMFMM", "FMMFF", "FMMFM", 
  # "FMMMF", "FMMMM", "MFFFF", "MFFFM", "MFFMF", "MFFMM", 
  # "MFMFF", "MFMFM", "MFMMF", "MFMMM", "MMFFF", "MMFFM", 
  # "MMFMF", "MMFMM", "MMMFF", "MMMFM", "MMMMF", "MMMMM", 
  # "FFFFFF", "FFFFFM", "FFFFMF", "FFFFMM", "FFFMFF", "FFFMFM", 
  # "FFFMMF", "FFFMMM", "FFMFFF", "FFMFFM", "FFMFMF", "FFMFMM", 
  # "FFMMFF", "FFMMFM", "FFMMMF", "FFMMMM", "FMFFFF", "FMFFFM", 
  # "FMFFMF", "FMFFMM", "FMFMFF", "FMFMFM", "FMFMMF", "FMFMMM", 
  # "FMMFFF", "FMMFFM", "FMMFMF", "FMMFMM", "FMMMFF", "FMMMFM", 
  # "FMMMMF", "FMMMMM", "MFFFFF", "MFFFFM", "MFFFMF", "MFFFMM", 
  # "MFFMFF", "MFFMFM", "MFFMMF", "MFFMMM", "MFMFFF", "MFMFFM", 
  # "MFMFMF", "MFMFMM", "MFMMFF", "MFMMFM", "MFMMMF", "MFMMMM", 
  # "MMFFFF", "MMFFFM", "MMFFMF", "MMFFMM", "MMFMFF", "MMFMFM", 
  # "MMFMMF", "MMFMMM", "MMMFFF", "MMMFFM", "MMMFMF", "MMMFMM", 
  # "MMMMFF", "MMMMFM", "MMMMMF", "MMMMMM")
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
    # allAncestors <- ancestors[ancestors[, 1] == pairs[i, 
    #                                                   1], 1:127] %in% ancestors[ancestors[, 1] == pairs[i, 
    #                                                                                                     2], 1:127]
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
    # indi1GGGP <- ancestors[ancestors[, 1] == pairs[i, 1], 
    #                        16:31]
    # indi2GGGP <- ancestors[ancestors[, 1] == pairs[i, 2], 
    #                        16:31]
    # indi1GGGGP <- ancestors[ancestors[, 1] == pairs[i, 1], 
    #                         32:63]
    # indi2GGGGP <- ancestors[ancestors[, 1] == pairs[i, 2], 
    #                         32:63]
    # indi1GGGGGP <- ancestors[ancestors[, 1] == pairs[i, 
    #                                                  1], 64:127]
    # indi2GGGGGP <- ancestors[ancestors[, 1] == pairs[i, 
    #                                                  2], 64:127]
    related[i] <- any(allAncestors)
    totalRelatives[i] <- sum(allAncestors)
    OneTwo[i] <- sum(c(indi1 %in% indi2Par, indi2 %in% indi1Par))
    OneThree[i] <- sum(c(indi1 %in% indi2GP, indi2 %in% 
                           indi1GP))
    OneFour[i] <- sum(c(indi1 %in% indi2GGP, indi2 %in% 
                          indi1GGP))
    # OneFive[i] <- sum(c(indi1 %in% indi2GGGP, indi2 %in% 
    #                       indi1GGGP))
    # OneSix[i] <- sum(c(indi1 %in% indi2GGGGP, indi2 %in% 
    #                      indi1GGGGP))
    # OneSeven[i] <- sum(c(indi1 %in% indi2GGGGGP, indi2 %in% 
    #                        indi1GGGGGP))
    TwoTwo[i] <- sum(c(indi1Par %in% indi2Par))
    TwoThree[i] <- sum(c(indi1Par %in% indi2GP, indi2Par %in% 
                           indi1GP))
    TwoFour[i] <- sum(c(indi1Par %in% indi2GGP, indi2Par %in% 
                          indi1GGP))
    # TwoFive[i] <- sum(c(indi1Par %in% indi2GGGP, indi2Par %in% 
    #                       indi1GGGP))
    # TwoSix[i] <- sum(c(indi1Par %in% indi2GGGGP, indi2Par %in% 
    #                      indi1GGGGP))
    # TwoSeven[i] <- sum(c(indi1Par %in% indi2GGGGGP, indi2Par %in% 
    #                        indi1GGGGGP))
    ThreeThree[i] <- sum(c(indi1GP %in% indi2GP))
    ThreeFour[i] <- sum(c(indi1GP %in% indi2GGP, indi2GP %in% 
                            indi1GGP))
    # ThreeFive[i] <- sum(c(indi1GP %in% indi2GGGP, indi2GP %in% 
    #                         indi1GGGP))
    # ThreeSix[i] <- sum(c(indi1GP %in% indi2GGGGP, indi2GP %in% 
    #                        indi1GGGGP))
    # ThreeSeven[i] <- sum(c(indi1GP %in% indi2GGGGGP, indi2GP %in% 
    #                          indi1GGGGGP))
    FourFour[i] <- sum(c(indi1GGP %in% indi2GGP))
    # FourFive[i] <- sum(c(indi1GGP %in% indi2GGGP, indi2GGP %in% 
    #                        indi1GGGP))
    # FourSix[i] <- sum(c(indi1GGP %in% indi2GGGGP, indi2GGP %in% 
    #                       indi1GGGGP))
    # FourSeven[i] <- sum(c(indi1GGP %in% indi2GGGGGP, indi2GGP %in% 
    #                         indi1GGGGGP))
    # FiveFive[i] <- sum(c(indi1GGGP %in% indi2GGGP))
    # FiveSix[i] <- sum(c(indi1GGGP %in% indi2GGGGP, indi2GGGP %in% 
    #                       indi1GGGGP))
    # FiveSeven[i] <- sum(c(indi1GGGP %in% indi2GGGGGP, indi2GGGP %in% 
    #                         indi1GGGGGP))
    # SixSix[i] <- sum(c(indi1GGGGP %in% indi2GGGGP))
    # SixSeven[i] <- sum(c(indi1GGGGP %in% indi2GGGGGP, indi2GGGGP %in% 
    #                        indi1GGGGGP))
    # SevenSeven[i] <- sum(c(indi1GGGGGP %in% indi2GGGGGP))
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
    mothers <- subset(indiv, indiv[, 2] == "F" &   # are they female
                        indiv[, 8] >= firstBreedFemale & # are they old enough to breed
                        is.na(indiv[, 6]))       # are they alive
    if (nrow(mothers) == 0) {
      warning("There are no carrying females in the population")
    }
    ## Subset potential fathers
    fathers <- subset(indiv, indiv[, 2] == "M" & 
                        indiv[, 8] >= firstBreedMale & # only include males that were mature a year ago
                        is.na(indiv[, 6])) # either alive 
    
    if (nrow(fathers) == 0) {
      warning("There were no mature males in the population one year ago.")
    }
  } else {
    ## Subset the birthing females
    mothers <- subset(indiv, indiv[, 2] == "F" &   # are they female
                        indiv[, 8] > firstBreedFemale & # are they old enough to breed
                        is.na(indiv[, 6]) &        # are they alive
                        indiv[, 10] == 1)          # are they pregnant?
    if (nrow(mothers) == 0) {
      warning("There are no carrying females in the population")
    }
    ## Subset the mating females
    mating_females <- subset(indiv, indiv[, 2] == "F" &   # are they female
                               indiv[, 8] >= firstBreedFemale & # are they old enough to breed
                               is.na(indiv[, 6]) &        # are they alive
                               indiv[, 10] == 0)          # are they pregnant?
    if (nrow(mating_females) == 0)  {
      warning("There are no mating females in the population")
    }
    ## Subset potential fathers = all mature males that were alive ONE YEAR AGO!
    fathers <- subset(indiv, indiv[, 2] == "M" & 
                        indiv[, 8] - 1 >= firstBreedMale & # only include males that were mature a year ago
                        is.na(indiv[, 6]) | indiv[, 6] == year) # either alive or were still alive one year ago
    
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
      indiv$Pregnant[indiv$Me %in% mothers$Me] <- FALSE # No longer pregnant
      indiv$Pregnant[indiv$Me %in% mating_females$Me] <- TRUE # Mating females now pregnant
    }
    clutch <- clutch[clutch > 0]
  }
  clutch[clutch > maxClutch] <- maxClutch
  
  sprog.m <- makeFounders(pop = 0) # this creates empty pop, but without column 10
  if (!no_gestation) {
    sprog.m$Pregnant <- character(0) # add column 10 for pregnancy
  }
  
  for (s in unique(mothers[, 7])) {
    mothersInStock <- mothers[mothers[, 7] == s, , drop = FALSE]
    clutchInStock <- clutch[mothers[, 7] == s]
    fathersInStock <- fathers[fathers[, 7] == s, , drop = FALSE]
    if (nrow(fathersInStock) == 0) {
      warning(paste("There were no mature males in stock ", 
                    s, ", so ", nrow(mothersInStock), " mature females did not produce offspring", 
                    sep = ""))
      sprog.stock <- makeFounders(pop = 0)
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
    sprog.stock[, 1] <- uuid(n = nrow(sprog.stock), drop_hyphens = TRUE)
    sprog.stock[, 2] <- sample(c("M", "F"), nrow(sprog.stock), 
                               TRUE, prob = osr)
    sprog.stock[, 5] <- c(rep(year, nrow(sprog.stock)))
    sprog.stock[, 6] <- c(rep(NA, nrow(sprog.stock)))
    sprog.stock[, 7] <- as.integer(rep(s), nrow(sprog.stock))
    sprog.stock[, 8] <- c(rep(0, nrow(sprog.stock)))
    sprog.stock[, 9] <- c(rep(NA, nrow(sprog.stock)))
    if (!no_gestation) {
      sprog.stock[, 10] <- c(rep(0, nrow(sprog.stock)))
    }
    sprog.m <- rbind(sprog.m, sprog.stock)
  }
  names(sprog.m) <- names(indiv)
  indiv <- rbind(indiv, sprog.m)
  
  if (!quiet) cat("There were", nrow(mothers), "breeding mothers this cycle.\n")
  return(indiv)
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
  
  sprog.m <- makeFounders(pop = 0) # this creates empty pop, but without column 10
  sprog.m$Recovery <- character(0) # add column 10 for recovery years remaining
  
  for (s in unique(mothers[, 7])) {
    mothersInStock <- mothers[mothers[, 7] == s, , drop = FALSE]
    clutchInStock <- clutch[mothers[, 7] == s]
    fathersInStock <- fathers[fathers[, 7] == s, , drop = FALSE]
    if (nrow(fathersInStock) == 0) {
      warning(paste("There were no mature males in stock ", 
                    s, ", so ", nrow(mothersInStock), " mature females did not produce offspring", 
                    sep = ""))
      sprog.stock <- makeFounders(pop = 0)
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
    sprog.stock[, 1] <- uuid(n = nrow(sprog.stock), drop_hyphens = TRUE)
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
