## ========================================================================== ##
## Some custom functions to add to fishSim to make it more similar to the     ##
## grey reef shark case study from Palmyra.                                   ##  
## ========================================================================== ##

findRelativesCustom <- function (indiv, sampled = TRUE, verbose = TRUE, 
                                 nCores = detectCores() - 1, delimitIndiv = TRUE) 
{
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
      ancestors[ancestors[, 1] == pairs[i,2], 1:15]
                                                                                                        
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



addPregnancy <- function(indiv, matingAges) {
  
  pregnant_females <- indiv$Sex == "F" & indiv$AgeLast %in% (matingAges)
  
  indiv$Pregnant <- pregnant_females
  
  return(indiv)
}


#' mateOrBirth
#'
#' @param indiv 
#' @param batchSize 
#' @param fecundityDist 
#' @param osr 
#' @param year 
#' @param firstBreedFemale 
#' @param firstBreedMale 
#' @param firstLitter 
#' @param type 
#' @param maxClutch 
#' @param singlePaternity 
#' @param exhaustFathers 
#' @param maturityCurve 
#' @param maleCurve 
#' @param femaleCurve 
#' @param quiet 
#'
#' @return
#' @export
#'
#' @examples
mateOrBirth <- function (indiv,
                         batchSize = 0.5, 
                         fecundityDist = "uniform", 
                         osr = c(0.5, 0.5), 
                         year = "-1", 
                         firstBreedFemale = 12, # only for females
                         firstBreedMale = 10,
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
  if (type == "ageSex") {
    
    mating_females <- mating_females[runif(nrow(mating_females)) < 
                                       femaleCurve[mating_females[, 8] + 1], , 
                       drop = FALSE]
    
    fathers <- fathers[runif(nrow(fathers)) < maleCurve[fathers[, 
                                                                8] + 1], , 
                       drop = FALSE]
    # if (fecundityDist == "poisson") {
    #   clutch <- rpois(n = nrow(mothers), lambda = batchSize)
    # }
    # if (fecundityDist == "truncPoisson") {
    #   clutch <- rTruncPoisson(n = nrow(mothers), T = batchSize)
    # }
    # if (fecundityDist == "binomial") {
    #   clutch <- rbinom(nrow(mothers), 1, prob = batchSize)
    # }
    if (fecundityDist == "uniform") {
      clutch <- sample(batchSize, size = nrow(mothers), replace = TRUE)
    }
    ## Here a line is added to set first litter at 1-2 pups, not 3-4. 
    clutch[mothers$AgeLast == firstBreedFemale + 1] <- sample(
      firstLitter,
      size = sum(mothers$AgeLast == firstBreedFemale + 1),
      replace = TRUE)
    
    mothers <- subset(mothers, clutch > 0)
    
    ## Add/remove pregnancy indicators
    indiv$Pregnant[indiv$Me %in% mothers$Me] <- FALSE # No longer pregnant
    indiv$Pregnant[indiv$Me %in% mating_females$Me] <- TRUE # Mating females now pregnant
    
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



#' addBirthRecovery()
#'
#' Add a 10th column to a starting population indiv object for birth recovery 
#' time left. Animals can only breed if this variable is 0. 
#' 
#' @param indiv The indiv object with the starting population
#' @param recoveryTime 
#' @param maturityAge 
#' @param random
#'
#' @return A indiv object with the 10th column added
#' @export
#'
#' @examples
#' # Five lighter shades
#' make_shades("goldenrod", 5)
#' # Five darker shades
#' make_shades("goldenrod", 5, lighter = FALSE)
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



#' mateWithRecovery()
#'
#' Add a mating cycle to a indiv object created by makeFounders(), where
#' females need to recover some time after breeding. Also added the breeding
#' option 'uniform'.
#' 
#' @param indiv The indiv object with most recent 
#' @param pause The number of shades to make
#' @param lighter Whether to make lighter (TRUE) or darker (FALSE) shades
#'
#' @return A vector of n colour hex codes
#' @export
#'
#' @examples
#' # Five lighter shades
#' make_shades("goldenrod", 5)
#' # Five darker shades
#' make_shades("goldenrod", 5, lighter = FALSE)
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



#' uniformCheckGrowthrate()
#'
#' A custom version of fishSim's check_growthrate() that allows for uniform
#' breeding and pauses after breeding cycles (the latter could also be 
#' achieved through a fixed biennial maturity curve).
#' 
#' @param indiv The indiv object with most recent 
#' @param pause The number of shades to make
#' @param lighter Whether to make lighter (TRUE) or darker (FALSE) shades
#'
#' @return A vector of n colour hex codes
#' @export
#'
#' @examples
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
  ageStockMort) # input for the mort() call
{
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

# female_curve <- c(rep(0, 19), 
#                   rep(c(1, 0), 50)) 
# 
# uniformCheckGrowthrate(
#   fecundityDist = "uniform", 
#   forceY1 = NA,      # year 1 mortality?
#   mateType = "ageSex", # mating type
#   mortType = "flat", # mortality constant or variable?
#   batchSize = c(3, 4, 5, 6), # expected litter size
#   firstBreed = 0, # first year of breeding (at least min maturity curve)
#   maxClutch = Inf, # max litter size
#   osr = c(0.5, 0.5), # sex ratio
#   maturityCurve = female_curve, # equivalent to femaleCurve, but used if mateType == "age"
#   femaleCurve = female_curve, # female maturity curve
#   maxAge = 60, # self explanatory
#   mortRate = 0.11, # rate of mortality
#   ageMort, # input for the mort() call
#   stockMort, # input for the mort() call
#   ageStockMort
# )
# 
# ## Not sure if this works, as I sometimes get negative growrates at values very close to positive growht reates.
# edit(PoNG)
# 
