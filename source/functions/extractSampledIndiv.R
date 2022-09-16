extractSampledIndiv <- function(indiv) {
  ## Extract the self captures
  self_captures <- indiv[indiv$no_samples > 1, ]

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