#Function

SearchDuplicates <- function(hilic, metadata, corr_cutoff = 0.9, rt_cutoff = 0.2, ppm_cutoff = 15, Adduct.Ru, nominal, condition_sets) {
  # some preliminary functions 
  ## function to flatten the correlation matrix
  flattenCorrMatrix <- function(cormat) {
    ut <- upper.tri(cormat)
    data.frame(peak1 = rownames(cormat)[row(cormat)[ut]], peak2 = rownames(cormat)[col(cormat)[ut]], corr_coeff =(cormat)[ut])
  }
  
  ## edit the hilic data to extract the mass or RT differences for each pair of peaks
  ### set factor_column to "mass" to get mass difference
  ### set factor_column to "RT" to get RT difference
  find_diff <- function(cordata, metadata, factor_column = "mass") {
    if(factor_column == "mass") by <- 2
    if(factor_column == "RT") by <- 3
    data_peak1 <- metadata[match(cordata$peak1, metadata[,1]), by]
    data_peak2 <- metadata[match(cordata$peak2, metadata[,1]), by]
    return(abs(data_peak1 - data_peak2))
  }
  
  # step 1 -  get correlation coefficient data
  cor.data <- flattenCorrMatrix(cormat = cor(t(hilic), method = "pearson"))
  ## data that pass the corr_cutoff
  cor.data <- cor.data[which(cor.data$corr_coeff >= corr_cutoff),]
  ## add mass_diff and RT_diff
  cor.data$RT_diff <- find_diff(cordata = cor.data, metadata = metadata, factor_column = "RT")
  cor.data$mass_diff <- find_diff(cordata = cor.data, metadata = metadata, factor_column = "mass")
  ## details about mode and mass data for ppm calculations
  cor.data$mass_peak1 <- metadata[match(cor.data$peak1, metadata[,1]), 2]
  cor.data$mass_peak2 <- metadata[match(cor.data$peak2, metadata[,1]), 2]
  cor.data$mode_peak1 <- metadata[match(cor.data$peak1, metadata[,1]),4]
  cor.data$mode_peak2 <- metadata[match(cor.data$peak2, metadata[,1]),4]  
  
  # step 2 - the RT cutoff
  cor.data <- cor.data[which(cor.data$RT_diff < rt_cutoff),]
  
  # step 3 - conditions for duplicate pairs
  ## condition set 1: no adducts/ru
  cond1.data <- cor.data
  ### ppm calculation
  cond1.data$ppm <- cond1.data$mass_diff*1000000/cond1.data$mass_peak1
  if(nominal == TRUE) cond1.data <- cond1.data[which(cond1.data$ppm < 15),]
  else {
    cond1.data <- cond1.data[which(cond1.data$mode_peak1 == cond1.data$mode_peak2),]
    cond1.data <- cond1.data[which(cond1.data$ppm < 15),]
  }
  cond1.data <- cond1.data[,-c(6:9)]#remove unnecessary columns
  cond1.data$RU_or_adduct_ID <- NA
  cond1.data$how_many_molecules <- 0
  
  # remove the peaks that passed the first condition from cor.data
  cor.data <- cor.data[-c(match(paste(cond1.data$peak1, cond1.data$peak2), paste(cor.data$peak1, cor.data$peak2))),]
  
  if(condition_sets == 1) {
    message(paste("You ran the SearchDuplicates function to return putative duplicates from the first set of condition only:\n - the correlation coefficient and RT difference cutoffs \n - condition set 1: Mass difference < ", ppm_cutoff ," ppm.\n", sep = ""))
    my_results_1set <- list(cond1.data)
    names(my_results_1set) <- c("condition1.results")
    return(my_results_1set)
  }
  
  ## condition set 2: one adduct/ru 
  cond2.data <- cor.data
  ### check if modes for all pairs are the same of different. 
  cond2.data$both_mode <- ifelse(cond2.data$mode_peak1 == cond2.data$mode_peak2, cond2.data$mode_peak1, "different")
  ### ppm calculation: show formula (this new version)
  cond2.data$ppm <- NA
  cond2.data$RU_or_adduct_ID <- NA
  cond2.data$how_many_molecules <- 1
  
  for(i in 1:nrow(cond2.data)){
    for(j in 1:nrow(Adduct.Ru)){
      # check if ionization modes are the same and masses are within  ppm_cutoff 
      if((cond2.data$both_mode[i] == Adduct.Ru[j,3]) & ((abs(Adduct.Ru[j,2] - cond2.data$mass_diff[i])*1000000/Adduct.Ru[j,2]) < ppm_cutoff)) {
        # ppm value
        cond2.data$ppm[i] <- abs(Adduct.Ru[j,2] - cond2.data$mass_diff[i])*1000000/Adduct.Ru[j,2]
        # the identity of the RU selected
        cond2.data$RU_or_adduct_ID[i] <- Adduct.Ru[j,1]
        next
      }
    }
  }

  cond2.data <- cond2.data[which(!is.na(cond2.data$ppm)),-c(6:10)]#remove unnecessary columns
  
  # remove peaks that passed the second condition from cor.data
  cor.data <- cor.data[-c(match(paste(cond2.data$peak1, cond2.data$peak2), paste(cor.data$peak1, cor.data$peak2))),]
  
  if(condition_sets == 2) {
    message(paste("You ran the SearchDuplicates function to return putative duplicates from two sets of condition :\n - the correlation coefficient and RT difference cutoffs \n - condition set 1:Mass difference  < ", ppm_cutoff ," ppm (cond1.data) \n - condition set 2: Difference in metabolite pairs is < ", ppm_cutoff ," ppm of one molecule of RU (cond2.data).\n", sep = ""))
    my_results_2sets <- list(cond1.data, cond2.data)
    names(my_results_2sets) <- c("condition1.results", "condition2.results")
    return(my_results_2sets)
  }
  
  ## 3. condition set 3: 2+ adducts/rus
  cond3.data <- cor.data
  ### check if modes for all pairs are the same of different. 
  cond3.data$both_mode <- ifelse(cond3.data$mode_peak1 == cond3.data$mode_peak2, cond3.data$mode_peak1, "different")
  
  ### maximum adduct possible
  max_molecules_possible <- vector( mode = "numeric", length = nrow(Adduct.Ru))
  names(max_molecules_possible) <- Adduct.Ru[,1]
  for(i in 1:length(max_molecules_possible)) {
    max_molecules_possible[i] <- ceiling(max(cond3.data$mass_diff)/Adduct.Ru[i,2]) #rounded up
  }
  
  ### ppm calculation: show formula (this new version)
  cond3.data$ppm <- NA
  cond3.data$RU_or_adduct_ID <- NA
  cond3.data$how_many_molecules <- NA
  
  for(i in 1:nrow(cond3.data)){
    for (j in 1:nrow(Adduct.Ru)) {
      for(k in 2:max_molecules_possible[j]){
        # if ionization modes are the same and mass difference is within ppm_cutoff of k molecules of RU
        if((cond3.data$both_mode[i] == Adduct.Ru[j,3]) & ((abs((k*Adduct.Ru[j,2]) - cond3.data$mass_diff[i])*1000000/(k*Adduct.Ru[j,2])) < ppm_cutoff)) {
          # ppm value
          cond3.data$ppm[i] <- abs((k*Adduct.Ru[j,2]) - cond3.data$mass_diff[i])*1000000/(k*Adduct.Ru[j,2])
          # the identity of the RU selected
          cond3.data$RU_or_adduct_ID[i] <- Adduct.Ru[j,1]
          # number of molecules of the RU selected
          cond3.data$how_many_molecules[i] <- k
          next
        }
      }
    }
  }
    
  cond3.data <- cond3.data[which(!is.na(cond3.data$ppm)),-c(6:10)]#remove unnecessary columns
  cor.data <- cor.data[-c(match(paste(cond3.data$peak1, cond3.data$peak2), paste(cor.data$peak1, cor.data$peak2))),]
  
  
  if(condition_sets == 3) {
    message(paste("You ran the SearchDuplicates function to return putative duplicates from three sets of condition :\n - after the correlation coefficient and RT difference cutoffs \n - condition set 1: Mass difference  < ", ppm_cutoff ," ppm (cond1.data) \n - condition set 2: Difference in metabolite pairs is < ", ppm_cutoff ," ppm of one molecule of RU (cond2.data) \n - condition set 3: Difference in metabolite pairs is within ", ppm_cutoff ," ppm of 2+ molecules of RU (cond3.data).\n", sep = ""))
    my_results_3sets <- list(cond1.data, cond2.data, cond3.data)
    names(my_results_3sets) <- c("condition1.results", "condition2.results", "condition3.results")
    return(my_results_3sets)
  }
}

