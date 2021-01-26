# Helper Function "find_diff" - calculates the absolute value of the difference in mass or in RT between each pairs of features in a flattened matrix.
# Input: 
  ## corMatrix: a flattened correlation matrix with 3 columns: 1) "peak1" : the first feature/peak in the pair; 2) "peak2": the second feature in the pair
  ## metadata: meta data on all of the HILIC features. Same as the second parameter in the function "SearchDuplicates"
  ## factor_column: the factor based on which the difference will be calculated. It can either be "mass" to get mass difference or "RT" to get RT difference
# Output: vector of numeric values representing the absolute value of the differences by factor_column
find_diff <- function(corMatrix, metadata, factor_column = "mass") {
  data_peak1 <- metadata[match(corMatrix[,"peak1"], metadata[,"peak"]), factor_column]
  data_peak2 <- metadata[match(corMatrix[,"peak2"], metadata[,"peak"]), factor_column]
  return(abs(data_peak1 - data_peak2))
}

# Main function
SearchDuplicates <- function(hilic, metadata, corr_cutoff = 0.9, rt_cutoff = 0.2, ppm_cutoff = 15, artifacts, neutral, condition_sets) {
  # step 1 -  get correlation coefficient data
  cor.data <- cor(t(hilic), method = "pearson")
    ## flattening the correlation matrix
  ut <- upper.tri(cor.data)
  cor.data <- data.frame(peak1 = rownames(cor.data)[row(cor.data)[ut]], peak2 = rownames(cor.data)[col(cor.data)[ut]], corr_coeff =(cor.data)[ut])
  rm(ut)
    ## data that pass the corr_cutoff
  cor.data <- cor.data[which(cor.data$corr_coeff >= corr_cutoff),]
    ## add mass_diff and RT_diff columns to the flattened matrix
  cor.data$RT_diff <- find_diff(corMatrix = cor.data, metadata = metadata, factor_column = "RT")
  cor.data$mass_diff <- find_diff(corMatrix = cor.data, metadata = metadata, factor_column = "mass")
    ## details about mode, mass, RT data for each individual peak for ppm calculations
  cor.data$mass_peak1 <- metadata[match(cor.data$peak1, metadata[,"peak"]), "mass"]
  cor.data$mass_peak2 <- metadata[match(cor.data$peak2, metadata[,"peak"]), "mass"]
  cor.data$RT_peak1 <- metadata[match(cor.data$peak1, metadata[,"peak"]), "RT"]
  cor.data$RT_peak2 <- metadata[match(cor.data$peak2, metadata[,"peak"]), "RT"]
  cor.data$mode_peak1 <- metadata[match(cor.data$peak1, metadata[,"peak"]), "mode"]
  cor.data$mode_peak2 <- metadata[match(cor.data$peak2, metadata[,"peak"]), "mode"]  
  
  # step 2 - pairs that pass the RT cutoff
  cor.data <- cor.data[which(cor.data$RT_diff < rt_cutoff),]
  
  # step 3 - condition sets for duplicate pairs
  ## condition set 1
  cond1.data <- cor.data[,-which(names(cor.data) %in% c("RT_peak1", "RT_peak2"))]
  ### ppm calculation
  cond1.data$ppm <- cond1.data$mass_diff*1000000/cond1.data$mass_peak1
  if(neutral == TRUE) {cond1.data <- cond1.data[which(cond1.data$ppm < 15),]
  } else {
    cond1.data <- cond1.data[which(cond1.data$mode_peak1 == cond1.data$mode_peak2),]
    cond1.data <- cond1.data[which(cond1.data$ppm < 15),]
  }
  cond1.data <- cond1.data[,-which(names(cond1.data) %in% c("mass_peak1", "mass_peak2", "mode_peak1", "mode_peak2"))] # remove columns I don't need anymore
  cond1.data$artifact_ID <- NA
  cond1.data$how_many_units <- 0
  
  ## remove the peaks that passed the first condition from cor.data
  cor.data <- cor.data[-c(match(paste(cond1.data$peak1, cond1.data$peak2), paste(cor.data$peak1, cor.data$peak2))),]
  
  if(condition_sets == 1) {
    message(paste("Success! You ran the SearchDuplicates function to return putative duplicates from the first set of condition only:\n - the correlation coefficient and RT difference cutoffs \n - condition set 1: Mass difference < ", ppm_cutoff ," ppm.\n", sep = ""))
    my_results_1set <- list(cond1.data, cor.data)
    names(my_results_1set) <- c("condition1_results", "residual_pais")
    return(my_results_1set)
  }
  
  ## condition set 2 
  cond2.data <- cor.data[,-which(names(cor.data) %in% c("RT_peak1", "RT_peak2"))]
  ### check if modes for all pairs are the same or different. 
  cond2.data$both_mode <- ifelse(cond2.data$mode_peak1 == cond2.data$mode_peak2, cond2.data$mode_peak1, "different")
  ### ppm calculation
  cond2.data$ppm <- NA
  cond2.data$artifact_ID <- NA
  cond2.data$how_many_units <- 1
  
  for(i in 1:nrow(cond2.data)){
    for(j in 1:nrow(artifacts)){
      # check if ionization modes are the same and masses are within  ppm_cutoff 
      if((cond2.data$both_mode[i] == artifacts$mode[j]) & ((abs(artifacts$mass[j] - cond2.data$mass_diff[i])*1000000/artifacts$mass[j]) < ppm_cutoff)) {
        # ppm value
        cond2.data$ppm[i] <- abs(artifacts$mass[j] - cond2.data$mass_diff[i])*1000000/artifacts$mass[j]
        # the identity of the same artifact different selected
        cond2.data$artifact_ID[i] <- artifacts$ID[j]
        next
      }
    }
  }

  cond2.data <- cond2.data[which(!is.na(cond2.data$ppm)),-which(names(cond2.data) %in% c("both_mode","mass_peak1", "mass_peak2", "mode_peak1", "mode_peak2"))]# remove columns I don't need anymore
  
  # remove peaks that passed the second condition from cor.data
  cor.data <- cor.data[-c(match(paste(cond2.data$peak1, cond2.data$peak2), paste(cor.data$peak1, cor.data$peak2))),]
  
  if(condition_sets == 2) {
    message(paste("Success! You ran the SearchDuplicates function to return putative duplicates from two sets of condition :\n - the correlation coefficient and RT difference cutoffs \n - condition set 1:Mass difference  < ", ppm_cutoff ," ppm (condition1_results) \n - condition set 2: Difference in metabolite pairs is < ", ppm_cutoff ," ppm of one unit of same artifact different (condition2_results).\n", sep = ""))
    my_results_2sets <- list(cond1.data, cond2.data, cor.data)
    names(my_results_2sets) <- c("condition1_results", "condition2_results", "residual_pairs")
    return(my_results_2sets)
  }
  
  ## condition set 3
  cond3.data <- cor.data[,-which(names(cor.data) %in% c("RT_peak1", "RT_peak2"))]
  ### check if modes for all pairs are the same or different. 
  cond3.data$both_mode <- ifelse(cond3.data$mode_peak1 == cond3.data$mode_peak2, cond3.data$mode_peak1, "different")
  
  ### maximum units of artifacts (contaminants/adducts/repeating units) possible
  max_units_possible <- vector( mode = "numeric", length = nrow(artifacts))
  names(max_units_possible) <- artifacts$ID
  for(i in 1:length(max_units_possible)) {
    max_units_possible[i] <- ceiling(max(cond3.data$mass_diff)/artifacts$mass[i]) #rounded up
  }
  
  ### ppm calculation
  cond3.data$ppm <- NA
  cond3.data$artifact_ID <- NA
  cond3.data$how_many_units <- NA
  
  for(i in 1:nrow(cond3.data)){
    for (j in 1:nrow(artifacts)) {
      for(k in 2:max_units_possible[j]){
        # if ionization modes are the same and mass difference is within ppm_cutoff of k units of same artifact different
        if((cond3.data$both_mode[i] == artifacts$mode[j]) & ((abs((k*artifacts$mass[j]) - cond3.data$mass_diff[i])*1000000/(k*artifacts$mass[j])) < ppm_cutoff)) {
          # ppm value
          cond3.data$ppm[i] <- abs((k*artifacts$mass[j]) - cond3.data$mass_diff[i])*1000000/(k*artifacts$mass[j])
          # the identity of the same artifact different selected
          cond3.data$artifact_ID[i] <- artifacts$ID[j]
          # number of units of the same artifact different selected
          cond3.data$how_many_units[i] <- k
          next
        }
      }
    }
  }
    
  cond3.data <- cond3.data[which(!is.na(cond3.data$ppm)),-which(names(cond3.data) %in% c("both_mode","mass_peak1", "mass_peak2", "mode_peak1", "mode_peak2"))]# remove columns I don't need anymore
  cor.data <- cor.data[-c(match(paste(cond3.data$peak1, cond3.data$peak2), paste(cor.data$peak1, cor.data$peak2))),]
  
  
  if(condition_sets == 3) {
    message(paste("Success! You ran the SearchDuplicates function to return putative duplicates from three sets of condition :\n - after the correlation coefficient and RT difference cutoffs \n - condition set 1: Mass difference  < ", ppm_cutoff ," ppm (condition1_results) \n - condition set 2: Difference in metabolite pairs is < ", ppm_cutoff ," ppm of one unit of same artifact different (condition2_results) \n - condition set 3: Difference in metabolite pairs is within ", ppm_cutoff ," ppm of 2+ units of same artifact different (condition3_results).\n", sep = ""))
    my_results_3sets <- list(cond1.data, cond2.data, cond3.data, cor.data)
    names(my_results_3sets) <- c("condition1_results", "condition2_results", "condition3_results", "residual_pairs")
    return(my_results_3sets)
  }
}

