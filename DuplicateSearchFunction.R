# CONFIGURING REPEATING UNITS AND ADDUCT DATA
# library(readxl)

# # One Unit
#artifacts_oneUnit <- data.frame(read_xlsx("Artifacts_data(repeating_units).xlsx", sheet = "OneUnit"))
# # Peaks like DMSO, MALDI, etc. were unlikely based on our experimental design. So, we can remove them.
# artifacts_oneUnit <- artifacts_oneUnit[-grep("DMSO", artifacts_oneUnit$Origin),]
# artifacts_oneUnit <- artifacts_oneUnit[-grep("MALDI", artifacts_oneUnit$Origin),]
# artifacts_oneUnit <- artifacts_oneUnit[-grep("cesium", artifacts_oneUnit$Origin),]
# artifacts_oneUnit <- artifacts_oneUnit[-grep("SDS", artifacts_oneUnit$Origin),]
# artifacts_oneUnit <- artifacts_oneUnit[,c(2, 1, 4)]
# colnames(artifacts_oneUnit) <- c("ID", "mass", "mode")
#
# # Multiple units
# artifacts_multUnit <- data.frame(read_xlsx("Artifacts_data(repeating_units).xlsx", sheet = "MultipleUnits"))
# # Peaks like DMSO, MALDI, etc. were unlikely based on our experimental design. So, we can remove them.
# artifacts_multUnit <- artifacts_multUnit[-grep("DMSO", artifacts_multUnit$Origin),]
# artifacts_multUnit <- artifacts_multUnit[-grep("MALDI", artifacts_multUnit$Origin),]
# artifacts_multUnit <- artifacts_multUnit[-grep("cesium", artifacts_multUnit$Origin),]
# artifacts_multUnit <- artifacts_multUnit[-grep("SDS", artifacts_multUnit$Origin),]
# artifacts_multUnit <- artifacts_multUnit[,c(2, 1, 4)]
# colnames(artifacts_multUnit) <- c("ID", "mass", "mode")
#
# # save these
# save(artifacts_multUnit, artifacts_oneUnit, file = "Articfacts_data.RData")


## function to flatten the correlation matrix
flattenCorrMatrix <- function(cormat) {
  ut <- upper.tri(cormat)
  data.frame(peak1 = rownames(cormat)[row(cormat)[ut]], peak2 = rownames(cormat)[col(cormat)[ut]], corr_coeff =(cormat)[ut])
}

## edit the hilic peak data to extract the mass or RT for each pair of peaks
find_diff <- function(corMatrix, metadata, factor_column = "mass") {
  data_peak1 <- metadata[match(corMatrix[,"peak1"], metadata[,"peak"]), factor_column]
  data_peak2 <- metadata[match(corMatrix[,"peak2"], metadata[,"peak"]), factor_column]
  return(abs(data_peak1 - data_peak2))
}

## function that calculates the ppm 
ppm_calc <- function(theoretical, actual){
  return(abs(theoretical-actual)*1000000/theoretical)
}

# function to replace the missing values by the half min abundance for that specific metabolite.
half_min_impute <- function(row){
  if(length(which(is.na(row))) > 0){
    min_abundance <- min(row, na.rm = TRUE)
    half_min <- min_abundance/2
    row[which(is.na(row))] <- half_min
    return(row)    
  } else {
    return(row)
  }
}


# The main function to ID duplicates and all the details

DuplicatesAction <- function(cor.data, salt_mass, salt_name, artifacts_oneUnit, artifacts_multUnit, ppm_cutoff, FUN, to_discard){
  # ----- STEP 1 ----- 
  step1.data <- cor.data
  salts_mass_list <- vector( mode = "numeric", length = 200) # create a vector with the masses of up to 200 units of salt
  for(i in 1:length(salts_mass_list)) {
  salts_mass_list[i] <- i*salt_mass
  names(salts_mass_list)[i] <- paste(i,"_units", sep = "")
  }
  ## search
  step1.data$ppm <- NA
  step1.data$how_many_units <- NA
  step1.data$salt_ID <- artifacts_oneUnit$ID[grep(salt_name, artifacts_oneUnit$ID)[1]]
  for(i in 1:nrow(step1.data)){
    for(j in 1:length(salts_mass_list)){
    if(FUN(theoretical = salts_mass_list[j], actual = step1.data$mass_diff[i]) < ppm_cutoff){
      step1.data$ppm[i] <- FUN(theoretical = salts_mass_list[j], actual = step1.data$mass_diff[i])
      step1.data$how_many_units[i] <- names(salts_mass_list)[j] # number of units of the salt selected
      next
      }
    }
  }
  step1.data <- step1.data[which(!is.na(step1.data$ppm)),]
  cor.data <- cor.data[-c(match(paste(step1.data$peak1, step1.data$peak2), paste(cor.data$peak1, cor.data$peak2))),]# remove these pairs from correlation data
  ## fate: find the pairs that match the sodium acetate clusters and remove all pairs in which those peaks are involved
  junk <- ifelse(step1.data$mass_peak1 > step1.data$mass_peak2, step1.data$peak1, step1.data$peak2)
  junk <- unique(junk)
  step1.data <- step1.data[,-c(6:11)]
  ### remove all pairs involving the highest mass peak from correlation data
  rows_to_remove <- c()
  for(i in 1:nrow(cor.data)){
    for(j in junk){
      if(cor.data$peak1[i] == j | cor.data$peak2[i] == j){
      rows_to_remove <- append(rows_to_remove, i)
      }
    }
  }
  rows_to_remove <- unique(rows_to_remove)
  cor.data <- cor.data[-rows_to_remove,]# remove from correlation data

  step1.removed <- nrow(initial)-nrow(cor.data)

  to_discard <- rbind(to_discard, data.frame(peak = junk, reason = "step 1 - sodium acetate clusters"))
  to_discard <- to_discard[-1,]
  rm(junk, rows_to_remove)
  
  # ----- STEP 2 -----
  ## peaks with masses very similar to each other (i.e. corresponding to no repeating units or other artifacts)
  step2a.data <- cor.data
  step2a.data$ppm <- FUN(theoretical = step2a.data$mass_peak1, actual = step2a.data$mass_peak2)
  step2a.data <- step2a.data[which(step2a.data$ppm < 15),-c(6:11)]
  step2a.data$artifact_ID <- NA
  step2a.data$how_many_units <- 0

  ## peaks with mass difference corresponding to one unit of artifact
  step2b.data <- cor.data
  step2b.data$both_mode <- ifelse(step2b.data$mode_peak1 == step2b.data$mode_peak2, step2b.data$mode_peak1, "different")
  step2b.data$ppm <- NA
  step2b.data$artifact_ID <- NA
  step2b.data$how_many_units <- "1 unit"

  for(i in 1:nrow(step2b.data)){
    for(j in 1:nrow(artifacts_oneUnit)){
    #if ionization modes are the same and masses are within ppm < ppm_cutoff 
      if((step2b.data$both_mode[i] == artifacts_oneUnit$mode[j]) & (FUN(theoretical = artifacts_oneUnit$mass[j], actual = step2b.data$mass_diff[i]) < ppm_cutoff)) {
        #ppm value
        step2b.data$ppm[i] <- FUN(theoretical = artifacts_oneUnit$mass[j], actual = step2b.data$mass_diff[i])
        #the identity of the artifact selected
        step2b.data$artifact_ID[i] <- artifacts_oneUnit$ID[j]
        next
      }
    }
  }
  step2b.data <- step2b.data[which(!is.na(step2b.data$ppm)),-c(6:12)]

  ## two or more units of artifacts
  step2c.data <- cor.data
  step2c.data$both_mode <- ifelse(step2c.data$mode_peak1 == step2c.data$mode_peak2, step2c.data$mode_peak1, "different")
  step2c.data$ppm <- NA
  step2c.data$artifact_ID <- NA
  step2c.data$how_many_units <- NA

  ### maximum adducts possible
  max_units_possible <- vector( mode = "numeric", length = nrow(artifacts_multUnit))
  names(max_units_possible) <- artifacts_multUnit$ID
  for(i in 1:length(max_units_possible)) {
    max_units_possible[i] <- ceiling(max(step2c.data$mass_diff)/artifacts_multUnit$mass[i])#rounded up
  }
  for(i in 1:nrow(step2c.data)){
    for (j in 1:nrow(artifacts_multUnit)) {
      for(k in 2:max_units_possible[j]){
        # if ionization modes are the same and mass difference is within 0.005
        if((step2c.data$both_mode[i] == artifacts_multUnit$mode[j]) & (FUN(theoretical = (k*artifacts_multUnit$mass[j]), actual = step2c.data$mass_diff[i]) < ppm_cutoff)) {
          #ppm value
          step2c.data$ppm[i] <- FUN(theoretical = (k*artifacts_multUnit$mass[j]), actual = step2c.data$mass_diff[i])
          #the identity of the artifact selected
          step2c.data$artifact_ID[i] <- artifacts_multUnit$ID[j]
          #number of units of the artifact selected
          step2c.data$how_many_units[i] <- paste(k, " units", sep="")
          next
        }
      }
    }
  }
  step2c.data <- step2c.data[which(!is.na(step2c.data$ppm)),-c(6:12)]
  # merging step2a, step2b, and step2c data
  step2.data <- rbind(step2a.data, step2b.data, step2c.data)
  cor.data <- cor.data[-c(match(paste(step2.data$peak1, step2.data$peak2), paste(cor.data$peak1, cor.data$peak2))),]# remove these pairs from correlation data

  ## fate: tag the lowest intensity peaks in the to_discard pile and update the cor.data 
  metab_lowest_median <- c()
  for(i in 1:nrow(step2.data)){
    lowest <- ifelse(median(unlist(hilic[match(step2.data$peak1[i], rownames(hilic)),])) < median(unlist(hilic[match(step2.data$peak2[i], rownames(hilic)),])), step2.data$peak1[i], step2.data$peak2[i])
    metab_lowest_median <- append(metab_lowest_median, lowest)
  }
  metab_lowest_median <- unique(metab_lowest_median)
  ### add these peaks to the to_discard bin
  to_discard <- rbind(to_discard, data.frame(peak = metab_lowest_median, reason = "step 2 - low intensity peak from no artifacts or non-salt adducts "))
  ### remove all pairs involving the lowest median peak from correlation data
  rows_to_remove <- c()
  for(i in 1:nrow(cor.data)){
    for(j in metab_lowest_median){
      if(cor.data$peak1[i] == j | cor.data$peak2[i] == j){
        rows_to_remove <- append(rows_to_remove, i)
      }
    }
  }
  rows_to_remove <- unique(rows_to_remove)
  cor.data <- cor.data[-rows_to_remove,]
  step2.removed <- nrow(initial)-step1.removed-nrow(cor.data)
  rm(rows_to_remove, metab_lowest_median, max_units_possible, lowest)
  
  
  # ----- STEP 3 -----
  step3.data <- cor.data
  count <- c()
  for(i in 1:nrow(step3.data)){
    if(step3.data$mass_peak1[i] < 1396 & step3.data$mass_peak2[i] < 1396 & step3.data$mass_diff[i] >= 0.999 
     & step3.data$mass_diff[i] <= 1.0063 & step3.data$RT_diff[i] <= 0.01){
      count <- append(count, i)
    }
  }
  step3.data <- step3.data[count,]
  cor.data <- cor.data[-c(match(paste(step3.data$peak1, step3.data$peak2), paste(cor.data$peak1, cor.data$peak2))),]

  ## fate: 
  metab_highest_mass <- c()
  for(i in 1:nrow(step3.data)){
    highest <- ifelse(step3.data$mass_peak1[i] > step3.data$mass_peak2[i], step3.data$peak1[i], step3.data$peak2[i])
    metab_highest_mass <- append(metab_highest_mass, highest)
  }

  metab_highest_mass <- unique(metab_highest_mass)
  ### add these peaks to the to_discard bin
  to_discard <- rbind(to_discard, data.frame(peak = metab_highest_mass, reason = "step 3 - isotope with the highest mass"))
  ### remove all pairs involving the highest mass peak from correlation data
  rows_to_remove <- c()
  for(i in 1:nrow(cor.data)){
    for(j in metab_highest_mass){
      if(cor.data$peak1[i] == j | cor.data$peak2[i] == j){
        rows_to_remove <- append(rows_to_remove, i)
      }
    }
  }
  rows_to_remove <- unique(rows_to_remove)
  cor.data <- cor.data[-rows_to_remove,]
  step3.data <- step3.data[,-c(6:11)]
  step3.removed <- nrow(initial)-step1.removed-step2.removed-nrow(cor.data)
  rm(count, highest, rows_to_remove, metab_highest_mass)
  

  # ----- STEP 4 -----
  step4.data <- cor.data
  step4.data$ratio <- NA
  step4.data$highest_mass <- NA
  step4.data$lowest_mass <- NA
  for(i in 1:nrow(step4.data)){
    step4.data$ratio[i] <- ifelse(step4.data$mass_peak1[i] < step4.data$mass_peak2[i], 
                                step4.data$mass_peak1[i]/step4.data$mass_peak2[i], 
                                step4.data$mass_peak2[i]/step4.data$mass_peak1[i])
    step4.data$highest_mass[i] <- ifelse(step4.data$mass_peak1[i] < step4.data$mass_peak2[i],
                                       "peak2", "peak1")
    step4.data$lowest_mass[i] <- ifelse(step4.data$highest[i] == "peak1", "peak2", "peak1")
  }
  count <- c()
  for(i in 1:nrow(step4.data)){
    if(((step4.data$RT_peak1[i]>=2.7) & (step4.data$RT_peak2[i]>=2.7) &
      (step4.data$RT_peak1[i]<=3.2) & (step4.data$RT_peak2[i]<=3.2)) |
     ((step4.data$RT_peak1[i]>=5.2) & (step4.data$RT_peak2[i]>=5.2) &
      (step4.data$RT_peak1[i]<=5.9) & (step4.data$RT_peak2[i]<=5.9))){ 
    if((step4.data$ratio[i]>=0.45 & step4.data$ratio[i]<=0.55) & 
       median(unlist(hilic[match(step4.data[i,step4.data$lowest_mass[i]], rownames(hilic)),])) >= median(unlist(hilic[match(step4.data[i,step4.data$highest_mass[i]], rownames(hilic)),]))){
        count <- append(count, i)
      }
    }
  }
  step4.data <- step4.data[count,-c(6:11)]
  cor.data <- cor.data[-c(match(paste(step4.data$peak1, step4.data$peak2), paste(cor.data$peak1, cor.data$peak2))),]

  ## fate: 
  metab_highest_mass <- c()
  for(i in 1:nrow(step4.data)){
    highest <- step4.data[i,step4.data$highest_mass[i]]
    metab_highest_mass <- append(metab_highest_mass, highest)
  }
  metab_highest_mass <- unique(metab_highest_mass)
  ### add these peaks to the to_discard bin
  to_discard <- rbind(to_discard, data.frame(peak = metab_highest_mass, reason = "step 4 - lipid dimer with the highest mass"))
  ### remove all pairs involving the lowest median peak from correlation data
  rows_to_remove <- c()
  for(i in 1:nrow(cor.data)){
    for(j in metab_highest_mass){
      if(cor.data$peak1[i] == j | cor.data$peak2[i] == j){
        rows_to_remove <- append(rows_to_remove, i)
      }
    }
  }
  rows_to_remove <- unique(rows_to_remove)
  cor.data <- cor.data[-rows_to_remove,]
  step4.data <- step4.data[,-c(7:8)]
  step4.removed <- nrow(initial)-step1.removed-step2.removed-step3.removed-nrow(cor.data)
  rm(rows_to_remove, metab_highest_mass,count,highest) 
  
  # ----- summary & return -----
  message(paste("Initially, we started with", nrow(initial), "pairs that had a correlation coefficient of >= to", corr_cutoff,"and a RT difference of <= to", rt_cutoff,".", sep = " "))
  message(paste("Step 1 removed", step1.removed ,"pairs of peaks.", sep = " "))
  message(paste("Step 2 removed", step2.removed ,"pairs of peaks.", sep = " "))
  message(paste("Step 3 removed", step3.removed,"pairs of peaks.", sep = " "))
  message(paste("Step 4 removed", step4.removed ,"pairs of peaks.", sep = " "))
  message(paste("At the end, we reduced this list to", nrow(cor.data) ,"pairs and got rid of", nrow(initial)-nrow(cor.data),"pairs total.", sep = " "))
  to_return <- list(step1.data, step2a.data, step2b.data, step2c.data, step2.data, step3.data, step4.data, cor.data, to_discard)
  names(to_return) <- c("step1.data", "step2a.data", "step2b.data", "step2c.data", "overall.step2.data", "step3.data", "step4.data", "remaining.data","to_remove_from_HILIC")
  return(to_return)
}


# EXAMPLE
# loading data 
load("../../Data/HILIC_processedData_WITH_duplicates.Rdata")
load("Articfacts_data.RData")
sample.metadata <- meta.data
rm(meta.data)
peak.metadata <- read.csv("../../Data/HILIC_peaks_metadata.csv")

# parameters and names
hilic <- processed.hilic.data
rm(processed.hilic.data)
corr_cutoff <- 0.9
rt_cutoff <- 0.2
ppm_cutoff <- 15
to_discard <- data.frame(peak = NA, reason = NA) #peak to discard from raw HILIC data

# Correlation and RT differences calculations
cor.data <- flattenCorrMatrix(cormat = cor(t(hilic), method = "pearson"))
## data that pass the corr_cutoff
cor.data <- cor.data[which(cor.data$corr_coeff >= corr_cutoff),]
## add mass_diff and RT_diff
cor.data$RT_diff <- find_diff(corMatrix = cor.data, metadata = peak.metadata, factor_column = "RT")
cor.data$mass_diff <- find_diff(corMatrix = cor.data, metadata = peak.metadata, factor_column = "mass")
## individual masses, RTs, and modes for ppm calculations
cor.data$mass_peak1 <- peak.metadata[match(cor.data$peak1, peak.metadata[,"peak"]), "mass"]
cor.data$mass_peak2 <- peak.metadata[match(cor.data$peak2, peak.metadata[,"peak"]), "mass"]
cor.data$RT_peak1 <- peak.metadata[match(cor.data$peak1, peak.metadata[,"peak"]), "RT"]
cor.data$RT_peak2 <- peak.metadata[match(cor.data$peak2, peak.metadata[,"peak"]), "RT"]
cor.data$mode_peak1 <- peak.metadata[match(cor.data$peak1, peak.metadata[,"peak"]), "mode"]
cor.data$mode_peak2 <- peak.metadata[match(cor.data$peak2, peak.metadata[,"peak"]), "mode"]
## data that pass the RT cutoff
cor.data <- cor.data[which(cor.data$RT_diff < rt_cutoff),]
initial <- cor.data
mass_NaCH3CO2 <- artifacts_oneUnit$mass[grep("NaCH3CO2", artifacts_oneUnit$ID)[1]]

method1.results <- DuplicatesAction(cor.data = cor.data,
             salt_mass = mass_NaCH3CO2,
             salt_name = "NaCH3CO2",
             artifacts_oneUnit = artifacts_oneUnit,
             artifacts_multUnit = artifacts_multUnit,
             ppm_cutoff = ppm_cutoff,
             FUN = ppm_calc,
             to_discard = to_discard
            )
write.csv(method1.results$to_remove_from_HILIC, "preprocessed_duplicates_new_version.csv")
