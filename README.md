# Duplicate-Search
A framework to search potential duplicates of metabolomics data and other artifacts such as false peaks. 

Hydrophilic Interaction Liquid Chromatography (HILIC) data are collected with information about feature mass and retention time (RT). However, it is not always obvious how to select which features are potential duplicates of each other, or those to which artifacts (e.g. salts) were added throughout the HILIC  process. 

This work is an attempt at a reproducible framework to identify potential metabolite duplicates and HILIC artifacts. The workflow of this particular quest consists of 3 sets of conditions (See section below).            

 ## Workflow 
 
<img width="925" alt="Screen Shot 2021-01-19 at 2 32 10 PM" src="https://user-images.githubusercontent.com/72724703/105083710-4a50bf00-5a63-11eb-8e7e-e65c9f7548bc.png">

In order to extract potential duplicate pairs and those related to specific HILIC articfacts, we will apply a set of condition to our data. 
The first line of conditions are correlation coefficient and RT conditions and these are employed in all 3 sets of conditions. Potential duplicate pairs must have:
      
      a. $correlation\ coefficient > corr_cutoff$
      b. $RT\ difference < rt_cutoff$
      c. Then, 3 sets of conditions are separately applied to the duplicate pairs data. In the following, these sets are applied:
            i. Condition set 1  : $ppm=\frac{(|mass1-mass2|)*1,000,000} {mass1} <= ppm_cutoff$
            ii. Condition set 2 : $ppm=\frac{(|mass_difference-mass_adduct|)*1,000,000} {mass_adduct} <= ppm_cutoff$
            iii. Condition set 3: $ppm=\frac{(|mass1_difference-(mass_adduct*k)|)*1,000,000} {(mass_adduct*k)} <= ppm_cutoff: n(>1)\ represents\ number\ of\ molecules\ of\ RUs.$
            
## Parameters for the function
  - hilic          = the data frame containing the HILIC data. Make sure the samples are in columns and the metabolites in rows. 
  - metadata       = more informations about the metabolites in the hilic data frame. This data frame must have four columns with the following names:
        - first column  : the metabolites or peak ID. Name of first column: "peak".
        - second column : the mass of the peak. Name of second column: "mass".
        - third column  : retention time (RT) of the peak. Name of third column: "RT".
        - fourth column : ionization mode at which the peak was collected. Name of fourth column: "mode".
   - corr_cutoff   = the correlation coefficient cutoff for potential duplicates. Set by default at 0.9.
   - rt_cutoff     = the retention time (RT) cutoff for potential duplicates. Set by default at 0.2.
   - ppm_cutoff    = the ppm cutoff for potential duplicates. Set by default at 15. The ppm cutoff depends on the condition_set described above. 
   - Artifact      = data frame with adduct or repeating unit(RU) information. The user can input their own data frame of artifacts making sure it has the three following columns with the corresponding names:
        - first column  : ID of RU or adduct 
        - second column : mass of the RU or adduct
        - third column  : ionization mode at which the RU or adduct was collected
   - nominal        = logical value. If TRUE, this means that peak masses are neutral. FALSE, otherwise.
   - condition_sets = 1, 2, or 3. 
        - if 1: condition set 1 only is performed. 
        - if 2: condition sets 1 and 2 are performed. 
        - if 3: condition sets 1, 2, and 3 are performed.
          

 
 ## Return
The SearchDuplicates function return a list depending on the condition_set paramater:
    - if 1: the function return a list of one data frame containing pairs that passed condition set 1
    - if 2: the function return a list of two data frames containing pairs that passed condition set 1 and condition set 2
    - if 3: the function return a list of three data frames containing pairs that passed condition set 1, 2, and 3.
    

