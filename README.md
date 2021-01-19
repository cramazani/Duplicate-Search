# Duplicate-Search
A framework to search potential duplicates of metabolomics data

This is an attempt reproducible framework at extracting metabolite duplicates based on 3 sets of conditions. The function parameters are explained below:

## Parameters for the function
  - hilic          = the data frame containing the metabolomics data. Make sure the samples are in columns and the metabolites in rows. 
  - metadata       = more informations about the metabolites in the hilic data frame. This data frame must have four columns:
        - first column  : the metabolites or peak ID
        - second column : the mass of the peak 
        - third column  : retention time (RT) of the peak 
        - fourth column : ionization mode at which the peak was collected
   - corr_cutoff   = the correlation coefficient cutoff for potential duplicates. Set by default at 0.9 
   - rt_cutoff     = the retention time (RT) cutoff for potential duplicates. Set by default at 0.2
   - ppm_cutoff    = the ppm cutoff for potential duplicates. Set by default at 15. The ppm cutoff depends on the condition/cutoff (see formulas in conditions section below below). 
   - Adduct.Ru   = data frame with adduct or repeating unit(RU) information. This data frame must have three columns:
        - first column  : ID of RU or adduct 
        - second column : mass of the RU or adduct
        - third column  : ionization mode at which the RU or adduct was collected
   - nominal        = logical value. If TRUE, this means that peak masses are neutral. FALSE, otherwise.
   - condition_sets = 1, 2, or 3. 
        - if 1: condition set 1 only is performed. 
        - if 2: condition sets 1 and 2 are performed. 
        - if 3: condition sets 1, 2, and 3 are performed.
          
## Conditions
In order to extract potential duplicate pairs, we will apply a set of condition to our data. 
The first line of conditions are correlation coefficient and RT conditions. Potential duplicate pairs must have:
      
      a. $correlation\ coefficient > 0.90$
      b. $RT\ difference < 0.2$
      c. Then, 3 sets of conditions are separately applied to the duplicate pairs data. In the following, these sets are applied:
            i. Condition set 1  : $ppm=\frac{(|mass1-mass2|)*1,000,000} {mass1} <= ppm_cutoff$
            ii. Condition set 2 : $ppm=\frac{(|mass_difference-mass_adduct|)*1,000,000} {mass_adduct} <= ppm_cutoff$
            iii. Condition set 3: $ppm=\frac{(|mass1_difference-(mass_adduct*k)|)*1,000,000} {(mass_adduct*k)} <= ppm_cutoff: n(>1)\ represents\ number\ of\ molecules\ of\ RUs.$

            
