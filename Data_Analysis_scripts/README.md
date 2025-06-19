### Reproduction of the results:
In the archive you will find 
- The main script `211.GeneRation_Analyses_code_Dryad.R` allowing for the reproduction of the data analysis workflow
- Two functions-carrying R files (`GeneRation_DataAnalysis_FUN.R`and `GeneRation_Fun_v1.R`)
- The folder `SIMULATIONS` containing raw data from 2 succesful simulations that were generated for the paper and used here as simple example for data treatment
- The folder `Extracted_data` used as destination folder to store extracted output from raw simulation data

The script `211.GeneRation_Analyses_code_Dryad.R` that is the main script to launch in R or Rstudio, containing blocks of code leading to
- the extraction of data from the simulations stored in the `SIMULATIONS` folder. This includes:
  - `storeFits`, a dataframe containing fitnesses of simulations, where each line corresponds to a simulation and each column to a generation
  - `storeJustLastFits`, a dataframe containing final fitnesses of simulations, where each line corresponds to a simulation with a column containing the final fitness
  - `storeLasts`, a dataframe containing the final row of 
