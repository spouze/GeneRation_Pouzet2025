# From SCRIPT 211 ## COPY AND BUILT ON SCRIPT 210. ## Based on SCRIPT 022

setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # Set working directory where the file is.
# getwd()
source("src/GeneRation_Fun_v1.R")
source("src/GeneRation_DataAnalysis_FUN.R")

#################################################################
### EXTRACT DATA FROM SIMULATIONS
#################################################################
# HERE, we will use 3 simulations per topologies, ie, 9 simulations.

simus_folder_path = "SIMULATIONS/" # ends with slash
destination_folder = "Extracted_data/"

## CHOOSE THE RIGHT SUFFIX TO SELECT A SUBSET OF SIMULATIONS: 
batch_name = "TUTO_SIMUS" # for the simus that have been generated as a tuto
#batch_name = "november"   # for the simus coming from the paper


#################################################################

burnin=F # is TRUE if two successive adaptations
all_batch_simus = list.files(paste(simus_folder_path), pattern = batch_name)
# Load the first one to estimate parameters
simu_id = all_batch_simus[1]
simu = LoadSimu(paste0(simus_folder_path, simu_id, "/"))
#fit1000 = simu[[1]]$MFit[1000]
generations = length(simu[[1]]$MFit)
if (burnin){generations = generations*2}
nb_simus_total = length(all_batch_simus)

# Start data collection  ########################################
# CREATE RECEPTACLES #############################
# 1. Create empty dataframe that will receive the fitness across generations
store_Fits_BATCH = data.frame(matrix(ncol=generations+2, nrow = nb_simus_total))
colnames(store_Fits_BATCH) = c("id", "topo", paste0("gen", 1:generations))
# 2. Create empty dataframe that will receive the last line of each simu w params
store_Last_BATCH = data.frame(matrix(ncol=263+10+10+2, nrow = nb_simus_total))
colnames(store_Last_BATCH) = c("id", "topo", colnames(simu[[1]]), paste0("FIT_OPT", 1:10), paste0("FIT_STR",1:10))
# 3. Create empty dataframe that will receive just the fitness of each simu
store_JustLastFit_BATCH = data.frame(matrix(ncol=1+2, nrow = nb_simus_total))
colnames(store_JustLastFit_BATCH) = c("id", "topo", "Fitness of pop")
# START COLLECT #############################
for (simu_id in all_batch_simus) {
  simuE1 = LoadSimu(paste0(simus_folder_path, simu_id, "/"))
  FitE1  = simuE1[[1]]$MFit
  i = match(simu_id, all_batch_simus)
  # 1. Save the last line in store_Last
  store_Last_BATCH[i,] = c(simu_id, 
                           substr(simu_id, 13,18), 
                           simuE1[[1]][generations,],
                           simuE1[[5]]$FITNESS_OPTIMUM,
                           simuE1[[5]]$FITNESS_STRENGTH
  )
  # 2. Save just the last fitness of the pop (at generation 1000, even if burnin)
  store_JustLastFit_BATCH[i,] = c(simu_id, 
                                  substr(simu_id, 13,18),
                                  simuE1[[1]][generations,]$MFit)
  # 3. Save the fitnesses in store_Fitnesses
  if (burnin) {
    simuE2 = LoadSimu(paste0(folder_name, "/", simu_id, "/", simu_id, ".burnin/"))
    FitE2  = simuE2[[1]]$MFit
    store_Fits_BATCH[i,] = c(simu_id, 
                             substr(simu_id, 13,18),
                             FitE1, FitE2) # IF BURNIN
  }else{
    store_Fits_BATCH[i,] = c(simu_id, 
                             substr(simu_id, 13,18), 
                             FitE1)
  }
  cat(paste0(i,".")) #Signal to user that this is done.
}

# SAVE ON DISK #############################
if(burnin){burninornot = "_burnin"}else{burninornot=""}
write.csv(store_Last_BATCH       , file = paste0(destination_folder, "storeLasts_",        batch_name, burninornot))
write.csv(store_Fits_BATCH       , file = paste0(destination_folder, "storeFits_",         batch_name, burninornot))
write.csv(store_JustLastFit_BATCH, file = paste0(destination_folder, "storeJustLastFits_", batch_name, burninornot))

cat("DATA EXTRACTION COMPLETED.\n")
