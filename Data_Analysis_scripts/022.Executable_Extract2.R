
#### INSTRUCTIONS ##############
# Lancer:
# Rscript 022.Executable_Extract_wBurnin.R -f 20240520_1x1500_burnin_AUTU1/SIMUS/
#                                          -p AUTU (poour éviter de prendre les error files)
#                                          -b TRUE (si burnin or not)
# The folder name is planned for to remove the "/SIMUS/ at the end


# 202405211030
source("GeneRation_Fun.R") # depends d'ou le script est lancé

cmd = commandArgs(trailingOnly=FALSE)

# path
my_path = cmd[grep(cmd, pattern="--file=")]
my_path = dirname(strsplit(my_path, split='=')[[1]][2])
#my_path = getwd()

#date=format(Sys.time(), "%m%d_%H%M%S_")

# -f fold = er_directory_target
which_f = which(cmd=="-f")
if (length(which_f)==0) {
  stop("You need to specify a folder name with -f folder_name/ and possibly a pattern with -p folder_pattern.\n")
  # Use cat() to display a console message.
} else {
  folder_name = cmd[which_f+1] 
  #folder_name = "20240416_test_AUTUMN/SIMUS/"
}

# -p pattern_to_check
which_p = which(cmd=="-p")
if (length(which_p)==0) {
  folder_pattern = ""
} else {
  folder_pattern = cmd[which_p+1]
}

# -b chech for Burnin
which_b = which(cmd=="-b")
input_b = toupper(cmd[which_b+1]) #toupper to make is case insensitive in a way
if (length(which_b)==0) {
  burnin = F
} else {
  if (input_b %in% c("T", "TRUE")){
    burnin = T
  } else if (input_b %in% c("F", "FALSE")){
    burnin = F
  } else {
    burnin = F
  }
}


batch_name = substr(folder_name, 1, nchar(folder_name)-7)
all_batch_simus = list.files(paste(my_path, folder_name, sep="/"), pattern = folder_pattern)

# Load the first one to estimate parameters
simu_id = all_batch_simus[1]
simu = LoadSimu(paste0(folder_name, "/", simu_id, "/"))
#fit1000 = simu[[1]]$MFit[1000]
generations = length(simu[[1]]$MFit)
if (burnin){generations = generations*2}
nb_simus_total = length(all_batch_simus)
cat(paste0(format(Sys.time(), "%Y%m%d %H:%M:%S - "), nb_simus_total, " simulations: \n"))

########################################
# Start data collection 
########################################

# CREATE RECEPTACLES
#############################
# 1. Create empty dataframe that will receive the fitness across generations
store_Fits_BATCH = data.frame(matrix(ncol=generations+2, nrow = nb_simus_total))
colnames(store_Fits_BATCH) = c("id", "topo", paste0("gen", 1:generations))
# 2. Create empty dataframe that will receive the last line of each simu w params
store_Last_BATCH = data.frame(matrix(ncol=263+10+10+2, nrow = nb_simus_total))
colnames(store_Last_BATCH) = c("id", "topo", colnames(simu[[1]]), paste0("FIT_OPT", 1:10), paste0("FIT_STR",1:10))
# 3. Create empty dataframe that will receive just the fiteness of each simu
store_JustLastFit_BATCH = data.frame(matrix(ncol=1+2, nrow = nb_simus_total))
colnames(store_JustLastFit_BATCH) = c("id", "topo", "Fitness of pop")

# START COLLECT
#############################
for (simu_id in all_batch_simus) {
  
  simuE1 = LoadSimu(paste0(folder_name, "/", simu_id, "/"))
  FitE1  = simuE1[[1]]$MFit
  
  i = match(simu_id, all_batch_simus)
  
  # 1. Save the last line in store_Last
  store_Last_BATCH[i,] = c(simu_id, 
                           substr(simu_id, 13,18), 
                           simuE1[[1]][1000,],
                           simuE1[[5]]$FITNESS_OPTIMUM,
                           simuE1[[5]]$FITNESS_STRENGTH
  )
  
  # 2. Save just the last fitness of the pop (at generation 1000, even if burnin)
  store_JustLastFit_BATCH[i,] = c(simu_id, 
                                  substr(simu_id, 13,18),
                                  simuE1[[1]][1000,]$MFit)
  
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

# SAVE ON DISK
#############################
cat("\nEXTRACTION COMPLETED\n")
if(burnin){burninornot = "_burnin"}else{burninornot=""}

write.csv(store_Last_BATCH       , file = paste0("storeLasts_",        batch_name, burninornot, "_", folder_pattern))
write.csv(store_Fits_BATCH       , file = paste0("storeFits_",         batch_name, burninornot, "_", folder_pattern))
write.csv(store_JustLastFit_BATCH, file = paste0("storeJustLastFits_", batch_name, burninornot, "_", folder_pattern))


######################################################################################
####################################################################################


# system(paste0("echo \"", paste(cmd, collapse = " "), "\" > output_file.txt"))
# system(paste0("echo \"", my_path, "\" >> output_file.txt"))
# system(paste0("echo \"", folder_name, " ", folder_pattern, "\" >> output_file.txt"))
#system(paste0("echo \"", "something here", "\" >> output_file.txt"))
cat(paste0("JOB COMPLETED ",format(Sys.time(), "%Y%m%d at %H:%M:%S"), "\n"))

# test avec burnin

# # test look at the files
# check = read.csv("storeFits_20240416_test_AUTUMN")[,-1]
# check = read.csv("storeJustLastFits_20240416_test_AUTUMN")[,-1]
# sum(check$Fitness.of.pop>.95)/nrow(check)
# check = read.csv("storeLasts_20240416_test_AUTUMN")[,-1]
# dim(check)



