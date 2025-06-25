# from SCRIPT 211

setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # Set working directory where the file is.
# getwd()
source("src/GeneRation_Fun_v1.R")
source("src/GeneRation_DataAnalysis_FUN.R")

#################################################################
### LOAD EXTRACTED DATA
#################################################################
destination_folder = "Extracted_data/"
folder = "Extracted_data/"
list.files(folder, pattern = "storeFits")

batch = "TUTO_SIMUS" # for the exemple for the script
# batch = "november" # on a few real simulations from the paper

three_files = list.files(folder, pattern = paste0(batch, "$"))

storeFits         = read.csv(paste0(folder, three_files[1]))
storeJustLastFits = read.csv(paste0(folder, three_files[2]))
storeLasts        = read.csv(paste0(folder, three_files[3]))

### SPLIT
storeFits_HIGHCO = storeFits[storeFits$topo=="HIGHCO",]
storeFits_RANDOM = storeFits[storeFits$topo=="RANDOM",]
storeFits_SCALEF = storeFits[storeFits$topo=="SCALEF",]

storeJustLastFits_HIGHCO = storeJustLastFits[storeFits$topo=="HIGHCO",]
storeJustLastFits_RANDOM = storeJustLastFits[storeFits$topo=="RANDOM",]
storeJustLastFits_SCALEF = storeJustLastFits[storeFits$topo=="SCALEF",]

storeLasts_HIGHCO = storeLasts[storeLasts$topo=="HIGHCO",]
storeLasts_RANDOM = storeLasts[storeLasts$topo=="RANDOM",]
storeLasts_SCALEF = storeLasts[storeLasts$topo=="SCALEF",]

#################################################################
### Successful simulations
#################################################################
## If I want a least 200...
# sum(storeJustLastFits_HIGHCO$Fitness.of.pop>0.9)
# sum(storeJustLastFits_RANDOM$Fitness.of.pop>0.9)
# sum(storeJustLastFits_SCALEF$Fitness.of.pop>0.9)
# which(storeJustLastFits_HIGHCO$Fitness.of.pop>0.9)

#success_threshold = 0.95
success_threshold = 0 # Set to zero here so that the code runs even with dummy simulations

dream_Fits_HIGHCO = storeFits_HIGHCO[which(storeJustLastFits_HIGHCO$Fitness.of.pop>success_threshold),]
dream_Fits_RANDOM = storeFits_RANDOM[which(storeJustLastFits_RANDOM$Fitness.of.pop>success_threshold),]
dream_Fits_SCALEF = storeFits_SCALEF[which(storeJustLastFits_SCALEF$Fitness.of.pop>success_threshold),]

dream_Lasts_HIGHCO = storeLasts_HIGHCO[which(storeJustLastFits_HIGHCO$Fitness.of.pop>success_threshold),]
dream_Lasts_RANDOM = storeLasts_RANDOM[which(storeJustLastFits_RANDOM$Fitness.of.pop>success_threshold),]
dream_Lasts_SCALEF = storeLasts_SCALEF[which(storeJustLastFits_SCALEF$Fitness.of.pop>success_threshold),]

# round(x, 2)
cat(paste0(
  "HIGHCO : ", nrow(dream_Lasts_HIGHCO), " sucessful ", "out of ", nrow(storeLasts_HIGHCO), ", ie ", round(nrow(dream_Lasts_HIGHCO)/nrow(storeLasts_HIGHCO)*100, 0), "%\n", 
  "RANDOM : ", nrow(dream_Lasts_RANDOM), " sucessful ", "out of ", nrow(storeLasts_RANDOM), ", ie ", round(nrow(dream_Lasts_RANDOM)/nrow(storeLasts_RANDOM)*100, 0), "%\n",
  "SCALEF : ", nrow(dream_Lasts_SCALEF), " sucessful ", "out of ", nrow(storeLasts_SCALEF), ", ie ", round(nrow(dream_Lasts_SCALEF)/nrow(storeLasts_SCALEF)*100, 0), "%\n"
))



#################################################################
#################################################################
#################################################################
### Carry out tests
#################################################################
#################################################################
#################################################################

# Let's take only the first 100 to start with
# n = 100
# dream_Lasts_HIGHCO = dream_Lasts_HIGHCO[c(1:n),]
# dream_Lasts_RANDOM = dream_Lasts_RANDOM[c(1:n),]
# dream_Lasts_SCALEF = dream_Lasts_SCALEF[c(1:n),]

# BASED ON SCRIPT 073
dream_SET = rbind(dream_Lasts_HIGHCO,
                 dream_Lasts_RANDOM,
                 dream_Lasts_SCALEF)
str(dream_SET)

{# PARAMETERS TO CARRY OUT THE TESTS
  simu_nb_start = 1  
  simu_nb_stop  = nrow(dream_SET)
  simu_seq = simu_nb_start:simu_nb_stop
  simu_nb = length(simu_seq)
  
  dup_muteff_to_test = c(1:5)
  del_muteff_to_test = c(1:5)
  reg_muteff_to_test = seq(-1,1,0.1) # c(0.5, -0.5) #seq is 21 seq(-1,1,0.1)
  cod_muteff_to_test = seq(-1,1,0.1) # c(0.5, -0.5)
  
  number_of_test_per_muttype = 20
  
  number_of_batches = 1 # nombre de fichiers
  start_batch = 1
  end_batch = number_of_batches
  batch_seq = start_batch:end_batch
  
  destination_folder = "Extracted_data/"
  check_folder_exists(destination_folder) # Check if the folder exists
  
} ######################################

## START TESTS ##################################################################
for (batch_x in batch_seq){ # start all batches
  # INIT :
  lines_number = number_of_test_per_muttype*(
    length(na.omit(reg_muteff_to_test))+
      length(na.omit(cod_muteff_to_test))+
      length(na.omit(dup_muteff_to_test))+
      length(na.omit(del_muteff_to_test))
  )*simu_nb
  
  allmutations = data.frame(matrix(ncol = length(cnames), nrow = 0))
  colnames(allmutations)=cnames
  very_init_time = Sys.time()
  counter=0
  
  for (selected_simu in simu_seq){ # start batch x
    counter=counter+1
    ind = IndivRecap(storeLasts_set = dream_SET, line = selected_simu, clean = F)
    cat(paste0(sprintf("%04d", counter), " / ", sprintf("%04d", simu_nb), " - "))
    init_time = Sys.time()
    for (test_number in 1:number_of_test_per_muttype){ # Choose the number of tests for each simu 
      
      # 1. REGULATORY MUTATIONS
      if (is.na(reg_muteff_to_test)){}else{
        allmutreg =  data.frame(matrix(ncol = length(cnames), nrow = 0))
        colnames(allmutreg)=cnames
        for (reg_muteff in reg_muteff_to_test){
          allmutreg = rbind(allmutreg, Mutate_in_new_environment(ind, "REG", reg_muteff = reg_muteff, id = counter))}}
      
      # 2. CODING MUTATIONS
      if (is.na(cod_muteff_to_test)){}else{
        allmutcod =  data.frame(matrix(ncol = length(cnames), nrow = 0))
        colnames(allmutcod)=cnames
        for (cod_muteff in cod_muteff_to_test){
          allmutcod = rbind(allmutcod, Mutate_in_new_environment(ind, "COD", cod_muteff = cod_muteff, id = counter))}}
      
      # 3. GENE DUPLICATIONS
      if (is.na(sum(dup_muteff_to_test))){ # on met la sum() pour pas que ca bug avec if(is.na(1:10)){} 20240417
        allmutdup =  data.frame(matrix(ncol = length(cnames), nrow = 0))
      }else{
        allmutdup =  data.frame(matrix(ncol = length(cnames), nrow = 0))
        colnames(allmutdup)=cnames
        for (dup_number in dup_muteff_to_test){ # because 10 genes
          allmutdup = rbind(allmutdup, Mutate_in_new_environment(ind, "DUP", dupdel_number = dup_number, id = counter))}}
      
      # 4. GENE DELETIONS
      if (is.na(sum(del_muteff_to_test))){
        allmutdel =  data.frame(matrix(ncol = length(cnames), nrow = 0))
      }else{
        allmutdel =  data.frame(matrix(ncol = length(cnames), nrow = 0))
        colnames(allmutdel)=cnames
        for (del_number in del_muteff_to_test){ # because 10 genes
          allmutdel = rbind(allmutdel, Mutate_in_new_environment(ind, "DEL", dupdel_number = del_number, id = counter))}}
      
      allmutations = rbind(allmutations, allmutreg, allmutcod, allmutdup, allmutdel)}
    final_time = Sys.time()
    cat(round(as.numeric(final_time - init_time, units = "secs"),3))
    cat(paste0(" secs.\n"))
  };{
    allmutations$id = 1:nrow(allmutations)
    
    very_final_time = Sys.time()
    cat(paste0("\nTotal time: ", 
               round(as.numeric(very_final_time - very_init_time, units = "mins"),3),
               " minutes.\n"))
  } # end batch x
  ## SAVING :
  write.csv(allmutations, 
            paste0(destination_folder, "allmutations_", format(very_init_time, "%Y%m%d%H%M%S"), "_",  sprintf("%03d",batch_x), "_", batch, ".csv"))
  cat("MUTATION TESTS COMPLETED.\n")
  
} # end all batches ## END OF "PERFORM TESTS"

####### IT'S OKAY IF THERE ARE MANY WARNINGS.
nrow(allmutations) / nrow(dream_SET)
