setwd("/Users/sylvain/Documents/DOCUMENTS/GeneRation/M2_PAPER")
source("../GeneRation_Fun.R")
source("R_SCRIPTS/Comfort_FUN.R")
source("R_SCRIPTS/128.New_Mutate_in_new_Env_4.R") # handles cistrans for DUPDEL
source("R_SCRIPTS/IndivRecap_fromSET_FUN.R")

# cnames
dream_SET = read.csv("Extracted_data/dream_SIMUSET_Lasts.csv", header = T)[,-1]

####################################################################
## PERFORM TESTS
####################################################################

###### PARAMS :

simu_nb_start = 1
simu_nb_stop  = nrow(dream_SET) # all simus = 3000

simu_seq = simu_nb_start:simu_nb_stop
simu_nb = length(simu_seq)

number_of_test_per_muttype = 1

reg_muteff_to_test = seq(-1, 1, 0.1)
cod_muteff_to_test = seq(-1, 1, 0.1)

dup_muteff_to_test = 0:10
del_muteff_to_test = 1:5

number_of_batches = 30 # de fichiers
start_batch = 1
end_batch = number_of_batches
batch_seq = start_batch:end_batch

destination_folder = "Extracted_data/" # avec "/" à la fin 
check_folder_exists(destination_folder) # Check if the folder exists
  
##########################################
{# PARAMETERS FOR FOR CIS/TRANS #########
  simu_nb_start = 1  
  simu_nb_stop  = nrow(dream_SET)
  simu_seq = simu_nb_start:simu_nb_stop
  simu_nb = length(simu_seq)
  dup_muteff_to_test = NA
  del_muteff_to_test = NA
  reg_muteff_to_test = 0.5 #seq is 21
  cod_muteff_to_test = 0.5
  number_of_test_per_muttype = 2
  
  number_of_batches = 20
  start_batch = 62
  end_batch = 200
  batch_seq = start_batch:end_batch
  
  destination_folder = "Extracted_data/allmutations_20240320_cistrans/"
  check_folder_exists(destination_folder) # Check if the folder exists
} ######################################
##########################################

##########################################
{# PARAMETERS FOR FOR FIG 06 ############ 20240417
  simu_nb_start = 1
  simu_nb_stop  = nrow(dream_SET)
  simu_seq      = simu_nb_start:simu_nb_stop
  simu_nb       = length(simu_seq)
  
  dup_muteff_to_test = 0:10
  del_muteff_to_test = 1:5
  reg_muteff_to_test = seq(-1, 1, 0.1)
  cod_muteff_to_test = seq(-1, 1, 0.1)
  
  number_of_test_per_muttype = 1
  
  batch_seq = c(22:(22+28))
  
  destination_folder = "Extracted_data/allmutations_20240417_fig06/"
  check_folder_exists(destination_folder)
  
} ######################################
##########################################

## START TESTS ##################################################################

for (batch_x in batch_seq){ # start all batches
  
  ###### INIT :
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
          allmutreg = rbind(allmutreg, Mutate_in_new_environment(ind, "REG", reg_muteff = reg_muteff, id = counter))
        }
      }
      
      # 2. CODING MUTATIONS
      if (is.na(cod_muteff_to_test)){}else{
        allmutcod =  data.frame(matrix(ncol = length(cnames), nrow = 0))
        colnames(allmutcod)=cnames
        for (cod_muteff in cod_muteff_to_test){
          allmutcod = rbind(allmutcod, Mutate_in_new_environment(ind, "COD", cod_muteff = cod_muteff, id = counter))
        }
      }

      
      # 3. GENE DUPLICATIONS
      if (is.na(sum(dup_muteff_to_test))){ # on met la sum() pour pas que ca bug avec if(is.na(1:10)){} 20240417
        allmutdup =  data.frame(matrix(ncol = length(cnames), nrow = 0))
      }else{
        allmutdup =  data.frame(matrix(ncol = length(cnames), nrow = 0))
        colnames(allmutdup)=cnames
        for (dup_number in dup_muteff_to_test){ # because 10 genes
          allmutdup = rbind(allmutdup, Mutate_in_new_environment(ind, "DUP", dupdel_number = dup_number, id = counter))
        }
      }
      
      # 4. GENE DELETIONS
      if (is.na(sum(del_muteff_to_test))){
        allmutdel =  data.frame(matrix(ncol = length(cnames), nrow = 0))
      }else{
        allmutdel =  data.frame(matrix(ncol = length(cnames), nrow = 0))
        colnames(allmutdel)=cnames
        for (del_number in del_muteff_to_test){ # because 10 genes
          allmutdel = rbind(allmutdel, Mutate_in_new_environment(ind, "DEL", dupdel_number = del_number, id = counter))
        }
      }
      
      
      allmutations = rbind(allmutations, allmutreg, allmutcod, allmutdup, allmutdel)
      
    }
    final_time = Sys.time()
    cat(round(as.numeric(final_time - init_time, units = "secs"),3))
    cat(paste0(" secs.\n"))
    
  };{
    allmutations$id = 1:nrow(allmutations)
    play()
    
    very_final_time = Sys.time()
    cat(paste0("\nTotal time: ", 
               round(as.numeric(very_final_time - very_init_time, units = "mins"),3),
               " minutes.\n"))
  } # end batch x
  
  ## SAVING :
  write.csv(allmutations, 
            paste0(destination_folder, "allmutations", sprintf("%03d",batch_x), ".csv"))

} # end all batches

####################################################################
## END OF "PERFORM TESTS"
####################################################################


###################################################
###################################
# CHECK 
allmutations=read.csv("Extracted_data/allmutations_20230314/allmutations01.csv")[,-1]

nrow(allmutations)
#allmutations$topology
nrow(allmutations[allmutations$notes=="NO MUTATION",]) # -> ou on veut les garder pour les barres à zero ! 
allmutations2 = allmutations[allmutations$notes!="NO MUTATION",]
nrow(allmutations2)


##########??
# Faudrait enlever celles ou ya "NO MUTATION" et ou c'est pas ZERO dans un des muteff.

subselect1 = allmutations[allmutations$fit_effect!="NEUTR" & allmutations$affected_genes==0 ,]
nrow(subselect1)
