setwd("Documents/DOCUMENTS/GeneRation/M2_PAPER/R_SCRIPTS/")
# NOTE: gerer launcher_temporary.txt et other launchers avec the ATOM logiciel. easier.

# General stuff:
#prefix = "system(\"Rscript "
#suffix = "\", intern = TRUE)" # Ca c'était avec read_line.R
prefix = "Rscript "
suffix = ""
path = "/shared/projects/evoplanet/Sylvain/"
exec = "Execution_files/"

# Names
folder = "20241017_controls/"
output_nametag = "a0.5"

#model = "Launch_BurninGeneRation#3.R"
model = "Launch_GeneRation.R"

# Parameters:
output_topo = "SCALEF"
#param = "quickparam.txt"
#param = "param_RANDOM.txt"
#param = "param_HIGHCO.txt"
#param = "param_SCALEF.txt"
#"param_AUTUMN.txt"
#param = "param_SCALEF_7_3.txt"
param = "param_SCALEF_a_0.5.txt"

# Number of simus
output_numbers = seq(1, 500, 1)
output_numbers
length(output_numbers)


## CA PART #########################################

for_sprintf = paste0("%0", nchar(as.character(length(output_numbers))), "d") #number of 000 to add
for (number in output_numbers){
  check = paste0(prefix, path, folder, exec, model, 
                 " -p ", path, folder, exec, "params/", param,
                 " -o ", output_topo, "_", sprintf(for_sprintf, number), "_", output_nametag, 
                 suffix)
  system(paste0("echo ", "\'", check, "\'", " >> launcher_temporary.txt"))
} ; play()




######################################### 
### testing and check
#system("echo 'coucou' >> launcher_temporary.txt") # add the carriage return automatically

# ## CHECK 
# example = "system('Rscript /shared/projects/evoplanet/Sylvain/20240518_slurmtest_AUTUMN/Execution_files/Launch_GeneRation.R -p /shared/projects/evoplanet/Sylvain/20240518_slurmtest_AUTUMN/Execution_files/params/param_AUTUMN.txt -o AUTUMN_051_cluster1stFullSimu', intern = TRUE)"
# example ; check
# ## Check character by character
# library(stringr)
# example_chars <- str_split(example, "", simplify = TRUE)
# check_chars <- str_split(check, "", simplify = TRUE)
# min_length <- min(length(example_chars), length(check_chars))
# comparison <- example_chars[1:min_length] == check_chars[1:min_length]
# print(comparison)
