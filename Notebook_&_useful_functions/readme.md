### Useful functions
Besides the functions necessary for the simuations, `GeneRation_Fun_v1.R`contains functions R functions for in-depth analyses of the simulations. Source the file in R to access the functions.

You will find here:
- the R script to use as template or exemple for further analyses
- a folder containing a simple default simulation (5 genes, 100 individuals, 100 generations)
- a folder containing an actual simulation used for the paper (10 genes, 1000 individuals, 1000 generations)

Using R studio (this is the content of the above "Tutorial" file):
```R
setwd("Documents/DOCUMENTS/GeneRation_Pouzet2024/Notebook_&_useful_functions/")

# Load functions
source("GeneRation_Fun_v1.R")
# Load required packages
library(rlist, RColorBrewer, igraph)

#########################################################
################################################
## Basic functions

# Load the simulation (default with 5 genes, 100 individuals, 100 generations)
simu=LoadSimu("1028_152647_simulation/")
# simu is a list containing 5 items:
# - Output Dataframe (generations 1 to 100)
nrow(simu[[1]]) # 100 generations
# - Initial population (gen 1) as list of individuals (list of lists)
length(simu[[2]]) # 100 individuals
simu[[2]][[1]] # individual #1 of intial population (low fitness)
# - Final population (gen 100) as list of individuals (list of lists)
length(simu[[3]]) # 100 individuals
simu[[3]][[1]] # individual #1 of intial population (low fitness)
# - Name of the simulation (string)
simu[[4]]
# - Parameters used (string)
simu[[5]]

# Check simu
graphall(simu[[1]])
dev.off()

# Focus on the final obtained matrix
# Starting with the initial avergae network:
W=WFinale(simu[[1]], gen = 1) # return the average final W interaction matrix from the output dataframe
WNetwork4(W, speak = T) # plot the network #error message is okay
# And the final aerage network:
W=WFinale(simu[[1]], gen = 100)
WNetwork4(W)
WTab3(W, cellnote = T)

# Look at the dynamics (developmental steps)
dev.off() # clean up plotting area first
WKinetics3(W, optima = simu[[5]]$FITNESS_OPTIMUM[1:5])


#########################################################
################################################
## Advanced functions

# Check matrix as points
W=WFinale(simu[[1]], gen = 100) # Final individual
WPlot(W) # Final average individual
# Check allele distribution in the population
WBoxplot(simu[[3]], stripchart = T) # all individuals of the final population

# Check Degree from matrix
WDegree(W)

# Single plotting functions
graphfit(simu[[1]])
graphexp(simu[[1]])
graphmtrxvalues(simu[[1]])
graphtransvector(simu[[1]])

# get average individual:
last_ind = IndivFinale(simu, gen = 100) # reconstruct the latest average individual from the output dataframe
# Get mean (or median) mtrix from the population (and NOT from the dataframe)
W = WFromPop(simu[[3]]) == IndivFinale(simu, gen = 100)$ind
# Those two previous operations are equivalent in general.


# integrate the activity (coding vector) to the matrix
W2 = WConvert(W, tfreg = "NOTunique")
# Plotting  the network from W2 (square matrix) doesnt work yet

# Clean up matrix from small values for display or other operations
W3 = WClean(W, threshold = 0.01)


#########################################################
################################################
## Mutations functions
# The scripts and function and hardcoded for Simulations with 10 genes.

# Get info from the simulation
# Directly use 022.Executable_Extract2.R ; launch with bash 
# Or by hand for one single simulation here: 
simu = LoadSimu("0614_014837_SCALEF_027_november/")
check = data.frame(matrix(ncol=263+10+10+2, nrow = 1))
colnames(check) = c("id", "topo", colnames(simu[[1]]), paste0("FIT_OPT", 1:10), paste0("FIT_STR",1:10))
check[1,] = c(simu[[4]], 
              substr(simu[[4]], 13,18), 
              simu[[1]][nrow(simu[[1]]),],
              simu[[5]]$FITNESS_OPTIMUM,
              simu[[5]]$FITNESS_STRENGTH)

# Need a specially formatted individual:
source("IndivRecap_fromSET_FUN.R") # in the "Data_Analysis_scripts" folder
new_ind = IndivRecap(check)

# Carry out a mutation test
source("128.New_Mutate_in_new_Env_4.R") # in the "Data_Analysis_scripts" folder
# Regulatory mutations (+0.5)
Mutate_in_new_environment(ind = new_ind, mut_type = "REG", reg_muteff = 0.5, id = 1)
# Coding mutations (-0.1)
Mutate_in_new_environment(ind = new_ind, mut_type = "COD", cod_muteff = -0.1, id = 2)
# Duplications (+1)
Mutate_in_new_environment(ind = new_ind, mut_type = "DUP", dupdel_number = 1, id = 3)
# Deletions (-2)
Mutate_in_new_environment(ind = new_ind, mut_type = "DEL", dupdel_number = 2, id = 4)
# The output includes:
# - an id (usually the counter)
# - the output name of the simu
# - the topology (here specifically formatted)
# - the mutation type ( REG COD DUP DELL)
# - the location of the mutation (ij)
# - the value of the mutation (init, final, delta)
# - the deplication/deletion number and genes concerned
# - Fitness info (population in first environement (E1), individual (ind) in E1, ind in second environment (E2), delta, effect, threshold)
# - Pleiotropy info (pleiotropy and threshold)
# - record of the mutation reg_muteff or cod_muteff
# - An impacted gene randomly picked: number and expression change
# - The max impacted gene : number and expression change
# - The min impacted gene : number and expression change
# - The cis-effect check (0 or 1).
```
