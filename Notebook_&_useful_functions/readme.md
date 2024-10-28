### Useful functions
Besides the functions necessary for the simuations, `GeneRation_Fun_v1.R`contains functions R functions for in-depth analyses of the simulations. Source the file in R to access the functions.

Using R studio: 
```R
# Load functions
source("GeneRation_Fun_v1.R")

########################
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


########################
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
IndivFinale(simu, gen = 100) # reconstruct the latest average individual from the output dataframe
# Get mean (or median) mtrix from the population (and NOT from the dataframe)
WFromPop(simu[[3]]) == IndivFinale(simu, gen = 100)$ind
# Those two previous operations are equivalent in general.


# integrate the activity (coding vector) to the matrix
W2 = WConvert(W, tfreg = "NOTunique")
# Plotting  the network from W2 (square matrix) doesnt work yet

# Clean up matrix from small values for display or other operations
W3 = WClean(W, threshold = 0.01)
```
