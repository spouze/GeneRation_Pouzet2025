### Useful functions
Besides the functions necessary for the simuations, `GeneRation_Fun_v1.R`contains functions R functions for in-depth analyses of the simulations. Source the file in R to access the functions.
- `LoadSimu()` loads a simulation, take a folder name as input: `simu=LoadSimu("1023_185639_simulation/")`. This will return a object `simu` as a list containing 5 objects:
  - Output Dataframe (generations 1 to 1000)
  - Initial population (gen 1) as list of individuals (list of lists)
  - Final population (gen 1000) as list of individuals (list of lists)
  - Name of the simulation (string)
  - Parameters used (string)
- `graphall(simu[[1]])` returns the various graphs automatically generated as png in the simulation file.
- `WFinale(simu[[1]])` return the average final W interaction matrix from the output dataframe
- `WKinetics3()` 

```R
# This is an example
mean(1,2)

```
