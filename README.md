# GeneRation software for Pouzet2024

This repository contains codes for the modified Wagner's model of Gene Regulatory Networks Evolution used in [Pouzet & Le Rouzic 2024](https://github.com/spouze/GeneRation_Pouzet2024/): "Gene network topology drives the mutational landscape of gene expression".

### Files
- `Launch_GeneRation.R` - The Launcher for adaptation in a single environment
- `Launch_BurninGeneRation.R` - The Launcher for two simulations: adaptation in two successive environements
- `GeneRation_Fun.R` - Functions, heart of the program
- Suite of R scripts for data analyses (see Folder)

### Packages required
All carried ou  on R version 4.04 (R Core Team 2021)
- `rlist` (v0.4.6.2) - save and load simulations
- `RColorBrewer` (v1.1.2) - color generation
- `igraph` (v1.2.10) - display networks

### Simple launch
```
./Launch_GeneRation.R
```
will launch a stand-alone default simulation
will return a dated folder MMDD_HHMMSS_simulation containing:
- `MMDD_HHMMSS_simulation.00.graphall.png` - Various graphs summarizing the evolution         
- `MMDD_HHMMSS_simulation.01.WInit.png` - Initial genotype (table and network)             
- `MMDD_HHMMSS_simulation.01.WInit.png` - Evolved genotype (table and network)             
- `MMDD_HHMMSS_simulation.03.WFinale_Kinetics.png` - Evolved devlopment and phenotype                 
- `MMDD_HHMMSS_simulation.finalpop.rds` - Evolved final population Rlist                   
- `MMDD_HHMMSS_simulation.initialpop.rds` - Initial final population Rlist                   
- `MMDD_HHMMSS_simulation.param` - Parameter file used for evolution                
- `MMDD_HHMMSS_simulation.table` - Generation by generation evolution table summary 



### Specific launch
```
./Launch_GeneRation.R -p param_file -o output_name -ipop MMDD_HHMMSS_simulation.finalpop.rds
```
- `-p` for the specific parameter file to be used
- `-o` for the specified output name
- `-ipop` for a specific .rds population to use as first generation for the subsequent evolution (often the final population of another simulation)

### Parameter file
- Use `RANDOM` for random optima

### Burnin

### License
<a rel="license" href="http://creativecommons.org/licenses/by-nc-sa/4.0/"><img alt="Licence Creative Commons" style="border-width:0" src="https://i.creativecommons.org/l/by-nc-sa/4.0/88x31.png" /></a><br />Ce(tte) œuvre est mise à disposition selon les termes de la <a rel="license" href="http://creativecommons.org/licenses/by-nc-sa/4.0/">Licence Creative Commons Attribution - Pas d’Utilisation Commerciale - Partage dans les Mêmes Conditions 4.0 International</a>.
