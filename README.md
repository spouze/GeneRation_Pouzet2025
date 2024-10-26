# GeneRation software for Pouzet2024

This repository contains codes for the modified Wagner's model of Gene Regulatory Networks Evolution used in [Pouzet & Le Rouzic 2024](https://github.com/spouze/GeneRation_Pouzet2024/): "Gene network topology drives the mutational landscape of gene expression".

### Files
- `Launch_GeneRation.R` - The Launcher for adaptation in a single environment
- `GeneRation_Fun_v1.R` - Functions, heart of the program
- Suite of R scripts for data analyses (`data_analysis_scripts` folder)
- Parameter files used for the paper (`parameter_files` folder)

### Packages required
This work was carried out on R version 4.04 (R Core Team 2021)
- `rlist` (v0.4.6.2) - save and load simulations
- `RColorBrewer` (v1.1.2) - color generation
- `igraph` (v1.2.10) - display networks

### Simple launch
with `Launch_GeneRation.R` and `GeneRation_Fun_v1.R` in the same folder.
```
./Launch_GeneRation.R
```
will launch a stand-alone default simulation. </br>
will return a dated folder `MMDD_HHMMSS_simulation` containing 8 files:
- `MMDD_HHMMSS_simulation.00.graphall.png` - Various graphs summarizing the evolution using the population mean: fitness, phenotype and genotype (regulatory and coding).
- `MMDD_HHMMSS_simulation.01.WInit.png` - Initial mean individual's genotype (table and network)
- `MMDD_HHMMSS_simulation.02.WFinal.png` - Evolved mean individual's genotype (table and network)
- `MMDD_HHMMSS_simulation.03.WFinale_Kinetics.png` - Evolved mean individual development and phenotype
- `MMDD_HHMMSS_simulation.finalpop.rds` - Final (evolved) population Rlist
- `MMDD_HHMMSS_simulation.initialpop.rds` - Initial population Rlist - clonal population in our simulations
- `MMDD_HHMMSS_simulation.param` - Parameter file used for evolution
- `MMDD_HHMMSS_simulation.table` - Generation by generation evolution table summary: id, generation, mean population (pop) phenotype, mean pop genotype, mean pop fitness

### Specific launch
Note that all launch parameters are independant and facultative.
```
./Launch_GeneRation.R -p param_file -o output_name -ipop MMDD_HHMMSS_simulation.finalpop.rds
```
- `-p` for the specific parameter file to be used - default is generated automatically, 5 genes, 2 under selection, 100 generations for a population of 100 individuals.
- `-o` for the specified output name - default is "simulation".
- `-ipop` for a specific .rds population to use as first generation for the subsequent evolution (often the final population of another simulation) - default is an automatically generated heterogeneous (non-clonal) population.
- `-fun` for a specific functions' file - default is `GeneRation_Fun_v1.R`.


### Other notes
- Parameter files detail: see "parameter_files" folder.
- Reproduction of the results / data pipeline : see "data_analysis_scripts" folder
- Launch Burnin - using `Launch_BurninGeneRation.R`: a first simulation is launched, then two expression optima are changed in a new param files (new second environment) and a new folder is created in the very folder containing the simulation for adaptation to the first environment.

---

### License
<a rel="license" href="http://creativecommons.org/licenses/by-nc-sa/4.0/"><img alt="Licence Creative Commons" style="border-width:0" src="https://i.creativecommons.org/l/by-nc-sa/4.0/88x31.png" /></a><br />Ce(tte) œuvre est mise à disposition selon les termes de la <a rel="license" href="http://creativecommons.org/licenses/by-nc-sa/4.0/">Licence Creative Commons Attribution - Pas d’Utilisation Commerciale - Partage dans les Mêmes Conditions 4.0 International</a>.
