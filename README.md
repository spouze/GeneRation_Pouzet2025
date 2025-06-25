# GeneRation software for Pouzet2025

This repository contains codes for the modified Wagner's model of Gene Regulatory Networks Evolution used in [Pouzet & Le Rouzic 2025 - Evolution](https://academic.oup.com/evolut/advance-article-abstract/doi/10.1093/evolut/qpaf068/8132774): "Gene network topology drives the mutational landscape of gene expression" (previously on [BioRxiv](https://www.biorxiv.org/content/10.1101/2024.11.28.625874v1)).


### Packages required
This work was carried out on R version 4.04 (R Core Team 2021)
- `rlist` (v0.4.6.2) - save and load simulations
- `RColorBrewer` (v1.1.2) - color generation
- `igraph` (v1.2.10) - display networks
- `MASS` (v7.3.53.1) - fits for degree distributions


### Reproduction of the results:
- Make sure R has the packages installed.
- Launch simulations: In bash, use `cd /SIMULATIONS` to move to the correct directory and launch `sh ../0-run_simulations.sh`, which will launch 9 small simulations labeled "TUTO_SIMUS"
- Extract simu features and isolate successful ones (`1-data_extract_from_simus.R`)
- Carry out mutation tests (`2-run_mutation_tests.R`)
- Plot figure 2: check and plot simu, fit and plot degree (`3-Plot_fig2.R`)
- Plot figure 3: compare adaptation profiles (`4-Plot_fig3.R`)
- Plot figure 4: Plot fitness effects (`5-Plot-fig4.R`)
- Plot figure 5: Plot cis-effects and pleiotropy (`6-Plot-fig5.R`)
- Plot figure 6: Enrichment analysis (`7-Plot-fig6.R`)

### Notes: 
- Figure 1 is the materials and methods figures, so not part of the results reproduction.
- Simulations are placed in the `/SIMUALTIONS` folder
- Results extracted from the simulations and mutation tests are placed in the `/Extracted_data` folder
- Generated figures are placed on the `/figures` folder
- The USED_DATA archive contains the data used in the paper, which can be used to recapitulate the figures of the paper based on the above scripts.
- Find more details on data structure and parameter files in the `/src` readme - It contains the 3 sets of parameter files used to obtain the 3 topologies described in the article, as well as a detail account of each parameter used in the model.
- **For a small tutorial and useful functions**, see the [Playground](Notebook_&_useful_functions) folder.</br>
  It contains a step by step template R script to get familiar with the treatment of the data generated from the simulations. Two simulation outputs are also included.
- Each subfolder has a corresponding README which will provide additional detail.
