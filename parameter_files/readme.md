## SCALEF Parameter file details: 

| Parameter               | Value                                                         | Details                                                   |
|-------------------------|---------------------------------------------------------------|-----------------------------------------------------------|
| SIMUL_GENER              | 1000                                                          | Number of simulation generations                           |
| INIT_PSIZE               | 1000                                                          | Population size                                    |
| GENET_NBLOC              | 10                                                            | Number of genes                             |
| INIT_BASAL               | 0.2                                                           | Basal initial condition value (constitutive expression)                              |
| INIT_CLONAL              | clonal                                                        | Indicates if the population is clonal or not (use "notclonal" for not clonal)               |
| INIT_CONDIAG             | 0                                                             | Initial and diagonal condition. If zero, will always remain zero, thus preventing autoregulation   |
| INIT_ALLELES             | 0.0 0.00002                                                   | Initial regulatory allele frequencies (normal law)      |
| INIT_TRANSALLELES        | 0.0 0.00002                                                   | Initial coding-allele frequencies (normal law)                          |
| TYPE_ALLELES             | zero                                                          | NOT USED                                        |
| GENET_MUTRATES           | 0.01                                                          | (cis)-regulatory mutation rate                                      |
| GENET_TRANSMUTRATES      | 0.01                                                          | Coding (activity) mutation rate                                         |
| GENET_MUTSD              | 0.5                                                           | Standard deviation for mutation size                     |
| FITNESS_OPTIMUM          | random                                                        | Optimum fitness values: "random" indicates that the optima (between 0 an 1) will be drawn randomly |
| FITNESS_STRENGTH         | 10 10 10 10 10 0 0 0 0 0                                      | Strength of fitness component (only the first 5 genes under selection)      |
| FITNESS_STABSTR          | 46000                                                         | Strength of stabilizing selection (prevents periodic/cycling networks)    |
| SIMUL_OUTPUT             | 1                                                             | Simulation output mode (if 1, writes the table file)           |
| DEV_TIMESTEPS            | 20                                                            | Developmental time steps to compute phenotype                |
| DEV_CALCSTEPS            | 2                                                             | Number of last developmental step to average to compute phenotype     |
| TF_REG                   | both                                                          | Transcription factor regulation type ("both" for acting both as activator and repressor. Use "unique" for being only activator or repressor. Careful, this changes the biological meaning.    |
| GENET_MUTTYPE            | individual                                                    | NOT USED                      |
| GENET_RECRATES           | 0.5                                                           | NOT USED - hardcoded to 0.5                                         |
| GENET_SELFING            | 0.0                                                           | NOT USED                                               |
| GENET_CLONAL             | 0.0                                                           | NOT USED                                            |
| GENET_PLOIDY             | 2                                                             | NOT USED - hardcoded to 2 (diploid)                                     |
| GENET_EPIGENET           | 0.0                                                           | NOT USED                                    |
| FITNESS_TYPE             | gaussian                                                      | NOT USED - hardcoded to gaussian                                      |
| INIT_CONNECT             | 1                                                             | NOT USED                                 |
| TYPE_SO                  | basal                                                         | NOT USED                             |
| FITNESS_STAB             | exponential_stab                                              | NOT USED                              |
| OUT_UNSTAB               | yes                                                           | NOT USED                          |
| OUT_GENO                 | yes                                                           | NOT USED                                  |
| OUT_CANAL_TESTS          | 0                                                             | NOT USED                                         |
| OUT_CANAL_MUTSD          | 0.5                                                           | NOT USED                   |
| OUT_HERIT_TESTS          | 0                                                             | NOT USED                                         |
| OUT_DIREPI_TESTS         | 0                                                             | NOT USED                                     |
| INIT_RECURRENCE          | 0                                                             | NOT USED                                         |
| FITNESS_FLUCT            | no_fluctuation                                                 | NOT USED                |
| SIMUL_MAXGEN             | 5000                                                          | NOT USED               |
| TYPE_ARCHI               | m2                                                            | NOT USED                     |
