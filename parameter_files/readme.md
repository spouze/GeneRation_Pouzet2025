## SCALEF Parameter file details: 
These are the parameters used to obtain the scale-free topology in the paper.

| Parameter               | Value                                                         | Details                                                   |
|-------------------------|---------------------------------------------------------------|-----------------------------------------------------------|
| SIMUL_GENER              | 1000                                                          | Number of simulation generations                           |
| INIT_PSIZE               | 1000                                                          | Population size                                    |
| GENET_NBLOC              | 10                                                            | Number of genes                             |
| INIT_BASAL               | 0.2                                                           | Default basal expression, when no regulation (constitutive expression)                              |
| INIT_CLONAL              | clonal                                                        | Indicates if the initial population is clonal or not (use "notclonal" for not clonal)               |
| INIT_CONDIAG             | 0                                                             | Initial and diagonal condition. If zero, will always remain zero, thus preventing autoregulation   |
| INIT_ALLELES             | 0.0 0.00002                                                   | Initial regulatory allele values (normal law: mean and sd)      |
| INIT_TRANSALLELES        | 0.0 0.00002                                                   | Initial coding-allele (normal law: mean and sd)                          |
| TYPE_ALLELES             | zero                                                          | _NOT USED_                                        |
| GENET_MUTRATES           | 0.01                                                          | (cis)-regulatory mutation rate                                      |
| GENET_TRANSMUTRATES      | 0.01                                                          | Coding (activity) mutation rate                                         |
| GENET_MUTSD              | 0.5                                                           | Standard deviation for mutation size                     |
| FITNESS_OPTIMUM          | random                                                        | Optimum fitness values: "random" indicates that the optima (between 0 an 1) will be drawn randomly. If written by hand, There should be as many values as there are genes, i.e., 10 values for 10 genes  - if not, the first value is used as default for missing values. |
| FITNESS_STRENGTH         | 10 10 10 10 10 0 0 0 0 0                                      | Strength of stabilizing selection (only the first 5 genes under selection (selection strenght is 10), the last 5 are free). There should be as many values as there are genes, i.e., 10 values for 10 genes - if not, the first value is used as default for missing values.      |
| FITNESS_STABSTR          | 46000                                                         | Strength of selection on stability (prevents periodic/cycling networks)    |
| SIMUL_OUTPUT             | 1                                                             | Frequency of the output files (if 1, writes at every generation)           |
| DEV_TIMESTEPS            | 20                                                            | Developmental time steps to compute phenotype                |
| DEV_CALCSTEPS            | 2                                                             | Number of last developmental step to average to compute phenotype     |
| TF_REG                   | both                                                          | Transcription factor regulation type ("both" for acting both as activator and repressor. Use "unique" for being only activator or repressor. Careful, this changes the biological meaning.    |
| GENET_MUTTYPE            | individual                                                    | _NOT USED_                      |
| GENET_RECRATES           | 0.5                                                           | _NOT USED_ - recombination rate hardcoded to 0.5                                         |
| GENET_SELFING            | 0.0                                                           | _NOT USED_                                               |
| GENET_CLONAL             | 0.0                                                           | _NOT USED_                                            |
| GENET_PLOIDY             | 2                                                             | _NOT USED_ - ploidy hardcoded to 2 (diploid)                                     |
| GENET_EPIGENET           | 0.0                                                           | _NOT USED_                                    |
| FITNESS_TYPE             | gaussian                                                      | _NOT USED_ - fitness type hardcoded to gaussian                                      |
| INIT_CONNECT             | 1                                                             | _NOT USED_                                 |
| TYPE_SO                  | basal                                                         | _NOT USED_                             |
| FITNESS_STAB             | exponential_stab                                              | _NOT USED_                              |
| OUT_UNSTAB               | yes                                                           | _NOT USED_                          |
| OUT_GENO                 | yes                                                           | _NOT USED_                                  |
| OUT_CANAL_TESTS          | 0                                                             | _NOT USED_                                         |
| OUT_CANAL_MUTSD          | 0.5                                                           | _NOT USED_                   |
| OUT_HERIT_TESTS          | 0                                                             | _NOT USED_                                         |
| OUT_DIREPI_TESTS         | 0                                                             | _NOT USED_                                     |
| INIT_RECURRENCE          | 0                                                             | _NOT USED_                                         |
| FITNESS_FLUCT            | no_fluctuation                                                 | _NOT USED_                |
| SIMUL_MAXGEN             | 5000                                                          | _NOT USED_               |
| TYPE_ARCHI               | m2                                                            | _NOT USED_                     |

NB: many parameters are not used because the parameter file originates from another version of the model that used more parameters.
