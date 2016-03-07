### Code utilized for generating table 1
This code is provided to reproduce the analysis for table 1 for BIOSTS-15205. For our differential expression analysis, we utilize the metadata and expression data from GSE61901

### Required R packages
 - utils (we utilized R 3.2.3's built-in version)
 - lme4 (we utilized version 1.1-10)
 - limma (we utilized version 3.26.5)
 - sva (we utilized version 3.18.0)
 
### Usage
Run the script using
```
$ Rscript table1_code.R
```
or simply `source('table1_code.R')` within R to download the data from GEO and run and output the analysis

### Additional information
The following methods utilize quantile Normalization, no averaging of technical replicates, and utilize CHIP as blocking variable for technical replicates

 * The `ranova` function utilizes a repeated-measures anova to take technical replicates into account as part of the differential expression analysis
 * The `lmem` function utilizes a linear mixed-effects model to take technical replicates into account as part of the differential expression analysis
 * The `limma_dupcorr` function utilizes LIMMA's duplicate correlation function to take technical replicates into account as part of the differential expression analysis


The following methods utilize quantile normalization, average technical replicates and apply ComBat including treatments as covariates

 * The `limma_noblocking` function utilizes LIMMA for differential expression analysis
 * The `twowayanova` function utilizes an ANOVA for differential expression analysis

The following methods utilize quantile normalization, average technical replicates and CHIP as blocking variable

 * The `limma_blocking` function utilizes LIMMA for differential expression analysis. We note that we are not utilizing this function and instead are citing Nygaard et al's result directly for Table 1. Nygaard et al's analysis utilizes the labels from GSE40566 and the expression data from GSE61901.
 * The `twowayanova_blocking` function utilizes ANOVA for differential expression analysis adjusting for CHIP