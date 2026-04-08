# PrediMeth
Git repository with scripts used to develop this project.

# Contents

## 1. Data preprocessing (Test cohort): *Stage1_DataPreprocessing.R*

### DESCRIPTION
In this script, the data from the discovery cohort is processed and prepared for the future EWAS. It involves phenotype and methylation data. The data from discovery cohort comes from IJC. 

### INPUT
* Phenotype data (diet adherence + covariables)
* Methylation data (beta values matrix already normalized)
* Annotation from EPIC version 2
* Blood cells counts

### OUTPUT 
* Data frame of phenotype (rows: samples, columns: variables)
* Matrix of beta values (rows: probe IDs, columns: samples)
* Matrix of M-values (rows: probe IDs, columns: samples)

2. EWAS of Mediterranean diet: *Stage1_FittingModel.R*

3. Results visualization: *Stage1_VisualizingEWASResults.R*

4. Annotation of top hits: *Stage1_AnnotatingEWASTopHits.R*
