# PrediMeth
Git repository with scripts used to develop this project.

# Contents

## 1. Data preprocessing (Test cohort): *Stage1_DataPreprocessing.R*

### Description
In this script, the data from the discovery cohort is processed and prepared for the future EWAS. It involves phenotype and methylation data. The data from discovery cohort comes from IJC. 

### Input
* Phenotype data (diet adherence + covariables)
* Methylation data (beta values matrix already normalized)
* Annotation from EPIC version 2
* Blood cells counts

### Output 
* Data frame of phenotype (rows: samples, columns: variables)
* Matrix of beta values (rows: probe IDs, columns: samples)
* Matrix of M-values (rows: probe IDs, columns: samples)

## 2. EWAS of Mediterranean diet: *Stage1_FittingModel.R*

### Description
In this script, the processed data is used to fit the best linear model (selected after sensitivity analyses) to each CpG, using limma R package. For this model, a topTable will be generated, including EPIC v2 annotation and calculated absolute delta beta (ADB). 
**Selected model**: Mvalues ~ diet + age + sex + smoking + BMI + cell counts
*Obs. diet may be binary (predimed_high) or multifactor (predimed_cat); some lines are commented depending on this.*

### Input 
* Data frame of phenotype (rows: samples, columns: variables)
* Matrix of beta values (rows: probe IDs, columns: samples)
* Matrix of M-values (rows: probe IDs, columns: samples)

### Output
* Results table (topTable + annotation + ADB)
* Bonferroni's threshold

## 3. Results visualization: *Stage1_VisualizingEWASResults.R*

### Description
In this script, EWAS results are visualized by plots, regarding different aspects. R packages like *ggplot2*, *qqman* and *EnhancedVolcano* will be used. 

### Input 
* EWAS results (topTable)

### Output 
* FDR-significant topTable subset
* CpGs in vectors (to use in enrichment)
* Volcano Plot
* Manhattan Plot
* QQplot (genomic inflation analysis)
* Genomic Categories Bar Plot

## 4. Annotation of top hits: *Stage1_AnnotatingEWASTopHits.R*

### Description
In this script, EWAS results will be annotated using missMethyl functions, to account for multi-probe bias, using some annotation databases: Gene Ontology (GO), Kyoto Encyclopedia of Genes and Genomes (KEGG), Reactome, and alternative gene sets. If not FDR-significant, results will be ranked by raw p-value, and showed by a dotplot or a barchart. 

### Input 
* EWAS results (CpG vectors)

### Output
* Gene Ontology (GO)
* Kyoto Encyclopedia of Genes and Genomes (KEGG)
* Reactome
* Alternative gene sets: MSigDB Hallmark, ImmuneSigDB, Wikipathways

## 5. Construction and calculation of Methylation Risk Scores: *Stage2_MRSconstruction.R*

### Description
In this script, CpGs are selected to construct different Methylation Risk Scores (MRS), which will be later calculated for each sample from tha validation cohort.  

### Input 
From discovery cohort: 
* Pheno data
* Mvalues
* Limma's topTables
From validation cohort: 
* Pheno data
* Betavalues

### Output
* MRSs calculated for validation cohort

## 6. Study of association between MRS and T2D: *Stage2_MRSandT2Dassociation.R*

### Description
In this script, we investigate if MRSs are associated with T2D in two different contexts. For prevalent T2D, we use a logistic regression and odds ratio. For incident T2D, we use Cox proportional hazards model and hazard ratio.   

### Input 
* Validation cohort (already preprocessed, with MRSs already calculated, as variables)

### Output
* Forest plot summarizing MRS and T2D associations
* Summary table (both analyses)
* Summary table of logistic regression
* Summary table of Cox proportional hazards
