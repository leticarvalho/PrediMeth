# PrediMeth
Git repository with scripts used to develop this project.

# Contents

## 1. Data preprocessing (Test cohort): *Stage1_DataPreprocessing.R*

### Description
In this script, the data from the discovery cohort is processed and prepared for the future EWAS. It involves phenotype and methylation data.  

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
In this script, the processed data is used to fit a linear model to each CpG, using limma R package. 
**Model**: Mvalues ~ MedDiet adherence + age + sex + BMI smoking + batch + cells
Distinct exposure variables to represent MedDiet adherence: 
* Continuous = predimed_score (0-14)
* Binary = predimed_high (high >= 9) - 553 (no) vs 421 (yes) 
* 3-factor categorical = predimed_cat (low <=4, high >=9): 38 (low) vs 515 (medium) vs 421 (high)
* Extreme = predimed_extreme (Low_Q1 <= 6, Mid_Q23 = 7-9, High_Q4 >= 10): 235 (Low_Q1) vs 262 (High_Q4)
For each model, a topTable will be generated, including EPIC v2 annotation and calculated absolute delta beta (ADB). Top hits will be FDR-adjusted.  

### Input 
* Data frame of phenotype (rows: samples, columns: variables)
* Matrix of beta values (rows: probe names, columns: samples)
* Matrix of M-values (rows: probe names, columns: samples)

### Output
* Lambdas comparison of all contrasts
For each contrast, a folder with: 
* Plot with beta values distribution (bimodal) 
* Results tables (topTable + annotation + ADB) (.Rdata)
* CpGs vectors to enrichment (.Rdata)
* Tables with fdr hits (.csv)
* Volcano plot
* Manhattan plot
* Genomic Categories barplot

## 3. DMR Analysis: *Stage1_AnalysisDMR.R*

### Description
In this script, EWAS results will be investigated for existence of Different Methylated Regions (DMRs), using DMRcate R package.

### Input 
Processed data (from discovery cohort): 
* phenotype data frame,
* beta values matrix
* m-values matrix

### Output 
For each contrast that has DMPs (DMRcate):
* DMPs 
* DMRs 

## 4. Data Preprocessing (Validation cohort): *Stage2_ValidationDataPreprocessing.R*

### Description
In this script, phenotype and methylation data from the validation cohort are preprocessed. 

### Input 
* Phenotype data (diet adherence + covariables)
* Blood cells counts

### Output
* Data frame of phenotype (rows: samples, columns: variables)
* Matrix of beta values (rows: probe names, columns: samples)
* Matrix of M-values (rows: probe names, columns: samples)

## 5. Construction and calculation of Methylation Risk Scores: *Stage2_MRSconstruction.R*

### Description
In this script, Methylation Risk Score (MRS) is constructed with different CpG selection strategies, and calculated for each sample from the validation cohort. 
CpG selection strategies:
1. MRS constructed with FDR hits from standard contrasts (HL + ML) using HL delta beta as weight (from discovery EWAS)
2. MRS constructed with top CpGs from continuous contrast, with pvalue < 1e-5, using Continuous delta beta per unit as weight (from discovery EWAS)
3. MRS constructed with top CpGs from continuous contrast, with pvalue < 1e-4, using Continuous delta beta per unit as weight (from discovery EWAS)
4. MRS constructed with an elastic net, with all CpGs, with alfa = 0.5, lambda min, 10-fold cv, residualizing for covariates  

### Input 
From discovery cohort: 
* Pheno table
* Betavalues
* Limma's topTable (continuous model)
* Limma's topTable (High vs low - 3 categories model)
From validation cohort:
* Pheno data
* Betavalues

### Output
* MRSs calculated for validation cohort (csv)
* MRSs names and weights (.R data)

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
