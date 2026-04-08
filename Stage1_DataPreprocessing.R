##########################################################
#
# Data Preprocessing - Stage 1 (Test Cohort)
#
##########################################################


# DESCRIPTION --------------------------------------------
# In this script, the data from the discovery cohort is 
# processed and prepared for the future EWAS. It involves
# phenotype and methylation data. The data from discovery 
# cohort comes from IJC. 

# INPUT --------------------------------------------------
# 	Phenotype data (diet adherence + covariables)
# 	Methylation data (beta values matrix already normalized)
# 	Annotation from EPIC version 2
# 	Blood cells counts

# OUTPUT -------------------------------------------------
# 	Data frame of phenotype (rows: samples, columns: variables)
# 	Matrix of beta values (rows: probe IDs, columns: samples)
#	  Matrix of M-values (rows: probe IDs, columns: samples)


##########################################################


#### LOAD LIBRARIES ####

library(dplyr)
library(data.table)
library(lumi)
library(ExperimentHub)
library(DMRcate)

#### RESOLVE PATHS ####

cat("Starting script Stage1_DataPreprocessing \n")
message("Starting script Stage1_DataPreprocessing : ", Sys.time() )

## Input paths -------------------------------------------

# Data 
methylation_path <- "/path_to_beta_values.R" # R object with beta values (already normalized)
phenotype_path <- "/path_to_phenotype_information.csv" # .csv file with phenotype information
annotation_path <- "/path_to_EPICv2_annotation.R" # R object with EPICv2 annotation
cellcounts_path <- "/path_to_cell_counts.csv" # .csv file with cell counts

# Project (PrediMeth)
predimeth_path <- "/path_to_project_folder/"
results_dir <- file.path(predimeth_path, "results")

## Output paths ------------------------------------------

results_folder <- file.path(results_dir, 'ProcessedData_Stage1')
dir.create(results_folder)

##########################################################

#### LOAD DATA ####

## Phenotype data ----------------------------------------

pheno <- fread(phenotype_path) # rows = sample ids, columns = phenotype variables

## Methylation data ---------------------------------------

# Beta values (already normalized)
Betavalues <- load(methylation_path) # rows = CpG probes names, columns = sample ids

# Annotation (EPIC v2)
load(annotation_path) # rows = CpG probes ids, columns = features 
EPICv2_annotation$Probe_name  <- stringr::str_remove(EPICv2_annotation$ProbeID, "_.*")
probes_Y <- subset(EPICv2_annotation, CpGchrm=="chrY") # probes from chrY
probes_X <- subset(EPICv2_annotation, CpGchrm=="chrX") # probes from chrX

# Blood cell counts (to after adjustment as covariates)
cell_counts <- data.table::fread(cellcounts_path) # rows = sample ids, columns = cell types

##########################################################

cat("Data loaded. \n\nStarting preprocessing of pheno data. \n")
message("Data loaded. \nStarting preprocessing of pheno data: ", Sys.time())

#### PREPROCESSING PHENOTYPE DATA ####

## Select variables of interest
pheno <- subset(pheno, select = c("sample_id", # identification
                                  "age", "sex", "smoking_habit", "bmi", # covariates
                                  "predimed_score", # PREDIMED score (adherence to MedDiet)
                                  "predimed_cat", # Adherence categories: 1 = low (0-4), 2 = medium (5-8), 3 = high (9-14)
                                  "predimed_high")) # High adherence (binary): 0 = no, 1 = yes
cell_counts <- cell_counts[, c("Bcell", "CD4T", "CD8T", "Mono", "Neu", "NK")]

## Adding info about cell counts
pheno <- merge(pheno, cell_counts, by = "sample_id")

## Setting rownames
rownames(pheno) <- pheno$sample_id

## Looking at dimensions
cat("	Dimensions after preprocessing of pheno data: ", 
    "\n		Bvalues = ", dim(Betavalues), 
    "\n		Pheno = ", dim(pheno), "\n")

##########################################################

cat("Preprocessing of pheno data completed. \n\nStarting preprocessing of methylation data: \n")
message("Preprocessing of pheno data completed. \nStarting preprocessing of methylation data: ", Sys.time())

#### PREPROCESSING METHYLATION DATA ####

## Filtering beta values to only baseline samples (there were some 2023 samples that have been removed from pheno data)
Betavalues <- subset(Betavalues, select = c(colnames(Betavalues) %in% rownames(pheno)))

## Removing missing values from meth set
Betavalues <- na.omit(Betavalues)
cat("	Bvalues dimensions without missing values: ", dim(Betavalues), "\n")

## Clamping beta-values = 0 (so they don't convert into infinite M-values)
cat("	Amount of Bvalues = 0: ", sum(Betavalues == 0), "\n")
cat("	Amount of Bvalues </= : 1e-6", sum(Betavalues <= 1e-6), "\n")
cat("	Amount of Bvalues = 1: ", sum(Betavalues == 1), "\n")
cat("	Amount of Bvalues >/= : 1 - 1e-6", sum(Betavalues >= (1 - 1e-6)), "\n")
Betavalues[Betavalues == 0] <- 1e-6
Betavalues[Betavalues == 1] <- 1 - 1e-6

## Changing names of probes (rownames), from probe names to probe IDs (EPIC v2)
cat("	Changing names of probes to ProbeIDs. \n")
id_map <- setNames(EPICv2_annotation$ProbeID, EPICv2_annotation$Probe_name)
new_names <- id_map[rownames(Betavalues)]
rownames(Betavalues) <- ifelse(is.na(new_names),
                                 rownames(Betavalues),
                                 new_names)

## Filtering probes with DMRcate - removes SNP-related probes and cross-hybridising probes
cat("Filtering probes with DMRcate. \n")
message("Filtering probes with DMRcate : ", Sys.time())
hub <- ExperimentHub()
setExperimentHubOption("CACHE", "/path")
Betavalues <- rmSNPandCH(Betavalues, 
							dist = 2, # maximum distance (from CpG to SNP/variant) of probes to be filtered out
							mafcut = 0.05, # minimum minor allele frequency of probes to be filtered out
							and = TRUE, 
							rmcrosshyb = TRUE, # filter cross-hybridized probes
							rmXY=FALSE) # filter sex probes
cat("Bvalues dimensions after SNP and CH related probes filtering: ", dim(Betavalues), "\n")

## Removing chrX or chrY associated probes from dataset (according to EPICv2 annotation)
Betavalues <- Betavalues[!(rownames(Betavalues) %in% probes_Y$ProbeID),] 
cat("Bvalues dimensions after chrY related probes filtering: ", dim(Betavalues), "\n")
Betavalues <- Betavalues[!(rownames(Betavalues) %in% probes_X$ProbeID),] 
cat("Bvalues dimensions after chrX related probes filtering: ", dim(Betavalues), "\n")

## Convert beta values to M-values
Mvalues <- minfi::logit2(Betavalues)
cat("Beta values converted to M-values. Dimensions: ", 
    "\nBvalues = ", dim(Betavalues), 
    "\nMvalues = ", dim(Mvalues), "\n")
    
## Analyse infinitive values (produced by Beta-Values = 0 or 1) 
inf_values <- !is.finite(Mvalues)
# per CpG
inf_per_cpg <- rowSums(!is.finite(Mvalues))
table(inf_per_cpg)
# per sample
inf_per_sample <- colSums(!is.finite(Mvalues))
table(inf_per_sample)
cat("Infinite Mvalues: ", sum(inf_values), 
    "\nPer CpG: ", table(inf_per_cpg), 
     "\nPer sample: ", table(inf_per_sample)) 
     
## Checking rownames 
cat("Checking methylation matrixes rownames: \n")
Betavalues[1:5, 1:5]
Mvalues[1:5, 1:5]

##########################################################

cat("Preprocessing of meth data completed. \nStarting update of pheno data: \n")
message("Preprocessing of meth data completed. Updating pheno data: ", Sys.time())

#### UPDATING PHENOTYPE DATA ####

## PrediMeth_all_merge has some samples that have been removed from beta values due to quality check
keep_samples <- pheno$sample_id %in% colnames(Betavalues)
pheno <- pheno[keep_samples,]

## Remove missing values from patients data (covariates)
pheno <- pheno[!pheno$predimed_high == "NULL", ] 
pheno <- pheno[!pheno$predimed_cat == "NULL", ] 
pheno <- pheno[!pheno$age == "NULL", ] 
pheno <- pheno[!pheno$sex == "NULL", ] 
pheno <- pheno[!pheno$bmi == "NULL", ] 
pheno <- pheno[!pheno$smoking_habit == "NULL", ] 
cat("Pheno data dimensions after removing missing values from covariables: ", dim(pheno), "\n")

## Setting rownames
rownames(pheno) <- pheno$sample_id

## Structure of variables
str(pheno)

##########################################################

cat("Update of pheno data completed. \nStarting update of meth data: \n")
message("Pheno data updated. Updating meth data: ", Sys.time())

#### UPDATING METHYLATION DATA ####

## M values
Mvalues <- subset(Mvalues, select = c(colnames(Mvalues) %in% pheno$sample_id)) 

## Betas values
Betavalues <- subset(Betavalues, select = c(colnames(Betavalues) %in% pheno$sample_id)) 

## Looking at dimensions
print("Dimensions after Updates of pheno and methylation data:")
cat("Dimensions after updates of pheno and methylation data: ", 
    "\nPhenoData (samples x variables) = ", dim(pheno),
    "\nBvalues (probes x samples) = ", dim(Betavalues), 
    "\nMvalues (probes x samples) = ", dim(Mvalues), "\n")

##########################################################

cat("Updates completed. \nSaving outputs: \n")
message("Meth data updated. Saving outputs: ", Sys.time())

#### SAVING OUTPUTS ####

save(pheno, Mvalues, Betavalues, file = file.path(results_folder, "processed_data.R"))

cat("Script completed. \n")
message("Script completed: ", Sys.time())

setExperimentHubOption("CACHE", NULL)  


