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
# 	Matrix of beta values (rows: probe IDs/names, columns: samples)
#	  Matrix of M-values (rows: probe IDs/names, columns: samples)


##########################################################


#### LOAD LIBRARIES ####

library(dplyr)
library(data.table)
library(lumi)
library(ExperimentHub)
library(DMRcate)

#### RESOLVE PATHS ####

## Define label ------------------------------------------

#label <- "_rmCH_XY_SNP" 
#~ label <- "_rmXY" 

label <- "_rmCH_XY_SNP_clamp"  # probe IDs
#~ label <- "_rmXY_clamp" # probe names 

cat("Starting script Stage1_DataPreprocessing with label ", label, " \n")
message("Starting script Stage1_DataPreprocessing with label ", label, " : ", Sys.time() )

## Input paths -------------------------------------------

# Data 
data_path <- "/imppc/labs/dnalab/studentdnabank/Data/"
methylation_dir <- file.path(data_path, 'omics', 'methylation', 'ijc')
questionary_dir <- file.path(data_path, 'questionnaire')

# Project (PrediMeth)
predimeth_path <- "/imppc/labs/dnalab/share/PrediMeth"
results_dir <- file.path(predimeth_path, "results")

## Output paths ------------------------------------------

results_folder <- file.path(results_dir, paste0('ProcessedData_Stage1', label))
dir.create(results_folder)

##########################################################

#### LOAD DATA ####

## Phenotype data ----------------------------------------

# Record linkage tables 
consulta <- fread(file.path(questionary_dir, "consulta_nivell2.csv"))
IJC_metadata <- fread(file.path(methylation_dir, "IJC_metadata_intern.csv"))

# Phenotype information
questionary_1 <- fread(file.path(questionary_dir, "questionari_nivell1.csv"))
measures_1 <- fread(file.path(questionary_dir, "mesures_nivell1.csv"))

## Methylation data ---------------------------------------

# Beta values (already normalized)
load("/imppc/labs/dnalab/studentdnabank/Data/omics/methylation/ijc/4.norm_pc_10_IJC_beta.Robj") 

# Annotation (EPIC v2)
load("/imppc/labs/dnalab/share/PrediMeth/0-data/EPICv2.annot.RData")
EPIC_MEsteller$Probe_name  <- stringr::str_remove(EPIC_MEsteller$ProbeID, "_.*")
probes_Y <- subset(EPIC_MEsteller, CpGchrm=="chrY") # probes from chrY
probes_X <- subset(EPIC_MEsteller, CpGchrm=="chrX") # probes from chrX

# Blood cell counts (to after adjustment as covariates)
cell_counts <- data.table::fread("/imppc/labs/dnalab/studentdnabank/Data/omics/methylation/ijc/cell_and_clocks/cell_counts.csv")
cell_counts$sample_id <- cell_counts$V1
cell_counts$V1 <- NULL

##########################################################

cat("Data loaded. \n\nStarting preprocessing of pheno data. \n")
message("Data loaded. \nStarting preprocessing of pheno data: ", Sys.time())

#### PREPROCESSING PHENOTYPE DATA ####

## Get entity_id for discovery samples
consulta$entity_id <- stringr::str_remove(consulta$entity_id, "=*")
consulta_merge <- merge(IJC_metadata, consulta, by = "entity_id")

## Filter baseline samples (exclude paired samples collected in 2023)
consulta_merge <- consulta_merge[consulta_merge$sample_year=="baseline",] 

## Select variables from other tables
consulta_merge <- subset(consulta_merge, 
                         select = c("entity_id", "SAMPLE_ID (xIJC)", "COUNTRY_BIRTH.y", 
                                    "METABOLOMED", "METABOLOMED_SAMPLE", "IA4T2D", "COMORBIDITY")) 
consulta_merge$sample_id <- consulta_merge$`SAMPLE_ID (xIJC)`
consulta_merge$`SAMPLE_ID (xIJC)` <- NULL

questionary_1_subset <- subset(questionary_1, 
                               select = c("entity_id", "EDAD_ANOS", "SEXO", 
                                          "smoking_habit", "predimed_score", "predimed_cat", "predimed_high",
                                          "sedentarisme", "mets_semana", "RURAL", 
                                          "REGION_RESIDENCE", "REGION_BIRTH", 
                                          "MUNICIPIO_RESIDENCIA", "MUNICIPIO_NACIMIENTO",
                                          "LABORAL_ESTADO", "SANIDAD", "INGRESOS", 
                                          "ESTUDIOS", "ESTADO_CIVIL", "ETNIA_PARTICIPANTE"))

measures_subset <- subset(measures_1, select = c("entity_id", "BMI", "CALC_AVG_PESO", 
                                               "CALC_AVG_HEIGHT", "CALC_AVG_CINTURA", 
                                               "CALC_AVG_CADERA", "WHR"))

## Information about patients (1 row per patient) - merged data
PrediMeth_all_merge <- merge(questionary_1_subset, measures_subset, by = "entity_id") 
PrediMeth_all_merge$entity_id <- stringr::str_remove(PrediMeth_all_merge$entity_id, "=*")
PrediMeth_all_merge <- merge(PrediMeth_all_merge, consulta_merge)

## Adding info about cell counts
PrediMeth_all_merge <- merge(PrediMeth_all_merge, cell_counts, by = "sample_id")

## Setting rownames
rownames(PrediMeth_all_merge) <- PrediMeth_all_merge$sample_id

## Looking at dimensions
cat("	Dimensions after preprocessing of pheno data: ", 
    "\n		Bvalues = ", dim(IJC_norm.beta), 
    "\n		Pheno = ", dim(PrediMeth_all_merge), "\n")

##########################################################

cat("Preprocessing of pheno data completed. \n\nStarting preprocessing of methylation data: \n")
message("Preprocessing of pheno data completed. \nStarting preprocessing of methylation data: ", Sys.time())

#### PREPROCESSING METHYLATION DATA ####

## Filtering beta values to only baseline samples (there were some 2023 samples that have been removed from pheno data)
IJC_norm.beta <- subset(IJC_norm.beta, select = c(colnames(IJC_norm.beta) %in% rownames(PrediMeth_all_merge)))

## Removing missing values from meth set
IJC_norm.beta <- na.omit(IJC_norm.beta)
cat("	Bvalues dimensions without missing values: ", dim(IJC_norm.beta), "\n")

## Clamping beta-values = 0 (so they don't convert into infinite M-values)
cat("	Amount of Bvalues = 0: ", sum(IJC_norm.beta == 0), "\n")
cat("	Amount of Bvalues </= : 1e-6", sum(IJC_norm.beta <= 1e-6), "\n")
cat("	Amount of Bvalues = 1: ", sum(IJC_norm.beta == 1), "\n")
cat("	Amount of Bvalues >/= : 1 - 1e-6", sum(IJC_norm.beta >= (1 - 1e-6)), "\n")
IJC_norm.beta[IJC_norm.beta == 0] <- 1e-6
IJC_norm.beta[IJC_norm.beta == 1] <- 1 - 1e-6

## Changing names of probes (rownames), from probe names to probe IDs (EPIC v2)
#~ cat("	Keeping names of probes as Probe_name. \n")
cat("	Changing names of probes to ProbeIDs. \n")
id_map <- setNames(EPIC_MEsteller$ProbeID, EPIC_MEsteller$Probe_name)
new_names <- id_map[rownames(IJC_norm.beta)]
rownames(IJC_norm.beta) <- ifelse(is.na(new_names),
                                 rownames(IJC_norm.beta),
                                 new_names)

#~ ## Filtering probes with DMRcate - removes SNP-related probes and cross-hybridising probes
cat("Filtering probes with DMRcate. \n")
message("Filtering probes with DMRcate : ", Sys.time())
hub <- ExperimentHub()
setExperimentHubOption("CACHE", "/home/labs/dnalab/studentdnabank/ExperimentHub")
IJC_norm.beta <- rmSNPandCH(IJC_norm.beta, 
							dist = 2, # maximum distance (from CpG to SNP/variant) of probes to be filtered out
							mafcut = 0.05, # minimum minor allele frequency of probes to be filtered out
							and = TRUE, 
							rmcrosshyb = TRUE, # filter cross-hybridized probes
							rmXY=FALSE) # filter sex probes
cat("Bvalues dimensions after SNP and CH related probes filtering: ", dim(IJC_norm.beta), "\n")

## Removing chrX or chrY associated probes from dataset (according to EPICv2 annotation)  ----------> CHANGE WHETHER IS PROBE NAME OR ID
IJC_norm.beta <- IJC_norm.beta[!(rownames(IJC_norm.beta) %in% probes_Y$ProbeID),] 
cat("Bvalues dimensions after chrY related probes filtering: ", dim(IJC_norm.beta), "\n")
IJC_norm.beta <- IJC_norm.beta[!(rownames(IJC_norm.beta) %in% probes_X$ProbeID),] 
cat("Bvalues dimensions after chrX related probes filtering: ", dim(IJC_norm.beta), "\n")
#~ IJC_norm.beta <- IJC_norm.beta[!(rownames(IJC_norm.beta) %in% probes_Y$Probe_name),] 
#~ cat("Bvalues dimensions after chrY related probes filtering: ", dim(IJC_norm.beta), "\n")
#~ IJC_norm.beta <- IJC_norm.beta[!(rownames(IJC_norm.beta) %in% probes_X$Probe_name),] 
#~ cat("Bvalues dimensions after chrX related probes filtering: ", dim(IJC_norm.beta), "\n")

## Convert beta values to M-values
IJC_m_values <- minfi::logit2(IJC_norm.beta)
cat("Beta values converted to M-values. Dimensions: ", 
    "\nBvalues = ", dim(IJC_norm.beta), 
    "\nMvalues = ", dim(IJC_m_values), "\n")
    
## Analyse infinitive values (produced by Beta-Values = 0 or 1) 
inf_values <- !is.finite(IJC_m_values)
# per CpG
inf_per_cpg <- rowSums(!is.finite(IJC_m_values))
table(inf_per_cpg)
# per sample
inf_per_sample <- colSums(!is.finite(IJC_m_values))
table(inf_per_sample)
cat("Infinite Mvalues: ", sum(inf_values), "\nPer CpG: ", table(inf_per_cpg), 
     "\nPer sample: ", table(inf_per_sample)) # 
     
## Checking rownames 
cat("Checking methylation matrixes rownames: \n")
IJC_norm.beta[1:5, 1:5]
IJC_m_values[1:5, 1:5]

##########################################################

cat("Preprocessing of meth data completed. \nStarting update of pheno data: \n")
message("Preprocessing of meth data completed. Updating pheno data: ", Sys.time())

#### UPDATING PHENOTYPE DATA ####

## PrediMeth_all_merge has some samples that have been removed from beta values 
## due to quality check - remove them
keep_samples <- PrediMeth_all_merge$sample_id %in% colnames(IJC_norm.beta)
PrediMeth_all_merge <- PrediMeth_all_merge[keep_samples,]

## Remove missing values from patients data (covariables)
PrediMeth_all_merge <- PrediMeth_all_merge[!PrediMeth_all_merge$predimed_high == "NULL", ] # predimeth_high
PrediMeth_all_merge <- PrediMeth_all_merge[!PrediMeth_all_merge$predimed_cat == "NULL", ] # predimeth_cat
PrediMeth_all_merge <- PrediMeth_all_merge[!PrediMeth_all_merge$EDAD_ANOS == "NULL", ] # age
PrediMeth_all_merge <- PrediMeth_all_merge[!PrediMeth_all_merge$SEXO == "NULL", ] # sex
PrediMeth_all_merge <- PrediMeth_all_merge[!PrediMeth_all_merge$BMI == "NULL", ] # BMI
PrediMeth_all_merge <- PrediMeth_all_merge[!PrediMeth_all_merge$smoking_habit == "NULL", ] # smoking
cat("Pheno data dimensions after removing missing values from covariables: ", dim(PrediMeth_all_merge), "\n")

## Setting rownames
rownames(PrediMeth_all_merge) <- PrediMeth_all_merge$sample_id

## Structure of variables
str(PrediMeth_all_merge)

##########################################################

cat("Update of pheno data completed. \nStarting update of meth data: \n")
message("Pheno data updated. Updating meth data: ", Sys.time())

#### UPDATING METHYLATION DATA ####

## M values
IJC_m_values <- subset(IJC_m_values, select = c(colnames(IJC_m_values) %in% PrediMeth_all_merge$sample_id)) 

## Betas values
IJC_norm.beta <- subset(IJC_norm.beta, select = c(colnames(IJC_norm.beta) %in% PrediMeth_all_merge$sample_id)) 

## Looking at dimensions
print("Dimensions after Updates of pheno and methylation data:")
cat("Dimensions after updates of pheno and methylation data: ", 
    "\nPhenoData (samples x variables) = ", dim(PrediMeth_all_merge),
    "\nBvalues (probes x samples) = ", dim(IJC_norm.beta), 
    "\nMvalues (probes x samples) = ", dim(IJC_m_values), "\n")

##########################################################

cat("Updates completed. \nSaving outputs: \n")
message("Meth data updated. Saving outputs: ", Sys.time())

#### SAVING OUTPUTS ####

save(PrediMeth_all_merge, IJC_m_values, IJC_norm.beta, file = file.path(results_folder, paste0("processed_data", label, ".R")))
#write.table(IJC_m_values, file = file.path(results_folder, paste0("IJC_m_values", label, ".txt")), sep = "\t", col.names = TRUE, row.names = TRUE)
#write.table(IJC_norm.beta, file = file.path(results_folder, paste0("IJC_norm.beta", label, ".txt")), sep = "\t", col.names = TRUE, row.names = TRUE)
#write.table(PrediMeth_all_merge, file = file.path(results_folder, paste0("PrediMeth_all_merge", label, ".txt")), col.names = TRUE, row.names = FALSE)

cat("Script completed. \n")
message("Script completed: ", Sys.time())

setExperimentHubOption("CACHE", NULL)  


