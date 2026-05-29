##########################################################
#
# Data Preprocessing - Stage 1 (Test Cohort)
#
##########################################################


# DESCRIPTION --------------------------------------------
# In this script, the data from the discovery cohort is 
# processed and prepared for the future EWAS. It involves
# phenotype and methylation data. 

# INPUT --------------------------------------------------
# 	Phenotype data (diet adherence + covariables)
# 	Methylation data (beta values matrix already normalized)
# 	Annotation from EPIC version 2
# 	Blood cells counts

# OUTPUT -------------------------------------------------
# 	Data frame of phenotype (rows: samples, columns: variables)
# 	Matrix of beta values (rows: probe names, columns: samples)
#	  Matrix of M-values (rows: probe names, columns: samples)

##########################################################

cat("\n\n\nStarting script - Data Preprocessing - Stage 1. \n")
message("Starting script: ", Sys.time() )

#### LOAD LIBRARIES ####

library(dplyr)
library(data.table)
library(lumi)

#### RESOLVE PATHS ####

## Input paths -------------------------------------------

methylation_dir <- ""
questionary_dir <- ""
path_to_normalized_betas <- "" # .Robj
path_to_EPIC2_annotation <- "" # .Rdata
path_to_disc_cell_counts <- "" # csv
path_to_CCI <- "" # csv of Charlson index

## Output paths ------------------------------------------

predimeth_path <- "" # project folder
results_dir <- file.path(predimeth_path, "results")
results_folder <- file.path(results_dir, 'OfficialAnalysis') 
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

# Annotation (EPIC v2)
load(path_to_EPIC2_annotation) # EPIC_MEsteller
EPIC_MEsteller$Probe_name  <- stringr::str_remove(EPIC_MEsteller$ProbeID, "_.*")
cat("Total annotation rows: ", nrow(EPIC_MEsteller), "\n")
cat("Unique ProbeIDs: ", dplyr::n_distinct(EPIC_MEsteller$ProbeID), "\n")                                                                                                                        
cat("Unique Probe_names: ", dplyr::n_distinct(EPIC_MEsteller$Probe_name), "\n")                                                                                                                     
dup_table <- table(EPIC_MEsteller$Probe_name) # how many Probe_names have > 1 ProbeID                                                                                                                                                             
cat("Probe_names with >1 ProbeID: ", sum(dup_table > 1), "\n") # should be zero                                                                                                                                              
cat("Max ProbeIDs per Probe_name: ", max(dup_table), "\n")
probes_Y <- subset(EPIC_MEsteller, CpGchrm=="chrY") # probes from chrY
probes_X <- subset(EPIC_MEsteller, CpGchrm=="chrX") # probes from chrX

# Beta values (already normalized)
load(path_to_normalized_betas) # IJC_norm.beta
in_anno <- rownames(IJC_norm.beta) %in% EPIC_MEsteller$Probe_name
cat("\nBeta-matrix rows in annotation: ", sum(in_anno), "/", nrow(IJC_norm.beta), "\n")    
unmatched <- setdiff(rownames(IJC_norm.beta), EPIC_MEsteller$Probe_name)  
cat("There are ", length(unmatched), " probes that are not in annotation and will be removed. \n")  # 825                                                                                                                                    
head(unmatched, 20)                                                                                                                                                                                           
table(substr(unmatched, 1, 3))  # rs / ch. / cg / nv- / ctl
# a chunk of rs (59 SNP probes for sample fingerprinting), some ch. (non-CpG), maybe a few legacy cg probes that were dropped from the v2 manifest.                                         
IJC_norm.beta <- IJC_norm.beta[rownames(IJC_norm.beta) %in% EPIC_MEsteller$Probe_name, ]
cat("Unmatched dropped:", length(unmatched), "| Remaining:", nrow(IJC_norm.beta), "\n") 

# Blood cell counts (to after adjustment as covariates)
cell_counts <- data.table::fread(path_to_disc_cell_counts)
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
consulta_merge$batch <- consulta_merge$`Plate name FINAL` # sample plate name (to after adjustment for batch effect)
consulta_merge <- subset(consulta_merge, 
                         select = c("entity_id", "SAMPLE_ID (xIJC)", "batch")) 
consulta_merge$sample_id <- consulta_merge$`SAMPLE_ID (xIJC)`
consulta_merge$`SAMPLE_ID (xIJC)` <- NULL

questionary_1_subset <- subset(questionary_1, 
                               select = c("entity_id", "EDAD_ANOS", "SEXO", 
                                          "smoking_habit", "predimed_score", "predimed_cat", "predimed_high",
                                          "mets_semana", "ETNIA_PARTICIPANTE",
                                          "PREDIMED_ACEITE_GRASA", "PREDIMED_ACEITE_CONSUMO", 
                                          "PREDIMED_VERDURAS", "PREDIMED_FRUTA", 
                                          "PREDIMED_CARNES_ROJAS", "PREDIMED_MANTEQUILLA", 
                                          "PREDIMED_BEBIDAS_CARBONATADAS", "PREDIMED_VINO", 
                                          "PREDIMED_LEGUMBRES", "PREDIMED_PESCADO", 
                                          "PREDIMED_REPOSTERIA", "PREDIMED_FRUTOS_SECOS", 
                                          "PREDIMED_CARNE", "PREDIMED_VEGETALES"))

measures_subset <- subset(measures_1, select = c("entity_id", "BMI"))

## Information about patients (1 row per patient) - merged data
PrediMeth_all_merge <- merge(questionary_1_subset, measures_subset, by = "entity_id") 
PrediMeth_all_merge$entity_id <- stringr::str_remove(PrediMeth_all_merge$entity_id, "=*")
PrediMeth_all_merge <- merge(PrediMeth_all_merge, consulta_merge, by = "entity_id")

## Adding info about cell counts
PrediMeth_all_merge <- merge(PrediMeth_all_merge, cell_counts, by = "sample_id")

## Remove missing values from patients data (covariables)
PrediMeth_all_merge <- PrediMeth_all_merge[!PrediMeth_all_merge$predimed_high == "NULL", ] # predimeth_high
PrediMeth_all_merge <- PrediMeth_all_merge[!PrediMeth_all_merge$predimed_cat == "NULL", ] # predimeth_cat
PrediMeth_all_merge <- PrediMeth_all_merge[!PrediMeth_all_merge$predimed_score == "NULL", ] # predimeth_score
PrediMeth_all_merge <- PrediMeth_all_merge[!PrediMeth_all_merge$EDAD_ANOS == "NULL", ] # age
PrediMeth_all_merge <- PrediMeth_all_merge[!PrediMeth_all_merge$SEXO == "NULL", ] # sex
PrediMeth_all_merge <- PrediMeth_all_merge[!PrediMeth_all_merge$BMI == "NULL", ] # BMI
PrediMeth_all_merge <- PrediMeth_all_merge[!PrediMeth_all_merge$smoking_habit == "NULL", ] # smoking
PrediMeth_all_merge <- PrediMeth_all_merge[!PrediMeth_all_merge$batch == "NULL", ] # batch
cat("Pheno data dimensions after removing missing values from covariables: ", dim(PrediMeth_all_merge), "\n")
str(PrediMeth_all_merge)

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
cat("	Dimensions after filtering not baseline samples: ", 
    "\n		Bvalues = ", dim(IJC_norm.beta), 
    "\n		Pheno = ", dim(PrediMeth_all_merge), "\n")

## Removing missing values from meth set
IJC_norm.beta <- na.omit(IJC_norm.beta) # excludes all rows with at least one NA value
cat("	Bvalues dimensions without missing values: ", dim(IJC_norm.beta), "\n") 

## Filtering recommended probes (according to annotation / Zhou's recommendation -> Mgeneral column)
cat("Filtering probes according to annotation recommendation. \n")
bad_cpgs <- unique(EPIC_MEsteller$Probe_name[EPIC_MEsteller$Mgeneral %in% TRUE])    
nBetasBeforeFilter <- nrow(IJC_norm.beta)                                                                                                                          
IJC_norm.beta <- IJC_norm.beta[!rownames(IJC_norm.beta) %in% bad_cpgs, ]                                                                                                                                          
cat("Masked CpGs removed:", nBetasBeforeFilter - nrow(IJC_norm.beta), "\n")                                                                                                                                      
cat("Bvalues dimensions after filtering: ", dim(IJC_norm.beta), "\n")

## Removing chrX or chrY associated probes from dataset (according to EPICv2 annotation) 
IJC_norm.beta <- IJC_norm.beta[!(rownames(IJC_norm.beta) %in% probes_Y$Probe_name),] 
cat("Bvalues dimensions after chrY related probes filtering: ", dim(IJC_norm.beta), "\n")
IJC_norm.beta <- IJC_norm.beta[!(rownames(IJC_norm.beta) %in% probes_X$Probe_name),] 
cat("Bvalues dimensions after chrX related probes filtering: ", dim(IJC_norm.beta), "\n")

#########################################

cat("Preprocessing of meth data completed. \nStarting update of pheno data: \n")
message("Preprocessing of meth data completed. Updating pheno data: ", Sys.time())

#### UPDATING PHENOTYPE DATA ####

## PrediMeth_all_merge has some samples that have been removed from beta values due to quality check - remove them
keep_samples <- PrediMeth_all_merge$sample_id %in% colnames(IJC_norm.beta)
PrediMeth_all_merge <- PrediMeth_all_merge[keep_samples,]

## Setting rownames
rownames(PrediMeth_all_merge) <- PrediMeth_all_merge$sample_id

##########################################################

cat("Update of methylation data completed. \nStarting alignment of data: \n")
message("Methylation data updated. Starting alignment of data: ", Sys.time())

#### ALIGNING SAMPLE ORDER #### 

## Renaming objects
beta_matrix <- IJC_norm.beta
pheno <- PrediMeth_all_merge
pheno <- as.data.frame(pheno)
rownames(pheno) <- pheno$sample_id

## Check if they match
cat("\nCheck if pheno rows math beta columns: ")
rownames(pheno)[1:5]
colnames(beta_matrix)[1:5]

## Reordering phenotype table row order to match Betavalues columns
common_samples <- intersect(colnames(beta_matrix), pheno$sample_id)
pheno <- pheno[match(common_samples, pheno$sample_id), , drop = FALSE]
rownames(pheno) <- pheno$sample_id
beta_matrix <- beta_matrix[, common_samples, drop = FALSE]

## Check
cat("Samples in pheno after alignment:", nrow(pheno), "\n")
rownames(pheno)[1:5]
cat("Samples in beta_matrix after alignment:", ncol(beta_matrix), "\n")
colnames(beta_matrix)[1:5]
stopifnot(all(rownames(pheno) %in% colnames(beta_matrix)))
stopifnot(identical(colnames(beta_matrix), rownames(pheno)))
cat("Sample alignment verified.\n")

##########################################################

cat("\nAlignment of data completed. Starting convertion of Beta values to Mvalues: \n")
message("\nAlignment of data completed. Starting convertion of Beta values to Mvalues: ", Sys.time())

#### CONVERTING BETAVALUES TO M-VALUES ####

cat("Conversion of Beta-values to M-values\n")

## Check extreme beta values
cat("beta = 0 :", sum(beta_matrix == 0, na.rm = TRUE), "\n")
cat("beta = 1 :", sum(beta_matrix == 1, na.rm = TRUE), "\n")
cpgs_beta_0 <- rownames(beta_matrix)[apply(beta_matrix == 0, 1, any)]
cat("\nCpGs with at least one value equal to zero : ", cpgs_beta_0)

## Avoid logarithmic division by zero
minBeta <- min(beta_matrix[beta_matrix > 0], na.rm = TRUE)
epsilon <- 1e-6
beta_matrix[beta_matrix == 0] <- epsilon
beta_matrix[beta_matrix == 1] <- 1 - epsilon
cat("\nMinimal Beta value (different of zero) =", minBeta, "\n\n")
cat("Zero beta values clamped to ", epsilon, "\n\n")

# Conversion of Betas to Mvalues
mval_matrix <- beta2m(beta_matrix)

## Checking
cat("\nBetavalues: \n")
beta_matrix[1:5,1:5]
cat("\nM-values: \n")
mval_matrix[1:5,1:5]

cat("Infinite Mvalues: ", sum(is.infinite(mval_matrix)), "\n") # Check infinite values = 0
cat("NaN values:", sum(is.nan(mval_matrix)), "\n") # Check NaN values = 0
cat("NA values:", sum(is.na(mval_matrix)), "\n") # Check NA values = 0

## Checking for zero variance CpGs 
probe_vars <- apply(mval_matrix, 1, var, na.rm = TRUE)
zero_var_probes <- names(probe_vars[probe_vars == 0])
cat("Zero-variance CpGs:", length(zero_var_probes), "\n") # Checked: 0

##########################################################

#### ADDING VARIABLES TO PHENOTYPE TABLE ####

cat("\nAdding some variables to pheno. \n")
message("\nAdding some variables to pheno: ", Sys.time())

# Extreme exposure variable
pheno$predimed_score <- as.numeric(pheno$predimed_score)
pheno$predimed_extreme <- factor(ifelse(pheno$predimed_score <= 6, "Low_Q1",
                                        ifelse(pheno$predimed_score >= 10, "High_Q4", "Mid_Q23")),
                                 levels = c("Low_Q1", "Mid_Q23", "High_Q4"))                                                                                                                                                                                                             
cat("Distribution of predimed_extreme (low <= 6, high >= 10:\n")                                                                                                                                                                    
print(table(pheno$predimed_extreme, useNA = "ifany"))    

## Load data from CCI (Charlson Index)

cci <- fread(path_to_CCI)
cci <- as.data.frame(cci)
dim(cci)

## Merge data

cci$entity_id <- stringr::str_remove(cci$entity_id, "=*")
cci <- cci[, c("entity_id", "ehr_charlson_score")]
pheno <- merge(pheno, cci, by = "entity_id", all.x = TRUE) # merge by entity_id, preserving all samples from pheno
table(pheno$ehr_charlson_score)
rownames(pheno) <- pheno$sample_id
dim(pheno)
str(pheno)

## Check alignment of data
cat("Check alingment after adding more phenotype variables: ")
cat("\nPheno table: \n")
rownames(pheno)[1:5]
cat("\nBetavalues: \n")
colnames(beta_matrix)[1:5]
cat("\nMvalues: \n")
colnames(mval_matrix)[1:5]

## Reordering phenotype table row order to match Betavalues columns
common_samples <- intersect(colnames(beta_matrix), pheno$sample_id)
pheno <- pheno[match(common_samples, pheno$sample_id), , drop = FALSE]
rownames(pheno) <- pheno$sample_id
beta_matrix <- beta_matrix[, common_samples, drop = FALSE]
mval_matrix <- mval_matrix[, common_samples, drop = FALSE]

## Check
cat("Samples in pheno after alignment:", nrow(pheno), "\n")
rownames(pheno)[1:5]
cat("Samples in beta_matrix after alignment:", ncol(beta_matrix), "\n")
colnames(beta_matrix)[1:5]
cat("Samples in mval_matrix after alignment:", ncol(mval_matrix), "\n")
colnames(mval_matrix)[1:5]
stopifnot(all(rownames(pheno) %in% colnames(beta_matrix)))
stopifnot(identical(colnames(beta_matrix), rownames(pheno)))
stopifnot(all(rownames(pheno) %in% colnames(mval_matrix)))
stopifnot(identical(colnames(mval_matrix), rownames(pheno)))
cat("Sample alignment verified.\n")

##########################################################

#### VARIABLES PREPARATION ####

## Preparing exposure variables
pheno$predimed_high <- factor(pheno$predimed_high, levels = c(0, 1), labels = c("no", "yes")) # binary: low (0-8), high (9-14)
pheno$predimed_cat <- factor(pheno$predimed_cat, levels = c(1,2,3), labels = c("low", "medium", "high")) # 3 cat original: low (0-4), medium (5-8), high (9-14)

## Preparing variables
pheno$age <- as.numeric(pheno$EDAD_ANOS)
pheno$sex <- factor(pheno$SEXO, levels = c(1,2), labels = c("men", "women"))
pheno$bmi <- as.numeric(pheno$BMI)
pheno$smoking_type <- factor(pheno$smoking_habit, levels = c(1,2,3), labels = c("smoker","exsmoker", "nonsmoker"))
pheno$batch <- as.factor(pheno$batch)

## 14-items PREDIMED
pheno$PREDIMED_ACEITE_GRASA <- as.factor(pheno$PREDIMED_ACEITE_GRASA)
pheno$PREDIMED_ACEITE_CONSUMO <- as.factor(pheno$PREDIMED_ACEITE_CONSUMO)
pheno$PREDIMED_VERDURAS <- as.factor(pheno$PREDIMED_VERDURAS)
pheno$PREDIMED_FRUTA <- as.factor(pheno$PREDIMED_FRUTA)
pheno$PREDIMED_CARNES_ROJAS <- as.factor(pheno$PREDIMED_CARNES_ROJAS)
pheno$PREDIMED_MANTEQUILLA <- as.factor(pheno$PREDIMED_MANTEQUILLA)
pheno$PREDIMED_BEBIDAS_CARBONATADAS <- as.factor(pheno$PREDIMED_BEBIDAS_CARBONATADAS)
pheno$PREDIMED_VINO <- as.factor(pheno$PREDIMED_VINO)
pheno$PREDIMED_LEGUMBRES <- as.factor(pheno$PREDIMED_LEGUMBRES)
pheno$PREDIMED_PESCADO <- as.factor(pheno$PREDIMED_PESCADO)
pheno$PREDIMED_REPOSTERIA <- as.factor(pheno$PREDIMED_REPOSTERIA)
pheno$PREDIMED_FRUTOS_SECOS <- as.factor(pheno$PREDIMED_FRUTOS_SECOS)
pheno$PREDIMED_CARNE <- as.factor(pheno$PREDIMED_CARNE)
pheno$PREDIMED_VEGETALES <- as.factor(pheno$PREDIMED_VEGETALES)

## Checking variables preparation
cat("Structure after variables preparation: ")
str(pheno)

## Checking missing values 
table(pheno$predimed_score, useNA = "ifany")                                                                                                                                                                  
table(pheno$smoking_type, useNA = "ifany")  

##########################################################

#### SAVING OUTPUTS - PROCESSED DATA ####

cat("\nSaving aligned processed data:  ") # variables already as corrected types
save(pheno, beta_matrix, mval_matrix, file = file.path(results_folder, "processed_data.R"))

cat("\nSaving aligned pheno data: pheno (.R and .txt)")
save(pheno, file = file.path(results_folder, "disc_pheno.R"))
write.csv(pheno, file = file.path(results_folder, "disc_pheno.csv"), row.names = FALSE)

cat("\nSaving version without samples identification.")
hidden_pheno_disc <- pheno
hidden_pheno_disc$entity_id <- NULL
hidden_pheno_disc$sample_id <- NULL
write.csv(hidden_pheno_disc, file = file.path(results_folder, "without_samples_disc_pheno.csv"), row.names = FALSE)

##########################################################

cat("Script completed. \n")
message("Script completed: ", Sys.time())

##########################################################


