##########################################################
#
# Data Preprocessing (Validation cohort) - Stage 2
#
##########################################################


# DESCRIPTION --------------------------------------------
# In this script, phenotype and methylation data from the 
# validation cohort are preprocessed. 

# INPUT --------------------------------------------------
# 	Phenotype data (diet adherence + covariables)
# 	Blood cells counts

# OUTPUT -------------------------------------------------
# 	Data frame of phenotype (rows: samples, columns: variables)
# 	Matrix of beta values (rows: probe names, columns: samples)
#	  Matrix of M-values (rows: probe names, columns: samples)


##########################################################

#### LOAD LIBRARIES ####

library(data.table)
library(dplyr)
library(lubridate)
library(ggplot2)
library(lumi)

#### RESOLVE PATHS ####

## Input paths -------------------------------------------

path_to_discovery_cohort <- "/imppc/labs/dnalab/share/PrediMeth/results/OfficialAnalysis/disc_pheno.R" # .Rdata
path_to_EPIC2_annotation <- "/imppc/labs/dnalab/share/PrediMeth/0-data/EPICv2.annot.RData" # .Rdata
path_to_valid_beta_values <- "/imppc/labs/dnalab/studentdnabank/Data/omics/methylation/cnio_ia4t2d/4.norm_pc_10_IA4T2D_beta.Robj" #.Robj

path_to_valid_cell_counts <- "/imppc/labs/dnalab/share/data_storage/gcat_epic/analysis/final_normalization/IA4T2D_samples/cells_and_clocks/cell_counts.csv"
path_to_metadata <- "/imppc/labs/dnalab/share/data_storage/IA4T2D_Methylation/IA4T2D_metadata.csv"
path_to_metadata2 <- "/imppc/labs/dnalab/share/data_storage/gcat_epic/data/metadata/20250305_IA4T2D_sampleSheet.csv"
path_to_t2d_diagnosis_final <- "/imppc/labs/dnalab/share/PrediMeth/0-data/ia4t2d_date_diagnosis_final.csv"

path_to_EHR_labs <- "/imppc/labs/dnalab/share/PrediMeth/0-data/2026/AP_Laboratoris.csv"
path_to_EHR_gluc <- "/imppc/labs/dnalab/share/PrediMeth/0-data/EHR_glucemia.csv"
path_to_EHR_glycosilated <- "/imppc/labs/dnalab/share/PrediMeth/0-data/EHR_glicosilada.csv"
path_to_EHR_BMI <- "/imppc/labs/dnalab/share/PrediMeth/0-data/2026/AP_IMC.csv"
path_to_cci_valid <- "/imppc/labs/dnalab/studentdnabank/Data/EHR/charlson_elixhauser_comorbidity_index_2025.csv"

data_path <- "/imppc/labs/dnalab/studentdnabank/Data/"
methylation_dir <- file.path(data_path, 'omics', 'methylation')
questionary_dir <- file.path(data_path, 'questionnaire')

## Output paths ------------------------------------------

predimeth_path <- "/imppc/labs/dnalab/share/PrediMeth"
results_dir <- file.path(predimeth_path, "results")
results_folder <- "/imppc/labs/dnalab/share/PrediMeth/results/OfficialAnalysis/"

##########################################################

#### LOAD DATA ####

## Validation cohort -------------------------------------

# Phenotype
consulta <- fread(file.path(questionary_dir, "consulta_nivell2.csv"))
questionary1 <- fread(file.path(questionary_dir, "questionari_nivell1.csv"))
t2d_incidence <- fread(path_to_t2d_diagnosis_final)
measures1 <- fread(file.path(questionary_dir, "mesures_nivell1.csv"))
metadata <- fread(path_to_metadata)
metadata2 <- fread(path_to_metadata2)
metadata2 <- as.data.frame(metadata2)
colnames(metadata2) <- metadata2[7,]
metadata2 <- metadata2[c(8:412),c(1,4)]
metadata2 <- rename(metadata2, sample_id = Sample_Name, batch = Sample_Plate)
metadata2$batch <- as.factor(metadata2$batch)

# EHR - laboratoris
EHR_laboratoris <- fread(path_to_EHR_labs)
EHRglucemia <- fread(path_to_EHR_gluc)
EHRglucemia$V1 <- NULL
EHRglicosilada <- fread(path_to_EHR_glycosilated)
EHRglicosilada$V1 <- NULL

# EHR - BMI
EHR_bmi <- fread(path_to_EHR_BMI)

# cell counts validation cohort
cell_counts <- fread(path_to_valid_cell_counts)
cell_counts$sample_id <- cell_counts$V1
cell_counts$V1 <- NULL

# Methylation data (Betavalues)
load(path_to_valid_beta_values)
Betavalues_valid <- IA4T2D_norm.beta

## Discovery cohort --------------------------------------
load(path_to_discovery_cohort)
discovery_cohort <- pheno
dim(discovery_cohort)
discovery_cohort$entity_id <- paste0("=", discovery_cohort$entity_id)

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

##########################################################

cat("Data loaded. \n\nStarting data preprocessing. \n")
message("Data loaded. \n\nStarting data preprocessing: ", Sys.time())

#### PREPROCESSING OF PHENOTYPE DATA ####

################################################

# Preprocessing data

# Merge record linkage table and validation samples table
valid.cohort <- merge(consulta, t2d_incidence, by = "entity_id")
dim(valid.cohort) # 424 samples
valid.cohort$individual_id <- valid.cohort$'PI-2024-06'

# Merge metadata table
valid.cohort <- merge(valid.cohort, metadata, by = "individual_id")
valid.cohort <- merge(valid.cohort, metadata2, by = "sample_id")
colnames(valid.cohort) # batch = sample plate name

# Exclude samples that are in discovery cohort
table(valid.cohort$entity_id %in% discovery_cohort$entity_id) # 5 samples in discovery cohort
valid.cohort <- valid.cohort[valid.cohort$entity_id %in% discovery_cohort$entity_id == FALSE,]
dim(valid.cohort) # 395 samples

# Check if they are baseline samples (extraction between 2014-2018)
min(valid.cohort$DATA_EXTRACCIO)
max(valid.cohort$DATA_EXTRACCIO)
dim(valid.cohort)

# Selecting variables of interest
names(valid.cohort)
valid.cohort <- valid.cohort[,c("entity_id", "individual_id", "sample_id", "etnia",
                                "DATE_FINAL", "METABOLOMED", "HbA1c", "DATA_EXTRACCIO", "batch")]

# Adding info about cell counts (only samples that passed quality check)
valid.cohort <- merge(valid.cohort, cell_counts, by = "sample_id")
dim(valid.cohort) # 388 samples

## Adding more variables 

# Selecting variables from questionary and measurements
questionary1_subset <- subset(questionary1, select = c("entity_id", "EDAD_ANOS", "SEXO", "smoking_habit", 
                                                       "predimed_high", "predimed_score", "predimed_cat"))
measures1_subset <- subset(measures1, select = c("entity_id", "BMI", "WHR"))

# Selecting interest variables from consulta
valid.cohort <- merge(valid.cohort, questionary1_subset, by = "entity_id")
valid.cohort <- merge(valid.cohort, measures1_subset, by = "entity_id")
names(valid.cohort)
dim(valid.cohort)

# Remove missing values from patients data (covariables) - no samples with missing data
valid.cohort <- valid.cohort[!valid.cohort$predimed_score == "NULL", ] # predimeth_score
valid.cohort <- valid.cohort[!valid.cohort$EDAD_ANOS == "NULL", ] # age
valid.cohort <- valid.cohort[!valid.cohort$SEXO == "NULL", ] # sex
valid.cohort <- valid.cohort[!valid.cohort$BMI == "NULL", ] # BMI
valid.cohort <- valid.cohort[!valid.cohort$smoking_habit == "NULL", ] # smoking
valid.cohort <- valid.cohort[!valid.cohort$batch == "NULL", ] # batch
cat("Pheno data dimensions after removing missing values from covariables: ", dim(valid.cohort), "\n")

# Preparing variables
valid.cohort$DATA_EXTRACCIO <- as.Date(valid.cohort$DATA_EXTRACCIO, format = "%Y-%m-%d") # blood extraction date
valid.cohort$DATE_FINAL <- as.Date(valid.cohort$DATE_FINAL, format = "%Y-%m-%d") # diagnosis date
table(valid.cohort$HbA1c) # only 57 have baseline HbA1c (388-331=57)
valid.cohort$HbA1c <- as.numeric(valid.cohort$HbA1c)
str(valid.cohort)

# Checking T2D groups: incident, prevalent, control
sum(is.na(valid.cohort$DATE_FINAL)) # 187 controls (sense diabetis) 
table(valid.cohort$DATE_FINAL > valid.cohort$DATA_EXTRACCIO) # 53 incidents, 148 prevalents
valid.cohort$T2D <- ifelse(is.na(valid.cohort$DATE_FINAL), "control", 
                           ifelse(valid.cohort$DATE_FINAL > valid.cohort$DATA_EXTRACCIO, 
                                  "incident", "prevalent"))
table(valid.cohort$T2D)
valid.cohort$T2D <- as.factor(valid.cohort$T2D)
str(valid.cohort)

# Prevalent T2D (binary variable)
valid.cohort$t2d_prevalent <- ifelse(valid.cohort$T2D == "prevalent", 1, 0)
valid.cohort$t2d_prevalent <- as.factor(valid.cohort$t2d_prevalent)

# Incident T2D (binary variable)
valid.cohort$t2d_incident <- ifelse(valid.cohort$T2D == "incident", 1, 0)
valid.cohort$t2d_incident <- as.factor(valid.cohort$t2d_incident)

# Time to event (in days)
# Censoring-time => last follow-up date in EHR (administrative censoring) = 2025-12-31
last_followup <- as.Date("2025-12-31", format = "%Y-%m-%d")
valid.cohort$TTE <- ifelse(valid.cohort$t2d_prevalent == 1, 0, # prevalents (0)
                           ifelse(valid.cohort$t2d_incident == 1,
                                  (valid.cohort$DATE_FINAL - valid.cohort$DATA_EXTRACCIO), # incidents
                                  (last_followup - valid.cohort$DATA_EXTRACCIO))) # controls
table(valid.cohort$TTE) # in days
valid.cohort$TTE <- as.numeric(valid.cohort$TTE) # time to event in days (for incidents)
valid.cohort$TTE <- valid.cohort$TTE/365.25 # Converting time to event to years 
summary(valid.cohort$TTE) # mean of 5.47 years, 
sd(valid.cohort$TTE) # sd = 5.11
str(valid.cohort)

# Preparing variables
str(valid.cohort)
valid.cohort$age <- as.numeric(valid.cohort$EDAD_ANOS)
valid.cohort$bmi <- as.numeric(valid.cohort$BMI)
valid.cohort$WHR <- as.numeric(valid.cohort$WHR)
valid.cohort$sex <- factor(valid.cohort$SEXO, levels = c(1,2), labels = c("men", "women"))
valid.cohort$smoking_habit <- factor(valid.cohort$smoking_habit, levels = c(1,2,3), labels = c("smoker","exsmoker", "nonsmoker"))
valid.cohort$mets_semana <- as.numeric(valid.cohort$mets_semana)
valid.cohort$sedentarisme <- as.numeric(valid.cohort$sedentarisme)
valid.cohort$predimed_high <- factor(valid.cohort$predimed_high, levels = c(0,1), labels = c("no", "yes"))
valid.cohort$predimed_cat <- factor(valid.cohort$predimed_cat, levels = c(1,2,3), labels = c("low", "medium", "high"))
valid.cohort$predimed_score <- as.numeric(valid.cohort$predimed_score)
valid.cohort$HbA1c_baseline <- valid.cohort$HbA1c
valid.cohort$HbA1c <- NULL
str(valid.cohort)


## Creating additional variables of interest

# age of onset
valid.cohort$age_diagnosis <- ifelse(valid.cohort$T2D=="incident", 
                                     (valid.cohort$age + time_length(interval(valid.cohort$DATA_EXTRACCIO, valid.cohort$DATE_FINAL), "years")),
                                     NA)
summary(valid.cohort$age_diagnosis) # 53 incidents, con respectiva edad del diagnostico
valid.cohort$OnsetBefore45 <- ifelse(valid.cohort$age_diagnosis < 45, TRUE, FALSE)

# incident cases
incidents_id <- valid.cohort$entity_id[valid.cohort$T2D=="incident"] # entity id dos incidents
length(incidents_id)
diagnosis_date <- valid.cohort[,c("entity_id", "DATE_FINAL")]
diagnosis_date$DATE_DIAGNOSIS <- diagnosis_date$DATE_FINAL
diagnosis_date$DATE_FINAL <- NULL
names(diagnosis_date)

# BMI at onset **************** see in EHR
test <- EHR_bmi %>% filter(entity_id %in% incidents_id) %>% merge(diagnosis_date, by = "entity_id") %>% 
  filter(Prova_data >= (DATE_DIAGNOSIS %m-% months(9)) & Prova_data <= (DATE_DIAGNOSIS %m+% months(3)))
length(unique(test$entity_id)) # 31 / 53 tienen BMI del onset
test <- test %>% mutate( days_from_BMIonset = as.numeric(Prova_data - DATE_DIAGNOSIS) ) %>%
  slice_min(abs(days_from_BMIonset), by = entity_id, n = 1, with_ties = FALSE)
dim(test)
test <- test %>% rename(BMIonset_date = Prova_data, BMIonset = Prova_resultat)
colnames(test)
test <- test[,c("entity_id", "BMIonset", "BMIonset_date")]
dim(valid.cohort)
valid.cohort <- merge(valid.cohort, test, by = "entity_id", all.x = TRUE)
valid.cohort$OnsetBMIless27 <- ifelse(valid.cohort$BMIonset < 27, TRUE, FALSE)
table(valid.cohort$OnsetBMIless27) # 8 true, 23 false
str(valid.cohort)


# HbA1c (EHR)

# glicada al baseline: 7/56 incidents, 47/162 prevalents, 10/201 controls
table(valid.cohort$T2D[!is.na(valid.cohort$HbA1c)]) 
valid.cohort$has_baseline_HbA1c <- ifelse(is.na(valid.cohort$HbA1c_baseline), FALSE, TRUE) 

# glicada al diagnostic dels incidents
dim(valid.cohort)
dim(EHRglicosilada)
glicosilada <- EHRglicosilada[EHRglicosilada$entity_id %in% incidents_id] # filter validation samples
dim(glicosilada)
length(unique(glicosilada$entity_id)) # 52 / 56 incidents tenen glicosilada al EHR
# filer glicosilada from diagnosis (6 months before - 1 months after) # 36/56
names(glicosilada)
glicosilada$DATE_HbA1c <- as.Date(glicosilada$DATA)
glicosilada$DATA <- NULL
glicosilada <- glicosilada %>% left_join(diagnosis_date, by = "entity_id")
test <- glicosilada %>% 
  filter(DATE_HbA1c >= (DATE_DIAGNOSIS %m-% months(6)) & DATE_HbA1c <= (DATE_DIAGNOSIS %m+% months(1)))
dim(test)
length(unique(test$entity_id))

test <- test %>% mutate( days_from_diagnosis_HbA1c = as.numeric(DATE_HbA1c - DATE_DIAGNOSIS) ) %>%
  slice_min(abs(days_from_diagnosis_HbA1c), by = entity_id, n = 1, with_ties = FALSE)

head(test)
test$HbA1c_diagnostic <- test$Resultat
test$HbA1c_diagnostic_unitat <- test$Unitats
names(test)
test <- test[,c("entity_id","DIABETES","PREDIABETES","DATE_HbA1c","DATE_DIAGNOSIS",
                "days_from_diagnosis_HbA1c","HbA1c_diagnostic","HbA1c_diagnostic_unitat")]

# merge with validation table
valid.cohort <- valid.cohort %>% left_join(test, by = "entity_id")
names(valid.cohort)
dim(valid.cohort)
valid.cohort$DIABETES_per_HbA1c <- valid.cohort$DIABETES
valid.cohort$PREDIABETES_per_HbA1c <- valid.cohort$PREDIABETES
valid.cohort$DIABETES <- NULL
valid.cohort$PREDIABETES <- NULL

# Blood glucose

# glucosa al diagnostic dels incidents
dim(valid.cohort)
dim(EHRglucemia)
incidents_id <- valid.cohort$entity_id[valid.cohort$T2D=="incident"] # entity id dos incidents
length(incidents_id)
glucemia <- EHRglucemia[EHRglucemia$entity_id %in% incidents_id] # filter validation samples
dim(glucemia)
length(unique(glucemia$entity_id)) # 56 / 56 incidents tenen glucemia al EHR
# filer glicosilada from diagnosis (6 months before - 1 months after) # 44/56
names(diagnosis_date)
names(glucemia)
glucemia$DATE_glucemia <- as.Date(glucemia$DATA)
glucemia$DATA <- NULL
glucemia <- glucemia %>% left_join(diagnosis_date, by = "entity_id")
test <- glucemia %>% 
  filter(DATE_glucemia >= (DATE_DIAGNOSIS %m-% months(6)) & DATE_glucemia <= (DATE_DIAGNOSIS %m+% months(1)))
dim(test)
length(unique(test$entity_id))

test <- test %>% mutate( days_from_diagnosis_gluc = as.numeric(DATE_glucemia - DATE_DIAGNOSIS) ) %>%
  slice_min(abs(days_from_diagnosis_gluc), by = entity_id, n = 1, with_ties = FALSE)
test$DATE_DIAGNOSIS <- NULL

head(test)
test$gluc_diagnostic <- test$Resultat
test$gluc_diagnostic_unitat <- test$Unitats
test$DIABETES_per_gluc <- test$DIABETES
test$PREDIABETES_per_gluc <- test$PREDIABETES

names(test)
test <- test[,c("entity_id","DIABETES_per_gluc","PREDIABETES_per_gluc","DATE_glucemia", "days_from_diagnosis_gluc","gluc_diagnostic","gluc_diagnostic_unitat")]


# merge with validation table
valid.cohort <- valid.cohort %>% left_join(test, by = "entity_id")
names(valid.cohort)
dim(valid.cohort)
str(valid.cohort)

#### ADDING VARIABLES TO PHENOTYPE TABLE ####

cat("\nAdding some variables to pheno. \n")
message("\nAdding some variables to pheno: ", Sys.time())

valid.cohort$predimed_score <- as.numeric(valid.cohort$predimed_score)
pheno$predimed_extreme <- factor(ifelse(pheno$predimed_score <= 6, "Low_Q1",
                                        ifelse(pheno$predimed_score >= 10, "High_Q4", "Mid_Q23")),
                                 levels = c("Low_Q1", "Mid_Q23", "High_Q4"))                                                                                                                                                                                                             
cat("Distribution of predimed_extreme (low <= 6, high >= 10:\n")                                                                                                                                                                    
print(table(pheno$predimed_extreme, useNA = "ifany"))    

## Load data from CCI (Charlson Index)

cci <- fread(path_to_cci_valid)
cci <- as.data.frame(cci)
dim(cci)

## Merge data

cci <- cci[, c("entity_id", "ehr_charlson_score")] # entity_id begins with =
valid.cohort <- merge(valid.cohort, cci, by = "entity_id", all.x = TRUE) # merge by entity_id, preserving all samples from valid.cohort
cat("\nDistribution Charlson Score: ")
table(valid.cohort$ehr_charlson_score)
rownames(valid.cohort) <- valid.cohort$sample_id
dim(valid.cohort)
str(valid.cohort)

####################################################

### PREPROCESSING OF METHYLATION DATA ###

## Preprocessing of Betavalues from validation cohort

## Checking rownames and colnames
rownames(Betavalues_valid)[1:5]
colnames(Betavalues_valid)[1:5]

## Filtering valid.cohort samples
Betavalues_valid <- subset(Betavalues_valid, select = c(colnames(Betavalues_valid) %in% rownames(valid.cohort)))
cat("Dimensions phenotype table (validation):", dim(valid.cohort), "\n")

## Filtering recommended probes (according to annotation -> Mgeneral column = )
cat("Filtering probes according to annotation recommendation. \n")
bad_cpgs <- unique(EPIC_MEsteller$Probe_name[EPIC_MEsteller$Mgeneral %in% TRUE])    
nBetasBeforeFilter <- nrow(Betavalues_valid)                                                                                                                          
Betavalues_valid <- Betavalues_valid[!rownames(Betavalues_valid) %in% bad_cpgs, ]                                                                                                                                          
cat("Masked CpGs removed:", nBetasBeforeFilter - nrow(Betavalues_valid), "\n")                                                                                                                                      
cat("Bvalues dimensions after filtering: ", dim(Betavalues_valid), "\n")

## Removing chrX or chrY associated probes from dataset (according to EPICv2 annotation) 
Betavalues_valid <- Betavalues_valid[!(rownames(Betavalues_valid) %in% probes_Y$Probe_name),] 
cat("Bvalues dimensions after chrY related probes filtering: ", dim(Betavalues_valid), "\n")
Betavalues_valid <- Betavalues_valid[!(rownames(Betavalues_valid) %in% probes_X$Probe_name),] 
cat("Bvalues dimensions after chrX related probes filtering: ", dim(Betavalues_valid), "\n")

## Alignment of Betavalues and valid.cohort ##

## Check if they match
rownames(valid.cohort)[1:5]
colnames(Betavalues_valid)[1:5]

## Reordering phenotype table row order to match Betavalues columns
common_samples <- intersect(colnames(Betavalues_valid), valid.cohort$sample_id)
valid.cohort <- valid.cohort[match(common_samples, valid.cohort$sample_id), , drop = FALSE]
rownames(valid.cohort) <- valid.cohort$sample_id
Betavalues_valid <- Betavalues_valid[, common_samples, drop = FALSE]

## Check
cat("Samples in pheno after alignment:", nrow(valid.cohort), "\n")
rownames(valid.cohort)[1:5]
cat("Samples in beta_matrix after alignment:", ncol(Betavalues_valid), "\n")
colnames(Betavalues_valid)[1:5]
stopifnot(all(rownames(valid.cohort) %in% colnames(Betavalues_valid)))
stopifnot(identical(colnames(Betavalues_valid), rownames(valid.cohort)))
cat("Sample alignment verified.\n")

## Prepare beta-values = 0 and beta-values = 1 (so they don't convert into infinite M-values)

cat("Conversion of Beta-values to M-values\n")
cat("beta = 0 :", sum(Betavalues_valid == 0, na.rm = TRUE), "\n")
cat("beta = 1 :", sum(Betavalues_valid == 1, na.rm = TRUE), "\n")

## 1) Avoid logarithmic division by zero
minBeta <- min(Betavalues_valid[Betavalues_valid > 0], na.rm = TRUE)
epsilon <- 1e-6
Betavalues_valid[Betavalues_valid == 0] <- epsilon
Betavalues_valid[Betavalues_valid == 1] <- 1 - epsilon
cat("\nMinimal Beta value (different of zero) =", minBeta, "\n\n")
cat("Zero beta values clamped to ", epsilon, "\n\n")

# Conversion of Betas to Mvalues
Mvalues_valid <- beta2m(Betavalues_valid)

## Check if infinitive values were produced 
inf_values <- !is.finite(Mvalues_valid)
cat("Infinite Mvalues: ", sum(inf_values), "\n") # Check infinite values
cat("NAs in Betavalues_valid:", sum(is.na(Betavalues_valid)), "\n") # Infinite values due to NAs
n_inf <- sum(is.infinite(Mvalues_valid), na.rm = TRUE) # true inf values (due to 0 or 1 betavalues)
cat("True Inf/-Inf values:", n_inf, "\n")  # should be 0 after clamping
n_na <- sum(is.na(Mvalues_valid)) # NA values
cat("NA values:", n_na, "\n")  # inherited from beta matrix
n_nan <- sum(is.nan(Mvalues_valid))
cat("NaN values:", n_nan, "\n") # from 0/0 or similar
cat("Total non-finite:", n_inf + n_na + n_nan, "\n")

## If n_inf is 0: NAs are expected missing values from probe-level QC failures in the validation cohort. 
## These will be handled in the MRS_calculator function, which skips missing values rather than propagating them into the score.

## Checking rownames 
cat("Checking methylation matrix rownames: \n")
Betavalues_valid[1:5, 1:5]
Mvalues_valid[1:5, 1:5]

##########################################################

#### SAVING OUTPUTS ####

cat("\nSaving aligned processed data: valid.cohort + Betavalues_valid + Mvalues_valid")
save(valid.cohort, Betavalues_valid, Mvalues_valid, file = file.path(results_folder, "processed_data_stage2.R"))

# Save valid.cohort dataframe
write.table(valid.cohort, file = file.path(results_folder, "validation_cohort.txt"), 
            sep = "\t", row.names = FALSE, col.names = TRUE)

##########################################################

cat("Script completed. \n")
message("Script completed: ", Sys.time())

##########################################################
