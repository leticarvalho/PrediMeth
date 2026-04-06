##########################################################
#
# Fitting Model - Stage 1
#
##########################################################


# DESCRIPTION --------------------------------------------
# In this script, the processed data is used to fit the 
# best linear model (selected after sensitivity analyses) 
# to each CpG, using limma R package. For this model, a
# topTable will be generated, including EPIC v2 annotation 
# and calculated absolute delta beta (ADB). 
# Selected model: Mvalues ~ diet + age + sex + smoking + BMI + cell counts
# Obs. diet may be binary (predimed_high) or multifactor (predimed_cat).

# INPUT --------------------------------------------------
# 	Data frame of phenotype (rows: samples, columns: variables)
# 	Matrix of beta values (rows: probe IDs, columns: samples)
#	  Matrix of M-values (rows: probe IDs, columns: samples)

# OUTPUT -------------------------------------------------
# 	Results table (topTable + annotation + ADB)
#   Bonferroni's threshold


##########################################################


#### LOAD LIBRARIES ####

library(dplyr)
library(data.table)
library(limma)
library(EnhancedVolcano)
library(qqman)

#### RESOLVE PATHS ####

## Define label ------------------------------------------

# Selected model: 
# βCpG ∼ High adherence to Med diet + Age + Sex + BMI + Smoking + Cell counts + SVs (full model) 

label <- "_3cat" # Probe ID
#label <- "_bin" # Probe ID

## Input paths -------------------------------------------

annotation_path <- "/path_to_EPICv2_annotation.R" # R object with EPICv2 annotation

predimeth_path <- "/path_to_project_folder/" # path to project folder
results_dir <- file.path(predimeth_path, "results") # path to results folder
processed_data <- file.path(results_dir, "ProcessedData_Stage1/processed_data.R")

## Output paths ------------------------------------------

results_folder <- file.path(results_dir, 'LimmaResults_Stage1')
dir.create(results_folder)

##########################################################

#### LOAD DATA ####

## Preprocessed data (Stage 1) ---------------------------

load(processed_data)

## Renaming objects
Mvalues <- IJC_m_values 
Betavalues <- IJC_norm.beta 
pheno <- PrediMeth_all_merge
rownames(pheno) <- pheno$sample_id 

## Looking at dimensions
print("Dimensions of loaded data:")
dim(IJC_norm.beta)
dim(IJC_m_values)
dim(PrediMeth_all_merge)

## Annotation (EPIC v2) ----------------------------------

load(annotation_path)
EPIC_MEsteller$Probe_name <- stringr::str_remove(EPIC_MEsteller$ProbeID, "_.*") 
print("Preparing annotation to add to topTable: \n")
EPIC_sub <- EPIC_MEsteller[match(rownames(Mvalues), EPIC_MEsteller$ProbeID),] # when rows = probe IDs

##########################################################

cat("Data loaded. \n\nStarting model design. \n")
message("Data loaded. \n\nStarting model design: ", Sys.time())

#### MODEL DESIGN ####

## Preparing variables

#~ predimed_type <- factor(pheno$predimed_high, levels = c(0, 1), labels = c("no", "yes"))
predimed_categ <- factor(pheno$predimed_cat, levels = c(1,2,3), labels = c("low", "medium", "high")) # "1=low [0-5), 2=medium [5-9), 3=high" [9-Inf

age <- as.numeric(pheno$EDAD_ANOS)
sex <- factor(pheno$SEXO, levels = c(1,2), labels = c("men", "women"))
bmi <- as.numeric(pheno$BMI)
smoking_type <- factor(pheno$smoking_habit, levels = c(1,2,3), labels = c("smoker","exsmoker", "nonsmoker"))

## Creating the design matrix -----------------------------

#~ design_matrix <- model.matrix(~0 + predimed_type + age + sex + bmi + smoking_type + Bcell + CD4T + CD8T + Mono + Neu + NK, data = pheno) 
#~ cat("\n Head of design matrix: \n")
#~ head(design_matrix)
#~ colnames(design_matrix)
# --------------
design_matrix <- model.matrix(~0 + predimed_categ + age + sex + bmi + smoking_type + Bcell + CD4T + CD8T + Mono + Neu + NK, data = pheno) 

cat("\n Head of design matrix: \n")
head(design_matrix)
colnames(design_matrix)

# Checking groups dimensions ---------------------------

#~ list_yes <- rownames(subset(PrediMeth_all_merge, predimed_type=="yes")) # yes
#~ cat("	Amount of high adherents to Med diet: ", length(list_yes), "\n")
#~ list_no <- rownames(subset(PrediMeth_all_merge, predimed_type=="no")) # no
#~ cat("	Amount of low adherents to Med diet: ", length(list_no), "\n")
# --------------
list_high <- rownames(subset(PrediMeth_all_merge, predimed_categ=="high")) # yes
cat("	Amount of high adherents to Med diet: ", length(list_high), "\n")
list_medium <- rownames(subset(PrediMeth_all_merge, predimed_categ=="medium")) # no
cat("	Amount of medium adherents to Med diet: ", length(list_medium), "\n")
list_low <- rownames(subset(PrediMeth_all_merge, predimed_categ=="low")) # no
cat("	Amount of low adherents to Med diet: ", length(list_low), "\n")

## Creating the contrast matrix ---------------------------

#~ contrast_matrix <- makeContrasts(predimed_typeyes_vs_predimed_typeno = predimed_typeyes-predimed_typeno, levels = design_matrix)
#~ contrast <- colnames(contrast_matrix) 
#~ cat(" Name of contrast: \n", contrast, "\n")
# --------------
levels(predimed_categ)
contrast_matrix <- makeContrasts(
  High_vs_Low = predimed_categhigh - predimed_categlow,
 	Medium_vs_Low = predimed_categmedium - predimed_categlow,
 	High_vs_Medium = predimed_categhigh - predimed_categmedium, 
 	levels = design_matrix)
contrasts <- colnames(contrast_matrix) 
cat(" Names of contrasts: \n", contrasts, "\n")

##########################################################

cat("Model design completed. \n\nStarting model fitting. \n")
message("Model design completed. \n\nStarting model fitting: ", Sys.time())

#### MODEL FITTING ####

## Checking dimensions
dim(Mvalues)
dim(pheno)

# Defining rownames 
rownames(pheno) <- pheno$sample_id

## Fitting the linear model to the data with limma's function
fit <- lmFit(Mvalues, design_matrix, na.action = na.exclude)

## Applying contrasts
fit2 <- contrasts.fit(fit, contrast_matrix)
summary(fit2$coefficients)

## Applying eBayes (trend = TRUE models mean-variance relationship) 
fit3 <- eBayes(fit2, trend = TRUE)
colnames(fit3$coefficients)
colnames(contrast_matrix)
dim(fit3$coefficients)
class(fit3)

##########################################################

cat("Model fitting completed. \n\nStarting results analysis. \n")
message("Model fitting completed. \n\nStarting results analysis: ", Sys.time())

#### RESULTS ANALYSIS ####

## Multiple testing
#~ results_fit <- decideTests(fit3)
#~ cat("Results summary: \n")
#~ summary(results_fit)

## Full Top Table with annotation  ---------------------------

#~ full_topTable <- limma::topTable(fit3, number = Inf,
#~                                  adjust.method = "BH", # Benjamini-Hochberg
#~                                  coef = contrast, sort.by = "p", 
#~                                  genelist = EPIC_sub)
High_vs_Low_topTable <- limma::topTable(fit3, number = Inf, adjust.method = "BH", coef = "High_vs_Low", sort.by = "p", genelist = EPIC_sub)
Medium_vs_Low_topTable <- limma::topTable(fit3, number = Inf, adjust.method = "BH", coef = "Medium_vs_Low", sort.by = "p", genelist = EPIC_sub)
High_vs_Medium_topTable <- limma::topTable(fit3, number = Inf, adjust.method = "BH", coef = "High_vs_Medium", sort.by = "p", genelist = EPIC_sub)

## Calculating and including ADB (Absolute Delta Beta) ---------------------------

#~ mean_cases <- rowMeans(Betavalues[, pheno$predimed_high == 1], na.rm = TRUE) # high diet adherence
#~ mean_controls <- rowMeans(Betavalues[, pheno$predimed_high == 0], na.rm=TRUE) # low diet adherence
#~ delta_beta <- mean_cases - mean_controls
#~ idx <- match(full_topTable$ProbeID, names(delta_beta))
#~ full_topTable$delta_beta <- delta_beta[idx]
#~ full_topTable$abs_delta_beta <- abs(full_topTable$delta_beta)
# -------------------
mean_low <- rowMeans(Betavalues[, pheno$predimed_cat == 1], na.rm=TRUE)
mean_medium <- rowMeans(Betavalues[, pheno$predimed_cat == 2], na.rm=TRUE)
mean_high <- rowMeans(Betavalues[, pheno$predimed_cat == 3], na.rm = TRUE) 
# HL
delta_beta_HL <- mean_high - mean_low
idx <- match(High_vs_Low_topTable$ProbeID, names(delta_beta_HL))
High_vs_Low_topTable$delta_beta <- delta_beta_HL[idx]
High_vs_Low_topTable$abs_delta_beta <- abs(High_vs_Low_topTable$delta_beta)
# ML
delta_beta_ML <- mean_medium - mean_low
idx <- match(Medium_vs_Low_topTable$ProbeID, names(delta_beta_ML))
Medium_vs_Low_topTable$delta_beta <- delta_beta_ML[idx]
Medium_vs_Low_topTable$abs_delta_beta <- abs(Medium_vs_Low_topTable$delta_beta)
# HM
delta_beta_HM <- mean_high - mean_medium
idx <- match(High_vs_Medium_topTable$ProbeID, names(delta_beta_HM))
High_vs_Medium_topTable$delta_beta <- delta_beta_HM[idx]
High_vs_Medium_topTable$abs_delta_beta <- abs(High_vs_Medium_topTable$delta_beta)


## Subset of FDR-significant hits 
#~ fdr_subset <- subset(full_topTable, abs(logFC)>0.26 & adj.P.Val<0.05)
fdr_HL <- subset(High_vs_Low_topTable, abs(logFC)>0.26 & adj.P.Val<0.05)
fdr_ML <- subset(Medium_vs_Low_topTable, abs(logFC)>0.26 & adj.P.Val<0.05)
fdr_HM <- subset(High_vs_Medium_topTable, abs(logFC)>0.26 & adj.P.Val<0.05)

## Bonferroni-threshold calculation
bonf_threshold <- 0.05 / dim(Mvalues)[1] # calculating Bonferroni's threshold

## Hits
#~ cat("FDR hits (binary):", dim(fdr_subset)[1], "\n")
cat("FDR hits: \n",
	"High_vs_Low: ", dim(fdr_HL)[1]," hits \n",
	"Medium_vs_Low: ", dim(fdr_ML)[1]," hits \n",
	"High_vs_Medium: ", dim(fdr_HM)[1]," hits \n")

##########################################################

cat("Results analysis completed. \n\nSaving outputs. \n")
message("Results analysis completed. \n\nSaving outputs: ", Sys.time())

#### SAVING OUTPUTS ####

## R object

#~ save(full_topTable, bonf_threshold,
#~      file = file.path(results_folder, paste0("limma_results", label, ".R")))
     
save(High_vs_Low_topTable,
	 Medium_vs_Low_topTable,
	 High_vs_Medium_topTable,
	 bonf_threshold,
     file = file.path(results_folder, paste0("limma_results", label, ".R")))

cat("Script completed. \n")
message("Script completed: ", Sys.time())

##########################################################
