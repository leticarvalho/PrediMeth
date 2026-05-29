##########################################################
#
# MRS construction and calculation (Stage 2)
#
##########################################################


# DESCRIPTION --------------------------------------------
# In this script, Methylation Risk Score (MRS) is constructed
# with different CpG selection strategies, and calculated
# for each sample from the validation cohort. 
# CpG selection strategies:
# 1) MRS constructed with FDR hits from standard contrasts (HL + ML)
#       using HL delta beta as weight (from discovery EWAS)
# 2) MRS constructed with top CpGs from continuous contrast,
#       with pvalue < 1e-5,
#       using Continuous delta beta per unit as weight (from discovery EWAS)
# 3) MRS constructed with top CpGs from continuous contrast,
#       with pvalue < 1e-4,
#       using Continuous delta beta per unit as weight (from discovery EWAS)
# 4) MRS constructed with an elastic net, with all CpGs,
#       with alfa = 0.5, lambda min, 10-fold cv, residualizing for covariates

# INPUT --------------------------------------------------
#   From discovery cohort: 
#       Pheno table
#       Betavalues
#       Limma's topTable (continuous model)
#       Limma's topTable (High vs low - 3 categories model)
#   From validation cohort: 
#       Pheno data
#       Betavalues

# OUTPUT -------------------------------------------------
# 	MRSs calculated for validation cohort (csv)
#   MRS names and weights (.R data)


##########################################################

#### LOAD LIBRARIES ####

library(dplyr)
library(data.table)
library(glmnet) # elastic net
library(ggplot2)

#### RESOLVE PATHS ####

## Input paths -------------------------------------------

path_to_processed_data_stage1 <- "" # Rdata
path_to_limma_results_continuous_stage_1 <- "" # Rdata
path_to_limma_results_3catoriginal_stage_1 <- ""  # Rdata
path_to_processed_data_stage2 <- ""  # Rdata

## Output paths ------------------------------------------

predimeth_path <- ""
results_dir <- file.path(predimeth_path, "results")
results_folder <- file.path(results_dir, "ValidationCohort")

##########################################################

#### LOAD DATA ####

## Discovery cohort --------------------------------------

# Processed data Stage 1 = pheno + beta_matrix + mval_matrix
load(path_to_processed_data_stage1) 
pheno_disc <- pheno
rownames(pheno_disc) <- pheno_disc$sample_id
Betavalues_disc <- beta_matrix

# Limma results from continuous model = full_topTable1, full_topTable2, full_topTable3,
load(path_to_limma_results_continuous_stage_1)

## Limma results from 3 categories original contrast, full model = High_vs_Low_topTable1, Medium_vs_Low_topTable1, High_vs_Medium_topTable1
load(path_to_limma_results_3catoriginal_stage_1)

## Validation cohort -------------------------------------

# Processed data Stage 2 = valid.cohort + Betavalues_valid + Mvalues_valid
load(path_to_processed_data_stage2)
valid.cohort <- as.data.frame(valid.cohort)
rownames(valid.cohort) <- valid.cohort$sample_id

##########################################################

cat("Data loaded. \n\nStarting data preparation. \n")
message("Data loaded. \n\nStarting data preparation: ", Sys.time())

## Checking validation cohort pheno data and Betavalues matrix

# Dimensions
cat("Dimensions of pheno validation table: ", dim(valid.cohort), "\n")
cat("Dimensions of Mvalues validation samples: ", dim(Mvalues_valid), "\n")

# Rownames
rownames(valid.cohort) <- valid.cohort$sample_id

# Checking alignment
cat("Checking alignment of pheno samples (rows) and Betavalues samples (columns): ", 
    "\n First pheno samples: ", rownames(valid.cohort)[1:5], "\n",
    "\n First Betavalues samples: ", colnames(Betavalues_valid)[1:5], "\n")
stopifnot(all(rownames(valid.cohort) %in% colnames(Betavalues_valid)))
stopifnot(identical(colnames(Betavalues_valid), rownames(valid.cohort)))
cat("Sample alignment verified.\n")


## Checking topTables

# Continuous model -->  weight = $delta_beta_per_unit
cont_topTable <- full_topTable1
cont_topTable <- as.data.frame(cont_topTable)
cat("Dimensions Continuous topTable: ", dim(cont_topTable), "\n")
rownames(cont_topTable) <- cont_topTable$Probe_name 

# HL FDR -->  weight = $delta_beta
HL_full <- as.data.frame(High_vs_Low_topTable1)
rownames(HL_full) <- HL_full$Probe_name
HL_lookup <- setNames(HL_full$delta_beta, HL_full$Probe_name)
fdr_HL1 <- as.data.frame(subset(High_vs_Low_topTable1, adj.P.Val<0.05))
fdr_ML1 <- as.data.frame(subset(Medium_vs_Low_topTable1, adj.P.Val<0.05))
rownames(fdr_HL1) <- fdr_HL1$Probe_name 
rownames(fdr_ML1) <- fdr_ML1$Probe_name 

##########################################################

cat("Data preparation concluded. \n\nStarting MRS construction. \n")
message("Data preparation concluded. \nStarting MRS construction: ", Sys.time())

#### MRS CONSTRUCTION ####

# MRS calculator (computes a weighted MRS for each sample)
# Parameters: 
#   Betavalues_valid = Betavalue matrix from valid cohort (CpG x samples)
#   cpgs = selected cpgs to form MRS
#   weights named numeric vector of weights (names = CpG names)
# Return: named numeric vector (one MRS per sample)
MRS_calculator <- function(Betavalues_valid, cpgs, weights) {
  cpgs_m <- intersect(cpgs, rownames(Betavalues_valid))
  missing <- length(cpgs) - length(cpgs_m)
  if(missing > 0) {warning(missing, " CpGs absent from Betavalues matrix and excluded.")}
  w <- weights[cpgs_m] 
  b <- Betavalues_valid[cpgs_m, , drop = FALSE]
  sapply(colnames(b), function(s) sum(w * b[, s], na.rm = TRUE))
}


## WEIGHTS DEFINITION ACCORDING TO STRATEGY #### 


# 1) MRS constructed with FDR hits from standard contrasts (HL + ML)
#     using delta-betas from HL topTable as weights

cpgs_1 <- unique(c(fdr_HL1$Probe_name, fdr_ML1$Probe_name))
weights_1 <- HL_lookup[cpgs_1]
weights_1 <- weights_1[!is.na(weights_1)]
cat("\nStrategy 1 - HL+ML FDR hits: N CpGs = ", length(cpgs_1))


# 2) MRS constructed with CpGs with pvalue < 1e-5 from the continuous contrast 
#     using delta_beta_per_unit from continuous topTable as weights

cont_pval5 <- cont_topTable[cont_topTable$P.Value < 1e-5,]
rownames(cont_pval5) <- cont_pval4$Probe_name

cpgs_2 <- cont_pval5$Probe_name
weights_2 <- setNames(cont_topTable[cpgs_2, "delta_beta_per_unit"], cpgs_2) 
cat("\nStrategy 2 - CpGs with pvalue < 1e-5 from the continuous contrast: N CpGs = ", length(cpgs_2))

# 3) MRS constructed with CpGs with pvalue < 1e-4 from the continuous contrast 
#     using delta_beta_per_unit from continuous topTable as weights

cont_pval4 <- cont_topTable[cont_topTable$P.Value < 1e-4,]
rownames(cont_pval4) <- cont_pval4$Probe_name

cpgs_3 <- cont_pval4$Probe_name
weights_3 <- setNames(cont_topTable[cpgs_3, "delta_beta_per_unit"], cpgs_3) 
cat("\nStrategy 3 - CpGs with pvalue < 1e-4 from the continuous contrast: N CpGs = ", length(cpgs_3))


# 4) MRS constructed with an elastic net, with all CpGs,
#       with alfa = 0.5, lambda min, 10-fold cv, residualizing for covariates

X <- t(Betavalues_disc) # transpose beta vlaues from discovery cohort
y <- pheno_disc$predimed_score[match(rownames(X), rownames(pheno_disc))] # outcome = predimed_score

cov_cols_disc <- c("age", "sex", "bmi", "smoking_type", "batch", "Bcell", "CD4T", "CD8T", "Mono", "NK")
covs_disc <- pheno_disc[match(rownames(X), rownames(pheno_disc)), cov_cols_disc]
y_resid <- residuals(lm(y ~ ., data = covs_disc)) # residualizing outcome (removing covariates effect)

keep <- complete.cases(X) & !is.na(y_resid) # excluding missing values
X_use <- X[keep, , drop = FALSE] # updating X
y_use <- y_resid[keep] # updating outcome y
cat("Samples used:", sum(keep), "/", length(y_resid), "\n")

set.seed(123) # setting random seed
cv_fit <- cv.glmnet(X_use, y_use, alpha = 0.5, family = "gaussian", standardize = TRUE, nfolds = 10)
coefs_min <- coef(cv_fit, s = "lambda.min")
nz_min <- coefs_min[coefs_min[, 1] != 0, , drop = FALSE] 
nz_min <- nz_min[rownames(nz_min) != "(Intercept)", , drop = FALSE]

weights_4 <- setNames(as.numeric(nz_min), rownames(nz_min))
cat("Non-zero CpGs (lambda.min):", length(weights_4), "\n")

##########################################

cat("MRS construction completed. \n\nStarting MRS calculation for the validation cohort. \n")
message("MRS construction completed. \nStarting MRS calculation for the validation cohort: ", Sys.time())

## Calculate MRSs from validation cohort (apply all weight sets and add to dataframe)

valid.cohort$MRS1 <- MRS_calculator(Betavalues_valid, names(weights_1), weights_1)
valid.cohort$MRS2 <- MRS_calculator(Betavalues_valid, names(weights_2), weights_2)
valid.cohort$MRS3 <- MRS_calculator(Betavalues_valid, names(weights_3), weights_3)
# valid.cohort$MRS4 <- MRS_calculator(Betavalues_valid, names(weights_4), weights_4) --> no CpGs

##########################################

cat("Calculation concluded. \n\nStarting tracking of non-missing CpGs. \n")
message("Calculation concluded. \nStarting tracking of non-missing CpGs: ", Sys.time())

## Track number of non-missing CpGs contributing to each MRS per sample
# Following the approach of Davyson et al. (AD MRS script - https://github.com/Elladavyson/AD_MRS/):

max_missing <- 0.10  # flag samples missing > 10% of MRS CpGs

count_contributing_cpgs <- function(Betavalues_valid, cpgs) {
  cpgs_m <- intersect(cpgs, rownames(Betavalues_valid))
  betaval_cpgs <- Betavalues_valid[cpgs_m, , drop = FALSE]
  count_used_cpgs <- colSums(!is.na(betaval_cpgs)) # for each sample, count how many CpGs are non-missing
  return(count_used_cpgs)
}

# Count for each strategy and add it in the pheno df
# Also report samples where MRS is based on < 90% of its CpGs

## S1
n_cpgs_S1 <- count_contributing_cpgs(Betavalues_valid, names(weights_1))
valid.cohort$n_cpgs_MRS1 <- n_cpgs_S1[rownames(valid.cohort)]
n_S1_expected <- length(weights_1)
low_coverage  <- sum(n_cpgs_S1 < (1 - max_missing) * n_S1_expected, na.rm = TRUE)
cat("\nSamples with <", (1 - max_missing) * 100, "% of MRS1 CpGs available:", low_coverage, "\n")
# Distribution plot of contributing CpGs for MRS1
cpg_contrib_df <- data.frame(sample = names(n_cpgs_S1), n_contrib = as.numeric(n_cpgs_S1), expected = n_S1_expected)
plot_contrib_S1 <- ggplot(cpg_contrib_df, aes(x = n_contrib)) + geom_histogram(bins = 20, fill = "#2166AC", alpha = 0.7) +
  geom_vline(xintercept = (1 - max_missing) * n_S1_expected, linetype = "dashed", colour = "red") +
  labs(title = "Number of CpGs contributing to MRS1 per sample", x = "Number of non-missing CpGs", y = "Number of samples") + theme_bw(base_size = 11)
ggsave(filename = file.path(results_folder, "contributing_cpgs_MRS1.png"), plot = plot_contrib_S1,
       width = 8, height = 5, units = "in", dpi = 300)

## S2
n_cpgs_S2 <- count_contributing_cpgs(Betavalues_valid, names(weights_2))
valid.cohort$n_cpgs_MRS2 <- n_cpgs_S2[rownames(valid.cohort)]
n_S2_expected <- length(weights_2)
print(n_S2_expected)
low_coverage  <- sum(n_cpgs_S2 < (1 - max_missing) * n_S2_expected, na.rm = TRUE)
cat("\nSamples with <", (1 - max_missing) * 100, "% of MRS2 CpGs available:", low_coverage, "\n")
# Distribution plot of contributing CpGs for MRS2
cpg_contrib_df <- data.frame(sample = names(n_cpgs_S2), n_contrib = as.numeric(n_cpgs_S2), expected = n_S2_expected)
plot_contrib_S2 <- ggplot(cpg_contrib_df, aes(x = n_contrib)) + geom_histogram(bins = 20, fill = "#2166AC", alpha = 0.7) +
  geom_vline(xintercept = (1 - max_missing) * n_S2_expected, linetype = "dashed", colour = "red") +
  labs(title = "Number of CpGs contributing to MRS2 per sample", x = "Number of non-missing CpGs", y = "Number of samples") + theme_bw(base_size = 11)
ggsave(filename = file.path(results_folder, "contributing_cpgs_MRS2.png"), plot = plot_contrib_S2,
       width = 8, height = 5, units = "in", dpi = 300)

## S3
n_cpgs_S3 <- count_contributing_cpgs(Betavalues_valid, names(weights_3))
valid.cohort$n_cpgs_MRS3 <- n_cpgs_S3[rownames(valid.cohort)]
n_S3_expected <- length(weights_3)
low_coverage  <- sum(n_cpgs_S3 < (1 - max_missing) * n_S3_expected, na.rm = TRUE)
cat("\nSamples with <", (1 - max_missing) * 100, "% of MRS3 CpGs available:", low_coverage, "\n")
# Distribution plot of contributing CpGs for MRS3
cpg_contrib_df <- data.frame(sample = names(n_cpgs_S3), n_contrib = as.numeric(n_cpgs_S3), expected = n_S3_expected)
plot_contrib_S3 <- ggplot(cpg_contrib_df, aes(x = n_contrib)) + geom_histogram(bins = 20, fill = "#2166AC", alpha = 0.7) +
  geom_vline(xintercept = (1 - max_missing) * n_S3_expected, linetype = "dashed", colour = "red") +
  labs(title = "Number of CpGs contributing to MRS3 per sample", x = "Number of non-missing CpGs", y = "Number of samples") + theme_bw(base_size = 11)
ggsave(filename = file.path(results_folder, "contributing_cpgs_MRS3.png"), plot = plot_contrib_S3,
       width = 8, height = 5, units = "in", dpi = 300)

# ## S4
# n_cpgs_S4 <- count_contributing_cpgs(Betavalues_valid, names(weights_4))
# valid.cohort$n_cpgs_MRS4 <- n_cpgs_S4[rownames(valid.cohort)]
# n_S4_expected <- length(weights_4)
# low_coverage  <- sum(n_cpgs_S4 < (1 - max_missing) * n_S4_expected, na.rm = TRUE)
# cat("\nSamples with <", (1 - max_missing) * 100, "% of MRS4 CpGs available:", low_coverage, "\n")
# # Distribution plot of contributing CpGs for MRS3
# cpg_contrib_df <- data.frame(sample = names(n_cpgs_S4), n_contrib = as.numeric(n_cpgs_S4), expected = n_S4_expected)
# plot_contrib_S4 <- ggplot(cpg_contrib_df, aes(x = n_contrib)) + geom_histogram(bins = 20, fill = "#2166AC", alpha = 0.7) +
#   geom_vline(xintercept = (1 - max_missing) * n_S3_expected, linetype = "dashed", colour = "red") +
#   labs(title = "Number of CpGs contributing to MRS4 per sample", x = "Number of non-missing CpGs", y = "Number of samples") + theme_bw(base_size = 11)
# ggsave(filename = file.path(results_folder, "contributing_cpgs_MRS4.png"), plot = plot_contrib_S4,
#        width = 8, height = 5, units = "in", dpi = 300)

##########################################

cat("Tracking concluded. \n\nStarting MRS standardization. \n")
message("Tracking concluded. \nStarting MRS standardization: ", Sys.time())


## Scale MRS (Z-score MRS columns)

# Transforming the score to have mean = 0 and SD = 1
# (For each MRS column, subtracting the validation cohort mean and dividing by the validation cohort SD)
# Important step, in order to use it in Cox model, so the HR is expressed as "per one standard deviation increase in MRS"

mrs_cols <- c("MRS1", "MRS2", "MRS3") # "MRS4" excluded for not having CpGs
mrs_cols <- intersect(mrs_cols, colnames(valid.cohort))
mrs_z_cols <- paste0(mrs_cols, "_z")

# Scale and save new value as a new column MRS_z
valid.cohort[mrs_z_cols] <- lapply(valid.cohort[mrs_cols], scale)

# Check Z-score scaling
round(colMeans(valid.cohort[mrs_z_cols], na.rm = TRUE), 6)  # should all be ~0
round(apply(valid.cohort[mrs_z_cols], 2, sd, na.rm = TRUE), 6)  # should all be ~1

##########################################

cat("Pipeline completed. \n\nSaving outputs. \n")
message("Pipeline completed. \n\nSaving outputs: ", Sys.time())

#### SAVING OUTPUTS ####

# Save 
write.csv(valid.cohort, file = file.path(results_folder, "validation_cohort_MRS.csv"), row.names = FALSE)
save(mrs_cols, mrs_z_cols, weights_1, weights_2, weights_3, weights_4, file = file.path(results_folder, "MRS_names_and_weights.R"))

cat("Script completed. \n")
message("Script completed: ", Sys.time())

##########################################################
