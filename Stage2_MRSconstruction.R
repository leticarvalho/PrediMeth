##########################################################
#
# MRS construction and calculation - Stage 2
#
##########################################################


# DESCRIPTION --------------------------------------------
# In this script, CpGs are selected to construct different 
# Methylation Risk Scores (MRS), which will be later 
# calculated for each sample from tha validation cohort. 

# INPUT --------------------------------------------------
#   From discovery cohort: 
#       Pheno data
#       Mvalues
#       Limma's topTables (HL and ML)
#   From validation cohort: 
#       Pheno data
#       Betavalues

# OUTPUT -------------------------------------------------
# 	MRSs calculated for validation cohort


##########################################################

#### LOAD LIBRARIES ####

library(dplyr)
library(data.table)
library(minfi)

#### RESOLVE PATHS ####

## Input paths -------------------------------------------

predimeth_path <- "project/path"
results_dir <- file.path(predimeth_path, "results")

## Output paths ------------------------------------------

results_folder <- file.path(results_dir, 'MRS')
dir.create(results_folder)

##########################################################

#### LOAD DATA ####

## Discovery cohort --------------------------------------

## Pheno data (as data frame) = pheno_disc

## Mvalues (as matrix) = Mvalues_disc

## Limma's topTables = High_vs_Low_topTable, Medium_vs_Low_topTable

## Validation cohort -------------------------------------

## Pheno data (as data frame) = valid.cohort

## Betavalues (as matrix) = Betavalues_valid

##########################################################

cat("Data loaded. \n\nStarting section 1. \n")
message("Data loaded. \n\nStarting section 1: ", Sys.time())

#### PREPARING DATA ####


## Checking discovery cohort pheno data and Mvalues matrix

# Dimensions
cat("Dimensions of pheno discovery table: ", dim(pheno_disc), "\n")
cat("Dimensions of Mvalues discovery samples: ", dim(Mvalues_disc), "\n")

# Rownames
rownames(pheno_disc) <- pheno_disc$sample_id

# Checking alignment
cat("Checking alignment of pheno samples (rows) and Mvalues samples (columns): ", 
    "\n First pheno samples: ", rownames(pheno_disc)[1:5], "\n",
    "\n First Mvalues samples: ", colnames(Mvalues_disc)[1:5], "\n")
stopifnot(all(rownames(pheno_disc) %in% colnames(Mvalues_disc)))
stopifnot(identical(colnames(Mvalues_disc), rownames(pheno_disc)))
cat("Sample alignment verified.\n")


## Checking topTables

# High vs Low
High_vs_Low_topTable <- as.data.frame(High_vs_Low_topTable)
cat("Dimensions HL topTable: ", dim(High_vs_Low_topTable), "\n")
rownames(High_vs_Low_topTable) <- High_vs_Low_topTable$Probe_name 

# Medium vs Low
Medium_vs_Low_topTable <- as.data.frame(Medium_vs_Low_topTable)
cat("Dimensions HL topTable: ", dim(Medium_vs_Low_topTable), "\n")
rownames(Medium_vs_Low_topTable) <- Medium_vs_Low_topTable$Probe_name 


## Subsets and hits

# High vs Low
fdr_HL <- High_vs_Low_topTable[, High_vs_Low_topTable$adj.P.val < 0.05]
rownames(fdr_HL) <- fdr_HL$Probe_name
cpg_HL <- fdr_HL$Probe_name

# Medium vs Low
fdr_ML <- Medium_vs_Low_topTable[, High_vs_Low_topTable$adj.P.val < 0.05]
rownames(fdr_ML) <- fdr_ML$Probe_name
cpg_ML <- fdr_ML$Probe_name 

# Shared hits
cpg_shared <- intersect(fdr_HL$Probe_name, fdr_ML$Probe_name) 

# Summarizing information
cat("\nShared hits between HL and ML:", length(cpg_shared), "\n")
cat("High-vs-Low hits :", length(cpg_HL),"\n")
cat("Medium-vs-Low hits :", length(cpg_ML),"\n")


##########################################################


cat("Data ready. \n\nStarting preprocessing of Betavalues from validation cohort. \n")
message("Data ready. \nStarting preprocessing of Betavalues from validation cohort: ", Sys.time())

## Preprocessing of Betavalues from validation cohort

## Checking rownames and colnames
rownames(Betavalues_valid)[1:5]
colnames(Betavalues_valid)[1:5]

## Filtering valid.cohort samples
Betavalues_valid <- subset(Betavalues_valid, select = c(colnames(Betavalues_valid) %in% rownames(valid.cohort)))

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
cat("Samples in pheno after subsetting:", nrow(valid.cohort), "\n")
rownames(valid.cohort)[1:5]
cat("Samples in beta_matrix after reordering:", ncol(Betavalues_valid), "\n")
colnames(Betavalues_valid)[1:5]
stopifnot(all(rownames(valid.cohort) %in% colnames(Betavalues_valid)))
stopifnot(identical(colnames(Betavalues_valid), rownames(valid.cohort)))
cat("Sample alignment verified.\n")

## Clamping beta-values = 0 and beta-values = 1 (so they don't convert into infinite M-values)
## (same approach as in discovery cohort)

cat("Beta clamping \n")
cat("beta = 0 cells:", sum(Betavalues_valid == 0, na.rm = TRUE), "\n")
cat("beta = 1 cells:", sum(Betavalues_valid == 1, na.rm = TRUE), "\n")

epsilon <- 0.001 # standard practice (Du et al. 2010; Maksimovic et al. 2017)
Betavalues_valid[Betavalues_valid < epsilon] <- epsilon
Betavalues_valid[Betavalues_valid > 1-epsilon] <- 1 - epsilon
n_low <- sum(Betavalues_valid < epsilon, na.rm = TRUE)                                                                                                                                                        
n_high <- sum(Betavalues_valid > 1 - epsilon, na.rm = TRUE)                                                                                                                                                        
cat("Clamped", n_low, "cells to", epsilon, "and", n_high, "cells to", 1 - epsilon, "(out of", length(Betavalues_valid), "total)\n")        

#### CONVERTING BETAVALUES TO M-VALUES ####

Mvalues_valid <- log2(Betavalues_valid / (1 - Betavalues_valid))
cat("M-value matrix:", nrow(Mvalues_valid), "probes x", ncol(Mvalues_valid), "samples\n\n")

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

## If n_inf is 0, everything is fine: clamping worked, and the NAs are expected 
## missing values from probe-level QC failures in the validation cohort. These 
## will be handled correctly by the na.rm = TRUE in the MRS_calculator function, 
## which skips missing values rather than propagating them into the score.

## Checking rownames 
cat("Checking methylation matrix rownames: \n")
Betavalues_valid[1:5, 1:5]
Mvalues_valid[1:5, 1:5]

##########################################

cat("Preprocessing concluded. \n\nStarting MRS construction. \n")
message("Preprocessing concluded. \nStarting MRS construction: ", Sys.time())

## MRS Construction

# MRS calculator (computes a weighted MRS for each sample)
# Parameters: 
#   Mvalues_valid = Mvalue matrix from valid cohort (CpG x samples)
#   cpgs = selected cpgs to form MRS
#   weights named numeric vector of weights (names = CpG names)
# Return: named numeric vector (one MRS per sample)
MRS_calculator <- function(Mvalues_valid, cpgs, weights) {
  cpgs_m <- intersect(cpgs, rownames(Mvalues_valid))
  missing <- length(cpgs) - length(cpgs_m)
  if(missing > 0) {warning(missing, " CpGs absent from Mvalue matrix and excluded.")}
  w <- weights[cpgs_m] 
  b <- Mvalues_valid[cpgs_m, , drop = FALSE] # prevents conversion to vector when single CpG
  # na.rm = TRUE: missing CpG values are skipped rather than propagating NA
  # This means each sample's score is based on however many CpGs are available
  sapply(colnames(b), function(s) sum(w * b[, s], na.rm = TRUE))
}

# Strategy 1 (S1): High vs Low HITS (with HL logFC as weights)

weights_S1 <- setNames(fdr_HL[cpg_HL, "logFC"], # vector of logFCs
                       cpg_HL) # names of the vector

# Strategy 2 (S2): Medium vs Low HITS (with ML logFC as weights)

weights_S2 <- setNames(fdr_ML[cpg_ML, "logFC"], # vector of logFCs
                       cpg_ML) # names of the vector

# Strategy 3 (S3) P-VALUE THRESHOLDS FROM HIGH VS LOW TOPTABLE

p_thresholds <- c(1e-7, 1e-6, 1e-5, 1e-4)

weights_S3_list <- lapply(p_thresholds, function(pt) {
  cpgs_pt <- rownames(High_vs_Low_topTable[High_vs_Low_topTable$P.Value < pt, ])
  setNames(High_vs_Low_topTable[cpgs_pt, "logFC"], cpgs_pt)
})
names(weights_S3_list) <- paste0("MRS3_p", gsub("-", "_", as.character(p_thresholds)))

cat("\nStrategy 3: CpG counts per p-value threshold:\n")
print(sapply(weights_S3_list, length))

message("Starting elastic net: ", Sys.time())

# 4) ELASTIC NET (using Mvalues from the discovery cohort)

# Preparing predimed_categ variable / predimed binary
predimed_type <- factor(pheno_disc$predimed_high, levels = c(0, 1), labels = c("no", "yes"))
#predimed_categ <- factor(pheno_disc$predimed_cat, levels = c(1, 2, 3), labels = c("low", "medium", "high"))
outcome_disc <- as.numeric(predimed_type) - 1L # because cv.glmnet expects a numeric 0/1/2 vector
#outcome_disc <- as.numeric(predimed_categ) - 1L # because cv.glmnet expects a numeric 0/1/2 vector
table(pheno_disc$predimed_high, outcome_disc) # check conversion
#table(pheno_disc$predimed_cat, outcome_disc) # check conversion
stopifnot(length(outcome_disc) == ncol(Mvalues_disc)) # check dimensions

cat("\nTraining outcome distribution (0=low, 1=medium, 2=high):\n")
print(table(outcome_disc))

# Pre-filter by variance across samples -> reduces computation without biasing toward EWAS result
cpg_sd <- apply(Mvalues_disc, 1, sd, na.rm = TRUE)
cpg_variable <- names(cpg_sd)[cpg_sd > quantile(cpg_sd, 0.75)] # keeps the top ~ 25% most variable CpGs
cat("CpGs passing variance filter:", length(cpg_variable), "\n")

# transpose to rows = samples, cols = CpGs, as glmnet expects
X_train <- t(Mvalues_disc[cpg_variable, ])
#X_train <- t(Mvalues_disc) 

set.seed(123)
cv_enet <- cv.glmnet(x=X_train, y=outcome_disc, family = "gaussian", alpha = 0.5, nfolds = 10, type.measure = "mse")
plot(cv_enet, main = "Elastic net CV: Strategy 4 (all EWAS-tested CpGs)")

coef_enet <- coef(cv_enet, s = cv_enet$lambda.1se)
cpg_enet <- rownames(coef_enet)[coef_enet[, 1] != 0 & rownames(coef_enet) != "(Intercept)"]

weights_S4 <- setNames(coef_enet[cpg_enet, 1], cpg_enet)

cat("\nElastic net selected", length(cpg_enet), "CpGs from",
    ncol(X_train), "candidates\n")

# how many selected CpGs overlap with EWAS hits
cat("Overlap with HL hits:", sum(cpg_enet %in% rownames(fdr_HL)), "\n") 
cat("Overlap with ML hits:", sum(cpg_enet %in% rownames(fdr_ML)), "\n")
cat("Novel CpGs (not in any hit list):", sum(!cpg_enet %in% union(rownames(fdr_HL), rownames(fdr_ML))), "\n")

##########################################

## Validation cohort processing - Update

# Checking rownames and colnames
cat("Rownames phenotype table (validation):", rownames(valid.cohort)[1:5], "\n")
cat("Colnames Mvalues table (validation):", colnames(Mvalues_valid)[1:5], "\n")

# Checking dimensions
cat("Dimensions phenotype table (validation):", dim(valid.cohort), "\n")
cat("Dimensions Mvalues table (validation):", dim(Mvalues_valid), "\n")

common_samples <- intersect(rownames(valid.cohort), colnames(Mvalues_valid))

cat("Samples in valid.cohort :", nrow(valid.cohort), "\n")
cat("Samples in Mvalues_valid :", ncol(Mvalues_valid), "\n")
cat("Common samples (will be used) :", length(common_samples),"\n")
cat("Dropped from valid.cohort :", nrow(valid.cohort) - length(common_samples), "\n")

# Subsetting validation cohort (excluding 7 samples that are not in Mvalues_valid)
valid.cohort <- as.data.frame(valid.cohort)
rownames(valid.cohort) <- valid.cohort$sample_id
print(rownames(valid.cohort)[!(rownames(valid.cohort) %in% common_samples)]) ## NEW
valid.cohort  <- valid.cohort[common_samples, , drop=FALSE]
cat("Samples in valid.cohort after subsetting:", nrow(valid.cohort), "\n")

# Check the first few of each side
cat("First 5 rownames valid.cohort :\n")
print(rownames(valid.cohort)[1:5])
cat("First 5 colnames Mvalues_valid:\n")
print(colnames(Mvalues_valid)[1:5])

# Check if they are the same set but in different order
cat("Same elements, different order?", 
    setequal(rownames(valid.cohort), colnames(Mvalues_valid)), "\n")

# Reordering Mvalues columns to match phenotype table row order
Mvalues_valid <- Mvalues_valid[, rownames(valid.cohort), drop = FALSE]
cat("Samples in Mvalues_valid after reordering:", ncol(Mvalues_valid), "\n")

# Final check
stopifnot(identical(rownames(valid.cohort), colnames(Mvalues_valid)))
cat("Sample alignment verified.\n")

# Checking dimensions
cat("Dimensions phenotype table (validation) after alignment:", dim(valid.cohort), "\n")
cat("Dimensions Mvalues table (validation) after alignment:", dim(Mvalues_valid), "\n")

##########################################

cat("MRS construction completed. \n\nStarting MRS calculation for the validation cohort. \n")
message("MRS construction completed. \nStarting MRS calculation for the validation cohort: ", Sys.time())

## Calculate MRSs from validation cohort (apply all weight sets and add to dataframe)

valid.cohort$MRS1 <- MRS_calculator(Mvalues_valid, names(weights_S1), weights_S1)
valid.cohort$MRS2 <- MRS_calculator(Mvalues_valid, names(weights_S2), weights_S2)

for (nm in names(weights_S3_list)) {
  w <- weights_S3_list[[nm]]
  valid.cohort[[nm]] <- MRS_calculator(Mvalues_valid, names(w), w)
}

#~ valid.cohort$MRS4 <- MRS_calculator(Mvalues_valid, names(weights_S4), weights_S4) # no CpGs selected by EN

##########################################

cat("Calculation concluded. \n\nStarting tracking of non-missing CpGs. \n")
message("Calculation concluded. \nStarting tracking of non-missing CpGs: ", Sys.time())

## Track number of non-missing CpGs contributing to each MRS per sample

# Following the approach of Davyson et al. (AD MRS script - https://github.com/Elladavyson/AD_MRS/):

max_missing <- 0.10  # flag samples missing > 10% of MRS CpGs

count_contributing_cpgs <- function(Mvalues_valid, cpgs) {
  cpgs_m <- intersect(cpgs, rownames(Mvalues_valid))
  mval_cpgs <- Mvalues_valid[cpgs_m, , drop = FALSE]
  count_used_cpgs <- colSums(!is.na(mval_cpgs)) # for each sample, count how many CpGs are non-missing
  return(count_used_cpgs)
}

# Count for each strategy and add it in the pheno df
# Also report samples where MRS is based on < 90% of its CpGs

## S1
n_cpgs_S1 <- count_contributing_cpgs(Mvalues_valid, names(weights_S1))
valid.cohort$n_cpgs_MRS1 <- n_cpgs_S1[rownames(valid.cohort)]
n_S1_expected <- length(weights_S1)
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
n_cpgs_S2 <- count_contributing_cpgs(Mvalues_valid, names(weights_S2))
valid.cohort$n_cpgs_MRS2 <- n_cpgs_S2[rownames(valid.cohort)]
n_S2_expected <- length(weights_S2)
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
plots_contrib_S3 <- list()  
for (mrs in names(weights_S3_list)) {
  w <- weights_S3_list[[mrs]]
  n_col <- paste0("n_cpgs_", mrs)
  counts <- count_contributing_cpgs(Mvalues_valid, names(w))
  valid.cohort[[n_col]] <- counts[rownames(valid.cohort)]
  n_expected <- length(w)
  low_coverage <- sum(counts < (1 - max_missing) * n_expected, na.rm = TRUE)
  percent_low <- round(low_coverage / ncol(Mvalues_valid) * 100, 1)
  cat("\nSamples with <", (1 - max_missing) * 100, "% of", mrs, "CpGs available:", low_coverage, "(", percent_low, "% of cohort)\n")
  # Distribution plot 
  cpg_contrib_df <- data.frame(sample = names(counts), n_contrib = as.numeric(counts))
  plots_contrib_S3[[mrs]] <- ggplot(cpg_contrib_df, aes(x = n_contrib)) +  
    geom_histogram(bins = 20, fill = "#1a9641", alpha = 0.7) +
    geom_vline(xintercept = (1 - max_missing) * n_expected, linetype = "dashed", colour = "red") +
    labs(title = paste0(mrs, ": contributing CpGs per sample (expected: ", n_expected, ")"),
         x = "Number of non-missing CpGs", y = "Number of samples") + theme_bw(base_size = 11)
}
n_s3_plots <- length(plots_contrib_S3)
png(filename = file.path(results_folder, "contributing_cpgs_MRS3.png"), 
    width = 2 * 5, height = ceiling(n_s3_plots / 2) * 4, units = "in", res = 300)
gridExtra::grid.arrange(grobs = plots_contrib_S3, ncol = 2, top = "S3: non-missing CpGs contributing to each threshold MRS")
dev.off()

#~ ## S4 # No CpGs selected
#~ n_cpgs_S4 <- count_contributing_cpgs(Mvalues_valid, names(weights_S4))
#~ valid.cohort$n_cpgs_MRS4 <- n_cpgs_S4[rownames(valid.cohort)]
#~ n_S4_expected <- length(weights_S4)
#~ low_coverage  <- sum(n_cpgs_S4 < (1 - max_missing) * n_S4_expected, na.rm = TRUE)
#~ cat("\nSamples with <", (1 - max_missing) * 100, "% of MRS4 CpGs available:", low_coverage, "\n")
#~ # Distribution plot of contributing CpGs for MRS4
#~ cpg_contrib_df <- data.frame(sample = names(n_cpgs_S4), n_contrib = as.numeric(n_cpgs_S4), expected = n_S4_expected)
#~ plot_contrib_S4 <- ggplot(cpg_contrib_df, aes(x = n_contrib)) + geom_histogram(bins = 20, fill = "#2166AC", alpha = 0.7) +
#~   geom_vline(xintercept = (1 - max_missing) * n_S4_expected, linetype = "dashed", colour = "red") +
#~   labs(title = "Number of CpGs contributing to MRS4 per sample", x = "Number of non-missing CpGs", y = "Number of samples") + theme_bw(base_size = 11)
#~ ggsave(filename = file.path(results_folder, "contributing_cpgs_MRS4.png"), plot = plot_contrib_S4,
#~        width = 8, height = 5, units = "in", dpi = 300)


##########################################

cat("Tracking concluded. \n\nStarting MRS standardization. \n")
message("Tracking concluded. \nStarting MRS standardization: ", Sys.time())

## Scale MRS (Z-score MRS columns)

# Transforming the score to have mean = 0 and SD = 1
# (For each MRS column, subtract the validation cohort mean and divide by the validation cohort SD)

# to use it in a logistic regression or Cox model, so 
# the OR and HR are expressed as "per one standard deviation increase in MRS"

mrs_cols <- c("MRS1", "MRS2", names(weights_S3_list))
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
write.table(valid.cohort, file = file.path(results_folder, "validation_cohort_MRS.txt"), sep = "\t", row.names = FALSE, col.names = TRUE)

cat("Script completed. \n")
message("Script completed: ", Sys.time())

##########################################################