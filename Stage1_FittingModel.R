##########################################################
#
# Fitting Model - Stage 1
#
##########################################################


# DESCRIPTION --------------------------------------------
# In this script, the processed data is used to fit a 
# linear model to each CpG, using limma R package. 
# Model: Mvalues ~ MedDiet adherence + age + sex + BMI smoking + batch + cells
# Distinct exposure variables to represent MedDiet adherence: 
#       Continuous = predimed_score (0-14)
#       Binary = predimed_high (high >= 9) - 553 (no) vs 421 (yes) 
#       3-factor categorical = predimed_cat (low <=4, high >=9): 38 (low) vs 515 (medium) vs 421 (high)
#       Extreme = predimed_extreme (Low_Q1 <= 6, Mid_Q23 = 7-9, High_Q4 >= 10): 235 (Low_Q1) vs 262 (High_Q4)
# For each model, a topTable will be generated, including 
# EPIC v2 annotation and calculated absolute delta beta (ADB). 
# Top hits will be FDR-adjusted, but top hits with significant
# nominal p-value will also be analysed.      

# INPUT --------------------------------------------------
# 	Data frame of phenotype (rows: samples, columns: variables)
# 	Matrix of beta values (rows: probes, columns: samples)
#	  Matrix of M-values (rows: probes, columns: samples)

# OUTPUT -------------------------------------------------
# 	For each contrast, a folder with: 
#       Plot with beta values distribution (bimodal) 
#       Results tables (topTable + annotation + ADB) (.Rdata)
#       CpGs vectors to enrichment (.Rdata)
#       Tables with fdr hits (.csv)
#       Volcano plot
#       Manhattan plot
#       Genomic Categories barplot
#   Lambdas comparison



##########################################################


#### LOAD LIBRARIES ####

library(dplyr)
library(data.table)
library(limma)
library(EnhancedVolcano)
library(ggplot2)
library(qqman)
library(RColorBrewer)

#### RESOLVE PATHS ####

## Input paths -------------------------------------------

path_to_EPIC2_annotation <- "/imppc/labs/dnalab/share/PrediMeth/0-data/EPICv2.annot.RData" # .Rdata
path_to_processed_data_stage1 <- "/imppc/labs/dnalab/share/PrediMeth/results/OfficialAnalysis/processed_data.R" # 

## Output paths -------------------------------------------

predimeth_path <- "/imppc/labs/dnalab/share/PrediMeth"
results_dir <- file.path(predimeth_path, "results")
results_folder <- file.path(results_dir, 'OfficialAnalysis')

##########################################################

#### LOAD DATA ####

## Preprocessed data (Stage 1) ---------------------------
load(path_to_processed_data_stage1) # pheno, beta_matrix, mval_matrix

## Annotation (EPIC v2) ----------------------------------
load(path_to_EPIC2_annotation) # EPIC_MEsteller
EPIC_MEsteller$Probe_name  <- stringr::str_remove(EPIC_MEsteller$ProbeID, "_.*")

##########################################################

cat("Data loaded. \nChecking information: \n")
message("Data loaded. Checking information: ", Sys.time())

#### CHECKING DATA STRUCTURE ####

## Renaming objects
Mvalues <- mval_matrix # rownames = Probe_name, colnames = IJC samples
Betavalues <- beta_matrix # rownames = Probe_name, colnames = IJC samples
pheno <- pheno
rownames(pheno) <- pheno$sample_id

# Checking alignment
print("Checking alignment: \n")
rownames(pheno)[1:5]
colnames(Mvalues)[1:5]
colnames(Betavalues)[1:5]

stopifnot(all(rownames(pheno) %in% colnames(Mvalues)))
stopifnot(all(rownames(pheno) %in% colnames(Betavalues)))
stopifnot(identical(colnames(Mvalues), rownames(pheno)))
stopifnot(identical(colnames(Betavalues), rownames(pheno)))
cat("Sample alignment verified.\n")

# Checking probe names order
print("Checking probe names order: \n")
colnames(Mvalues)[1:5]
colnames(Betavalues)[1:5]

## Checking missing values in Mvalues
cat("NA cells in Mvalues: ", sum(is.na(Mvalues)), "\n")  
rows_with_NA <- rowSums(is.na(Mvalues)) > 0                                                                                                                                                                   
cat("Probes with any NA:", sum(rows_with_NA), "/", nrow(Mvalues), "\n")  

## Checking pheno structure
str(pheno)

## Checking distribution of beta values in discovery cohort

cases <- pheno$sample_id[pheno$predimed_high == "yes"]
controls <- pheno$sample_id[pheno$predimed_high == "no"]
beta_matrix_noID <- beta_matrix
colnames(beta_matrix_noID)[colnames(beta_matrix_noID) %in% cases] <- "CASE"
colnames(beta_matrix_noID)[colnames(beta_matrix_noID) %in% controls] <- "CONTROL"
case_cols <- which(colnames(beta_matrix_noID) == "CASE")
control_cols <- which(colnames(beta_matrix_noID) == "CONTROL")

png(filename=file.path(results_folder, "methylation_profile_discovery_cohort.png"), width=1600, height=950)

plot(density(beta_matrix_noID[, case_cols[1]]),
     col="yellowgreen",
     xlab="Beta value",
     ylim=c(0,6),
     main="PrediMeth samples",
     lwd=2)

invisible(sapply(case_cols[-1], function(x)
  lines(density(beta_matrix_noID[,x]), col="yellowgreen", lwd=1)
))

invisible(sapply(control_cols, function(x)
  lines(density(beta_matrix_noID[,x]), col="cornflowerblue", lwd=1)
))

legend("topleft",
       c("CASES", "CONTROLS"),
       text.col=c("yellowgreen", "cornflowerblue"))

dev.off()

rm(beta_matrix_noID)

##########################################################

cat("Data checked. \n\nPreparing functions. \n")
message("Data checked. \n\nPreparing functions: ", Sys.time())

#### PREPARING FUNCTIONS ####

## To evaluate p-value inflation
measureInflation <- function(topTable, label, output_folder) {
  pvals  <- topTable$P.Value
  pvals  <- pvals[!is.na(pvals)]                                                                                                                                                                              
  chisq  <- qchisq(1 - pvals, df = 1)
  lambda <- median(chisq) / qchisq(0.5, df = 1)                                                                                                                                                               
  
  observed <- -log10(sort(pvals))                                                                                                                                                                             
  expected <- -log10(ppoints(length(pvals)))
  
  # QQ plot
  pdf(file = file.path(output_folder, paste0("pvalue_inflation_", label, ".pdf")), width  = 6, height = 6)                                                                                                                                                                                
  plot(expected, observed, pch = 20, cex = 0.6, xlab = "Expected -log10(p)", ylab = "Observed -log10(p)", 
       main = sprintf("QQ plot - %s Model \n(lambda = %.3f, N = %d probes)", label, lambda, length(pvals)))                                                                                                                                                                                     
  abline(0, 1, col = "red")
  dev.off()                                                                                                                                                                                                   
  
  invisible(lambda)                                                                                                                                                                                           
}

## To add delta beta in continuous model (predimed_score)
addDeltaBeta_Cont <- function(fit_beta, sd_score, topTable) {
  d_unit <- fit_beta$coefficients[, "predimed_score"] # expected change in beta per 1-point increase in predimed_score
  idx <- match(topTable$Probe_name, names(d_unit)) 
  topTable$delta_beta_per_unit <- d_unit[idx] 
  topTable$delta_beta_per_SD <- topTable$delta_beta_per_unit * sd_score 
  topTable$abs_delta_beta_per_SD <- abs(topTable$delta_beta_per_SD) 
  topTable                                                                                                                                                                                                    
} 

## To add delta beta in category models (predimed_score)
addDeltaBeta_Cat <- function(topTable, delta_beta) {
  idx <- match(topTable$Probe_name, names(delta_beta))
  topTable$delta_beta <- delta_beta[idx]
  topTable$abs_delta_beta <- abs(topTable$delta_beta)
  topTable
}

## To create a Volcano plot
volcano_plot <- function(topTable, contrast, output_folder) {
  plot <- EnhancedVolcano(topTable, x = "logFC", y = "adj.P.Val",
                          lab = "", pCutoff = 0.05, FCcutoff = log2(1.2),
                          pointSize = 3) + 
    geom_hline(yintercept = -log10(bonf_threshold), linetype = "dashed", color = "coral", linewidth = 0.5) + 
    ggplot2::labs(title = paste0("Contrast: ", contrast))
  ggsave(file.path(output_folder, paste0("Volcano_Plot_", contrast, ".pdf")), plot = plot, width = 8, height = 7)
}

## To create a Manhattan plot
MAplot <- function(topTable, contrast) {
  ma_data <- subset(topTable, select = c("CpGchrm", "CpGbeg", "adj.P.Val", "Probe_name"))
  ma_data$CpGchrm <- stringr::str_remove(ma_data$CpGchrm, "^chr")
  ma_data <- ma_data[ma_data$CpGchrm != "M",]
  table(ma_data$CpGchrm)
  ma_data$CpGchrm <- as.numeric(ma_data$CpGchrm)
  ma_data <- na.omit(ma_data)
  
  ma_plot <- manhattan(ma_data, chr = "CpGchrm", bp = "CpGbeg", 
                       p = "adj.P.Val", snp = "Probe_name",
                       genomewideline = -log10(5.722212e-08),
                       suggestiveline = -log10(0.05),
                       col = c("steelblue", "coral"), cex = 0.6,
                       main = paste0("EWAS Manhattan Plot - ", contrast))
  
  return(ma_plot)
}

## To create a bar plot to show genomic categories distribution
genomic_categories <- function(signif_subset, contrast, significance_type, output_folder) {
  df_summary <- signif_subset %>% group_by(RelationToGeneCategory) %>% summarise(count = n()) %>% arrange(desc(count))
  plot <- ggplot(df_summary, aes(x = reorder(RelationToGeneCategory, -count), y = count)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    theme_bw() +
    labs(x = "Genomic Category",
         y = "Number of CpGs",
         title = paste0(contrast, " ", significance_type, " significant CpGs by Genomic Category")) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(file.path(output_folder, paste0("GenomicCategoriesPlot_", contrast, ".pdf")), plot = plot)
}

#### BONFERRONI THRESHOLD ####

bonf_threshold <- 0.05 / dim(Mvalues)[1] # calculating Bonferroni's threshold
cat("Bonferroni threshold: ", bonf_threshold, "\n")

#### VECTOR WITH ALL CPGS IN THE ANALYSIS ####

all_probes_names <- rownames(Mvalues)
cat("Number of all CpGs used in the analysis: ", length(all_probes_names), "\n")
all_probes_ids <- EPIC_MEsteller[match(rownames(Mvalues), EPIC_MEsteller$Probe_name), ]$ProbeID

#### ALIGN ANNOTATION ####

genelist_aligned <- EPIC_MEsteller[match(rownames(Mvalues), EPIC_MEsteller$Probe_name), ]                                                                                                                                                      
stopifnot(identical(genelist_aligned$Probe_name, rownames(Mvalues)))

##########################################################

cat("Variables ready. \n\n Starting continuous models. \n")
message("Variables ready. \n\n Starting continuous models: ", Sys.time())

#### CONTINUOUS MODEL ####

## Output folder -----------------------------
continuous_folder <- file.path(results_folder, 'ContinuousContrast')
dir.create(continuous_folder, recursive = TRUE, showWarnings = FALSE)

## Check distribution of exposure variable -----------------
cat("Distribution of exposure variable: predimed_score")
table(pheno$predimed_score)

## Model designs -----------------------------

dm1 <- model.matrix(~ predimed_score + age + sex + bmi + smoking_type + batch + Bcell + CD4T + CD8T + Mono + NK, data = pheno) 
cat("\n Head of design matrix 1: \n")
head(dm1)

# Checking missing values in design matrix  -----------------------------
cat("NA in design matrix 1: ", sum(is.na(dm1)), "\n") # MUST be 0   

# Checking dimensions  -----------------------------
stopifnot(ncol(Mvalues) == nrow(dm1))                                                                                                                                                                         
stopifnot(identical(colnames(Mvalues), rownames(dm1))) 

## Model fitting -----------------------------

fit1 <- lmFit(Mvalues, dm1, na.action = na.exclude)
fit1 <- eBayes(fit1, trend = TRUE) ## Applying eBayes (trend = TRUE models mean-variance relationship)
colnames(fit1$coefficients)

## TopTable -----------------------------
full_topTable1 <- limma::topTable(fit1, number = Inf, adjust.method = "BH", coef = "predimed_score", genelist = genelist_aligned)

## Checking p-value inflation  -----------------------------
lambda1cont <- measureInflation(full_topTable1, "Continuous", continuous_folder)
cat("Lambda Model 1 (Full) = ", lambda1cont, "\n")

## Checking number of FDR-significant hits ------------- 
fdr_hits1 <- subset(full_topTable1, adj.P.Val<0.05)
cat("FDR hits - Model 1 (full) - Continuous: ", dim(fdr_hits1)[1], "\n")
bonf_hits1 <- subset(full_topTable1, adj.P.Val<bonf_threshold)
cat("Bonferroni hits - Full Model - Continuous: ", dim(bonf_hits1)[1], "\n")


## Calculating ADB (Absolute Delta Beta) and adding it to topTable -------------
fit_beta1 <- lmFit(Betavalues, dm1)
sd_score <- sd(pheno$predimed_score, na.rm = TRUE)  
full_topTable1 <- addDeltaBeta_Cont(fit_beta1, sd_score, full_topTable1)

## Updating subset of FDR-significant hits ------------- 
fdr_hits1 <- subset(full_topTable1, adj.P.Val<0.05)

## Subsetting sugestive hits (pvalue < 1e-4) ------------- 
suggestive_hits1 <- subset(full_topTable1, P.Value < 1e-4)
cat("CpGs with p-value < 1e-4 - Full model - Continuous: ", dim(suggestive_hits1)[1], "\n")

## Subsetting nominal hits (nominal pvalue < 0.05) ------------- 
nominal_hits1 <- subset(full_topTable1, P.Value < 0.05)
cat("CpGs with p-value < 0.05 - Full model - Continuous: ", dim(nominal_hits1)[1], "\n")

## CpGs (vectors) to posterior annotation
fdr_names1 <- fdr_hits1$ProbeID
suggestive_names1 <- suggestive_hits1$ProbeID
nominal_names1 <- nominal_hits1$ProbeID

## Saving outputs -----------------------------

cat("Saving outputs - Continuous: \n")
message("Saving outputs - Continuous: ", Sys.time())

save(full_topTable1, file = file.path(continuous_folder, "limma_results_continuous.R"))

save(fdr_names1, suggestive_names1, nominal_names1, all_probes_ids, file = file.path(continuous_folder, "hits_vectors_continuous.R"))

write.csv(fdr_hits1, file = file.path(continuous_folder, "fdr_hits_fullModel.csv"), row.names = FALSE)
write.csv(suggestive_hits1, file = file.path(continuous_folder, "suggestive_hits_fullModel.csv"), row.names = FALSE)
write.csv(nominal_hits1, file = file.path(continuous_folder, "nominal_hits_fullModel.csv"), row.names = FALSE)

## Volcano plot ---------------------------------
volcano_plot(full_topTable1, "continuous", continuous_folder)

## Manhattan plot ---------------------------------
png(file = file.path(continuous_folder, "MAplot_continuous.png"), width = 1600, height = 900, res = 150)
MAplot(full_topTable1, "Continuous") 
dev.off()

## Genomic categories barplot ---------------------
genomic_categories(suggestive_hits1, "Continuous", "p-value", continuous_folder)

##########################################################

cat("Continuous models pipeline completed. \n\n Starting binary models. \n")
message("Continuous models pipeline completed. \n\n Starting binary models: ", Sys.time())

#### BINARY MODEL ####

## Output folder -----------------------------
binary_folder <- file.path(results_folder, 'BinaryContrast')
dir.create(binary_folder, recursive = TRUE, showWarnings = FALSE)

## Check distribution of exposure variable -----------------
cat("Distribution of exposure variable: predimed_high")
table(pheno$predimed_high)

## Model designs -----------------------------

dm1 <- model.matrix(~0 + predimed_high + age + sex + bmi + smoking_type + batch + Bcell + CD4T + CD8T + Mono + NK, data = pheno) 
cat("\n Head of design matrix 1: \n")
head(dm1)

# Checking missing values in design matrix  -----------------------------
cat("NA in design matrix 1: ", sum(is.na(dm1)), "\n") # MUST be 0   

# Checking dimensions  -----------------------------
stopifnot(ncol(Mvalues) == nrow(dm1))                                                                                                                                                                         
stopifnot(identical(colnames(Mvalues), rownames(dm1))) 

## Contrast designs -----------------------------
cm1 <- makeContrasts(predimed_highyes_vs_predimed_highno = predimed_highyes-predimed_highno, levels = dm1)
contrast1 <- colnames(cm1) 
cat(" Name of contrast 1: \n", contrast1, "\n")

## Model fitting -----------------------------
fit1 <- lmFit(Mvalues, dm1, na.action = na.exclude)
fit1 <- contrasts.fit(fit1, cm1)
fit1 <- eBayes(fit1, trend = TRUE) 
colnames(fit1$coefficients)

## TopTable -----------------------------
full_topTable1 <- limma::topTable(fit1, number = Inf, adjust.method = "BH", coef = contrast1, genelist = genelist_aligned)

## Checking p-value inflation  -----------------------------
lambda1bin <- measureInflation(na.omit(full_topTable1), "Binary", binary_folder)
cat("Lambda Model 1 (Full) = ", lambda1bin, "\n")

## Calculating ADB (Absolute Delta Beta) and adding to topTable -------------
mean_cases <- rowMeans(Betavalues[, pheno$predimed_high == "yes"], na.rm = TRUE) # high diet adherence
mean_controls <- rowMeans(Betavalues[, pheno$predimed_high == "no"], na.rm=TRUE) # low diet adherence
delta_beta <- mean_cases - mean_controls
full_topTable1 <- addDeltaBeta_Cat(full_topTable1, delta_beta)

## Subset of FDR-significant hits ------------- 
fdr_hits1 <- subset(full_topTable1, adj.P.Val<0.05)
cat("FDR hits - Full Model - Binary: ", dim(fdr_hits1)[1], "\n")
bonf_hits1 <- subset(full_topTable1, adj.P.Val<bonf_threshold)
cat("Bonferroni hits - Full Model - Binary: ", dim(bonf_hits1)[1], "\n")

## Subsetting sugestive hits (pval < 1e-5) ------------- 
suggestive_hits1 <- subset(full_topTable1, P.Value < 1e-4)
cat("CpGs with p-value < 1e-4 - Full model - Binary: ", dim(suggestive_hits1)[1], "\n")

## Subsetting nominal hits (nominal pvalue < 0.05) ------------- 
nominal_hits1 <- subset(full_topTable1, P.Value < 0.05)
cat("CpGs with p-value < 0.05 - Full model - Binary: ", dim(nominal_hits1)[1], "\n")

## CpGs (vectors) to posterior annotation
fdr_names1 <- fdr_hits1$ProbeID
suggestive_names1 <- suggestive_hits1$ProbeID
nominal_names1 <- nominal_hits1$ProbeID

## Saving outputs -----------------------------

cat("Saving outputs - binary: \n")
message("Saving outputs - binary: ", Sys.time())

save(full_topTable1, file = file.path(binary_folder, "limma_results_binary.R"))

save(fdr_names1, suggestive_names1, nominal_names1, all_probes_ids, file = file.path(binary_folder, "hits_vectors_binary.R"))

write.csv(fdr_hits1, file = file.path(binary_folder, "fdr_hits_fullModel.csv"), row.names = FALSE)
write.csv(suggestive_hits1, file = file.path(binary_folder, "suggestive_hits_fullModel.csv"), row.names = FALSE)
write.csv(nominal_hits1, file = file.path(binary_folder, "nominal_hits_fullModel.csv"), row.names = FALSE)

## Volcano plot ---------------------------------
volcano_plot(full_topTable1, "binary", binary_folder)

## Manhattan plot ---------------------------------
png(file = file.path(binary_folder, "MAplot_binary.png"), width = 1600, height = 900, res = 150)
MAplot(full_topTable1, "Binary") 
dev.off()

## Genomic categories barplot ---------------------
genomic_categories(suggestive_hits1, "Binary", "p-value", binary_folder)

##########################################################

cat("Binary models pipeline completed. \n\n Starting Original 3-categorical models. \n")
message("Binary models pipeline completed. \n\n Starting Original 3-categorical models: ", Sys.time())

#### ORIGINAL 3-CATEGORICAL MODELS ####

## Output folder -----------------------------
cat3orig_folder <- file.path(results_folder, 'Original3CategoricalContrast')
dir.create(cat3orig_folder, recursive = TRUE, showWarnings = FALSE)

## Check distribution of exposure variable -----------------
cat("Distribution of exposure variable: predimed_cat")
table(pheno$predimed_cat)

## Model designs -----------------------------

dm1 <- model.matrix(~0 + predimed_cat + age + sex + bmi + smoking_type + batch + Bcell + CD4T + CD8T + Mono + NK, data = pheno) 
cat("\n Head of design matrix 1: \n")
head(dm1)

# Checking missing values in design matrix  -----------------------------
cat("NA in design matrix 1: ", sum(is.na(dm1)), "\n") # MUST be 0   

# Checking dimensions  -----------------------------
stopifnot(ncol(Mvalues) == nrow(dm1))                                                                                                                                                                         
stopifnot(identical(colnames(Mvalues), rownames(dm1))) 

## Contrast designs -----------------------------
cm1 <-makeContrasts(High_vs_Low = predimed_cathigh - predimed_catlow,
                    Medium_vs_Low = predimed_catmedium - predimed_catlow,
                    High_vs_Medium = predimed_cathigh - predimed_catmedium,
                    levels = dm1)
contrasts1 <- colnames(cm1) 
cat(" Names of contrasts 1: \n", contrasts1, "\n")

## Model fitting -----------------------------
fit1 <- lmFit(Mvalues, dm1, na.action = na.exclude)
fit1 <- contrasts.fit(fit1, cm1)
fit1 <- eBayes(fit1, trend = TRUE) 
colnames(fit1$coefficients)

## TopTable -----------------------------
High_vs_Low_topTable1 <- limma::topTable(fit1, number = Inf, adjust.method = "BH", coef = "High_vs_Low", genelist = genelist_aligned)
Medium_vs_Low_topTable1 <- limma::topTable(fit1, number = Inf, adjust.method = "BH", coef = "Medium_vs_Low", genelist = genelist_aligned)
High_vs_Medium_topTable1 <- limma::topTable(fit1, number = Inf, adjust.method = "BH", coef = "High_vs_Medium", genelist = genelist_aligned)

## Checking p-value inflation  -----------------------------
lambda1HL <- measureInflation(High_vs_Low_topTable1, "High_vs_Low", cat3orig_folder)
lambda1ML <- measureInflation(Medium_vs_Low_topTable1, "Medium_vs_Low", cat3orig_folder)
lambda1HM <- measureInflation(High_vs_Medium_topTable1, "High_vs_Medium", cat3orig_folder)
cat("Lambdas Model 1 (Full adjustment): ", 
    "\nHigh vs Low = ", lambda1HL, "\n",
    "\nMedium vs Low = ", lambda1ML, "\n",
    "\nHigh vs Medium = ", lambda1HM, "\n")

## Calculating ADB (Absolute Delta Beta) and adding to topTable -------------

mean_low <- rowMeans(Betavalues[, pheno$predimed_cat == "low"], na.rm=TRUE)
mean_medium <- rowMeans(Betavalues[, pheno$predimed_cat == "medium"], na.rm=TRUE)
mean_high <- rowMeans(Betavalues[, pheno$predimed_cat == "high"], na.rm = TRUE) 

delta_beta_HL <- mean_high - mean_low
High_vs_Low_topTable1 <- addDeltaBeta_Cat(High_vs_Low_topTable1, delta_beta_HL)

delta_beta_ML <- mean_medium - mean_low
Medium_vs_Low_topTable1 <- addDeltaBeta_Cat(Medium_vs_Low_topTable1, delta_beta_ML)

delta_beta_HM <- mean_high - mean_medium
High_vs_Medium_topTable1 <- addDeltaBeta_Cat(High_vs_Medium_topTable1, delta_beta_HM)

## Subset of FDR-significant hits ------------- 

fdr_HL1 <- subset(High_vs_Low_topTable1, adj.P.Val<0.05)
fdr_ML1 <- subset(Medium_vs_Low_topTable1, adj.P.Val<0.05)
fdr_HM1 <- subset(High_vs_Medium_topTable1, adj.P.Val<0.05)
cat("FDR hits - Model 1 (full) - 3 Categories (original): \n",
    "High_vs_Low: ", dim(fdr_HL1)[1]," hits \n",
    "Medium_vs_Low: ", dim(fdr_ML1)[1]," hits \n",
    "High_vs_Medium: ", dim(fdr_HM1)[1]," hits \n")

## CpGs (vectors) to posterior annotation
fdr_HL1_names <- fdr_HL1$ProbeID
fdr_ML1_names <- fdr_ML1$ProbeID
fdr_HM1_names <- fdr_HM1$ProbeID

## Subset of Bonferroni significant hits -------
bonf_HL1 <- subset(High_vs_Low_topTable1, P.Value < bonf_threshold)
bonf_ML1 <- subset(Medium_vs_Low_topTable1, P.Value < bonf_threshold)
cat("Bonferroni hits - Model 1 (full) - 3 Categories (original): \n",
    "High_vs_Low: ", dim(bonf_HL1)[1]," hits \n",
    "Medium_vs_Low: ", dim(bonf_ML1)[1]," hits \n")

## Saving outputs -----------------------------

cat("Saving outputs: \n")
message("Saving outputs: ", Sys.time())

save(High_vs_Low_topTable1,
     Medium_vs_Low_topTable1,
     High_vs_Medium_topTable1,
     file = file.path(cat3orig_folder, "limma_results_model1_cat3original.R"))

save(fdr_HL1_names, fdr_ML1_names, fdr_HM1_names, 
     all_probes_ids, file = file.path(cat3orig_folder, "hits_vectors_cat3original.R"))

write.csv(fdr_HL1, file = file.path(cat3orig_folder, "fdr_hits_HL_fullModel.csv"), row.names = FALSE)
write.csv(fdr_ML1, file = file.path(cat3orig_folder, "fdr_hits_ML_fullModel.csv"), row.names = FALSE)
write.csv(fdr_HM1, file = file.path(cat3orig_folder, "fdr_hits_HM_fullModel.csv"), row.names = FALSE)

## Volcano plots ---------------------------------
volcano_plot(High_vs_Low_topTable1, "High_vs_Low", cat3orig_folder)
volcano_plot(Medium_vs_Low_topTable1, "Medium_vs_Low", cat3orig_folder)
volcano_plot(High_vs_Medium_topTable1, "High_vs_Medium", cat3orig_folder)

## Manhattan plot ---------------------------------
png(file = file.path(cat3orig_folder, "MAplot_High_vs_Low.png"), width = 1600, height = 900, res = 150)
MAplot(High_vs_Low_topTable1, "High vs Low") 
dev.off()
png(file = file.path(cat3orig_folder, "MAplot_Medium_vs_Low.png"), width = 1600, height = 900, res = 150)
MAplot(Medium_vs_Low_topTable1, "Medium vs Low") 
dev.off()
png(file = file.path(cat3orig_folder, "MAplot_High_vs_Medium.png"), width = 1600, height = 900, res = 150)
MAplot(High_vs_Medium_topTable1, "High vs Medium") 
dev.off()

## Genomic categories barplot ---------------------
genomic_categories(fdr_HL1, "High_vs_Low", "FDR", cat3orig_folder)
genomic_categories(fdr_ML1, "Medium_vs_Low", "FDR", cat3orig_folder)
genomic_categories(fdr_HM1, "High_vs_Medium", "FDR", cat3orig_folder)

##########################################################

cat("Original 3-categorical models pipeline completed. \n\n Starting extreme model pipeline. \n")
message("Original 3-categorical models pipeline completed. \n\n Starting extreme model pipeline: ", Sys.time())

#### EXTREME (FIRST AND FOURTH "QUARTILES") ####

#### PERCENTILE-EXTREME CONTRAST ####                                                                                                                                                                         

# Low_Q1 <= 6 (n = 235 ~ 24,1%) , Mid_Q23 = 7-9, High_Q4 >= 10 (n = 262 ~ 26,9%)                                                                                                                                                       

## Output folder -----------------------------                                                                                                                                                                
extreme_folder <- file.path(results_folder, 'ExtremeContrast')                                                                                                                                      
dir.create(extreme_folder, recursive = TRUE, showWarnings = FALSE)                                                                                                                                            

## Model designs -----------------------------                                                                                                                                                                
dm1 <- model.matrix(~0 + predimed_extreme + age + sex + bmi + smoking_type + batch + Bcell + CD4T + CD8T + Mono + NK, data = pheno)                                                                                   
stopifnot(sum(is.na(dm1)) == 0)                                                                                                                                                         
stopifnot(identical(colnames(Mvalues), rownames(dm1)))

## Contrasts -----------------------------          
cm1 <- makeContrasts(                                                                                                                                                                                         
  HighQ4_vs_LowQ1   = predimed_extremeHigh_Q4 - predimed_extremeLow_Q1,
  MidQ23_vs_LowQ1   = predimed_extremeMid_Q23 - predimed_extremeLow_Q1,                                                                                                                                       
  HighQ4_vs_MidQ23  = predimed_extremeHigh_Q4 - predimed_extremeMid_Q23,                                                                                                                                      
  levels = dm1)                                                                                                                                                                                                             

## Fit -----------------------------                                                                                                                                                                          
fit1 <- lmFit(Mvalues, dm1, na.action = na.exclude)
fit1 <- contrasts.fit(fit1, cm1)
fit1 <- eBayes(fit1, trend = TRUE)                                                                      

## TopTables -----------------------------                                                                                                                                                                    
HighQ4_vs_LowQ1_topTable1  <- limma::topTable(fit1, number = Inf, adjust.method = "BH", coef = "HighQ4_vs_LowQ1",  genelist = genelist_aligned)                                                               
MidQ23_vs_LowQ1_topTable1  <- limma::topTable(fit1, number = Inf, adjust.method = "BH", coef = "MidQ23_vs_LowQ1",  genelist = genelist_aligned)                                                               
HighQ4_vs_MidQ23_topTable1 <- limma::topTable(fit1, number = Inf, adjust.method = "BH", coef = "HighQ4_vs_MidQ23", genelist = genelist_aligned)

## Lambdas -----------------------------                                                                                                                                                                      
lambda1_HQ4LQ1 <- measureInflation(HighQ4_vs_LowQ1_topTable1, "model1_HQ4_LQ1", extreme_folder)                                                                                                               
lambda1_MidLQ1 <- measureInflation(MidQ23_vs_LowQ1_topTable1, "model1_MidQ23_LQ1", extreme_folder)                                                                                                            
lambda1_HQ4Mid <- measureInflation(HighQ4_vs_MidQ23_topTable1, "model1_HQ4_MidQ23", extreme_folder)                                                                                                           

## Delta beta (HighQ4 vs LowQ1, the primary contrast) -----------------------------                                                                                                                           
mean_lowQ1  <- rowMeans(Betavalues[, pheno$predimed_extreme == "Low_Q1"],  na.rm = TRUE)                                                                                                                      
mean_midQ23 <- rowMeans(Betavalues[, pheno$predimed_extreme == "Mid_Q23"], na.rm = TRUE)                                                                                                                      
mean_highQ4 <- rowMeans(Betavalues[, pheno$predimed_extreme == "High_Q4"], na.rm = TRUE)                                                                                                                      
delta_beta_HQ4LQ1 <- mean_highQ4 - mean_lowQ1                                                                                                                                                                 
HighQ4_vs_LowQ1_topTable1 <- addDeltaBeta_Cat(HighQ4_vs_LowQ1_topTable1, delta_beta_HQ4LQ1)                                                                                                                   

## FDR hits (genome-wide) -----------------------------                                                                                                                                                       
fdr_HQ4LQ1_1 <- subset(HighQ4_vs_LowQ1_topTable1, adj.P.Val < 0.05)                                                                                                                                          
cat("FDR hits (genome-wide) - extreme HighQ4_vs_LowQ1:\n",                                                                                                                                                    
    "  Model 1 (full):", nrow(fdr_HQ4LQ1_1), "\n")     

## Suggestive hits
suggestive_HQ4LQ1_1 <- subset(HighQ4_vs_LowQ1_topTable1, P.Value < 1e-4)                                                                                                                                          
cat("CpGs with nominal p-value < 1e-4 - extreme HighQ4_vs_LowQ1 - Full Model: ", nrow(suggestive_HQ4LQ1_1), "\n")    

## Subsetting nominal hits (nominal pvalue < 0.05) ------------- 
nominal_hits1 <- subset(HighQ4_vs_LowQ1_topTable1, P.Value < 0.05)
cat("CpGs with p-value < 0.05 - Full model - Extreme: ", dim(nominal_hits1)[1], "\n")

## CpGs (vectors) to posterior annotation
fdr_names1 <- fdr_HQ4LQ1_1$ProbeID
suggestive_names1 <- suggestive_HQ4LQ1_1$ProbeID
nominal_names1 <- nominal_hits1$ProbeID

## Saving outputs -----------------------------                                                                                                                                                     
save(HighQ4_vs_LowQ1_topTable1, MidQ23_vs_LowQ1_topTable1, HighQ4_vs_MidQ23_topTable1,                                                                                                                        
     file = file.path(extreme_folder, "limma_results_model1_extreme.R"))                                                                                                                                      

save(fdr_names1, suggestive_names1, nominal_names1, all_probes_ids, 
     file = file.path(extreme_folder, "hits_vectors_extreme.R"))

write.csv(fdr_HQ4LQ1_1, file = file.path(extreme_folder, "fdr_hits_HQ4LQ1_fullModel.csv"), row.names = FALSE)                                                                             
write.csv(suggestive_HQ4LQ1_1, file = file.path(extreme_folder, "suggestive_hits_HQ4LQ1_fullModel.csv"), row.names = FALSE)                                                                             
write.csv(nominal_hits1, file = file.path(extreme_folder, "nominal_hits_HQ4LQ1_fullModel.csv"), row.names = FALSE)                                                                             

## Volcano plots ---------------------------------
volcano_plot(HighQ4_vs_LowQ1_topTable1, "HighQ4_vs_LowQ1", extreme_folder)

## Manhattan plot ---------------------------------
png(file = file.path(extreme_folder, "MAplot_extreme.png"), width = 1600, height = 900, res = 150)
MAplot(HighQ4_vs_LowQ1_topTable1, "Extreme (HighQ4 vs LowQ1)") 
dev.off()

## Genomic categories barplot ---------------------
genomic_categories(suggestive_HQ4LQ1_1, "Extreme", "p-value", extreme_folder)

##########################################################

#### LAMBDAS COMPARISON ####

lambdas <- data.frame(                                                                                                                                                                                        
  analysis = c("continuous_full",                                                                                                                                    
               "binary_full", 
               "cat3_HL_full", "cat3_ML_full", "cat3_HM_full",
               "extreme_HQ4LQ1_full", "extreme_MidLQ1_full", "extreme_HQ4Mid_full"),                                                                                                                           
  lambda = c(lambda1cont,                                                                                                                                                       
             lambda1bin,              
             lambda1HL, lambda1ML, lambda1HM,   
             lambda1_HQ4LQ1, lambda1_MidLQ1, lambda1_HQ4Mid))
write.csv(lambdas, file.path(results_folder, "lambdas.csv"), row.names = FALSE)

##########################################################

cat("Script completed. \n")
message("Script completed: ", Sys.time())

##########################################################
