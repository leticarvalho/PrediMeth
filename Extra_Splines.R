##########################################################
#
# Splines Exploratory Analysis - Extra
#
##########################################################


# DESCRIPTION --------------------------------------------
# In this script, an exploratory analysis of data is 
# carried, using natural splines (df = 3), in order to 
# try to capture non-linearity in a possible association
# of predimed score and differential methylated CpGs.

# INPUT --------------------------------------------------
#   Annotation (EPIC v2)
#   Processed data from discovery cohort (Betas, Mvalues, pheno)
#   TopTable from continuos contrast
#   TopTable from HL contrast

# OUTPUT -------------------------------------------------
# 	Table with splines hits (possibly curved association)
#   Plots of FDR significant curved CpGs, of predimed score 
#       in x-axis and predicted betas in y-axis (adjusted to mean covariates)

##########################################################

cat("\n\n\nStarting script. \n")
message("Starting script: ", Sys.time() )

#### LOAD LIBRARIES ####

library(dplyr)
library(ggplot2)
library(ggrepel)
library(qvalue)       
library(data.table)
library(limma)
library(splines)   

#### RESOLVE PATHS ####

## Input paths -------------------------------------------

path_to_EPIC2_annotation <- "/imppc/labs/dnalab/share/PrediMeth/0-data/EPICv2.annot.RData" # .Rdata
path_to_processed_data_stage1 <- "/imppc/labs/dnalab/share/PrediMeth/results/OfficialAnalysis/processed_data.R" # official (clamp6)
path_to_continuous_limma_results <- "/imppc/labs/dnalab/share/PrediMeth/results/OfficialAnalysis/ContinuousContrast/limma_results_continuous.R"

## Output paths ------------------------------------------

results_folder <- "/imppc/labs/dnalab/share/PrediMeth/results/OfficialAnalysis/"
splines_folder <- file.path(results_folder, 'SplinesTest')
dir.create(splines_folder, recursive = TRUE, showWarnings = FALSE)

##########################################################

#### LOAD DATA ####

## Preprocessed data (Stage 1) ---------------------------
load(path_to_processed_data_stage1) # pheno, beta_matrix, mval_matrix

## Continuous limma results
load(path_to_continuous_limma_results)

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
cont_topTable1_saved <- full_topTable1
HL_topTable <- High_vs_Low_topTable1
str(pheno)

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

#### VECTOR WITH ALL CPGS IN THE ANALYSIS ####

all_probes_names <- rownames(Mvalues)
cat("Number of all CpGs used in the analysis: ", length(all_probes_names), "\n")
all_probes_ids <- EPIC_MEsteller[match(rownames(Mvalues), EPIC_MEsteller$Probe_name), ]$ProbeID

#### ALIGN ANNOTATION ####

genelist_aligned <- EPIC_MEsteller[match(rownames(Mvalues), EPIC_MEsteller$Probe_name), ]                                                                                                                                                      
stopifnot(identical(genelist_aligned$Probe_name, rownames(Mvalues)))


#### PREPARING FUNCTIONS ####

## To evaluate p-value inflation
measureInflation <- function(topTable, label, output_folder) {
  pvals <- topTable$P.Value
  pvals <- pvals[!is.na(pvals)]                                                                                                                                                                              
  chisq <- qchisq(1 - pvals, df = 1)
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

## To compare linear and spline results for each Cpg with an ANOVA on beta-values (for interpretability), 
## and check whether it suggests curvature 
visualize_hit <- function(probe, Betavalues, pheno) {                                                                                                                                                         
  d <- data.frame(                                                                                                                                                                                            
    beta = Betavalues[probe, ],
    score = pheno$predimed_score,                                                                                                                                                                             
    age = pheno$age, sex = pheno$sex, bmi = pheno$bmi, smk = pheno$smoking_type, batch = pheno$batch,
    Bcell = pheno$Bcell, CD4T = pheno$CD4T, CD8T = pheno$CD8T, Mono = pheno$Mono, NK = pheno$NK)                                                                                                                                                                                                           
  d <- na.omit(d) # NA is not expected, but just to be sure
  fit_lin <- lm(beta ~ score + age + sex + bmi + smk + batch + Bcell + CD4T + CD8T + Mono + NK, data = d)                                                                                                                                          
  fit_spline <- lm(beta ~ splines::ns(score, df = 3) + age + sex + bmi + smk + batch + Bcell + CD4T + CD8T + Mono + NK, data = d)                                                                                                                                          
  
  # does curvature adds anything beyond a straight line? test it with ANOVA (true if pvalue < 0.05)
  anova <- anova(fit_lin, fit_spline)
  p_nonlin <- anova(fit_lin, fit_spline)[2, "Pr(>F)"]
  cat("\n", probe, " - results anova(lin, spline)", anova, "\np_nonlin = ", p_nonlin,"\n")
  
  # results
  list(fit_lin = fit_lin, fit_spline = fit_spline, p_nonlinearity = p_nonlin)   
} 

## To plot the predicted curve 
plot_hit <- function(probe, fit_spline, pheno) {                                                                                                                                                              
  grid <- data.frame(
    # fix all covariates at representative values ("reference individual"): 
    score = seq(min(pheno$predimed_score), max(pheno$predimed_score), length = 100),
    age = median(pheno$age), # median age
    sex = factor("women", levels = levels(pheno$sex)), # predominant sex
    bmi = median(pheno$bmi), # median BMI
    smk = factor("nonsmoker", levels = levels(pheno$smoking_type)), # nonsmoker
    batch = names(sort(table(pheno$batch), decreasing = TRUE))[1], # most common batch
    Bcell = mean(pheno$Bcell), CD4T = mean(pheno$CD4T), CD8T  = mean(pheno$CD8T),  Mono = mean(pheno$Mono), NK = mean(pheno$NK)) # mean cell counts
  # predict beta value and CI
  pred <- predict(fit_spline, newdata = grid, se.fit = TRUE)
  grid$fit <- pred$fit                                                                                                                                                                                        
  grid$lo <- pred$fit - 1.96 * pred$se.fit                                                                                                                                                                   
  grid$hi <- pred$fit + 1.96 * pred$se.fit  
  # draw plot with predicted curve
  ggplot2::ggplot(grid, aes(score, fit)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.2) + geom_line() +   
    scale_x_continuous( breaks = 0:14, limits = c(0, 14)) +
    labs(title = probe, x = "predimed_score", y = "predicted beta value (covariates fixed)") + theme_minimal()                                                                                                                                                          
}

## To create a Manhattan plot
MAplot <- function(topTable, contrast) {
  ma_data <- subset(topTable, select = c("CpGchrm", "CpGbeg", "adj.P.Val", "Probe_name"))
  #~ 	###ma_data$CpGchrm <- stringr::str_remove(ma_data$CpGchrm, "chr*")
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
    theme_bw() + labs(x = "Genomic Category", y = "Number of CpGs",
                      title = paste0(contrast, " ", significance_type, " significant CpGs by Genomic Category")) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(file.path(output_folder, paste0("GenomicCategoriesPlot_", contrast, ".pdf")), plot = plot)
}


##########################################################

cat("Starting splines test to catch curvature. \n")
message("Starting splines test to catch curvature: ", Sys.time())

#### CHECK FOR NON-LINEARITY - SPLINES (df = 3) ####

# Due to suspicion that some CpGs may have a threshold shape (effect at HL and ML contrast, but none in HM contrast) 
# Natural splines (Restricted cubic splines) were used to test for any non-linear effect of diet, 
# modeling curvature with 3 extra df instead of 13 df (14 scores - 1);
# if the splines outperform continuous, it suggests curvature; 
# if they don't, the relationship is essentially linear and continuous was already optimal -> no need to insist on different categorization

# Model 
dm1 <- model.matrix(~ ns(predimed_score, df = 3) + age + sex + bmi + smoking_type + batch + Bcell + CD4T + CD8T + Mono + NK, data = pheno)                                                                                                                                                                                                             
fit1 <- lmFit(Mvalues, dm1)
fit1 <- eBayes(fit1, trend = TRUE)  

# Testing the joint significance of the spline ("any effect of score")                                                                                                                                     
cols <- grep("ns\\(predimed_score", colnames(dm1))
top1 <- topTable(fit1, coef = cols, number = Inf, adjust.method = "BH", genelist = genelist_aligned)  
# computes a moderated F-statistic that jointly tests all three spline coefficients simultaneously (H0 = Coef1 = Coef2 = Coef3 = 0)
# to reject H0 means that some function of predimed_score (linear or curved) is associated with methylation at the CpG
lambda_spline <- measureInflation(top1, "spline_omnibus", splines_folder)                                                                                                                                     
cat("Lambda (spline) =", lambda_spline, "\n")

# FDR significant CpG sites
tophits <- rownames(top1)[top1$adj.P.Val < 0.05]
cat("Non-zero relationship in", length(tophits), "probes (possibly curved) \n")    
tophits

# Plotting predicted curves
nonlin <- numeric(length(tophits))
names(nonlin) <- tophits
splines_curves_folder <- file.path(splines_folder, 'Curves')
dir.create(splines_curves_folder)

# For loop: for each of the top hits, avaliate nonlinearity and, if so, construct the predicted curve
for (p in seq_along(tophits)) {                                                                                                                                                                               
  probe <- tophits[p]
  results <- visualize_hit(probe, Betavalues, pheno)                                                                                                                                                          
  nonlin[p] <- results$p_nonlinearity
  if (results$p_nonlinearity < 0.05) {
    png(file = file.path(splines_curves_folder, paste0(probe, "_SplineCurve.png")), width = 1600, height = 900, res = 150)
    print(plot_hit(probe, results$fit_spline, pheno))                                                                                                                                                         
    dev.off()
  }                                                                                                                                                                                                           
}               

## Organizing splines hits topTable 
Splines_topTable1 <- as.data.frame(top1)
spline_fdr_hits1 <- Splines_topTable1[Splines_topTable1$adj.P.Val < 0.05,]
spline_fdr_hits1 <- rename(spline_fdr_hits1, "spline_P" = P.Value, "spline_FDR" = adj.P.Val)

## Adding continuous findings to spline topTable
spline_fdr_hits1$linear_P <- cont_topTable1_saved[tophits, "P.Value"]
spline_fdr_hits1$linear_FDR <- cont_topTable1_saved[tophits, "adj.P.Val"]
spline_fdr_hits1$delta_beta_per_SD <- cont_topTable1_saved[tophits, "delta_beta_per_SD"]
spline_fdr_hits1$delta_beta_per_unit <- cont_topTable1_saved[tophits, "delta_beta_per_unit"]

## Organizing non-linear findings 
spline_fdr_hits1$p_nonlinearity <- nonlin[tophits]
spline_fdr_hits1$relationship_type <- ifelse(spline_fdr_hits1$p_nonlinearity < 0.05, "curved", "linear")

## Saving outputs
save(Splines_topTable1, file = file.path(splines_folder, "Splines_topTables.R"))
write.csv(spline_fdr_hits1, file = file.path(splines_folder, "spline_fdr_hits_full.csv"))

## CpGs (vectors) to posterior annotation
spline_fdr_cpgs <- spline_fdr_hits1$ProbeID
save(spline_fdr_cpgs, all_probes_ids, file = file.path(splines_folder, "hits_vectors_spline.R"))


##########################################################

cat("Script completed. \n")
message("Script completed: ", Sys.time())

##########################################################
