##########################################################
#
# Studying MRS and T2D associations
#
##########################################################


# DESCRIPTION --------------------------------------------
# In this script, we investigate if MRSs are associated 
# with T2D in two different contexts. For prevalent T2D,
# we use a logistic regression and odds ratio. For 
# incident T2D, we use Cox proportional hazards model and 
# hazard ratio. 

# INPUT --------------------------------------------------
# 	Validation cohort (already preprocessed)

# OUTPUT -------------------------------------------------
# 	Forest plot summarizing MRS and T2D associations
#   file 2




##########################################################

#### LOAD LIBRARIES ####


library(dplyr)
library(data.table)
library(minfi)
library(glmnet)
library(gridExtra)
library(pROC)
library(tibble)
library(broom)
library(tidyr)
library(ggplot2)

#### RESOLVE PATHS ####

## Input paths -------------------------------------------

predimeth_path <- "/imppc/labs/dnalab/share/PrediMeth"
results_dir <- file.path(predimeth_path, "results")

## Output paths ------------------------------------------

results_folder <- file.path(results_dir, 'MRS_T2D_association')
dir.create(results_folder)

##########################################################

#### LOAD DATA ####

## Data 1 ------------------------------------------------

valid.cohort <- fread("path/to/valid/cohort.txt")
valid.cohort <- as.data.frame(valid.cohort)
rownames(valid.cohort) <- valid.cohort$sample_id

##########################################################

cat("Data loaded. \n\nStarting subsets definition. \n")
message("Data loaded. \n\nStarting subsets definition: ", Sys.time())

#### DIVIDING VALIDATION COHORT IN SUBSETS ####

## Preparing variables
valid.cohort$age <- as.numeric(valid.cohort$EDAD_ANOS)
valid.cohort$sex <- factor(valid.cohort$SEXO, levels = c(1,2), labels = c("men", "women"))
valid.cohort$bmi <- as.numeric(valid.cohort$BMI)
valid.cohort$smoking_type <- factor(valid.cohort$smoking_habit, levels = c(1,2,3), labels = c("smoker","exsmoker", "nonsmoker"))

## Define analyses subsets (prevalent and incident T2D)

# pheno_prev: all validation samples
##### cases: t2d_prevalent == 1  (T2D at blood extraction)
##### controls: t2d_prevalent == 0 (control at blood extraction)
pheno_prev <- valid.cohort

# pheno_inc: samples without T2D at blood extraction, excluding prevalent cases (t2d_prevalent == 0)
##### cases: t2d_incident == 1  
##### controls: t2d_incident == 0 
pheno_inc <- subset(valid.cohort, t2d_prevalent == 0)

##########################################################

cat("Subsets defined. \n\nStarting logistic regression. \n")
message("Subsets defined. \nStarting logistic regression: ", Sys.time())

## Logistic regression - prevalent T2D 

# Empty lists to store results
logistic_models <- list() # stores each fitted object
logistic_table <- list() # stores a summary for the MRS term (for later forest plot / summary table)

# Covariates
covariates <- c("age", "sex", "bmi", "smoking_type", "Bcell", "CD4T", "CD8T", "Mono", "NK")
cov_formula <- paste(covariates, collapse = " + ")
print(cov_formula)

# Function to fit all calculated MRS_z scores (logistic regression)
logistic_MRS <- function(mrs) {
  form <- as.formula(paste0("t2d_prevalent ~ ", mrs, " +" , cov_formula))
  fit <- glm(form, data = pheno_prev, family = binomial())
  result <- broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE) 
  result <- filter(result, term == mrs) %>% mutate(MRS = mrs)
  return(result)
}

# Apply function and creat result table
logistic_results <- lapply(mrs_z_cols, logistic_MRS)
logistic_table <- bind_rows(logistic_results)

# See results
cat("\n Logistic regression: Prevalent T2D (OR per SD) \n")
print(logistic_table[, c("MRS","estimate","conf.low","conf.high","p.value")])

write.table(logistic_table, file = file.path(results_folder, "MRS_logistic_prevalent.txt"), sep = "\t", row.names = FALSE, col.names = TRUE)

# Model evaluation: function to calculate AUC for each model
calc_AUC <- function(mrs) {
  roc <- pROC::roc(pheno_prev$t2d_prevalent, pheno_prev[[mrs]], quiet = TRUE) # calculate ROC object
  auc <- pROC::auc(roc)
  auc <- as.numeric(auc)
  return(auc)
}
auc_prev <- sapply(mrs_z_cols, calc_AUC)
names(auc_prev) <- mrs_cols # names without "_z"
cat("\nAUC (prevalent T2D):\n")
print(round(sort(auc_prev, decreasing = TRUE), 3))

##########################################

## Cox proportional hazards - incident T2D

# Function to fit all calculated MRS_z scores (cox proportional hazards)
cox_MRS <- function(mrs) {
  form <- as.formula(paste("Surv(TTE, t2d_incident) ~", mrs, "+", cov_formula))
  fit <- coxph(form, data = pheno_inc)
  result <- broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE) 
  result <- filter(result, term == mrs) %>% mutate(MRS = mrs)
  return(result)
}

# Apply function and create result table
cox_results <- lapply(mrs_z_cols, cox_MRS)
cox_table <- bind_rows(cox_results)

# See results
cat("\n Cox proportional hazards: Incident T2D (HR per SD) \n")
print(cox_table[, c("MRS","estimate","conf.low","conf.high","p.value")])

write.table(cox_table, file = file.path(results_folder, "MRS_cox_incident.txt"), sep = "\t", row.names = FALSE, col.names = TRUE)

# Model evaluation: function to calculate Cindex for each model
calc_cindex <- function(mrs) {
  form <- as.formula(paste("Surv(TTE, t2d_incident) ~", mrs, "+", cov_formula))
  cindex <- summary(coxph(form, data = pheno_inc))$concordance["C"] # C-index
  return(cindex)
}
cindex_inc <- sapply(mrs_z_cols, calc_cindex)
names(cindex_inc) <- mrs_cols # names without "_z"
cat("\nC-index (incident T2D):\n")
print(round(sort(cindex_inc, decreasing = TRUE), 3))

##########################################

## Proportional hazards check

# Choosing the MRS with the best C-index
best_mrs_z <- mrs_z_cols[which.max(cindex_inc)]

# Fitting the model 
form <- as.formula(paste("Surv(TTE, t2d_incident) ~", best_mrs_z, "+", cov_formula))
model <- coxph(form, data = pheno_inc)
prop_haz_test <- cox.zph(model) # Schoenfeld residuals test (SRT)
# SRT: for each covariate it tests whether the residuals are correlated with time; 
# if they are, the effect of that variable is changing over time, violating the proportional hazards assumption.
# Wanted: non-significant p-values (p > 0.05) for all terms, especially for your MRS

# Seeing results
cat("\nSchoenfeld test (", best_mrs_z, ") \n", sep = "")
print(prop_haz_test)
plot(prop_haz_test, var = best_mrs_z,
     main = paste("Schoenfeld residuals -", best_mrs_z))

##########################################

## Summary table (prevalence and incidence) ---------------------------------

n_cpgs <- c(MRS1 = length(weights_S1),
            MRS2 = length(weights_S2),
            sapply(weights_S3_list, length))[mrs_cols]

# Mean contributing CpGs per sample 
n_contrib_cols <- paste0("n_cpgs_", mrs_cols)
mean_contrib <- sapply(n_contrib_cols, function(col) {
  if (col %in% colnames(valid.cohort)) {
    mean(valid.cohort[[col]], na.rm = TRUE)
  } else {
    NA_real_ # without it, sapply might return a logical vector instead of a numeric
  }
})
names(mean_contrib) <- mrs_cols
mean_contrib_pct <- round(mean_contrib / n_cpgs * 100, 1)

summary_df <- data.frame(MRS = mrs_cols, n_CpGs = n_cpgs, mean_contrib_pct = paste0(mean_contrib_pct, "%"),
                         AUC_prev = round(auc_prev[mrs_cols], 3),
                         OR_per_SD = round(logistic_table$estimate, 2),
                         OR_95CI = paste0("[", round(logistic_table$conf.low,  2), "-", round(logistic_table$conf.high, 2), "]"),
                         OR_p = signif(logistic_table$p.value, 3),
                         Cindex_inc = round(cindex_inc[mrs_cols], 3),
                         HR_per_SD = round(cox_table$estimate, 2),
                         HR_95CI = paste0("[", round(cox_table$conf.low,  2), "-", round(cox_table$conf.high, 2), "]"),
                         HR_p = signif(cox_table$p.value, 3),
                         row.names = NULL)
cat("\nSummary table: \n")
print(summary_df)
write.table(summary_df, file = file.path(results_folder, "MRS_summary.txt"), sep = "\t", row.names = FALSE, col.names = TRUE)


## Forest plots (prevalence and incidence) ---------------------------------

forest_plot_OR <- logistic_table %>% mutate(MRS = factor(MRS, levels = rev(MRS))) %>%
  ggplot(aes(x = estimate, xmin = conf.low, xmax = conf.high, y = MRS)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey60") +
  geom_errorbarh(height = 0.25, colour = "#2166AC", linewidth = 0.7) +
  geom_point(size = 3, colour = "#2166AC") +
  scale_x_log10() +
  labs(title = "Prevalent T2D", x = "OR per SD (95% CI)", y = NULL) +
  theme_bw(base_size = 11)

forest_plot_HR <- cox_table %>% mutate(MRS = factor(MRS, levels = rev(MRS))) %>%
  ggplot(aes(x = estimate, xmin = conf.low, xmax = conf.high, y = MRS)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey60") +
  geom_errorbarh(height = 0.25, colour = "#D6604D", linewidth = 0.7) +
  geom_point(size = 3, colour = "#D6604D") +
  scale_x_log10() +
  labs(title = "Incident T2D", x = "HR per SD (95% CI)", y = NULL) +
  theme_bw(base_size = 11)

png(file.path(results_folder, "forestplots_plot.png"), width = 800, height = 600)
gridExtra::grid.arrange(forest_plot_OR, forest_plot_HR, ncol = 2,
                        top = "MRS associations with Type 2 Diabetes (validation cohort)")
dev.off()


##########################################################

cat("Script completed. \n")
message("Script completed: ", Sys.time())

##########################################################