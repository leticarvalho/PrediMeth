##########################################################
#
# DMR analysis - Stage 1
#
##########################################################


# DESCRIPTION --------------------------------------------
# In this script, EWAS results will be investigated for
# existence of Different Methylated Regions (DMRs).

# INPUT --------------------------------------------------
# 	Processed data (from discovery cohort): # pheno, beta_matrix, mval_matrix

# OUTPUT -------------------------------------------------
# 	For each contrast that has DMPs (DMRcate)
#       DMPs 
#       DMRs 

##########################################################

cat("\n\n\nStarting script - Data Preprocessing - Stage 1. \n")
message("Starting script: ", Sys.time() )

#### HUB CONFIGURATION - TO USE DMRCATE ####

scratch <- Sys.getenv("SCRATCH")
if (scratch == "") scratch <- "/tmp"
cache_dir <- file.path(scratch, "R_cache")
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(HOME = scratch)
Sys.setenv(XDG_CACHE_HOME = cache_dir)

eh_cache  <- file.path(cache_dir, "ExperimentHub")
ah_cache  <- file.path(cache_dir, "AnnotationHub")
bfc_cache <- file.path(cache_dir, "BiocFileCache")
dir.create(eh_cache,  recursive = TRUE, showWarnings = FALSE)
dir.create(ah_cache,  recursive = TRUE, showWarnings = FALSE)
dir.create(bfc_cache, recursive = TRUE, showWarnings = FALSE)

library(ExperimentHub)
library(AnnotationHub)
setExperimentHubOption("CACHE", eh_cache)
setAnnotationHubOption("CACHE", ah_cache)

#### LOAD LIBRARIES ####

library(Biobase)              
library(SummarizedExperiment) 
library(minfi)                 
library(DMRcate)             
library(GenomicRanges)
library(limma)
library(missMethyl)
library(dplyr)
library(data.table)
library(edgeR)
library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)

loadNamespace("minfi")
cat("All packages loaded.\n")

#### RESOLVE PATHS ####

## Input paths -------------------------------------------

path_to_processed_data_stage1 <- "" # R.data

## Output paths ------------------------------------------

predimeth_path <- ""
results_dir <- file.path(predimeth_path, "results")
results_folder <- file.path(results_dir, 'OfficialAnalysis', 'DMR' )
dir.create(results_folder)

##########################################################

#### LOAD DATA ####

## Preprocessed data (Stage 1) ---------------------------
load(path_to_processed_data_stage1) # pheno, beta_matrix, mval_matrix
Mvalues <- mval_matrix 
Betavalues <- beta_matrix

## Annotation to hg38 -------------------------------------
ann.hg38 <- getAnnotation(IlluminaHumanMethylationEPICv2anno.20a1.hg38)

##########################################################

#### DATA PREPARATION ####

## Alignment
stopifnot(all(rownames(pheno) %in% colnames(Mvalues)))
stopifnot(identical(colnames(Mvalues), rownames(pheno)))
cat("Sample alignment verified.\n")

## Change rownames of Mvalues (adding EPIC v2 suffix)
cat("Adding EPIC V2 suffix to row names, according to hg38 annotation: ")
base_ids    <- sub("_.*", "", rownames(ann.hg38))

matched_idx <- match(rownames(Mvalues), base_ids)
rownames(Mvalues) <- rownames(ann.hg38)[matched_idx]
cat(rownames(Mvalues)[1:5], "\n")

matched_idx <- match(rownames(Betavalues), base_ids)
rownames(Betavalues) <- rownames(ann.hg38)[matched_idx]
cat(rownames(Betavalues)[1:5], "\n")

## Checking structure of pheno table
str(pheno)

##########################################################

cat("Data loaded \nStarting DMR analysis with DMRcate: \n")
message("Data loaded. Starting DMR analysis with DMRcate: ", Sys.time())

#### DMR ANALYSES WITH DMRCATE ####
setClassUnion("ExpData", c("matrix", "SummarizedExperiment"))

## Continuous ##

dm1_cont <- model.matrix(~ predimed_score + age + sex + bmi + smoking_type + batch + Bcell + CD4T + CD8T + Mono + NK, data = pheno) 


# DMP 
cpg_anno <- cpg.annotate(
  datatype = "array",
  object = Mvalues,
  what = "M",
  annotation = ann.hg38,
  arraytype = "EPICv2",
  analysis.type = "differential",
  design = dm1_cont,
  contrasts = FALSE,
  fdr = 0.05,  
  coef = "predimed_score",
  epicv2Filter = "mean")

cat("\n Continuous CpG-annotated object summary:\n")
print(cpg_anno)

########################

## Binary ##

# Design matrix
dm1_bin <- model.matrix(~0 + predimed_high + age + sex + bmi + smoking_type + batch + Bcell + CD4T + CD8T + Mono + NK, data = pheno) 

# Contrast matrix
cm1_bin <- makeContrasts(predimed_highyes_vs_predimed_highno = predimed_highyes-predimed_highno, levels = dm1_bin)
contrasts <- colnames(cm1_bin)
print(contrasts)

# DMP 
cpg_anno <- cpg.annotate(
  datatype = "array",
  object = Mvalues,
  what = "M",
  annotation = ann.hg38,
  arraytype = "EPICv2",
  analysis.type = "differential",
  design = dm1_bin,
  contrasts = TRUE,
  cont.matrix = cm1_bin,
  fdr = 0.05,  # more permissive
  coef = contrasts[1],
  epicv2Filter = "mean")

cat("\n Binary CpG-annotated object summary:\n")
print(cpg_anno)

########################

## 3-categorical contrast ##

# Design matrix
dm1_3cat <- model.matrix(~0 + predimed_cat + age + sex + bmi + smoking_type + batch + Bcell + CD4T + CD8T + Mono + NK, data = pheno) 

# Contrast matrix
cm1_3cat <-makeContrasts(High_vs_Low = predimed_cathigh - predimed_catlow,
                         Medium_vs_Low = predimed_catmedium - predimed_catlow,
                         High_vs_Medium = predimed_cathigh - predimed_catmedium,
                         levels = dm1_3cat)
contrasts <- colnames(cm1_3cat)
print(contrasts)

# DMP High vs Low
cpg_annoHL <- cpg.annotate(
  datatype = "array",
  object = Mvalues,
  what = "M",
  annotation = ann.hg38,
  arraytype = "EPICv2",
  analysis.type = "differential",
  design = dm1_3cat,
  contrasts = TRUE,
  cont.matrix = cm1_3cat,
  fdr = 0.05,
  coef = "High_vs_Low",
  epicv2Filter = "mean"
)

cat("\n High vs Low CpG-annotated object summary:\n")
print(cpg_annoHL)

# DMP Medium vs Low
cpg_annoML <- cpg.annotate(
  datatype = "array",
  object = Mvalues,
  what = "M",
  annotation = ann.hg38,
  arraytype = "EPICv2",
  analysis.type = "differential",
  design = dm1_3cat,
  contrasts = TRUE,
  cont.matrix = cm1_3cat,
  fdr = 0.05, 
  coef = "Medium_vs_Low",
  epicv2Filter = "mean"
)

cat("\n Medium vs Low CpG-annotated object summary:\n")
print(cpg_annoML)

save(cpg_annoHL,
     file = file.path(results_folder, "DMP_HL.R"))
save(cpg_annoML,
     file = file.path(results_folder, "DMP_ML.R"))

####

## DMR High vs Low ##

dmr_resultsHL <- dmrcate(cpg_annoHL,
                         lambda = 1000, # bandwidth in bp 
                         C = 2, # scales lambda; effective window = C * lambda
                         min.cpgs = 2) # minimum probes per DMR

cat("\n High vs Low DMR results summary:\n")
print(dmr_resultsHL)

## Extract DMR data frame 
dmr_rangesHL <- extractRanges(dmr_resultsHL, genome = "hg38")
dmr_dfHL <- as.data.frame(dmr_rangesHL) %>% arrange(Stouffer)   # sort by combined p-value

cat("\nTop 20 DMRs HL:\n")
print(head(dmr_dfHL, 20))

## DMRcate built-in DMR plot (top 3 DMRs)
if (nrow(dmr_dfHL) >= 1) {

  n_plot <- min(3, nrow(dmr_dfHL))

  # Open PDF ONCE
  pdf(file.path(results_folder, "DMRcate_top_DMRsHL.pdf"),
      width = 10, height = 6)

  for (i in seq_len(n_plot)) {

    tryCatch({

      DMR.plot(
        ranges     = dmr_rangesHL,
        dmr        = i,
        CpGs       = Betavalues,
        what       = "Beta",
        arraytype  = "EPICv2",
        genome     = "hg38",
        phen.col   = ifelse(pheno$predimed_cat == "high",
                            "green",
                            ifelse(pheno$predimed_cat == "medium",
                                   "blue", "red"))
      )

      title(main = paste0(
        "DMR #", i, " – ",
        dmr_dfHL$seqnames[i], ":",
        dmr_dfHL$start[i], "-",
        dmr_dfHL$end[i]
      ))

    }, error = function(e) {
      message("Could not plot DMR #", i, ": ", conditionMessage(e))
    })
  }

  # Close PDF ONCE
  dev.off()

  cat("DMR HL plots saved to: DMRcate_top_DMRsHL.pdf\n")
}


####

## DMR Medium vs Low ##

dmr_resultsML <- dmrcate(cpg_annoML, lambda = 1000, C = 2, min.cpgs = 2)

cat("\n Medium vs Low DMR results summary:\n")
print(dmr_resultsML)

## Extract DMR data frame 
dmr_rangesML <- extractRanges(dmr_resultsML, genome = "hg38")
dmr_dfML <- as.data.frame(dmr_rangesML) %>% arrange(Stouffer)   # sort by combined p-value

cat("\nTop 20 DMRs ML:\n")
print(head(dmr_dfML, 20))

## DMRcate built-in DMR plot (top 3 DMRs)
if (nrow(dmr_dfML) >= 1) {

  n_plot <- min(3, nrow(dmr_dfML))

  # Open PDF ONCE
  pdf(file.path(results_folder, "DMRcate_top_DMRsML.pdf"),
      width = 10, height = 6)

  for (i in seq_len(n_plot)) {

    tryCatch({

      DMR.plot(
        ranges     = dmr_rangesML,
        dmr        = i,
        CpGs       = Betavalues,
        what       = "Beta",
        arraytype  = "EPICv2",
        genome     = "hg38",
        phen.col   = ifelse(pheno$predimed_cat == "high",
                            "green",
                            ifelse(pheno$predimed_cat == "medium",
                                   "blue", "red"))
      )

      title(main = paste0(
        "DMR #", i, " – ",
        dmr_dfHL$seqnames[i], ":",
        dmr_dfHL$start[i], "-",
        dmr_dfHL$end[i]
      ))

    }, error = function(e) {
      message("Could not plot DMR #", i, ": ", conditionMessage(e))
    })
  }

  # Close PDF ONCE
  dev.off()

  cat("DMR ML plots saved to: DMRcate_top_DMRsML.pdf\n")
}

###########

# Enrichment with goregion() - HL

sig_maskHL <- dmr_dfHL$Stouffer < 0.05
sig_dmr_rangesHL <- dmr_rangesHL[sig_maskHL]
all_dmr_rangesHL <- dmr_rangesHL 
cat("\nSignificant DMRs HL for enrichment:", sum(sig_maskHL), "\n")

save(dmr_resultsHL, sig_dmr_rangesHL, all_dmr_rangesHL,
     file = file.path(results_folder, "DMR_HL.R"))

if (sum(sig_maskHL) == 0) {
  warning("No significant DMRs HL at Stouffer p < 0.05. ",
          "Consider relaxing the threshold or checking your data.")
} else {
  ## GO ##
  go_results <- goregion(regions = sig_dmr_rangesHL, all.cpg  = NULL, collection = "GO",
                         array.type = "EPICv2", anno = ann.hg38, plot.bias = FALSE)
  
  ## Filter for significance
  go_dmr_sig <- go_results %>% filter(FDR < 0.05) %>% arrange(FDR)
  cat("\nSignificant GO HL terms (FDR < 0.05):\n")
  print(head(go_dmr_sig[, c("ONTOLOGY", "TERM", "N", "DE", "FDR")], 20))
  
  ########################
  
  ## KEGG ##
  kegg_results <- goregion(regions = sig_dmr_rangesHL, all.cpg = NULL, collection = "KEGG",
                           array.type = "EPICv2", anno = ann.hg38, plot.bias = FALSE)
  ## Filter for significance
  kegg_dmr_sig <- kegg_results %>% filter(FDR < 0.05) %>% arrange(FDR)
  cat("\nSignificant KEGG HL pathways (FDR < 0.05):\n")
  print(head(kegg_dmr_sig[, c("TERM", "N", "DE", "FDR")], 20))
}


     
# Enrichment with goregion() - ML

sig_maskML <- dmr_dfML$Stouffer < 0.05
sig_dmr_rangesML <- dmr_rangesML[sig_maskML]
all_dmr_rangesML <- dmr_rangesML       
cat("\nSignificant DMRs ML for enrichment:", sum(sig_maskML), "\n")

save(dmr_resultsML, sig_dmr_rangesML, all_dmr_rangesML,
     file = file.path(results_folder, "DMR_ML.R"))

if (sum(sig_maskML) == 0) {
  warning("No significant DMRs ML at Stouffer p < 0.05. ",
          "Consider relaxing the threshold or checking your data.")
} else {
  ## GO ##
  go_results <- goregion(regions = sig_dmr_rangesML, all.cpg  = NULL, collection = "GO",
                         array.type = "EPICv2", anno = ann.hg38, plot.bias = FALSE)
  
  ## Filter for significance
  go_dmr_sig <- go_results %>% filter(FDR < 0.05) %>% arrange(FDR)
  cat("\nSignificant GO ML terms (FDR < 0.05):\n")
  print(head(go_dmr_sig[, c("ONTOLOGY", "TERM", "N", "DE", "FDR")], 20))
  
  ########################
  
  ## KEGG ##
  kegg_results <- goregion(regions = sig_dmr_rangesML, all.cpg = NULL, collection = "KEGG",
                           array.type = "EPICv2", anno = ann.hg38, plot.bias = FALSE)
  ## Filter for significance
  kegg_dmr_sig <- kegg_results %>% filter(FDR < 0.05) %>% arrange(FDR)
  cat("\nSignificant KEGG ML pathways (FDR < 0.05):\n")
  print(head(kegg_dmr_sig[, c("TERM", "N", "DE", "FDR")], 20))
}

##########################################################

cat("Script completed. \n")
message("Script completed: ", Sys.time())

##########################################################


