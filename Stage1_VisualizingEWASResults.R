##########################################################
#
# Visualizing EWAS Results - Stage 1
#
##########################################################


# DESCRIPTION --------------------------------------------
# In this script, EWAS results are visualized by plots, 
# regarding different aspects. R packages like ggplot2,
# qqman and EnhancedVolcano will be used. 

# INPUT --------------------------------------------------
# 	EWAS results (topTable) 

# OUTPUT -------------------------------------------------
# 	FDR-significant topTable subset
# 	CpGs in vectors (to use in enrichment)
# 	Volcano Plot
#   Manhattan Plot
# 	QQplot (genomic inflation analysis)
# 	Genomic Categories Bar Plot


##########################################################

#### LOAD LIBRARIES ####

library(ggplot2)
library(EnhancedVolcano)
library(qqman)
library(RColorBrewer)
library(dplyr)


#### RESOLVE PATHS ####

## Define label ------------------------------------------

label <- "_3cat" # Probe ID
#label <- "_bin" # Probe ID

## Input paths -------------------------------------------

# Project (PrediMeth)
predimeth_path <- "/path_to_project_folder/" # path to project folder
results_dir <- file.path(predimeth_path, "results") # path to results folder

# Data 

limma_results <- file.path(results_dir, "LimmaResults_Stage1/limma_results_3cat.R") 
# contains: High_vs_Low_topTable, Medium_vs_Low_topTable, High_vs_Medium_topTable, bonf_threshold

#~ limma_results <- file.path(results_dir, "LimmaResults_Stage1/limma_results_bin.R")
# contains: full_topTable, bonf_threshold

## Output paths ------------------------------------------

results_folder <- file.path(results_dir, paste0('VisualizingResults_Stage1', label))
dir.create(results_folder)

##########################################################

#### LOAD DATA ####

## Data ------------------------------------------------

load(limma_results)

##########################################################

cat("Data loaded. \n\nStarting results analysis. \n")
message("Data loaded. \n\nStarting results analysis: ", Sys.time())

#### SECTION 1 ####

## Subset of top hits 

# High_vs_Low Adherence
fdr_HL <- subset(High_vs_Low_topTable, abs(logFC)>0.26 & adj.P.Val<0.05)
cat("FDR HL hits : ", dim(fdr_HL)[1], "\n")
write.table(fdr_HL, file = file.path(results_folder, "High_vs_Low.txt"), sep = "\t", col.names = TRUE, row.names = FALSE)
bonf_HL <- subset(High_vs_Low_topTable, abs(logFC)>0.26 & adj.P.Val<bonf_threshold) 
cat("Bonf HL hits : ", dim(bonf_HL)[1], "\n")

# Medium_vs_Low Adherence
fdr_ML <- subset(Medium_vs_Low_topTable, abs(logFC)>0.26 & adj.P.Val<0.05)
cat("FDR ML hits : ", dim(fdr_ML)[1], "\n")
write.table(fdr_HL, file = file.path(results_folder, "Medium_vs_Low.txt"), sep = "\t", col.names = TRUE, row.names = FALSE)
bonf_ML <- subset(Medium_vs_Low_topTable, abs(logFC)>0.26 & adj.P.Val<bonf_threshold) 
cat("Bonf ML hits : ", dim(bonf_ML)[1], "\n")

# High_vs_Medium Adherence
fdr_HM <- subset(High_vs_Medium_topTable, abs(logFC)>0.26 & adj.P.Val<0.05) 
cat("FDR HM hits : ", dim(fdr_HM)[1], "\n")
bonf_HM <- subset(High_vs_Medium_topTable, abs(logFC)>0.26 & adj.P.Val<bonf_threshold) 
cat("Bonf HM hits : ", dim(bonf_HM)[1], "\n")

## Generating vectors with CpGs IDs and names (for later enrichment)

# All the CpGs probe IDs tested 
all_probes_ids_HL <- High_vs_Low_topTable$ProbeID
length(all_probes_ids_HL)
all_probes_ids_ML <- Medium_vs_Low_topTable$ProbeID
length(all_probes_ids_ML)
identical(all_probes_ids_HL, all_probes_ids_ML) # FALSE: the order is not the same,
table(all_probes_ids_HL %in% all_probes_ids_ML) # TRUE: but all probes are in both

all_probes_ids <- all_probes_ids_HL
all_probes_names <- High_vs_Low_topTable$Probe_name
cat("All Cpgs: ", length(all_probes_ids), "\n")

# High vs Low FDR hits 
fdr_HL_ids <- fdr_HL$ProbeID
fdr_HL_names <- fdr_HL$Probe_name

# Medium vs Low FDR hits 
fdr_ML_ids <- fdr_ML$ProbeID
fdr_ML_names <- fdr_ML$Probe_name

# Saving

save(all_probes_ids, all_probes_names,
     fdr_HL_ids, fdr_HL_names, 
     fdr_ML_ids, fdr_ML_names,
     file = file.path(results_folder, "hits_in_vectors.R"))

##########################################################

cat("Results analysis completed. \n\nVisualizing results. \n")
message("Results analysis completed. \n\nVisualizing results: ", Sys.time())

#### VOLCANO PLOT ####

volcano_plot_HL <- EnhancedVolcano(High_vs_Low_topTable, x = "logFC", y = "adj.P.Val",
                                lab = "", pCutoff = 0.05, FCcutoff = log2(1.2),
                                pointSize = 3) + ggplot2::labs(title = "High vs Low Adherence to Mediterranean diet") 

volcano_plot_ML <- EnhancedVolcano(Medium_vs_Low_topTable, x = "logFC", y = "adj.P.Val",
                                   lab = "", pCutoff = 0.05, FCcutoff = log2(1.2),
                                   pointSize = 3) + ggplot2::labs(title = "Medium vs Low Adherence to Mediterranean diet")

# to add bonf line: geom_hline(yintercept = -log10(bonf_threshold), linetype = "dashed", color = "coral", linewidth = 0.5)

ggsave(file.path(results_folder, "Volcano_Plot_HL.pdf"), plot = volcano_plot_HL, width = 8, height = 7)
ggsave(file.path(results_folder, "Volcano_Plot_ML.pdf"), plot = volcano_plot_ML, width = 8, height = 7)


#### MANHATTAN PLOT ####

MAplot <- function(topTable) {
  ma_data <- subset(topTable, select = c("CpGchrm", "CpGbeg", "adj.P.Val", "Probe_name"))
	ma_data$CpGchrm <- stringr::str_remove(ma_data$CpGchrm, "^chr")
 	ma_data <- ma_data[ma_data$CpGchrm != "M",]
 	table(ma_data$CpGchrm)
 	ma_data$CpGchrm <- as.numeric(ma_data$CpGchrm)
 	ma_data <- na.omit(ma_data)
	
 	ma_plot <- manhattan(ma_data, chr = "CpGchrm", bp = "CpGbeg", 
 						 p = "adj.P.Val", snp = "Probe_name",
 						 genomewideline = -log10(bonf_threshold),
 						 suggestiveline = -log10(0.05),
 						 col = c("steelblue", "coral"), cex = 0.6,
 						 main = "EWAS Manhattan Plot")
	
 	return(ma_plot)
}

png(file = file.path(results_folder, "High_vs_Low_MAplot.png"), width = 1600, height = 900, res = 150)
MAplot(High_vs_Low_topTable) 
dev.off()

png(file = file.path(results_folder, "Medium_vs_Low_MAplot.png"), width = 1600, height = 900, res = 150)
MAplot(Medium_vs_Low_topTable)
dev.off()


#### QQPLOT AND INFLATION ####


## High vs Low adherence

pvals <- High_vs_Low_topTable$P.Value
table(is.na(pvals))
length(pvals)
chisq <- qchisq(1 - pvals, df = 1)
lambda <- median(chisq) / qchisq(0.5, df = 1)
cat("Lambda High vs Low = ", lambda, "\n")

observed <- -log10(sort(pvals))
expected <- -log10(ppoints(length(pvals)))

pdf(file = file.path(results_folder, "qqplot_pvalues_HL.pdf"), width = 6, height = 6)
plot(expected, observed, pch = 20, cex = 0.6, 
     xlab = "Expected -log10(p)", ylab = "Observed -log10(p)")
abline(0, 1, col = "red")
dev.off()

## Medium vs Low adherence

pvals <- Medium_vs_Low_topTable$P.Value
table(is.na(pvals))
length(pvals)
chisq <- qchisq(1 - pvals, df = 1)
lambda <- median(chisq) / qchisq(0.5, df = 1)
cat("Lambda Medium vs Low = ", lambda, "\n")

observed <- -log10(sort(pvals))
expected <- -log10(ppoints(length(pvals)))

pdf(file = file.path(results_folder, "qqplot_pvalues_ML.pdf"), width = 6, height = 6)
plot(expected, observed, pch = 20, cex = 0.6, 
     xlab = "Expected -log10(p)", ylab = "Observed -log10(p)")
abline(0, 1, col = "red")
dev.off()


#### GENOMIC CATEGORIES OF HITS ####

df_summaryHL <- High_vs_Low_topTable %>%
		group_by(RelationToGeneCategory) %>%
			summarise(count = n()) %>%
				arrange(desc(count))

plotHL <- ggplot(df_summaryHL, aes(x = reorder(RelationToGeneCategory, -count), y = count)) +
			geom_bar(stat = "identity", fill = "steelblue") +
			theme_bw() +
			labs(x = "Genomic Category",
				 y = "Number of CpGs",
				 title = "High vs Low FDR-significant CpGs by Genomic Category") +
			theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(results_folder, "FDRhits_by_GenomicCategory_HL.pdf"), plot = plotHL)

df_summaryML <- Medium_vs_Low_topTable %>%
		group_by(RelationToGeneCategory) %>%
			summarise(count = n()) %>%
				arrange(desc(count))

plotML <- ggplot(df_summaryML, aes(x = reorder(RelationToGeneCategory, -count), y = count)) +
			geom_bar(stat = "identity", fill = "steelblue") +
			theme_bw() +
			labs(x = "Genomic Category",
				 y = "Number of CpGs",
				 title = "Medium vs Low FDR-significant CpGs by Genomic Category") +
			theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(results_folder, "FDRhits_by_GenomicCategory_ML.pdf"), plot = plotML)

##########################################################

cat("Script completed. \n")
message("Script completed: ", Sys.time())

##########################################################
