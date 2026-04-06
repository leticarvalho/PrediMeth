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

#~ label <- "_3cat_rmXY_M4" # 3 categories (low, medium, high adherence), rmXY, model 4 (full)
label <- "_3cat_rmCH_XY_SNP_M4" # 3 categories (low, medium, high adherence), rmCH_XY_SNP, model 4 (full)


## Input paths -------------------------------------------

predimeth_path <- "/imppc/labs/dnalab/share/PrediMeth"
results_dir <- file.path(predimeth_path, "results")

## Output paths ------------------------------------------

results_folder <- file.path(results_dir, paste0('VisualizingResults_Stage1', label))
#~ dir.create(results_folder)

##########################################################

#### LOAD DATA ####

## Data 1 ------------------------------------------------

#~ load("/imppc/labs/dnalab/share/PrediMeth/results/LimmaResults_Stage1_3cat_rmXY_M4/limma_results_3cat_rmXY_M4.R")
load("/imppc/labs/dnalab/share/PrediMeth/results/LimmaResults_Stage1_3cat_rmCH_XY_SNP_M4/limma_results_3cat_rmCH_XY_SNP_M4.R")


##########################################################

#~ cat("Data loaded. \n\nStarting results analysis. \n")
#~ message("Data loaded. \n\nStarting results analysis: ", Sys.time())

#~ #### SECTION 1 ####

#~ ## Subset of top hits 

#~ # High_vs_Low Adherence
#~ fdr_HL <- subset(High_vs_Low_topTable, abs(logFC)>0.26 & adj.P.Val<0.05)
#~ cat("FDR HL hits : ", dim(fdr_HL)[1], "\n")
#~ write.table(fdr_HL, file = file.path(results_folder, "High_vs_Low.txt"), sep = "\t", col.names = TRUE, row.names = FALSE)
#~ bonf_HL <- subset(High_vs_Low_topTable, abs(logFC)>0.26 & adj.P.Val<bonf_threshold) # zero
#~ cat("Bonf HL hits : ", dim(bonf_HL)[1], "\n")

#~ # Medium_vs_Low Adherence
#~ fdr_ML <- subset(Medium_vs_Low_topTable, abs(logFC)>0.26 & adj.P.Val<0.05)
#~ cat("FDR ML hits : ", dim(fdr_ML)[1], "\n")
#~ write.table(fdr_HL, file = file.path(results_folder, "Medium_vs_Low.txt"), sep = "\t", col.names = TRUE, row.names = FALSE)
#~ bonf_ML <- subset(Medium_vs_Low_topTable, abs(logFC)>0.26 & adj.P.Val<bonf_threshold) # zero
#~ cat("Bonf ML hits : ", dim(bonf_ML)[1], "\n")

#~ # High_vs_Medium Adherence
#~ fdr_HM <- subset(High_vs_Medium_topTable, abs(logFC)>0.26 & adj.P.Val<0.05) # zero
#~ cat("FDR HM hits : ", dim(fdr_HM)[1], "\n")
#~ bonf_HM <- subset(High_vs_Medium_topTable, abs(logFC)>0.26 & adj.P.Val<bonf_threshold) # zero
#~ cat("Bonf HM hits : ", dim(bonf_HM)[1], "\n")

#~ ## Generating vectors with CpGs IDs and names (for later enrichment)

#~ # All the CpGs probe IDs tested 
#~ all_probes_ids_HL <- High_vs_Low_topTable$ProbeID
#~ length(all_probes_ids_HL)
#~ all_probes_ids_ML <- Medium_vs_Low_topTable$ProbeID
#~ length(all_probes_ids_ML)
#~ identical(all_probes_ids_HL, all_probes_ids_ML) # FALSE: the order is not the same,
#~ table(all_probes_ids_HL %in% all_probes_ids_ML) # TRUE: but all probes are in both

#~ all_probes_ids <- all_probes_ids_HL
#~ all_probes_names <- High_vs_Low_topTable$Probe_name
#~ cat("All Cpgs: ", length(all_probes_ids), "\n")

#~ # High vs Low FDR hits 
#~ fdr_HL_ids <- fdr_HL$ProbeID
#~ fdr_HL_names <- fdr_HL$Probe_name

#~ # Medium vs Low FDR hits 
#~ fdr_ML_ids <- fdr_ML$ProbeID
#~ fdr_ML_names <- fdr_ML$Probe_name

#~ # Saving

#~ save(all_probes_ids, all_probes_names,
#~      fdr_HL_ids, fdr_HL_names, 
#~      fdr_ML_ids, fdr_ML_names,
#~      file = file.path(results_folder, "hits_in_vectors.R"))


########### Comparing results between rmXY and rmCH_XY_SNP

#~ load("~/Desktop/TFM/results/Categories_Results_rmXY/top_hits_vectors.R")
#~ rmxyhits_HL <- data.table::fread("~/Desktop/TFM/results/Categories_Results_rmXY/High_vs_Low.txt")
#~ rmxyhits_ML <- data.table::fread("~/Desktop/TFM/results/Categories_Results_rmXY/Medium_vs_Low.txt")

#~ rmSNPhits_HL <- fdr_HL
#~ rmSNPhits_ML <- fdr_ML

# High vs Low

#~ bothHL <- intersect(rmSNPhits_HL$ProbeID, fdr_HL_ids) # 49 

#~ onlySNP_HL <- setdiff(rmSNPhits_HL$ProbeID, fdr_HL_ids)
#~ onlySNP_HL 
#~ # cg17539962_BC21, chr1, C1orf159, 5'UTR, mRNA, Shelf, 15_Quies, LINE, logFC 0.288, delta beta NA (pq??)
#~ rmSNPhits_HL[rmSNPhits_HL$ProbeID==onlySNP_HL,]

#~ onlyXY_HL <- setdiff(fdr_HL_ids, rmSNPhits_HL$ProbeID)
#~ onlyXY_HL
#~ # cg04833433_BC21, chr3, ENSG00000288703, Body, lncRNA, OpenSea, 15_Quies, rmsk1=DNA, logFC 0.26, delta beta 0.036
#~ # cg15386636_TC21, chr10, CPEB3;NHP2P1, 5'UTR, mRNA, OpenSea, 15_Quies, logFC 0.397, delta beta 0.033
#~ rmxyhits_HL[rmxyhits_HL$ProbeID==onlyXY_HL,]

#~ # Medium vs Low

#~ bothML <- intersect(rmSNPhits_ML$ProbeID, fdr_ML_ids) # 31 

#~ onlySNP_ML <- setdiff(rmSNPhits_ML$ProbeID, fdr_ML_ids)
#~ onlySNP_ML # ZERO

#~ onlyXY_ML <- setdiff(fdr_ML_ids, rmSNPhits_ML$ProbeID)
#~ onlyXY_ML 
#~ # cg15386636_TC21, chr10, CPEB3;NHP2P1, 5'UTR, mRNA, OpenSea, 15_Quies, logFC 0.391, delta beta 0.032
#~ rmxyhits_ML[rmxyhits_ML$ProbeID==onlyXY_ML,]


#~ ####### Comparing results between High vs Low, and Medium vs Low

#~ bothXY <- intersect(rmxyhits_HL$ProbeID, rmxyhits_ML$ProbeID) # 26

#~ bothSNP <- intersect(rmSNPhits_HL$ProbeID, rmSNPhits_ML$ProbeID) # 26

#~ both <- intersect(bothSNP, bothXY) # 25 CpGs are in all results (HL and ML from rmXY and rmCH_XY_SNP)

#~ both_topTable <- rmxyhits_HL[rmxyhits_HL$ProbeID %in% both,]


#~ ### medium menys columnas

#~ MLresum <- rmxyhits_ML[,c("CpGchrm", "genesUniq", "logFC", "delta_beta")]
#~ names(rmxyhits_ML)

##########################################################

#~ cat("Results analysis completed. \n\nVisualizing results. \n")
#~ message("Results analysis completed. \n\nVisualizing results: ", Sys.time())

#~ #### VOLCANO PLOT ####

#####to add bonf line: geom_hline(yintercept = -log10(bonf_threshold), linetype = "dashed", color = "coral", linewidth = 0.5)

#~ volcano_plot_HL <- EnhancedVolcano(High_vs_Low_topTable, x = "logFC", y = "adj.P.Val",
#~                                 lab = "", pCutoff = 0.05, FCcutoff = log2(1.2),
#~                                 pointSize = 3) + ggplot2::labs(title = "High vs Low Adherence to Mediterranean diet") 

#~ volcano_plot_ML <- EnhancedVolcano(Medium_vs_Low_topTable, x = "logFC", y = "adj.P.Val",
#~                                    lab = "", pCutoff = 0.05, FCcutoff = log2(1.2),
#~                                    pointSize = 3) + ggplot2::labs(title = "Medium vs Low Adherence to Mediterranean diet")

#~ ggsave(file.path(results_folder, "Volcano_Plot_HL.pdf"), plot = volcano_plot_HL, width = 8, height = 7)
#~ ggsave(file.path(results_folder, "Volcano_Plot_ML.pdf"), plot = volcano_plot_ML, width = 8, height = 7)

#~ #### MANHATTAN PLOT ####

#~ MAplot <- function(topTable) {
#~ 	ma_data <- subset(topTable, select = c("CpGchrm", "CpGbeg", "adj.P.Val", "Probe_name"))
	###ma_data$CpGchrm <- stringr::str_remove(ma_data$CpGchrm, "chr*")
#~ 	ma_data$CpGchrm <- stringr::str_remove(ma_data$CpGchrm, "^chr")
#~ 	ma_data <- ma_data[ma_data$CpGchrm != "M",]
#~ 	table(ma_data$CpGchrm)
#~ 	ma_data$CpGchrm <- as.numeric(ma_data$CpGchrm)
#~ 	ma_data <- na.omit(ma_data)
	
#~ 	ma_plot <- manhattan(ma_data, chr = "CpGchrm", bp = "CpGbeg", 
#~ 						 p = "adj.P.Val", snp = "Probe_name",
#~ 						 genomewideline = -log10(bonf_threshold),
#~ 						 suggestiveline = -log10(0.05),
#~ 						 col = c("steelblue", "coral"), cex = 0.6,
#~ 						 main = "EWAS Manhattan Plot")
	
#~ 	return(ma_plot)
#~ }

#~ png(file = file.path(results_folder, "High_vs_Low_MAplot.png"), width = 1600, height = 900, res = 150)
#~ MAplot(High_vs_Low_topTable) 
#~ dev.off()

#~ png(file = file.path(results_folder, "Medium_vs_Low_MAplot.png"), width = 1600, height = 900, res = 150)
#~ MAplot(Medium_vs_Low_topTable)
#~ dev.off()


#### QQPLOT AND INFLATION ####

#~ measureInflation <- function(topTable, label) {
#~ 	pvals <- topTable$P.Value
#~ 	table(is.na(pvals))
#~ 	length(pvals)
#~ 	chisq <- qchisq(1 - pvals, df = 1)
#~ 	lambda <- median(chisq) / qchisq(0.5, df = 1)
#~ 	lambda
	
#~ 	observed <- -log10(sort(pvals))
#~ 	expected <- -log10(ppoints(length(pvals)))
	
#~ 	pdf(file = file.path(results_folder, paste0("qqplot_pvalues_", label, ".pdf")), width = 6, height = 6)
#~ 	plot(expected, observed, pch = 20, cex = 0.6, 
#~ 		 xlab = "Expected -log10(p)", ylab = "Observed -log10(p)")
#~     abline(0, 1, col = "red")
#~     dev.off()
#~ }

#~ measureInflation(High_vs_Low_topTable, "HL") # High vs Low
#~ measureInflation(Medium_vs_Low_topTable, "ML") # Medium vs Low

#~ pvals <- High_vs_Low_topTable$P.Value
#~ table(is.na(pvals))
#~ length(pvals)
#~ chisq <- qchisq(1 - pvals, df = 1)
#~ lambda <- median(chisq) / qchisq(0.5, df = 1)
#~ cat("Lambda High vs Low = ", lambda, "\n")

#~ pvals <- Medium_vs_Low_topTable$P.Value
#~ table(is.na(pvals))
#~ length(pvals)
#~ chisq <- qchisq(1 - pvals, df = 1)
#~ lambda <- median(chisq) / qchisq(0.5, df = 1)
#~ cat("Lambda Medium vs Low = ", lambda, "\n")


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

ggsave(file.path(results_folder, paste0("FDRhits_by_GenomicCategory_ML.pdf")), plot = plotML)

##########################################################

cat("Script completed. \n")
message("Script completed: ", Sys.time())

##########################################################
