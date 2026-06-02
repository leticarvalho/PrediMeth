##########################################################
#
# Annotating EWAS Top Hits - Stage 3
#
##########################################################


# DESCRIPTION --------------------------------------------
# In this script, EWAS results will be annotated using
# missMethyl functions, to account for multi-probe bias, 
# using some annotation databases: Gene Ontology (GO), 
# Kyoto Encyclopedia of Genes and Genomes (KEGG), 
# Reactome and Wikipathways. 
# If not FDR-significant, results will be ranked by raw
# p-value, and showed by a dotplot or a barchart. 

# INPUT --------------------------------------------------
# 	EWAS results (topTable)

# OUTPUT -------------------------------------------------
# 	Gene Ontology (GO)
#   Kyoto Encyclopedia of Genes and Genomes (KEGG)
# 	Reactome
# 	Wikipathways


##########################################################

#### LOAD LIBRARIES ####

library(missMethyl) # gometh(), gsameth()
library(msigdbr) # Reactome and Wikipathways
library(clusterProfiler) # plots
library(enrichplot) # plots
library(data.table) # read files
library(ggplot2) # plots

#### RESOLVE PATHS ####

## Input paths -------------------------------------------

path_to_continuous_hits_in_vector <- "" # Rdata
path_to_3catoriginal_hits_in_vector <- ""  # Rdata
path_to_spline_hits_in_vector <- ""  # Rdata

## Output paths ------------------------------------------

results_folder <- ""
dir.create(results_folder)

##########################################################

#### LOAD DATA ####

# Categorical contrast
load(path_to_3catoriginal_hits_in_vector) # fdr_HL1_names, fdr_ML1_names, fdr_HM1_names, all_probes_ids

# Spline analysis
load(path_to_spline_hits_in_vector) # spline_fdr_cpgs, all_probes_ids

# Continuous contrast
load(path_to_continuous_hits_in_vector) # fdr_names1, all_probes_ids

##########################################################

cat("Data loaded. \n\nStarting Over-Representation Analyses (ORA). \n")
message("Data loaded. \n\nStarting Over-Representation Analyses (ORA): ", Sys.time())

#### VARIABLES #### 

## Choose CpGs that are going to pass through ORA 
sig_cpgs <- fdr_ML1_names # example: a vector with ML FDR hits 
all_cpgs <- all_probes_ids # all probes tested 

## Check if they are in the right format -> they have to be on EPIC v2 format (ending in suffix: "ProbeID")
sig_cpgs[1:5]
all_cpgs[1:5]

#### GET ANNOTATION ####

## Annotation from hg38 version of human genome
ann.hg38 <- getAnnotation(IlluminaHumanMethylationEPICv2anno.20a1.hg38)

##########################################################

cat("ORA with GO. \n")
message("ORA with GO: ", Sys.time())

#### ORA with GO ####

## ORA using gometh() for bias-correction
v2.gometh_GO <- gometh(sig.cpg = sig_cpgs, all.cpg = all_cpgs, collection = "GO", 
                       array.type = "EPIC_V2", sig.genes = TRUE)

## Select significant GO terms (by FDR or nominal p-value)
go_sig <- v2.gometh_GO[v2.gometh_GO$P.DE<0.05,]

## Create a object "enrichResult" in order to visualize results
er <- new("enrichResult", 
          result = data.frame(ID = rownames(go_sig), 
                              Description = go_sig$TERM,
                              GeneRatio = paste(go_sig$DE, go_sig$N, sep = "/"),
                              BgRatio = paste(go_sig$N,  nrow(v2.gometh_GO), sep = "/"),
                              pvalue = go_sig$P.DE,
                              p.adjust = go_sig$P.DE, # change it for FDR 
                              qvalue = go_sig$P.DE, # changeit for FDR 
                              geneID = go_sig$SigGenesInSet,
                              Count = go_sig$DE),
          pvalueCutoff = 0.05, pAdjustMethod = "BH",
          organism = "Homo sapiens", keytype = "ENTREZID", ontology = "BP",
          gene = sig_cpgs, # vector with significant CpGs
          universe = all_cpgs, # all CpGs tested
          geneSets = list(), readable = FALSE)

## Simplify redundant GO terms
er_simple <- simplify(er, cutoff = 0.7, # lower = more aggressive pruning
                      by = "p.adjust", select_fun = min, measure = "Wang")

## Visualise by dot plot
dotplot(er_simple, showCategory = 16, x = "GeneRatio", 
        color = "pvalue", font.size = 7, label_format = 30) +
  scale_color_gradient(low = "#E8593C", high = "#3B8BD4", name = "Nominal p-value") +
  ggtitle("Simplified GO ORA \nMedium vs Low FDR hits") +
  theme_minimal(base_size = 11)

## Construct GO-bias plot to check multi-probe bias correction
gometh_go_biasplot <- gometh(sig.cpg = sig_cpgs, all.cpg = all_probes_ids, collection = "GO",
                             array.type = "EPIC_V2", sig.genes = TRUE, 
                             anno = ann.hg38, plot.bias = TRUE)

## Optative: ORAs by regions 

# Promoter 
v2.gometh_GO_promoter <- gometh(sig.cpg = fdr_HL_ids, all.cpg = all_probes_ids, 
                                collection = "GO", array.type = "EPIC_V2",
                                genomic.features = c("TSS200", "TSS1500", "5'UTR", "1stExon"), # promoter regions
                                sig.genes = TRUE)

# Body
v2.gometh_GO_body <- gometh(sig.cpg = fdr_HL_ids, all.cpg = all_probes_ids, 
                            collection = "GO", array.type = "EPIC_V2",
                            genomic.features = c("Body"), # body region
                            sig.genes = TRUE)

##########################################################

cat("ORA with KEGG \n")
message("ORA with KEGG: ", Sys.time())

#### ORA with KEGG ####

## ORA using gometh() for bias-correction
v2.gometh_KEGG <- gometh(sig.cpg = sig_cpgs, all.cpg = all_cpgs, collection = "KEGG",
                         array.type = "EPIC_V2", sig.genes = TRUE)

## Select significant pathways (by FDR or nominal p-value)
kegg_sig <- v2.gometh_KEGG[v2.gometh_KEGG$P.DE<0.05,]

## Create a object "enrichResult" in order to visualize results 
er_kegg <- new("enrichResult", 
               result = data.frame(ID = rownames(kegg_sig),
                                   Description = kegg_sig$Description,  
                                   GeneRatio = paste(kegg_sig$DE, kegg_sig$N, sep = "/"),
                                   BgRatio = paste(kegg_sig$N, nrow(v2.gometh_KEGG), sep = "/"),
                                   pvalue = kegg_sig$P.DE,
                                   p.adjust = kegg_sig$P.DE, # change it for FDR 
                                   qvalue = kegg_sig$P.DE, # change it for FDR 
                                   geneID = "",
                                   Count = kegg_sig$DE),
               pvalueCutoff = 0.05, pAdjustMethod = "none",
               organism = "hsa", keytype = "ENTREZID", ontology = "KEGG",
               gene = sig_cpgs, # vector with significant CpGs
               universe = all_cpgs, # all CpGs tested
               geneSets = list(), readable = FALSE)

## Visualise by bar plot
barplot(er_kegg, showCategory = 20, # how many categories to display
        x = "GeneRatio", font.size = 10, title = "KEGG ORA: Medium vs Low") 

## Construct KEGG-bias plot to check multi-probe bias correction
gometh_kegg_biasplot <- gometh(sig.cpg = sig_cpgs, all.cpg = all_expanded, collection = "KEGG",
                               array.type = "EPIC_V2", sig.genes = TRUE, 
                               anno = ann.hg38, plot.bias = TRUE)

## Optative: ORAs by regions 

# Promoter 
v2.gometh_KEGG_promoter <- gometh(sig.cpg = fdr_HL_ids, all.cpg = all_probes_ids, 
                                  collection = "KEGG", array.type = "EPIC_V2",
                                  genomic.features = c("TSS200", "TSS1500", "5'UTR", "1stExon"), # promoter
                                  sig.genes = TRUE)

# Body
v2.gometh_KEGG_body <- gometh(sig.cpg = fdr_HL_ids, all.cpg = all_probes_ids, 
                              collection = "KEGG", array.type = "EPIC_V2",
                              genomic.features = c("Body"), # body
                              sig.genes = TRUE)

##########################################################

cat("ORA with Reactome \n")
message("ORA with Reactome: ", Sys.time())

#### ORA with Reactome ####

## Via msigdbr package (category C2, subcategory CP:REACTOME)
msig_reactome_raw <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME") %>%
  dplyr::select(gs_name, entrez_gene) %>% group_by(gs_name) %>%
  summarise(genes = list(as.character(entrez_gene)), .groups = "drop")
reactome_all <- setNames(msig_reactome_raw$genes, msig_reactome_raw$gs_name)

## ORA using gsameth() for bias-correction
gsa_reactome_all <- gsameth(sig.cpg = sig_cpgs, # significant CpGs
                            all.cpg = all_cpgs, # all CpGs tested
                            collection = reactome_all, # reactome
                            array.type = "EPIC_V2", anno = ann.hg38,
                            prior.prob = TRUE, fract.counts = TRUE) %>%
  tibble::rownames_to_column("GeneSet") %>%
  mutate(FDR = p.adjust(P.DE, method = "BH")) %>% arrange(P.DE)

## Select significant pathways (by FDR or nominal p-value)
reactome_sig <- gsa_reactome_all[gsa_reactome_all$P.DE < 0.05, ]

## Clean pathway names
reactome_sig$GeneSet_clean <- gsub("^REACTOME_", "", reactome_sig$GeneSet)
reactome_sig$GeneSet_clean <- gsub("_", " ", reactome_sig$GeneSet_clean)
reactome_sig$GeneSet_clean <- tolower(reactome_sig$GeneSet_clean)
reactome_sig$GeneSet_clean <- stringr::str_to_sentence(reactome_sig$GeneSet_clean)

## Create a object "enrichResult" in order to visualize results
er_reactome <- new("enrichResult",
                   result = data.frame(
                     ID = reactome_sig$GeneSet,
                     Description = reactome_sig$GeneSet_clean,
                     GeneRatio = paste(reactome_sig$DE, reactome_sig$N, sep = "/"),
                     BgRatio = paste(reactome_sig$N, nrow(gsa_reactome_all), sep = "/"),
                     pvalue = reactome_sig$P.DE,
                     p.adjust = reactome_sig$P.DE, qvalue = reactome_sig$P.DE, # change them for FDR 
                     geneID = "", Count = reactome_sig$DE),
                   pvalueCutoff = 0.05, pAdjustMethod = "BH",
                   organism = "Homo sapiens", keytype = "ENTREZID", ontology = "Reactome",
                   gene = sig_cpgs,
                   universe = all_cpgs,
                   geneSets = list(), readable = FALSE)

## Visualise by dot plot
dotplot(er_reactome, showCategory = 20, x = "GeneRatio", color = "pvalue",
        font.size = 10, label_format = 45) +
  scale_color_gradient(low = "#E8593C", high = "#3B8BD4", name = "Raw p-value") +
  ggtitle("ORA - Reactome: \nMedium vs Low") +
  theme_minimal(base_size = 11)

##########################################################

cat("ORA with Wikipathways \n")
message("ORA with Wikipathways: ", Sys.time())

#### ORA with Wikipathways ####

## Via msigdbr package (category C2, subcategory CP:WIKIPATHWAYS)
msig_wikipath <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:WIKIPATHWAYS") %>%
  dplyr::select(gs_name, entrez_gene) %>%
  group_by(gs_name) %>% summarise(genes = list(as.character(entrez_gene)), .groups = "drop")
wikipathways <- setNames(msig_wikipath$genes, msig_wikipath$gs_name)

## ORA using gsameth() for bias-correction
gsa_wikipathways <- gsameth(sig.cpg = sig_cpgs, # significant CpGs
                            all.cpg = all_cpgs, # CpGs tested
                            collection = wikipathways, # wikipathways
                            array.type = "EPIC_V2", anno = ann.hg38,
                            prior.prob = TRUE, fract.counts = TRUE) %>%
  tibble::rownames_to_column("GeneSet") %>% 
  mutate(FDR = p.adjust(P.DE, method = "BH")) %>% arrange(P.DE)

## Select significant pathways (by FDR or nominal p-value)
wikipathways_sig <- gsa_wikipathways[gsa_wikipathways$P.DE < 0.05, ]

## Clean pathway names
wikipathways_sig$GeneSet_clean <- gsub("^WP_", "", wikipathways_sig$GeneSet)
wikipathways_sig$GeneSet_clean <- gsub("_", " ", wikipathways_sig$GeneSet_clean)
wikipathways_sig$GeneSet_clean <- tolower(wikipathways_sig$GeneSet_clean)
wikipathways_sig$GeneSet_clean <- stringr::str_to_sentence(wikipathways_sig$GeneSet_clean)

## Create a object "enrichResult" in order to visualize results
er_wikipathways <- new("enrichResult",
                       result = data.frame(ID = wikipathways_sig$GeneSet, 
                                           Description = wikipathways_sig$GeneSet_clean,
                                           GeneRatio = paste(wikipathways_sig$DE, wikipathways_sig$N, sep = "/"),
                                           BgRatio = paste(wikipathways_sig$N, nrow(gsa_wikipathways), sep = "/"),
                                           pvalue = wikipathways_sig$P.DE,
                                           p.adjust = wikipathways_sig$P.DE, qvalue = wikipathways_sig$P.DE, # change them for FDR 
                                           geneID = "", Count = wikipathways_sig$DE),
                       pvalueCutoff  = 0.05, pAdjustMethod = "BH",
                       organism = "Homo sapiens", keytype = "ENTREZID", ontology = "Reactome",
                       gene = sig_cpgs, # significant CpGs
                       universe = all_cpgs, # all Cpgs
                       geneSets = list(), readable = FALSE)

## Visualise by dot plot
dotplot(er_wikipathways, showCategory = 20, x = "GeneRatio", color = "pvalue",
        font.size = 10, label_format = 45) +
  scale_color_gradient(low = "#E8593C", high = "#3B8BD4", name = "Raw p-value") +
  ggtitle("ORA - Wikipathways: \nMedium vs Low") +
  theme_minimal(base_size = 11)

##########################################################

cat("Saving outputs. \n")
message("Saving outputs: ", Sys.time())

save(v2.gometh_GO, v2.gometh_KEGG, gsa_reactome_all, gsa_wikipathways, 
     file = file.path(results_folder, "ORAnalyses.R"))

cat("Script completed. \n")
message("Script completed: ", Sys.time())

##########################################################
