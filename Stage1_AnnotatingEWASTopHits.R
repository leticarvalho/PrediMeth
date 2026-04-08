##########################################################
#
# Annotating EWAS Top Hits - Stage 1
#
##########################################################


# DESCRIPTION --------------------------------------------
# In this script, EWAS results will be annotated using
# missMethyl functions, to account for multi-probe bias, 
# using some annotation databases: Gene Ontology (GO), 
# Kyoto Encyclopedia of Genes and Genomes (KEGG), 
# Reactome, and alternative gene sets. 
# If not FDR-significant, results will be ranked by raw
# p-value, and showed by a dotplot or a barchart. 

# INPUT --------------------------------------------------
# 	EWAS results (CpG vectors)

# OUTPUT -------------------------------------------------
# 	Gene Ontology (GO)
#   Kyoto Encyclopedia of Genes and Genomes (KEGG)
# 	Reactome
# 	Alternative gene sets: 
#       MSigDB Hallmark, ImmuneSigDB, Wikipathways

# Obs.: For enrichment analyses with alternative gene sets, 
#       some help for code was used, from Claude. 


##########################################################

#### LOAD LIBRARIES ####

library(missMethyl)
library(msigdbr)
library(clusterProfiler)
library(enrichplot)

#### RESOLVE PATHS ####

## Define label ------------------------------------------

label <- "_3cat" # Probe ID
#label <- "_bin" # Probe ID

## Paths -------------------------------------------

predimeth_path <- "/path_to_project_folder/" # path to project folder
results_dir <- file.path(predimeth_path, "results") # path to results folder

results_folder <- file.path(results_dir, paste0('VisualizingResults_Stage1', label))

CpGvectors <- file.path(results_folder, "hits_in_vectors.R")

##########################################################

#### LOAD DATA ####

## Data 1 ------------------------------------------------

load(CpGvectors)

# Annotation  ------------------------------------------------

ann.hg38 <- getAnnotation(IlluminaHumanMethylationEPICv2anno.20a1.hg38)

##########################################################

cat("Data loaded. \n\nStarting section 1. \n")
message("Data loaded. \n\nStarting section 1: ", Sys.time())


### Enrichment Analysis with Gene Ontology ###

v2.gometh_GO <- gometh(sig.cpg = fdr_HL_ids, all.cpg = all_probes_ids, collection = "GO",
                       array.type = "EPIC_V2", sig.genes = TRUE)

# no FDR significant
go_sig <- v2.gometh_GO[v2.gometh_GO$P.DE<0.05,]

# Creating a enrichResult object (in order to do a dot plot later)
er <- new("enrichResult",
          result = data.frame(ID = rownames(go_sig),
                              Description = go_sig$TERM,
                              GeneRatio = paste(go_sig$DE, go_sig$N, sep = "/"),
                              BgRatio = paste(go_sig$N,  nrow(v2.gometh_GO), sep = "/"),
                              pvalue = go_sig$P.DE,
                              p.adjust = go_sig$P.DE, # because there are no FDR-significant hit
                              qvalue = go_sig$P.DE,
                              geneID = go_sig$SigGenesInSet,
                              Count = go_sig$DE),
          pvalueCutoff  = 0.05, pAdjustMethod = "BH", organism = "Homo sapiens",
          keytype = "ENTREZID", ontology = "BP", # BP, MF or CC
          gene = fdr_HL_ids, universe = all_probes_ids,
          geneSets = list(), readable = FALSE)

# Removing semantically redundant GO terms
er_simple <- simplify(er, cutoff = 0.7, by = "p.adjust", select_fun  = min, measure = "Wang")

# Creating a Dotplot
dotplot(er_simple, showCategory = 20, x = "GeneRatio", color = "p.adjust", 
        font.size = 8, label_format = 30) +
  scale_color_gradient(low = "#E8593C", high = "#3B8BD4", name = "Raw p-value") +
  ggtitle("Simplified GO enrichment (gometh)") + theme_minimal(base_size = 11)


# GO-bias plot #

gometh_go_biasplot <- gometh(sig.cpg = fdr_HL_ids, all.cpg = all_probes_ids, collection = "GO",
                             array.type = "EPIC_V2", sig.genes = TRUE, 
                             anno = ann.hg38, plot.bias = TRUE)

ggsave(file.path(results_folder, "GO_biasplot.pdf"), plot = gometh_go_biasplot, width = 8, height = 7)

# Only Promoter #

v2.gometh_GO_promoter <- gometh(sig.cpg = fdr_HL_ids, all.cpg = all_probes_ids, collection = "GO",
                                array.type = "EPIC_V2",
                                genomic.features = c("TSS200", "TSS1500", "5'UTR", "1stExon"),
                                sig.genes = TRUE)

# Only Body #

v2.gometh_GO_body <- gometh(sig.cpg = fdr_HL_ids, all.cpg = all_probes_ids, collection = "GO",
                            array.type = "EPIC_V2",
                            genomic.features = c("Body"),
                            sig.genes = TRUE)

# Only Hypomethylated #

dim(fdr_HL[fdr_HL$logFC<0,]) # FDR
v2.hypo <- fdr_HL$ProbeID[fdr_HL$logFC < 0]

v2.gometh_GO_hypo <- gometh(sig.cpg = v2.hypo, all.cpg = all_probes_ids, collection = "GO",
                            array.type = "EPIC_V2", sig.genes = TRUE)

# Only Hypermethylated #

dim(fdr_HL[fdr_HL$logFC>0,]) # FDR
v2.hyper <- fdr_HL$ProbeID[fdr_HL$logFC > 0]

v2.gometh_GO_hyper <- gometh(sig.cpg = v2.hyper, all.cpg = all_probes_ids, collection = "GO",
                            array.type = "EPIC_V2", sig.genes = TRUE)





### Enrichment Analysis with KEGG ###

v2.gometh_KEGG <- gometh(sig.cpg = fdr_HL_ids, all.cpg = all_probes_ids, collection = "KEGG",
                       array.type = "EPIC_V2", sig.genes = TRUE)
# no FDR significant
kegg_sig <- v2.gometh_KEGG[v2.gometh_KEGG$P.DE<0.05,]

# Creating a enrichResult object (in order to do a barplot later)
er_kegg <- new("enrichResult",
               result = data.frame(ID = rownames(kegg_sig), 
                                   Description = kegg_sig$Description,
                                   GeneRatio = paste(kegg_sig$DE, kegg_sig$N, sep = "/"),
                                   BgRatio = paste(kegg_sig$N, nrow(v2.gometh_KEGG), sep = "/"),
                                   pvalue = kegg_sig$P.DE, 
                                   p.adjust = kegg_sig$P.DE, # because there are no FDR-significant hit
                                   qvalue = kegg_sig$P.DE,
                                   geneID = "", 
                                   Count = kegg_sig$DE),
               pvalueCutoff  = 0.05, pAdjustMethod = "none", organism = "hsa",
               keytype = "ENTREZID", ontology = "KEGG",
               gene = fdr_HL_ids, universe = all_probes_ids,
               geneSets = list(), readable = FALSE)

# Barplot
barplot(er_kegg, showCategory = 20, x = "GeneRatio", font.size = 10) + theme_minimal()

# KEGG-bias plot #

gometh_kegg_biasplot <- gometh(sig.cpg = fdr_HL_ids, all.cpg = all_probes_ids, collection = "KEGG",
                             array.type = "EPIC_V2", sig.genes = TRUE, 
                             anno = ann.hg38, plot.bias = TRUE)

ggsave(file.path(results_folder, "KEGG_biasplot.pdf"), plot = gometh_kegg_biasplot, width = 8, height = 7)

# Only Promoter # 

v2.gometh_KEGG_promoter <- gometh(sig.cpg = fdr_HL_ids, all.cpg = all_probes_ids, 
                                  collection = "KEGG", array.type = "EPIC_V2",
                                  genomic.features = c("TSS200", "TSS1500", "5'UTR", "1stExon"),
                                  sig.genes = TRUE)

# Only Body # 

v2.gometh_KEGG_body <- gometh(sig.cpg = fdr_HL_ids, all.cpg = all_probes_ids, 
                              collection = "KEGG",
                            array.type = "EPIC_V2",
                            genomic.features = c("Body"),
                            sig.genes = TRUE)

# Only Hypomethylated # 

v2.gometh_KEGG_hypo <- gometh(sig.cpg = v2.hypo, all.cpg = all_probes_ids, 
                              collection = "KEGG",
                            array.type = "EPIC_V2", sig.genes = TRUE)

# Only Hypermethylated #

v2.gometh_KEGG_hyper <- gometh(sig.cpg = v2.hyper, all.cpg = all_probes_ids, 
                              collection = "KEGG",
                              array.type = "EPIC_V2", sig.genes = TRUE)



### Enrichment Analysis with MSigDB Hallmark gene sets (bias-corrected via gsameth) ###

msig_h <- msigdbr(species = "Homo sapiens", category = "H") %>%
  dplyr::select(gs_name, entrez_gene) %>%
  group_by(gs_name) %>%
  summarise(genes = list(as.character(entrez_gene)), .groups = "drop")

hallmark_sets <- setNames(msig_h$genes, msig_h$gs_name)
  
gsa_hallmarks <- gsameth(sig.cpg = fdr_HL_ids, all.cpg = all_probes_ids,
                         collection = hallmark_sets, array.type = "EPIC_V2",
                         anno = ann.hg38, prior.prob = TRUE, fract.counts = TRUE) %>%
  tibble::rownames_to_column("GeneSet") %>%
  mutate(FDR = p.adjust(P.DE, method = "BH")) %>%
  arrange(P.DE)


### Enrichment with Reactome (via msigdbr (C2 / CP:REACTOME)) ###

msig_reactome_raw <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME") %>%
  dplyr::select(gs_name, entrez_gene) %>% group_by(gs_name) %>% 
  summarise(genes = list(as.character(entrez_gene)), .groups = "drop")

reactome_all <- setNames(msig_reactome_raw$genes, msig_reactome_raw$gs_name)

gsa_reactome_all <- gsameth(sig.cpg = fdr_HL_ids, all.cpg = all_probes_ids,
                            collection = reactome_all, array.type = "EPIC_V2",
                            anno = ann.hg38, prior.prob = TRUE, fract.counts = TRUE) %>%
  tibble::rownames_to_column("GeneSet") %>%
  mutate(FDR = p.adjust(P.DE, method = "BH")) %>%
  arrange(P.DE)

# Filter significant sets
reactome_sig <- gsa_reactome_all[gsa_reactome_all$P.DE < 0.05, ]

# Clean up pathway names (MSigDB uses REACTOME_PATHWAY_NAME format)
reactome_sig$GeneSet_clean <- gsub("^REACTOME_", "", reactome_sig$GeneSet)
reactome_sig$GeneSet_clean <- gsub("_", " ", reactome_sig$GeneSet_clean)
reactome_sig$GeneSet_clean <- tolower(reactome_sig$GeneSet_clean)
reactome_sig$GeneSet_clean <- stringr::str_to_sentence(reactome_sig$GeneSet_clean)

# Coerce to enrichResult
er_reactome <- new("enrichResult",
                   result = data.frame(ID = reactome_sig$GeneSet,
                                       Description = reactome_sig$GeneSet_clean,
                                       GeneRatio = paste(reactome_sig$DE, reactome_sig$N, sep = "/"),
                                       BgRatio = paste(reactome_sig$N, nrow(gsa_reactome_all), sep = "/"),
                                       pvalue = reactome_sig$P.DE,
                                       p.adjust = reactome_sig$P.DE, # because there are no FDR-significant hit
                                       qvalue = reactome_sig$P.DE,
                                       geneID = "",
                                       Count = reactome_sig$DE),
                   pvalueCutoff  = 0.05, pAdjustMethod = "BH", organism = "Homo sapiens", 
                   keytype = "ENTREZID", ontology = "Reactome",
                   gene = fdr_HL_ids, universe = all_probes_ids,
                   geneSets = list(), readable = FALSE)

# Dot plot
dotplot(er_reactome, showCategory = 20, x = "GeneRatio", color = "p.adjust", 
        font.size = 10, label_format = 45) +
  scale_color_gradient(low = "#E8593C", high = "#3B8BD4", name = "Raw p-value") +
  ggtitle("Reactome enrichment (gsameth, bias-corrected)") +
  theme_minimal(base_size = 11)
  

### Enrichment Analysis with ImmuneSigDB gene sets from MSigDB (bias-corrected via gsameth) ###

msig_immune <- msigdbr(species = "Homo sapiens", category = "C7", subcategory = "IMMUNESIGDB") %>%
  dplyr::select(gs_name, entrez_gene) %>%
  group_by(gs_name) %>%
  summarise(genes = list(as.character(entrez_gene)), .groups = "drop")

immunesig_sets <- setNames(msig_immune$genes, msig_immune$gs_name)
  
gsa_immunesig <- gsameth(sig.cpg = fdr_HL_ids, all.cpg = all_probes_ids, 
                         collection = immunesig_sets, array.type = "EPIC_V2",
                         anno = ann.hg38, prior.prob = TRUE, fract.counts = TRUE) %>%
  tibble::rownames_to_column("GeneSet") %>%
  mutate(FDR = p.adjust(P.DE, method = "BH")) %>%
  arrange(P.DE)

### Enrichment with Wikipathways (via msigdbr (C2 / CP:WIKIPATHWAYS)) ###

msig_wikipath <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:WIKIPATHWAYS") %>%
  dplyr::select(gs_name, entrez_gene) %>% group_by(gs_name) %>%
  summarise(genes = list(as.character(entrez_gene)), .groups = "drop")

wikipathways <- setNames(msig_wikipath$genes, msig_wikipath$gs_name)

gsa_wikipathways <- gsameth(sig.cpg = fdr_HL_ids, all.cpg = all_probes_ids,
                            collection = wikipathways, array.type = "EPIC_V2",
                            anno = ann.hg38, prior.prob = TRUE, fract.counts = TRUE) %>%
  tibble::rownames_to_column("GeneSet") %>%
  mutate(FDR = p.adjust(P.DE, method = "BH")) %>%
  arrange(P.DE)

##########################################################

#### SAVING OUTPUTS ####

## R object

save(v2.gometh_GO, v2.gometh_GO_promoter, v2.gometh_GO_body, v2.gometh_GO_hypo, v2.gometh_GO_hyper,
     file = file.path(results_folder, paste0("enrich_GO", label, ".R")))

save(v2.gometh_KEGG, v2.gometh_KEGG_promoter, v2.gometh_KEGG_body, v2.gometh_KEGG_hypo, v2.gometh_KEGG_hyper,
     file = file.path(results_folder, paste0("enrich_KEGG", label, ".R")))

save(gsa_hallmarks, gsa_reactome_all, gsa_immunesig, gsa_wikipathways,
     file = file.path(results_folder, paste0("enrich_others", label, ".R")))

##########################################################

cat("Script completed. \n")
message("Script completed: ", Sys.time())

##########################################################
