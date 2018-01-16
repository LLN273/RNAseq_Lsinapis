##################################################################################################################################################################### 
# 
# Script used in:
#
# "Gene expression profiling across ontogenetic stages in wood white (Leptidea sinapis) reveals pathways linked to butterfly diapause regulation."
# Luis Leal, Venkat Talla, Thomas Källman, Magne Friberg, Christer Wiklund, Vlad Dincă, Roger Vila, Niclas Backström
# Mol. Eco. (2018)
#
#
##################################################################################################################################################################### 
#  
# Luis Leal
# Uppsala University, Uppsala, Sweden, 2017
#
##################################################################################################################################################################### 


################################################################################################### 
#
# R/3.4.0 script used to carry out gene enrichment analysis using the TopGO R Bioconductor package 
#
# Dependencies: DESeq2_OUTPUT_sva.csv >> DESeq2 results file
#               OUTPUT_trinity_to_GO.txt >> two-column file listing each Trinity gene and its associated GO terms 
#               
#
###################################################################################################



########## Libraries

##Install TopGO (you only need to do this once)
#source("http://bioconductor.org/biocLite.R")
#biocLite("topGO")


#load libraries
library(topGO)





########## Clear all states and remove all plots
rm(list=ls(all=TRUE))
if(!is.null(dev.list())) dev.off()





########## input files and paths

### Path to folder containing DESeq2 results file
base_deseq <- "/DESeq2"
list.files(base_deseq)            

### Name of DESeq2 results file
###
###             IMPORTANT>> MUST EDIT DESeq2 file and add "gene_ID"	to header
###
deseq_file <- "DESeq2_OUTPUT_sva.csv"


### Import DESeq2 results file
deseq_data <- read.table(file.path(base_deseq, deseq_file), header = TRUE)
head(deseq_data, 10)

### Path to folder containing Trinotate export file (two-column file listing each Trinity gene and its associated GO terms)
base_trinotate <- "/Trinotate"
list.files(base_trinotate)            

### Two-column file listing trinity gene names and associated Go terms (IN CSV FORMAT) 
##    NOTE:  headers should be as follows:
##           "gene_id" <tab> "gene_ontology"
##
trinotate_file <- "OUTPUT_trinity_to_GO.txt"

### import Trinotate export file
annot <- read.table(file.path(base_trinotate, trinotate_file), header = TRUE)
head(annot, 50)




########## Define gene universe and the set of interesting genes

## genes of interest (p-adj<0.5)
goi_aux <- deseq_data[complete.cases(deseq_data), ]           # filter-out genes with p-adj=N/A
goi_aux <- goi_aux[goi_aux[, "padj"] < 0.05,]                 # keep genes with p-adj < 0.05
goi_aux <- goi_aux[abs(goi_aux[, "log2FoldChange"]) > 1,]     # keep genes with log2 fold change > 1
goi_aux <- goi_aux[abs(goi_aux[, "baseMean"]) > 10,]          # keep genes with baseMean > 10
myInterestingGenes <-  as.vector(goi_aux$gene_ID)             # genes of interest
head(myInterestingGenes,10)
Numb_genes_of_interest <- length(myInterestingGenes)
Numb_genes_of_interest

## gene universe:  genes with (baseMean >= 6)  +  myInterestingGenes
guniv_aux <- deseq_data[deseq_data[, "baseMean"] >= 6,]  # keep genes with baseMean >= 6
guniv_names <- as.vector(guniv_aux$gene_ID)              # gene names
head(guniv_names,10)
geneNames <- c(myInterestingGenes, guniv_names)          # combine the two genes lists 
geneNames <-unique(geneNames)                            # remove duplicate entries
Numb_background_genes <- length(geneNames)
Numb_background_genes

## Create a factor - geneList - which indicates for every gene in our universe 
#  (union of background and DE-genes), whether it is differentially expressed or not.
geneList <- factor(as.integer(geneNames %in% myInterestingGenes))
names(geneList) <- geneNames
str(geneList)




########## Establish association between genes and GO terms
GO_list <- as.list(as.character(annot[['gene_id']])) # generates list; element names are gene IDs
GO_list <- as.list(setNames(as.character(annot$gene_ontology), as.character(annot$gene_id))) # adds Gene Ontology data to list
GO_list <- lapply(GO_list, function(x) unlist(strsplit(x, split='\\`'))) # split single GO terms string into a character vector, one element per term
GO_list <- lapply(GO_list, function(x) substr(x, 1, 10)) # strips away GO descriptions, leaving just ID number for each element
geneID2GO <- GO_list
str(head(geneID2GO,50))
## make full list of transcript names, geneNames
#geneNames <- names(GO_list)




########## build the topGOdata object
## This  object  will  contain  all  gene identifiers  and  their  scores,  the  GO  annotations,  
##    the  GO  hierarchical  structure  and  all  other  information needed to perform the desired
##    enrichment analysis.


GOdata <- new("topGOdata", 
              description = "DESeq2+sva+Trinotate-invertebrates session",
              ontology = "BP", 
              allGenes = geneList, 
              geneSel = topDiffGenes,
              nodeSize = 1,
              annot = annFUN.gene2GO, 
              gene2GO = geneID2GO)

## notes:
## description : character string containing a short description of the study [optional].
## ontology : character string specifying the ontology of interest 
##            >> "BP" (Biological Process), "MF" (Molecular Function) or "CC" (Cellular Component).
## allGenes : named vector of type numeric or factor.  The names attribute contains the genes identifiers.  
##            The genes listed in this object define the gene universe.
## geneSel : list of interesting genes (eg, genes with padj<0.05). 
## nodeSize : an integer larger or equal to 1.  This parameter is used to prune the GO hierarchy 
##            from the terms which have less than nodeSize annotated genes (after the true path rule 
##            is applied)
## annot : annotationFun is a function  which  maps  genes  identifiers  to  GO  terms. The annotation
##         functions take three arguments.  One of those arguments is specifying where the mappings 
##         can be found, and needs to be provided by the user. There are several possibilies. Here we use: 
##         >> annFUN.gene2GO: this function is used when the annotations are provided as a gene-to-GOs mapping.
## gene2GO : mapping provided by the user





## Number of genes annotated to the GO terms
numGenes(GOdata)

## basic stats of GOdata object:
termStat(GOdata)




####### Perform enrichment tests

## Fisher (classic) >> high false positive rate
resultFisher.classic <- runTest(GOdata, algorithm = "classic", statistic = "fisher")

## Fisher (elim) [elimination algorithm] >> low # false positives; risk of discarding some relevant nodes]
resultFisher.elim <- runTest(GOdata, algorithm = "elim", statistic = "fisher")

## Fisher (weight) [not as strict as elim, and thus a higher # of false positives]
resultFisher.weight <- runTest(GOdata, algorithm = "weight", statistic = "fisher")

## Fisher (weight01)  [mixture between the elim and the weight algorithms   >>> TopGO default ]
resultFisher.weight01 <- runTest(GOdata, algorithm = "weight01", statistic = "fisher")

## Fisher (parentchild)
resultFisher.parentChild <- runTest(GOdata, algorithm = "parentChild", statistic = "fisher")




####### Analysis of results
allRes <- GenTable( GOdata, 
                    classic_Fisher = resultFisher.classic, 
                    elim_Fisher = resultFisher.elim,
                    weight_Fisher = resultFisher.weight,
                    weight01_Fisher = resultFisher.weight01,
                    parentChild_Fisher = resultFisher.parentChild,
                    orderBy = "elim_Fisher",
                    ranksOf = "classic_Fisher",
                    topNodes = 1000)


allRes




####### save results to file
#write.table(allRes, file = "TopGO_OUTPUT.csv", sep="\t")




