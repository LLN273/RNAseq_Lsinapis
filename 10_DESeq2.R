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
# R/3.4.0 script used to carry out diferential expression analysis using the DESeq2 R Bioconductor package 
#
# Dependencies: tximport_OUTPUT_CONSOLIDATED.csv >> consolidated tximport-processed kallisto file with counts per gene (all samples)
#               samples_ID.txt >> auxiliary table describing the experimental design 
#               
#
###################################################################################################



##Install DESeq2
#source("http://bioconductor.org/biocLite.R")
#biocLite("DESeq2")

##Install sva
#biocLite("sva")

## Install lpsymphony package (required for IHW package)
#biocLite("lpsymphony")

## Install IHW package (Independent Hypothesis Weighting)
#biocLite("IHW")

## Install pcaExplorer (used for PCA analysis)
#biocLite("pcaExplorer")

## Install ReportingTools package
#biocLite("ReportingTools")


#####

#load libraries
library(DESeq2)
library(sva)


#####

#Clear all states and remove all plots
rm(list=ls(all=TRUE))
if(!is.null(dev.list())) dev.off()



########## input files and paths

### Path to folder containing ximport-processed kallisto files in count format (non-normalized);
base_dir <- "/tximport"
list.files(base_dir)            

### name of tximport-processed kallisto file with counts per gene (all samples)
txi_file <- "tximport_OUTPUT_CONSOLIDATED.csv" 

### import tximport-processed kallisto file with counts per gene
txi <- read.table(file.path(base_dir, txi_file), header = TRUE)
head(txi)


### Specify location of auxiliary table describing the experimental design
table_dir <- "/tximport"

### Specify name of auxiliary table describing the experimental design
auxTable_filename <- "samples_ID.txt"

### Read auxiliary table
samples_table <- read.table(file.path(table_dir, auxTable_filename), header = TRUE)
samples_table





######## Experiment number

experimentNumber <- 1         # 01_SDIII_LDIII  >> 1
                              # 02_SDV_LDV      >> 2
                              # 03_SDP_LDP      >> 3
                              # 04_SDIII_SDV    >> 4
                              # 05_SDV_SDP      >> 5
                              # 07_LDIII_LDV    >> 7
                              # 08_LDV_LDP      >> 8
                              # 09_LDP_LDA      >> 9



####################################################### Process input data

######## Round estimated counts to closest integer 
txi_r <- round(txi)
head(txi_r)



######## Construct a DESeqDataSet object

### if assessing differences across DAY-LENGTH
dds <- DESeqDataSetFromMatrix(countData = txi_r,
                              colData = samples_table,
                              design = ~ daylength)               
dds                                   

### if assessing differences across STAGE
#dds <- DESeqDataSetFromMatrix(countData = txi_r,
#                              colData = samples_table,
#                              design = ~ stage)               
#dds 







######## Filtering: remove genes with baseMean (average across all samples) lower than 1

dds <- dds[ rowMeans(counts(dds)) > 1, ]
dds







######## Select samples for DGE downstream analysis 


if ( experimentNumber == 1) {
    ##
    ## 01_SDIII_LDIII
    dds <- dds[,c("P5052_101_S39","P5052_109_S40","P5052_117_S42","P5052_201_S47","P5052_209_S50","P5052_217_S54")]
    dds
    #
} else if ( experimentNumber == 2) {
    #
    ## 02_SDV_LDV
    dds <- dds[,c("P5052_102_S27","P5052_103_S28","P5052_110_S41","P5052_111_S31","P5052_118_S43","P5052_119_S44","P5052_202_S48","P5052_203_S49","P5052_210_S51","P5052_211_S52","P5052_218_S55","P5052_219_S56")]
    dds
    #
} else if ( experimentNumber == 3) {
    #
    ## 03_SDP_LDP
    dds <- dds[,c("P5052_104_S29","P5052_105_S30","P5052_112_S32","P5052_113_S33","P5052_120_S45","P5052_121_S46","P5052_204_S34","P5052_205_S35","P5052_212_S36","P5052_213_S53","P5052_220_S57","P5052_221_S58")]
    dds
    #
} else if ( experimentNumber == 4) {    
    #
    ## 04_SDIII_SDV
    dds <- dds[,c("P5052_101_S39","P5052_102_S27","P5052_103_S28","P5052_109_S40","P5052_110_S41","P5052_111_S31","P5052_117_S42","P5052_118_S43","P5052_119_S44")]
    dds
    #
} else if ( experimentNumber == 5) {     
    #
    ## 05_SDV_SDP
    dds <- dds[,c("P5052_102_S27","P5052_103_S28","P5052_104_S29","P5052_105_S30","P5052_110_S41","P5052_111_S31","P5052_112_S32","P5052_113_S33","P5052_118_S43","P5052_119_S44","P5052_120_S45","P5052_121_S46")]
    dds
    #
} else if ( experimentNumber == 7) { 
    #
    ## 07_LDIII_LDV
    dds <- dds[,c("P5052_201_S47","P5052_202_S48","P5052_203_S49","P5052_209_S50","P5052_210_S51","P5052_211_S52","P5052_217_S54","P5052_218_S55","P5052_219_S56")]
    dds
    #
} else if ( experimentNumber == 8) { 
    #
    ## 08_LDV_LDP
    dds <- dds[,c("P5052_202_S48","P5052_203_S49","P5052_204_S34","P5052_205_S35","P5052_210_S51","P5052_211_S52","P5052_212_S36","P5052_213_S53","P5052_218_S55","P5052_219_S56","P5052_220_S57","P5052_221_S58")]
    dds
    #
} else if ( experimentNumber == 9) { 
    #
    ## 09_LDP_LDA
    dds <- dds[,c("P5052_204_S34","P5052_205_S35","P5052_212_S36","P5052_213_S53","P5052_220_S57","P5052_221_S58","P5052_226_S59","P5052_227_S38","P5052_233_S60","P5052_234_S61","P5052_235_S62")]
    dds
    #
}






######## drop levels

dds$family <- droplevels(dds$family)
dds$family
dds$stage <- droplevels(dds$stage)
dds$stage
dds$sex <- droplevels(dds$sex)
dds$sex
dds$replicate <- droplevels(dds$replicate)
dds$replicate
dds$daylength <- droplevels(dds$daylength)
dds$daylength
dds






######## change variable names after filtering
dds_f <- dds






######## MODEL BATCH EFFECTS: Remove hidden batch effects using the sva package

dds_sva <- dds_f
dds_sva <- estimateSizeFactors(dds_sva)                   # normalization
dat  <- counts(dds_sva, normalized = TRUE)
idx  <- rowMeans(dat) > 1                                 #filter after normalization
dat  <- dat[idx, ]

# declare model with variable of interest
mod  <- model.matrix(~ daylength, colData(dds_sva))  # when looking for DGE across ligth treatments (day-length)
#mod  <- model.matrix(~ stage, colData(dds_sva))     # when looking for DGE across stages

#declare known adjustment variables, or let sva discover batch effects (creates new surrogate variables)
mod0 <- model.matrix(~ 1, colData(dds_sva))


# estimate number of surrogate variables
svseq <- svaseq(dat, mod, mod0)    ##   If  the sva function  is  called  without  the n.sv argument  specified, the  number  
                                   ##   of factors will be estimated for you
svseq$sv



## Use sva to remove any effect of surrogate variables
## (adjust according to number of surrogate variables)
ddssva_f <- dds_f
ddssva_f$SV1 <- svseq$sv[,1]
ddssva_f$SV2 <- svseq$sv[,2]
#ddssva_f$SV3 <- svseq$sv[,3]
#ddssva_f$SV4 <- svseq$sv[,4]
design(ddssva_f) <- ~ SV1 + SV2 + daylength                  # 01_SDIII_LDIII
#design(ddssva_f) <- ~ SV1 + SV2 + SV3 + SV4 + daylength      # 02_SDV_LDV 
#design(ddssva_f) <- ~ SV1 + SV2 + SV3 + daylength            # 03_SDP_LDP 
#design(ddssva_f) <- ~ SV1 + SV2 + SV3 + stage                # 04_SDIII_SDV, 05_SDV_SDP, 06_LDV_LDP, 07_LDIII_LDV, 09_LDP_LDA







######## Set reference condition
dds_f$day_length <- relevel(dds_f$daylength, ref="SD")
#dds_f$stage <- relevel(dds_f$stage, ref="III")
#dds_f$stage <- relevel(dds_f$stage, ref="V")
#dds_f$stage <- relevel(dds_f$stage, ref="P")

ddssva_f$day_length <- relevel(ddssva_f$daylength, ref="SD")
#ddssva_f$stage <- relevel(ddssva_f$stage, ref="III")
#ddssva_f$stage <- relevel(ddssva_f$stage, ref="V")
#ddssva_f$stage <- relevel(ddssva_f$stage, ref="P")








###################################### Perform differential gene expression analysis

dds_f <- DESeq(dds_f)
res <- results(dds_f)
res

ddssva_f <- DESeq(ddssva_f)
resSVA <- results(ddssva_f)
resSVA




##### check conditions being compared
resultsNames(dds_f)
resultsNames(ddssva_f)

mcols(res, use.names = TRUE)
mcols(resSVA, use.names = TRUE)



##### summarize results for alpha=0.1 (associated to DESeq's 'independent filtering')
##use when summarizing DGE results FDR threshold of 10%
#sum(res$padj < 0.1, na.rm=TRUE)
#sum(res$padj < 0.05, na.rm=TRUE)

#sum(resSVA$padj < 0.1, na.rm=TRUE)
#sum(resSVA$padj < 0.05, na.rm=TRUE)




##### summarize results for alpha=0.05
## use when summarizing DGE results FDR threshold of 5%

res05 <- results(dds_f, alpha=.05)
#table(res05$padj < .05) 
#summary(res05)
sum(res05$padj < 0.1, na.rm=TRUE)
sum(res05$padj < 0.05, na.rm=TRUE)

resSVA05 <- results(ddssva_f, alpha=.05)
#table(resSVA05$padj < .05) 
#summary(resSVA05)
sum(resSVA05$padj < 0.1, na.rm=TRUE)
sum(resSVA05$padj < 0.05, na.rm=TRUE)





##### Order results table by the smallest adjusted p value (alpha=0.1)
#resOrdered <- res[order(res$padj),]
#resOrdered

#resOrderedSVA <- resSVA[order(resSVA$padj),]
#resOrderedSVA





####### Order results table by the smallest adjusted p value (alpha=0.05)
resOrdered <- res05[order(res05$padj),]
resOrdered

resOrderedSVA <- resSVA05[order(resSVA05$padj),]
resOrderedSVA





##### Write results of DGE analysis to file in CSV format

#write.table(resOrderedSVA, file = "DESeq2_OUTPUT_sva.csv", sep="\t")
#write.table(resOrdered, file = "DESeq2_OUTPUT_BASIC-ANALYSIS.csv", sep="\t")






