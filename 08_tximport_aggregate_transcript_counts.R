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
# R/3.4.0 script used to convert Kallisto output files from count-per-transcript to count-per-gene
#
# Dependencies: abundance.tsv >> kallisto output file
#               samples_ID_tximport.txt >> basic identification of each sample (experimental design)
#               gene_names.csv >> associates each transcript to its parental gene
#
###################################################################################################



###Install tximport
#source("https://bioconductor.org/biocLite.R")
#biocLite("tximport")

### Install HDF5 interface to R package
#biocLite("rhdf5")

#### Install readr package
#install.packages("readr")



###load tximport and rhdf5 libraries
library(tximport)
library(rhdf5)
library(readr)


###Clear all states
rm(list=ls(all=TRUE))



##### PATHS

### locate the directory containing the kallisto files; list files in directory
base_dir <- "/kallisto_OUT"
list.files(base_dir)

### Specify where auxiliary table describing the experimental design is stored
table_dir <- "/tximport"

### Specify name of auxiliary table describing the experimental design
auxTable_filename <- "samples_ID_tximport.txt"

### Specify output folder
out_dir <- "/tximportOUT"




##### Import auxiliary files

### Read auxiliary table; show table
samples <- read.table(file.path(table_dir, auxTable_filename), header = TRUE)
samples

### Read file that associates transcripts with gene IDs (necessary for gene-level summarization)
tx2gene <- read.csv(file.path(table_dir, "gene_names.csv"), header = TRUE)
head(tx2gene)




##### Import Kallisto file
files <- file.path(base_dir, "results", samples$experiment, "abundance.tsv")
names(files) <-t(samples[1])
files



##### Run tximport
### Note: do NOT normalize by library size (DESeq2 will perform this step)
txi <- tximport(files, type = "kallisto", tx2gene = tx2gene)
names(txi)
head(txi$counts)

# Write data to file in CSV format
write.table(txi$counts, file.path(out_dir, "tximport_OUTPUT.csv"), sep="\t")


