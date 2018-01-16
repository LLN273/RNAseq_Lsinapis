####################################################################################################################
# 
# Pipeline used in:
#
# "Gene expression profiling across ontogenetic stages in wood white (Leptidea sinapis) reveals pathways linked to butterfly diapause regulation."
# Luis Leal, Venkat Talla, Thomas Källman, Magne Friberg, Christer Wiklund, Vlad Dincă, Roger Vila, Niclas Backström
# Mol. Eco. (2018)
#
#
########################################################################################################################################
#  
# Luis Leal
# Uppsala University, Uppsala, Sweden, 2017
#
##################################################################################################################################################################### 







##################################################################################################################################################################### 
#
# 1. QUALITY CONTROL		
#	
# -Filter libraries using TrimGalore		
# -Mask libraries using Fastq Masker		
# -Filter libraries using prinseq		
# -Filter libraries using condetri		
# -Screen libraries using fastQ Screen		
# -Reconcile paired-end libraries after using fastQ Screen (removes non-paired reads)		
# Note: Index fasta files as required.		
 
01_quality_control.sh





##################################################################################################################################################################### 
#
# 2. DE NOVO TRANSCRIPTOME ASSEMBLY USING TRINITY
#
# # Note: Index fasta files as required.
		
02_trinity_transcriptome_assembly.sh




	
		
##################################################################################################################################################################### 
#
# 3. FUNCTIONAL ANNOTATION BASED ON TRINOTATE PROTOCOL 		


# a. Identify each trinity transcript using modified Trinotate protocol

03_trinotate_functional_annotation.sh


# b. Build local UniProt databases

04_build_local_UniProt_db.sh


# c. Find GO terms associated to each Trinity gene

05_associate_GO_terms_to_genes.sh





##################################################################################################################################################################### 
#
# 4. TRANSCRIPT QUANTIFICATION USING KALLISTO
# Note: Index files as required.

06_kallisto_transcript_quantification.sh






##################################################################################################################################################################### 
#
# 5. GENE COUNT CONSOLIDATION USING TXIMPORT		
		

# a. Convert counts-per-transcript to counts-per-gene using tximport		

#		a1. Associate Trinity transcript name to its gene name

07_associate_transcript_to_gene.py


# 		a2. Convert Kallisto output files from count-per-transcript to count-per-gene

08_tximport_aggregate_transcript_counts.R


# b. Consolidate annotation file & gene counts		

09_tximport_consolidateGenes.py





##################################################################################################################################################################### 
#
# 6. DIFFERENTIAL GENE EXPRESSION ANALYSIS USING DESEQ2		
		
10_DESeq2.R


	
		
		
##################################################################################################################################################################### 
#
# 7. GENE ONTOLOGY ENRICHMENT ANALYSIS USING TOPGO		

	
# a. Create two-column file that lists each Trinity gene and its associated GO terms

11_gene_to_GO.py


# b. Perform gene ontology enrichment analysis 

12_TopGO.R



##################################################################################################################################################################### 
##################################################################################################################################################################### 




