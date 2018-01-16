##################################################################################################################################################################### 
# 
# Pipeline used in:
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




# ********************************************************************************************* 
#
# Functional annotation using modified Trinotate, v.​ ​3.0.2​​ ​(Haas,​ ​2017a) protocol
# For pipeline details, see https://trinotate.github.io/
# Dependencies: Build_Trinotate_Boilerplate_SQLite_db_MODIFIED_INVERTEBRATES.pl
# 
# ********************************************************************************************* 








###################### PART A: Obtain Sequence Databases Required & Run Sequence Analyses



##### PART 1: Prepare 'Sequence Databases Required' and obtain 'Files needed for execution'


### Download invertebrates UniProt libraries (sprot and trembl)
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/uniprot_sprot_invertebrates.dat.gz
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/uniprot_trembl_invertebrates.dat.gz



### Run modified Trinotate script that downloads other relevant libraries and initializes sqlite database
./Build_Trinotate_Boilerplate_SQLite_db_MODIFIED_INVERTEBRATES.pl Trinotate-invertebrates       



### Use TransDecoder,​ ​v.​ ​3.0.1 to get the longest-ORF peptide candidates generated from the Trinity Assembly.
TransDecoder.LongOrfs -t Trinity.fasta
TransDecoder.Predict -t Trinity.fasta



### Prepare the protein database for blast (v. 2.6.0) searches:
makeblastdb -in uniprot_sprot.pep -dbtype prot



### Uncompress and prepare the Pfam database for use with hmmscan (hmmer, v. 3.1b2):
gunzip -c Pfam-A.hmm.gz > Pfam-A.hmm
hmmpress Pfam-A.hmm







##### PART 2. Capturing BLAST Homologies


### search Trinity transcripts (blast, v. 2.6.0)
blastx -query Trinity.fasta -db uniprot_sprot.pep -num_threads 16 -max_target_seqs 1 -outfmt 6 > blastx.outfmt6



### search Transdecoder-predicted proteins (blast, v. 2.6.0)
blastp -query Trinity.fasta.transdecoder.pep -db uniprot_sprot.pep -num_threads 16 -max_target_seqs 1 -outfmt 6 > blastp.outfmt6



### Blast Trinity transcripts against custom databases

## Download proteomes​​ from Uniprot database

#Bombyx_mori
ftp://​ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000005204_7091.fasta.gz​

#Danaus_plexippus​
ftp://​ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000007151_13037.fasta.gz​

#Papilio_machaon​
ftp://​ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000053240_76193.fasta.gz​

#Drosophila_melanogaster
ftp://​ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000000803_7227.fasta.gz

#Drosophila_melanogaster (ensembl)
wget ftp://ftp.ensembl.org/pub/release-89/fasta/drosophila_melanogaster/pep/Drosophila_melanogaster.BDGP6.pep.all.fa.gz

#Heliconius_melpomene (lepbase)
#Note: this file includes some troublesome characters (ä and ü). Copy file, edit and replace ü for u, and ä for a. Save as  Heliconius_melpomene_melpomene_Hmel2_-_proteins_MODIFIED.fa
wget http://download.lepbase.org/v4/sequence/Heliconius_melpomene_melpomene_Hmel2_-_proteins.fa.gz


## Blast
## Note: must create protein blast database (makeblastdb) for each reference genome before running the following commands

##drosophila_melanogaster (Ensembl) >> not used during downstream analysis
blastx -query Trinity.fasta -db Drosophila_melanogaster.BDGP6.pep.all.fa -num_threads 16 -max_target_seqs 1 -outfmt 6 > blastx-Drosophila_melanogaster_Ensembl.outfmt6
blastp -query Trinity.fasta.transdecoder.pep -db Drosophila_melanogaster.BDGP6.pep.all.fa -num_threads 16 -max_target_seqs 1 -outfmt 6 > blastp-Drosophila_melanogaster_Ensembl.outfmt6

#drosophila_melanogaste
blastx -query Trinity.fasta -db Drosophila_melanogaster_UP000000803_7227.fasta -num_threads 16 -max_target_seqs 1 -outfmt 6 > blastx-Drosophila_melanogaster_Uniprot.outfmt6
blastp -query Trinity.fasta.transdecoder.pep -db Drosophila_melanogaster_UP000000803_7227.fasta -num_threads 16 -max_target_seqs 1 -outfmt 6 > blastp-Drosophila_melanogaster_Uniprot.outfmt6

##Heliconius_melpomene (lepbase) >> not used during downstream analysis
blastx -query Trinity.fasta -db Heliconius_melpomene_melpomene_Hmel2_-_proteins_MODIFIED.fa -num_threads 16 -max_target_seqs 1 -outfmt 6 > blastx-Heliconius_melpomene_Ensembl.outfmt6
blastp -query Trinity.fasta.transdecoder.pep -db Heliconius_melpomene_melpomene_Hmel2_-_proteins_MODIFIED.fa -num_threads 16 -max_target_seqs 1 -outfmt 6 > blastp-Heliconius_melpomene_Ensembl.outfmt6

#Bombyx_mori
blastx -query Trinity.fasta -db Bombyx_mori_UP000005204_7091.fasta -num_threads 16 -max_target_seqs 1 -outfmt 6 > blastx-Bombyx_mori_Uniprot.outfmt6
blastp -query Trinity.fasta.transdecoder.pep -db Bombyx_mori_UP000005204_7091.fasta -num_threads 16 -max_target_seqs 1 -outfmt 6 > blastp-Bombyx_mori_Uniprot.outfmt6

#Danaus_plexippus
blastx -query Trinity.fasta -db Danaus_plexippus_UP000007151_13037.fasta -num_threads 16 -max_target_seqs 1 -outfmt 6 > blastx-Danaus_plexippus_Uniprot.outfmt6
blastp -query Trinity.fasta.transdecoder.pep -db Danaus_plexippus_UP000007151_13037.fasta -num_threads 16 -max_target_seqs 1 -outfmt 6 > blastp-Danaus_plexippus_Uniprot.outfmt6

#Papilio_machaon
blastx -query Trinity.fasta -db Papilio_machaon_UP000053240_76193.fasta -num_threads 16 -max_target_seqs 1 -outfmt 6 > blastx-Papilio_machaon_Uniprot.outfmt6
blastp -query Trinity.fasta.transdecoder.pep -db Papilio_machaon_UP000053240_76193.fasta -num_threads 16 -max_target_seqs 1 -outfmt 6 > blastp-Papilio_machaon_Uniprot.outfmt6





##### PART 3. Run HMMER, v. 3.1b2, to identify protein domains
hmmscan --cpu 1 --domtblout TrinotatePFAM.out Pfam-A.hmm Trinity.fasta.transdecoder.pep > pfam.log





##### PART 4. Run signalP, v. 4.1c, to predict presence and location of signal peptides cleavage sites
signalp -f short -n signalp.out Trinity.fasta.transdecoder.pep





##### PART 5. Run tmHMM, v. 2.0c, to predict transmembrane regions
tmhmm --short Trinity.fasta.transdecoder.pep > tmhmm.out




##### PART 6. Run RNAMMER, v. 1.2, to identify rRNA transcripts
/Trinotate-3.0.2/util/rnammer_support/RnammerTranscriptome.pl --transcriptome Trinity.fasta --path_to_rnammer /sw/apps/bioinfo/rnammer/1.2/rnammer








 
###################### PART B: Loading Above Results into a Trinotate SQLite Database (sqlite, v.3.16.2)



##### PART 1. Load transcripts and coding regions

/trinity/2.3.2/util/support_scripts/get_Trinity_gene_to_trans_map.pl Trinity.fasta > Trinity.fasta.gene_trans_map


/Trinotate-3.0.2/Trinotate Trinotate-invertebrates.sqlite init --gene_trans_map Trinity.fasta.gene_trans_map --transcript_fasta Trinity.fasta --transdecoder_pep Trinity.fasta.transdecoder.pep





##### PART 2. Loading BLAST homologies

### A: Load main swissprot library

## load protein hits
/Trinotate-3.0.2/Trinotate Trinotate-invertebrates.sqlite LOAD_swissprot_blastp blastp.outfmt6

## load transcript hits
/Trinotate-3.0.2/Trinotate Trinotate-invertebrates.sqlite LOAD_swissprot_blastx blastx.outfmt6



### Load hits based on custom libraries 

## Drosophila_melanogaster_Ensembl >> not used during downstream analysis
/Trinotate-3.0.2/Trinotate Trinotate-invertebrates.sqlite LOAD_custom_blast --outfmt6 blastp-Drosophila_melanogaster_Ensembl.outfmt6 --prog blastp --dbtype Drosophila_Ensembl_blastp
/Trinotate-3.0.2/Trinotate Trinotate-invertebrates.sqlite LOAD_custom_blast --outfmt6 blastx-Drosophila_melanogaster_Ensembl.outfmt6 --prog blastx --dbtype Drosophila_Ensembl_blastx


## Drosophila_melanogaster
/Trinotate-3.0.2/Trinotate Trinotate-invertebrates.sqlite LOAD_custom_blast --outfmt6 blastp-Drosophila_melanogaster_Uniprot.outfmt6 --prog blastp --dbtype Drosophila_Uniprot_blastp
/Trinotate-3.0.2/Trinotate Trinotate-invertebrates.sqlite LOAD_custom_blast --outfmt6 blastx-Drosophila_melanogaster_Uniprot.outfmt6 --prog blastx --dbtype Drosophila_Uniprot_blastx


## Heliconius_melpomene_Ensembl >> not used during downstream analysis
/Trinotate-3.0.2/Trinotate Trinotate-invertebrates.sqlite LOAD_custom_blast --outfmt6 blastp-Heliconius_melpomene_Ensembl.outfmt6 --prog blastp --dbtype Heliconius_blastp
/Trinotate-3.0.2/Trinotate Trinotate-invertebrates.sqlite LOAD_custom_blast --outfmt6 blastx-Heliconius_melpomene_Ensembl.outfmt6 --prog blastx --dbtype Heliconius_blastx


## Bombyx_mori
/Trinotate-3.0.2/Trinotate Trinotate-invertebrates.sqlite LOAD_custom_blast --outfmt6 blastp-Bombyx_mori_Uniprot.outfmt6 --prog blastp --dbtype Bombyx_blastp
/Trinotate-3.0.2/Trinotate Trinotate-invertebrates.sqlite LOAD_custom_blast --outfmt6 blastx-Bombyx_mori_Uniprot.outfmt6 --prog blastx --dbtype Bombyx_blastx


## Danaus_plexippus
/Trinotate-3.0.2/Trinotate Trinotate-invertebrates.sqlite LOAD_custom_blast --outfmt6 blastp-Danaus_plexippus_Uniprot.outfmt6 --prog blastp --dbtype Danaus_blastp
/Trinotate-3.0.2/Trinotate Trinotate-invertebrates.sqlite LOAD_custom_blast --outfmt6 blastx-Danaus_plexippus_Uniprot.outfmt6 --prog blastx --dbtype Danaus_blastx


## Papilio_machaon
/Trinotate-3.0.2/Trinotate Trinotate-invertebrates.sqlite LOAD_custom_blast --outfmt6 blastp-Papilio_machaon_Uniprot.outfmt6 --prog blastp --dbtype Papilio_blastp
/Trinotate-3.0.2/Trinotate Trinotate-invertebrates.sqlite LOAD_custom_blast --outfmt6 blastx-Papilio_machaon_Uniprot.outfmt6 --prog blastx --dbtype Papilio_blastx





##### PART 3. Load Pfam domain entries
/Trinotate-3.0.2/Trinotate Trinotate-invertebrates.sqlite LOAD_pfam TrinotatePFAM.out





##### PART 4. Load transmembrane domains
/Trinotate-3.0.2/Trinotate Trinotate-invertebrates.sqlite LOAD_tmhmm tmhmm.out





##### PART 5. Load signal peptide predictions
/Trinotate-3.0.2/Trinotate Trinotate-invertebrates.sqlite LOAD_signalp signalp.out








###################### PART C:  Combine search results and create annotation report 

/Trinotate-3.0.2/Trinotate Trinotate-invertebrates.sqlite report > trinotate_annotation_report-invertebrates.xls





