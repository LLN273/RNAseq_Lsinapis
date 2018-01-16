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
# Build local UniProt databases
# Dependencies: glean_uniprot_spot_DB.py
# 
# ********************************************************************************************* 


# if you haven't done it yet, download the invertebrates UniProt libraries (sprot and trembl)
# wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/uniprot_sprot_invertebrates.dat.gz
# wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/uniprot_trembl_invertebrates.dat.gz

# Unzip both files and run the following commands:


python3 glean_uniprot_spot_DB.py uniprot_sprot_invertebrates.dat Insecta
mv sprot_OUT.txt sprot_OUT_Insecta.txt

python3 glean_uniprot_spot_DB.py uniprot_trembl_invertebrates.dat Insecta
mv sprot_OUT.txt trembl_OUT_Insecta.txt


