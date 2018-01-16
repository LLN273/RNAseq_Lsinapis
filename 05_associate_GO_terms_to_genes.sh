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
# Associate GO terms to each Trinity gene
# Dependencies: GO_Associations.py
# 
# ********************************************************************************************* 





# Trinotate output file
# Note: trinotate_annotation_report-invertebrates.xls must be converted to csv format prior to running the following commands.
TRINOUT=trinotate_annotation_report-invertebrates.csv

# uniprot-KB db
UNIPROT_sprot=sprot_OUT_Insecta.txt

# uniprot-trembl db
UNIPROT_trembl=trembl_OUT_Insecta.txt

# Results file
OUTFILE=OUT_GO_Associations.csv





python3 GO_Associations.py $TRINOUT $UNIPROT_sprot $UNIPSPROT_trembl $OUTFILE






