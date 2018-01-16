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
# Transcriptome assembly using Trinity, v. 2.3.2 (Grabherr et al. 2011)
# Dependencies: sample_list.txt (adjust as required)
# 
# ********************************************************************************************* 

Trinity --seqType fq \
        --samples_file $SAMPLE_LIST \
        --SS_lib_type RF \
        --min_contig_length 126 \
        --CPU 8 \
        --max_memory 50G \
        --output ./

 
