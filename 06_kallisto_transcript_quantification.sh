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
# Transcript quantification using kallisto, ​v.​ ​0.43.0 ​(Bray​ ​​et​ ​al.​,​ ​2016)
# 
# ********************************************************************************************* 

kallisto quant \
            --index $TRANSDIR/index_kallisto/index.idx \
            -o ./kallisto \
            --bootstrap-samples 100 \
            --rf-stranded \
            --threads 1 \
            $READ1_fqs_filtered_adjusted \
            $READ2_fqs_filtered_adjusted
