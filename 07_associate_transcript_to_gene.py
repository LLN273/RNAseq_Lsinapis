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
# Written in Python 3
#
##################################################################################################################################################################### 




##################################################################################################################################################################### 
#
# Script used to associate Trinity transcript name to its gene name, based on Kallisto output file. 
# Input: kallisto's abundance.tsv file
# (Trinity transcript name format assumed, eg TRINITY_DN0_c0_g1_i1, where 'i1' is the isoform number).
#
##################################################################################################################################################################### 




print('\n Collecting TPM data ...')

######################################################### USAGE

error_message1 = ' \n \
 USAGE: $python3 transcript_to_gene.py <inputfile.txt> <outputfile.txt>'


# Example:
# python3 transcript_to_gene.py /kallisto/abundance.tsv gene_names.txt

######################################################### LOAD STANDARD MODULES

import sys
import re                                                 # 're' module provides regular expression matching operations
import ast                                                # required to convert reference tree from string to list
import os                                                 # required to create new folders
import time



######################################################## OPEN INPUT FILES

try:
    inputfile1 = open(sys.argv[1], 'r')                  # open Kallisto output file (transcript names located on first column)
except:
    print('\n Error: input files missing.')
    print(error_message1) 
    exit() 

try:
    outputfilename = sys.argv[2]      					     # output filename
except:
    print('\n Error: input value missing.')
    print(error_message1) 
    exit() 


    
######################################################## GLEAN DATA; SAVE DATA TO FILE 



### read input file, line by line 

headerFlag = 0

transcript_list=list()
gene_list=list()
for line in inputfile1:                             #get path to each subfolder
    line = line.rstrip()                            #remove break line mark from line
    if headerFlag == 0 :
        headerFlag = 1 
    else :
        auxID_marker = line.find('\t') 					#find first tab
        line_aux_1 = line[:(auxID_marker)]              #get transcript name
        transcript_list.append(line_aux_1)              #save transcript name to a list
        auxID_marker2 = line.rfind('_')                 #locate mark before isoform number
        line_aux_2 = line_aux_1[:(auxID_marker2)] 	    #get gene name
        gene_list.append(line_aux_2)                    #save gene name to a list



### save transcript name and gene name to file

outfile1 = open(outputfilename, 'w')


headers = 'TXNAME' + ',' + 'GENEID' + '\n'
outfile1.write(headers)
     
for k in range(len(transcript_list)) :   
    aux_summary = transcript_list[k] + ',' + gene_list[k] + '\n'   
    outfile1.write(aux_summary)
        
outfile1.close()									#close output file    

###



