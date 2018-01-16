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
# Script used to parse paired-end files after 'fastQ Screen' filtering.
# Goes through both files and removes non-paired reads.
#
##################################################################################################################################################################### 




print('\n Parsing started ...')



######################################################### USAGE

error_message1 = ' \n \
 USAGE: $python3 adjust_paired_files.py <file1.fq> <file2.fq> '


######################################################### LOAD STANDARD MODULES

import sys
import re                                                 # 're' module provides regular expression matching operations
import ast                                                # required to convert reference tree from string to list
import os                                                 # required to create new folders
import time



######################################################## OPEN INPUT FILES

try:
    inputfile1_R = open(sys.argv[1], 'r')                  # open first fq file (reverse reads)
except:
    print('\n Error: input files missing.')
    print(error_message1) 
    exit() 

try:
    inputfile2_F = open(sys.argv[2], 'r')                  # open second fq file (forward reads)
except:
    print('\n Error: input files missing.')
    print(error_message1) 
    exit() 


######################################################## PUT READS FROM FILE2 INTO DICTIONARY FORMAT

#Header-sequence dictionary
FILE2_F_dict_SEQ = dict()

#Header-quality dictionary
FILE2_F_dict_QUAL = dict()

counterAux1 = 1

for line in inputfile2_F:
    if counterAux1 == 1:                      #isolate header
        marker1 = line.find(' ')              
        headerREF = line[:(marker1)]
    if counterAux1 == 2:
        seqREF = line
        FILE2_F_dict_SEQ[headerREF] = seqREF
        #print('seqREF:',seqREF)
    if counterAux1 == 4:
        qualREF = line
        #print('qualREF:',qualREF)
        FILE2_F_dict_QUAL[headerREF] = qualREF
    counterAux1 = counterAux1 + 1
    if counterAux1 == 5: counterAux1 = 1

        


# number of items in dictionaries
print('number of items in SEQ F-dictionary:', len(FILE2_F_dict_SEQ))
print('number of items in QUAL F-dictionary:', len(FILE2_F_dict_QUAL))









######################################################## CHECK WHETHER READS IN FILE1_R ARE PRESENT IN FILE2_F;
#                                                        SAVE PAIRED READS TO TWO NEW FILES


#output files

straux1=str(inputfile1_R)
nameseek = straux1.find("name='")
straux1 = straux1[(nameseek + 6):]
nameseek = straux1.find("' mode=")
straux1 = straux1[:nameseek]
#print('straux1',straux1)

straux2=str(inputfile2_F)
nameseek = straux2.find("name='")
straux2 = straux2[(nameseek + 6):]
nameseek = straux2.find("' mode=")
straux2 = straux2[:nameseek]
#print('straux2',straux2)

outFileName1_R = 'adjusted-' + straux1
outFileName2_F = 'adjusted-' + straux2

outfile1_R = open(outFileName1_R, 'w')                           
outfile2_F = open(outFileName2_F, 'w')


#match files

counterAux2 = 1
counterReadsRev = 0
 
for line in inputfile1_R:                     #for each read in FILE1_R
    counterReadsRev += 1
    if counterAux2 == 1:                      #isolate header
        marker1 = line.find(' ')           
        #print('marker1:', marker1)
        headerREF = line[:(marker1)]
    if counterAux2 == 2:
        seqREF = line
        #print('seqREF:',seqREF)
    if counterAux2 == 4:
        qualREF = line
        #print('qualREF:',qualREF)
    counterAux2 = counterAux2 + 1
    if counterAux2 == 5: 
        counterAux2 = 1
        seekDictSEQ = FILE2_F_dict_SEQ.get(headerREF, 'MinkoSays:Hobiron!')     #search header in dictionary (File2_F)
        #print(seekDictSEQ)
        if seekDictSEQ != 'MinkoSays:Hobiron!':                                #if read is present in dictionary
            #print('score!')
            headN=headerREF
            headmarker = headN.find('.')
            headN = headN[(headmarker + 1):]
            #save record to file1_R
            aux1_R = headerREF + ' ' + '1' + '\n'
            aux2_R = seqREF
            aux3_R = '+' + '\n'
            aux4_R = qualREF
            outfile1_R.write(aux1_R)
            outfile1_R.write(aux2_R)
            outfile1_R.write(aux3_R)
            outfile1_R.write(aux4_R)
            #save record to file2_F
            seekDictQUAL = FILE2_F_dict_QUAL.get(headerREF, 'MinkoSaysHobiron')
            aux1_F = headerREF + ' ' + '2' + '\n'
            aux2_F = seekDictSEQ
            aux3_F = '+' + '\n'
            aux4_F = seekDictQUAL
            outfile2_F.write(aux1_F)
            outfile2_F.write(aux2_F)
            outfile2_F.write(aux3_F)
            outfile2_F.write(aux4_F)
        


print('Number of reads in input FILE1_R:', int(counterReadsRev/4))
print()



outfile1_R.close()
outfile2_F.close()



