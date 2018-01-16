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
# Script used to create two-column file that lists each Trinity gene and its associated GO terms
#
# Input files:	  -> trinotate-invertebrates-CONSOLIDATED.csv (consolidated trinotate file produced by 09_tximport_consolidateGenes.py)
#                 -> text file containing list of Trinity genes
#
# Output files:	  -> file listing genes and their associated GO terms ('OUTPUT_trinity_to_GO.txt')
#                 -> file listing GO codes and their associated genes ('OUTPUT_GO_to_trinity.txt')
#                 -> file associating transcripts to gene names ('OUTPUT_trinity_to_gene.txt')
#
##################################################################################################################################################################### 




print('\n Gleaning GO terms ...')

######################################################### USAGE

error_message1 = ' \n \
 USAGE: $python3 gene_to_GO.py <annotation_file.csv> <gene_list.txt>'


# Example:
# python3 gene_to_GO.py trinotate-invertebrates-CONSOLIDATED.csv genelist_test1.txt

######################################################### LOAD STANDARD MODULES

import sys
import re                                                 # 're' module provides regular expression matching operations
import ast                                                # required to convert reference tree from string to list
import os                                                 # required to create new folders
import time
import csv                                                #required to open csv file



######################################################## OPEN INPUT FILES	

GENE_NAME_DIC = dict()
GENE_ANNOT_DIC = dict()
try:   
    with open(sys.argv[1], newline='') as csvfile1:               # open consolidated Trinotate file
        inputfile_1 = csv.reader(csvfile1, delimiter='\t')
        for row in inputfile_1:
            trinity_aux = row[0]
            gene_name_aux = row[4]
            annot_aux = row[17]
            GENE_NAME_DIC[trinity_aux] = gene_name_aux            #saves entry to dictionary (removing quotes)
            GENE_ANNOT_DIC[trinity_aux] = annot_aux               #saves entry to dictionary (removing quotes)
except:
    print('\n Error[1]: input file missing.')
    print(error_message1) 
    exit() 





try:
    inputfile_2 = open(sys.argv[2], 'r')      			          # open file containing list of genes of interest
except:
    print('\n Error[2]: input file missing.')
    print(error_message1) 
    exit() 


######################################################## PROCESS DATA




### store entries from file 2 as a list 
goi_list = list()
for line in inputfile_2:                          
    line = line.rstrip()                            
    goi_list.append(line)                           







### Get GO terms associated to genes in gene list
GO_list = list()
for gene in goi_list :
        aux_1 = GENE_ANNOT_DIC.get(gene ,0)
        GO_list.append(aux_1)
        #print(gene, aux_1)


#print(GO_list)


### Get gene name associated to genes in gene list
geneName_list = list()
for gene in goi_list :
        aux_1 = GENE_NAME_DIC.get(gene ,0)
        geneName_list.append(aux_1)
        #print(gene, aux_1)






### get genes associated to each GO terms (for genes of interest only)

GO_DIC = dict()           #dictionary: Go term and associated genes
GO_ind = list()           #list of all GO terms
counter_1 = -1

for line in GO_list:
    #print(line)
    aux_list1 = list()
    counter_1 += 1
    GO_aux = ""
    if line != '.' :                                          					# if GO term != '.' (empty)
#        print(line + '\n')
        for m in re.finditer('GO:', line):                    					# find all GO terms in description
            aux_list1.append(m.start())
#            print(aux_list1)     
        for n in range(len(aux_list1)) :                                		# for all GO terms in description, get GO term (based on index)
#            print(line[aux_list1[n]])
            if n != (len(aux_list1)-1) :                                
                GO_aux = line[aux_list1[n]:(aux_list1[(n+1)]-1)]
            else :
                GO_aux = line[aux_list1[n]:]
            if GO_aux != "" :                                           		# if GO term is not empty	
                GO_ind.append(GO_aux)
                if GO_aux in GO_DIC :                                           # if GO term already in dictionary
                    entry_aux = GO_DIC.get(GO_aux ,0) 
                    if goi_list[counter_1] not in entry_aux :                   # make sure gene has not been associated yet to GO term in dictionary
                        entry_GO_aux = entry_aux + "," +  goi_list[counter_1]
                        GO_DIC[GO_aux] = entry_GO_aux
                else :
                    GO_DIC[GO_aux] = goi_list[counter_1]
            


GO_ind = sorted(list(set(GO_ind)))     #remove duplicates and sort list of GO terms alphabetically

#print dictionary
#for term in GO_ind:
#    print(term, GO_DIC[term])




######################################################## SAVE RESULTS TO FILE

outFileName1 = 'OUTPUT_trinity_to_GO.txt'
outFileName2 = 'OUTPUT_GO_to_trinity.txt'
outFileName3 = 'OUTPUT_trinity_to_gene.txt'

outfile1 = open(outFileName1, 'w') 
outfile2 = open(outFileName2, 'w') 
outfile3 = open(outFileName3, 'w') 

for n in range(len(goi_list)) :
    aux_save = goi_list[n] + '\t' + GO_list[n] + '\n'
    outfile1.write(aux_save)


for term in GO_ind :
    aux_save = term[0:10] + '\t' + GO_DIC[term] + '\n'
    outfile2.write(aux_save)

for n in range(len(goi_list)) :
    aux_save = goi_list[n] + '\t' + geneName_list[n] + '\n'
    outfile3.write(aux_save)

outfile1.close()
outfile2.close()
outfile3.close()



