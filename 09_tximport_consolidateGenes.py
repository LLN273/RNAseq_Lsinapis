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
# Script used to consolidate Trinity entries and counts associated to the same gene
#
# Input files: 1 -> tximport output file (tximport_OUTPUT.csv)
#              2 -> Annotation file (OUT_GO_Associations.csv)  
#
# Output: 2 files -> updated tximport output file (consolidated counts)
#                 -> updated annotation file (consolidated GO terms)
#
##################################################################################################################################################################### 




print('\n Gleaning data ...')

######################################################### USAGE

error_message1 = ' \n \
 USAGE: $python3 tximport_consolidateGenes.py <tximport_OUTPUT.csv> <annotation_file.csv> \n'


######################################################### LOAD STANDARD MODULES

import sys
import re                                                 # 're' module provides regular expression matching operations
import ast                                                # required to convert reference tree from string to list
import os                                                 # required to create new folders
import time
import csv                                                #required to open csv file

from pdb import set_trace as bp                           # when using break points during debugging







######################################################## OPEN INPUT FILES	


try:
    tximport_list = list()                                    #open tximport output file as csv
    with open(sys.argv[1], newline='') as csvfile1:         
        inputfile_1 = csv.reader(csvfile1, delimiter='\t')
        for row in inputfile_1:
            tximport_list.append(row)
#            print(', '.join(row))
except:
    print('\n Error[1]: input file missing.')
    print(error_message1) 
    exit() 



try:   
    annotation_list = list()                                   #open annotation file as csv
    with open(sys.argv[2], newline='') as csvfile2:         
        inputfile_2 = csv.reader(csvfile2, delimiter='\t')
        for row in inputfile_2:
            annotation_list.append(row)
except:
    print('\n Error[2]: input file missing.')
    print(error_message1) 
    exit() 









####################################################### FUNCTION: readGENE

def readGene(geneName_long):
    if geneName_long != "." :
        mark_g = geneName_long.find('|')
        if mark_g != -1 :
            aux_g = geneName_long[(mark_g+1):]
            mark_g2 = aux_g.find('|')
            if mark_g2 != -1 :
                geneName_long = aux_g[:mark_g2]

    return geneName_long
    
####### END









####################################################### FUNCTION: updateList_v2

def updateList_v2(newS,referenceS) :
    if newS == '.' or newS == '' : 
        aux_output = referenceS
    elif referenceS == '.' or referenceS == '' :
        aux_output = newS
    elif newS in referenceS :
        aux_output = referenceS
    elif referenceS in newS :
        aux_output = newS
    else :
        aux_output = referenceS + "`" + newS

    return aux_output

####### END








####################################################### FUNCTION: cleanGO


def  cleanGO(entry_aux) :
    aux_list1 = list()
    aux_list_terms = list()
    if entry_aux != "." and entry_aux != "":
        for m in re.finditer('GO:', entry_aux):                    			    # find all GO terms in description
            aux_list1.append(m.start())
        for n in range(len(aux_list1)) :                                		# for all GO terms in description, get GO term (based on index)
            GO_aux3 = entry_aux[aux_list1[n]:(aux_list1[(n)]+10)]
            #print(GO_aux3)
            aux_list_terms.append(GO_aux3)
        aux_list_terms = sorted(list(set(aux_list_terms)))                      # remove duplicates and sort list of GO terms alphabetically
        Go_clean = "^NA`".join(aux_list_terms)                                  # convert list to string (Trinotate format)
    else :
        Go_clean = entry_aux
        
    return Go_clean

####### END








####################################################### FUNCTION: entryCompare


def entryCompare(entry_1, entry_2) :
    dif_counter_1 = 0
    dif_counter_2 = 0
    #
    if entry_1 == entry_2 :                                                                # the two lists of gene-hits are identical
        compare_out = 1
    else:
        entry_list_1 = entry_1.split('`')
        entry_list_2 = entry_2.split('`')
        #
        for gene in entry_list_1 :
            if gene not in entry_list_2 : dif_counter_1 += 1
        #
        for gene in entry_list_2 :
            if gene not in entry_list_1 : dif_counter_2 += 1
        #
        if (dif_counter_1 == 0) and (len(entry_list_1) >= 2) :                              #largest list contains shortest list and smallest list is at least 2 elements long
            compare_out = 2          
        elif (dif_counter_2 == 0) and (len(entry_list_2) >= 2) :                            #largest list contains shortest list and smallest list is at least 2 elements long
            compare_out = 2    
        else : 
            compare_out = 0                                                                 #the two lists are different
    
    return compare_out


####### END










####################################################### FUNCTION: updateCounts



def updateCounts(countsList_1, countsList_2) :
    outList = list()
    for k in range(len(countsList_1)) :
        aux_1 = float(countsList_1[k]) + float(countsList_2[k])
        outList.append(aux_1)

    return outList


####### END










####################################################### FUNCTION: updateAnnotation


def updateAnnotation(annList_1, annList_2) :
    outList = list()
    for k in range(len(annList_1)) :
        if annList_1[k] == annList_2[k] :
            outList.append(annList_1[k])
        else :
            aux_1 = updateList_v2(annList_1[k],annList_2[k])
            outList.append(aux_1)

    return outList


####### END










####################################################### FUNCTION: checkCounts


def checkCounts(countsList) :
    counter_aux = 0
    for entry in countsList :
        if float(entry) > 0.5 : counter_aux += 1
    if counter_aux > 1 :
        outFlag = 1
    else :
        outFlag = 0

    return outFlag
 
####### END








######################################################## PROCESS DATA



##### read tximport file
tximport_TRINITY_name_list = list()
tximport_counts_list = list()
Flag_header = 1

for line in tximport_list :
    if Flag_header == 0 :
        tximport_TRINITY_name_list.append(line[0])          # store Trinity gene name
        tximport_counts_list.append(line[1:])               # store counts associated to each gene (for all samples)
    else :
        tximportHeader = line
        Flag_header = 0


##### read annoation file
annoation_ALL = list()
Flag_header = 1

for line in annotation_list :
    if Flag_header == 0 :
        annoation_ALL.append(line)
    else :
        Flag_header = 0





##### Find Trinity entries associated to the same gene; adjust tximport and annotation file accordingly.
##### This is done by sorting gene names (UniProt entries) alphabetically and then comparing each entry with the previous one


### sort annotation file alphabetically

annoation_SORTED_list = list()
annoation_SORTED_list = sorted(annoation_ALL, key=lambda entry:entry[1])              #  sort annotation alphabetically by UniProt's entry name
#for line in annoation_SORTED_list : print('\n', line)

annoation_TRINITY_name_list = list()
annoation_entry_name_list = list()
annotation_geneInfo_list = list()
for entry in annoation_SORTED_list :
        annoation_TRINITY_name_list.append(entry[0])               # store Trinity gene name
        annoation_entry_name_list.append(entry[1])                 # store UniProt entry name associated to each Trinity entry
        annotation_geneInfo_list.append(entry[1:])                 # store info associated to each gene (protein name, Go terms, etc)



### check for presence of several Trinity entries associated to the same gene

Annotation_DIC = dict()                           # dictionary: keeps info about unique Trinity genes (gene duplicates are consolidated into just one entry) 
Counts_DIC = dict()                               # dictionary: consolidated counts
Removed_DIC = dict()                              # dictionary: duplicate entries

geneRepeat_counter = 0 
geneReject_counter = 0 
counter_MAIN = 0
counter_b = 0  

# First dictionary entry
Annotation_DIC[annoation_TRINITY_name_list[0]] =  annotation_geneInfo_list[0]
aux_i = tximport_TRINITY_name_list.index(annoation_TRINITY_name_list[0])                         #index location of Trinity entry in tximport lists
Counts_DIC[tximport_TRINITY_name_list[aux_i]] =  tximport_counts_list[aux_i]
                        
store_Annot = annoation_entry_name_list[0]
store_Trinity = annoation_TRINITY_name_list[0]

for k in range(len(annoation_TRINITY_name_list)) :
     #
     aux_i = tximport_TRINITY_name_list.index(annoation_TRINITY_name_list[k])                 #index location of Trinity entry in tximport lists   
     #
     if (k > 0) :
         #
         counter_MAIN += 1
         counter_b += 1
         if  counter_b == 1000 :
             print(counter_MAIN)                                                              #counter used to keep track of things during long runs
             counter_b = 0
         #
         flagCounts = checkCounts(tximport_counts_list[aux_i])                                    # checks whether there are counts observed for more than one sample
         if (flagCounts == 0) and (annoation_entry_name_list[k] == '.') :                      # filter out non-annotated genes with counts for only one sample
              geneReject_counter += 1
              print(annoation_TRINITY_name_list[k])
              #
              store_Annot = store_Annot
              store_Trinity = store_Trinity
              #
              Removed_DIC[annoation_TRINITY_name_list[k]] =  annotation_geneInfo_list[k]
              #
         elif (annoation_entry_name_list[k] == '.') or (store_Annot == '.') :                   #no annotation info (or comparing with gene with no annotation info)
              Annotation_DIC[annoation_TRINITY_name_list[k]] =  annotation_geneInfo_list[k]
              Counts_DIC[tximport_TRINITY_name_list[aux_i]] =  tximport_counts_list[aux_i]
              #
              store_Annot = annoation_entry_name_list[k]
              store_Trinity = annoation_TRINITY_name_list[k]
         else :
             Flag_sameness = entryCompare(annoation_entry_name_list[k], store_Annot)              #compares list of gene-hits associated to two Trinity entries
             #print(annoation_TRINITY_name_list[k],Flag_sameness)
             if Flag_sameness == 0 :                                                              #gene-hit lists are different >> different gene
                 #
                 Annotation_DIC[annoation_TRINITY_name_list[k]] =  annotation_geneInfo_list[k]
                 Counts_DIC[tximport_TRINITY_name_list[aux_i]] =  tximport_counts_list[aux_i]
                 #
                 store_Annot = annoation_entry_name_list[k]
                 store_Trinity = annoation_TRINITY_name_list[k]
                 #
             elif Flag_sameness == 1 :                                                            #gene-hit lists are exactly the same >> same gene >> must update counts
                 geneRepeat_counter += 1
                 #
                 ind_ref = tximport_TRINITY_name_list.index(store_Trinity) 
                 counts_ref_list = Counts_DIC[tximport_TRINITY_name_list[ind_ref]]
                 newCounts = updateCounts(tximport_counts_list[aux_i], counts_ref_list)
                 Counts_DIC[tximport_TRINITY_name_list[ind_ref]] =  newCounts
                 #print('\n', '\n', tximport_counts_list[aux_i])
                 #print('\n', counts_ref_list)
                 #print('\n', newCounts)
                 #
                 store_Annot = store_Annot
                 store_Trinity = store_Trinity
                 #
                 Removed_DIC[annoation_TRINITY_name_list[k]] =  annotation_geneInfo_list[k]       
             else :                                                                               # gene-hit lists almost the same >> same gene  >> must update counts and annotation info
                 geneRepeat_counter += 1
                 #print(annoation_TRINITY_name_list[k])
                 #
                 ind_ref = tximport_TRINITY_name_list.index(store_Trinity) 
                 counts_ref_list = Counts_DIC[tximport_TRINITY_name_list[ind_ref]]
                 newCounts = updateCounts(tximport_counts_list[aux_i], counts_ref_list)           #computes new counts
                 Counts_DIC[tximport_TRINITY_name_list[ind_ref]] =  newCounts                     #updates counts in dictionary
                 #
                 annot_Ref_list = Annotation_DIC[store_Trinity]
                 newAnnot_list = updateAnnotation(annotation_geneInfo_list[k], annot_Ref_list)
                 Annotation_DIC[store_Trinity] = newAnnot_list                                    #updates annotation info in dictionary
                 #
                 store_Annot = newAnnot_list[0]
                 store_Trinity = store_Trinity  
                 #
                 Removed_DIC[annoation_TRINITY_name_list[k]] =  annotation_geneInfo_list[k] 
                 





### Number or sequences removed
print('\n', 'Number of sequences removed:', geneRepeat_counter)
print('\n', 'Number of sequences rejected:', geneReject_counter)





##### List of genes in final dictionary
trinity_list_final = list()
for entry in Annotation_DIC :
    trinity_list_final.append(entry)
    #print(entry)
   
trinity_list_final = sorted(trinity_list_final)   ## sort alphabetically

             







##### remove duplicate GO entries from Trinity_DIC dictionary
for gene in trinity_list_final :
    entry_aux_list = Annotation_DIC.get(gene ,0)      
    if type(entry_aux_list) == tuple :  entry_aux_list = list(entry_aux_list)        #converts tuple to list
    if entry_aux_list != "0" and entry_aux_list != 0: 
        #
        entry_goBP = entry_aux_list[7]               ##GO-BP
        new_goBP = cleanGO(entry_goBP)
        #
        entry_goMF = entry_aux_list[8]               ##GO-MF
        new_goMF = cleanGO(entry_goMF)
        #
        entry_goCC = entry_aux_list[9]               ##GO-CC
        new_goCC = cleanGO(entry_goCC)
        #
        entry_gopfam = entry_aux_list[10]            ##GO-Pfam
        new_gopfam = cleanGO(entry_gopfam)
        #
        allGo = new_goBP + "^NA`" + new_goMF + "^NA`" + new_goCC + "^NA`" + new_gopfam		## All GO terms
        new_allGo = cleanGO(allGo)
        #
        entry_aux_list = entry_aux_list[:12]         #remove previous GO_clean entries; add the new ones
        entry_aux_list.append(new_goBP)
        entry_aux_list.append(new_goMF)
        entry_aux_list.append(new_goCC)
        entry_aux_list.append(new_gopfam)
        entry_aux_list.append(new_allGo)
        #
        Annotation_DIC[gene] = entry_aux_list         # new columns containing clean GO terms add to each entry
    


#print dictionaries
#for gene in trinity_list_final:
#    aux_dic_list2 = Annotation_DIC.get(gene ,0)
#    print()
#    print(gene, aux_dic_list2)


#for gene in trinity_list_final:
#    aux_dic_list2 = Counts_DIC.get(gene ,0)
#    print()
#    print(gene, aux_dic_list2)


#for gene in Removed_DIC:
#    aux_dic_list2 = Removed_DIC.get(gene ,0)
#    print()
#    print(gene, aux_dic_list2)


#bp()



######################################################## SAVE RESULTS TO FILE


## Consolidated annotation file
outFileName1 = 'trinotate-invertebrates-CONSOLIDATED.txt'
outfile1 = open(outFileName1, 'w') 

aux_save = '"Trinity_gene_ID"' + '\t' + '"Entry_name"' + '\t' + '"Accession_code"' + '\t' + '"Protein_name"' + '\t' + '"Gene_name"' + '\t' + '"Organism_name"' + '\t' + '"KEGG"' + '\t' + '"FlyBase"' + '\t' + '"GO_BioProcess"' + '\t' + '"GO_MolFunction"' + '\t' + '"GO_CelComponent"' + '\t' + '"GO_Pfam"' + '\t' + '"KB_or_trembl"' + '\t' + '"GO_BioProcess_clean"' + '\t' + '"GO_MolFunction_clean"' + '\t' + '"GO_CelComponent_clean"' + '\t' + '"GO_Pfam_clean"' + '\t' + '"All_GO_terms_clean"' + '\n'          #write header

outfile1.write(aux_save)

for gene in trinity_list_final :
    entry_aux_list = Annotation_DIC.get(gene ,0)
    #print('\n', entry_aux_list)
    if entry_aux_list != "" and entry_aux_list != "." and entry_aux_list != 0 :
        aux_save = '"' + gene + '"'
        for k in range(len(entry_aux_list)) :
            #print(entry_aux_list[k])
            if entry_aux_list[k] == '' : entry_aux_list[k] = '.'
            aux_save = aux_save + '\t' + '"' + entry_aux_list[k] + '"'
    else :
        aux_save = '"' + gene + '"'
        for k in range(0,17) :
            aux_save = aux_save + '\t' + '"."'
    aux_save = aux_save + '\n'        
    outfile1.write(aux_save)




## Genes excluded from consolidated annotation file
trinity_list_rejected = list()
for entry in Removed_DIC :
    trinity_list_rejected.append(entry)
    #print(entry)
   
trinity_list_rejected = sorted(trinity_list_rejected)   ## sort alphabetically

outFileName2 = 'OUT_rejectedTrinityGenes.txt'
outfile2 = open(outFileName2, 'w') 

aux_save = '"Trinity_gene_ID"' + '\t' + '"Entry_name"' + '\t' + '"Accession_code"' + '\t' + '"Protein_name"' + '\t' + '"Gene_name"' + '\t' + '"Organism_name"' + '\t' + '"KEGG"' + '\t' + '"FlyBase"' + '\t' + '"GO_BioProcess"' + '\t' + '"GO_MolFunction"' + '\t' + '"GO_CelComponent"' + '\t' + '"GO_Pfam"' + '\t' + '"KB_or_trembl"' + '\t' + '"GO_BioProcess_clean"' + '\t' + '"GO_MolFunction_clean"' + '\t' + '"GO_CelComponent_clean"' + '\t' + '"GO_Pfam_clean"' + '\t' + '"All_GO_terms_clean"' + '\n'          #write header

outfile2.write(aux_save)

for gene in trinity_list_rejected :
    entry_aux_list = Removed_DIC.get(gene ,0)
    #print('\n', entry_aux_list)
    if entry_aux_list != "" and entry_aux_list != "." and entry_aux_list != 0 :
        aux_save = '"' + gene + '"'
        for k in range(len(entry_aux_list)) :
            #print(entry_aux_list[k])
            if entry_aux_list[k] == '' : entry_aux_list[k] = '.'
            aux_save = aux_save + '\t' + '"' + entry_aux_list[k] + '"'
    else :
        aux_save = '"' + gene + '"'
        for k in range(0,17) :
            aux_save = aux_save + '\t' + '"."'
    aux_save = aux_save + '\n'        
    outfile2.write(aux_save)



## Consolidated counts

tximportHeader

outFileName3 = 'tximport-CONSOLIDATED.txt'
outfile3 = open(outFileName3, 'w') 

aux_save = ''
for k in range(len(tximportHeader)) :
    if k == 0 :
        aux_save = '"' + str(tximportHeader[k]) + '"'
    else :
        aux_save = aux_save + '\t' + '"' + str(tximportHeader[k]) + '"'
aux_save = aux_save + '\n'        
outfile3.write(aux_save)                  #write header

for gene in trinity_list_final :
    entry_aux_list = Counts_DIC.get(gene ,0)
    aux_save = '"' + gene + '"'
    for k in range(len(entry_aux_list)) :
        aux_save = aux_save + '\t' + str(entry_aux_list[k])
    aux_save = aux_save + '\n'        
    outfile3.write(aux_save)





outfile1.close()
outfile2.close()
outfile3.close()

print()
print('\n Done!')



