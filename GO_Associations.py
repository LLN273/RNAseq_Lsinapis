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
# Script used to associate GO terms to each transcript/gene
#
##################################################################################################################################################################### 

##################################################################################################################################################################### 
#
#
# Input files: 1 -> Trinotate annotation report, in csv format
#              2 -> sprot_OUT_Insecta.txt file (based on uniprot_sprot_invertebrates.dat.gz (Knowledge Base) reference file, filtered for Insecta only)
#              3 -> trembl_OUT_Insecta.txt file (based on uniprot_trembl_invertebrates.dat.gz (trembl) reference file, filtered for Insecta only)   
#
# Output: 1 files -> file associating each Trinity gene entry to a known gene and its associated GO terms
#           
#
# Strategy:    
# 1. Associate Drosophila and lepidopteran hits to GO-BioProcess terms, using UniProt-sprot db
# 2. If 1. doesn't work, use UniProt-trembl
# 3. If 1-2 don't work, just get the gene name, as well as GO-MolFunction and GO-CellComponent (as well as Pfam GO terms)
# 4. Otherwise, file Trinity gene as a new gene
#
#
#
#
###

print('\n Gleaning data ...')

######################################################### USAGE

error_message1 = ' \n \
 USAGE: $python3 GO_Associations.py <trinotate_annotation_report-invertebrates.csv> <sprot_OUT_Insecta.txt> <sprot_OUT_trembl_Insecta.txt> <results_OUT.csv> \n'


######################################################### LOAD STANDARD MODULES

import sys
import re                                                 # 're' module provides regular expression matching operations
import ast                                                # required to convert reference tree from string to list
import os                                                 # required to create new folders
import time
import csv                                                #required to open csv file

from pdb import set_trace as bp                           # when using break points during debugging >> useage: bp()

######################################################## OPEN INPUT FILES	


try:
    trinotate_list = list()                                  #open Trinotate output file as csv
    with open(sys.argv[1], newline='') as csvfile1:         
        inputfile_1 = csv.reader(csvfile1, delimiter='\t')
        for row in inputfile_1:
            trinotate_list.append(row)
#            print(', '.join(row))
except:
    print('\n Error[1]: input file missing.')
    print(error_message1) 
    exit() 



try:   
    sprot_inv_list = list()                                    #open sprot_OUT_invertebrates file as csv
    with open(sys.argv[2], newline='') as csvfile2:         
        inputfile_4 = csv.reader(csvfile2, delimiter='\t')
        for row in inputfile_4:
            sprot_inv_list.append(row)
except:
    print('\n Error[2]: input file missing.')
    print(error_message1) 
    exit() 



try:   
    sprot_trembl_list = list()                                    #open trembl_OUT_invertebrates (Knowledge Base) file as csv
    with open(sys.argv[3], newline='') as csvfile3:         
        inputfile_5 = csv.reader(csvfile3, delimiter='\t')
        for row in inputfile_5:
            sprot_trembl_list.append(row)
except:
    print('\n Error[3]: input file missing.')
    print(error_message1) 
    exit() 



try:
    outFileName_user = sys.argv[4]                                         # read output file name
    if outFileName_user == '' :
        print('\n Error[4]: input file missing.')
        print(error_message1) 
        exit() 
except:
    print('\n Error[4]: input file missing.')
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




####################################################### FUNCTION: updateList

def updateList(referenceS, newS) :
    if newS == '.' or newS == '' : 
        aux_output = referenceS
    elif referenceS == '.' or referenceS == '' :
        aux_output = newS
    elif newS not in referenceS :
        aux_output = referenceS + "`" + newS
    else:
        aux_output = referenceS

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







######################################################## PROCESS DATA




## read sprot_OUT_invertebrates file
## note: each Entry_name (UniProt-KB) entry may be associated to more than one Accession code (UniProt)
sprotKB_Entry_name_list = list()
sprotKB_Accession_list = list()
sprotKB_ProteinName_list = list()
sprotKB_GeneName_list = list()
sprotKB_Organism = list()
sprotKB_KEGG_list = list()
sprotKB_FlyBase_ID_list = list()
sprotKB_GO_BioProcess_list = list()
sprotKB_GO_MolFunction_list = list()
sprotKB_GO_CelComponent_list = list()
sprotKB_Pfam_list = list()
for line in sprot_inv_list :
    aux_A = line[1]                        # reads Accession code(s)
    aux_A_list = aux_A.split(';')          # gets individual Accession codes (if there is more than one)
    aux_en = line[0]                       # reads Entry_name
    aux_pn = line[2]                       # reads Protein Name
    aux_gn = line[4]                       # reads Gene Name
    aux_o = line[5]                        # reads Organism name
    aux_k = line[7]                        # reads KEGG
    aux_fb = line[8]                       # reads FlyBase ID
    aux_goBP = line[11]                    # reads GO_BioProcess codes
    aux_goMF = line[12]                    # reads GO_MolFunction codes
    aux_goCC = line[13]                    # reads GO_CelComponent codes
    aux_pf = line[14]                      # reads pfam code(s)
    #
    aux_2 = aux_gn[:5]                     #if there is a clear gene name, (line starts with 'Name=', get it)
    if aux_2 == "Name=" :
        aux_3 = aux_gn[5:]
        mark_gn = aux_3.find(' ')
        if mark_gn == -1 :
            aux_gn = aux_3
        else : 
             aux_gn = aux_3[:(mark_gn-1)]
    #
    for code in aux_A_list : 
        sprotKB_Accession_list.append(code)
        sprotKB_Entry_name_list.append(aux_en)
        sprotKB_ProteinName_list.append(aux_pn)
        sprotKB_GeneName_list.append(aux_gn)
        sprotKB_Organism.append(aux_o)
        sprotKB_KEGG_list.append(aux_k)
        sprotKB_FlyBase_ID_list.append(aux_fb)
        sprotKB_GO_BioProcess_list.append(aux_goBP)
        sprotKB_GO_MolFunction_list.append(aux_goMF)
        sprotKB_GO_CelComponent_list.append(aux_goCC)
        sprotKB_Pfam_list.append(aux_pf)
        #print(aux_en, code,aux_pn,aux_gn,aux_o,aux_k,aux_fb,aux_goBP,aux_pf)



## read trembl_invertebrates file
## note: each Entry_name (UniProt-trembl) entry may be associated to more than one Accession code (UniProt)
sprot_trembl_Entry_name_list = list()
sprot_trembl_Accession_list = list()
sprot_trembl_ProteinName_list = list()
sprot_trembl_GeneName_list = list()
sprot_trembl_Organism = list()
sprot_trembl_KEGG_list = list()
sprot_trembl_FlyBase_ID_list = list()
sprot_trembl_GO_BioProcess_list = list()
sprot_trembl_GO_MolFunction_list = list()
sprot_trembl_GO_CelComponent_list = list()
sprot_trembl_Pfam_list = list()
for line in sprot_trembl_list :
    aux_A = line[1]                        # reads Accession code(s)
    aux_A_list = aux_A.split(';')          # gets individual Accession codes (if there is more than one)
    aux_en = line[0]                       # reads Entry_name
    aux_pn = line[2]                       # reads Protein Name
    aux_gn = line[4]                       # reads Gene Name
    aux_o = line[5]                        # reads Organism name
    aux_k = line[7]                        # reads KEGG
    aux_fb = line[8]                       # reads FlyBase ID
    aux_goBP = line[11]                    # reads GO_BioProcess codes
    aux_goMF = line[12]                    # reads GO_MolFunction codes
    aux_goCC = line[13]                    # reads GO_CelComponent codes
    aux_pf = line[14]                      # reads pfam code(s)
    #
    aux_2 = aux_gn[:5]                     #if there is a clear gene name, (line starts with 'Name=', get it)
    if aux_2 == "Name=" :
        aux_3 = aux_gn[5:]
        mark_gn = aux_3.find(' ')
        if mark_gn == -1 :
            aux_gn = aux_3
        else : 
             aux_gn = aux_3[:(mark_gn-1)]
    #
    for code in aux_A_list : 
        sprot_trembl_Accession_list.append(code)
        sprot_trembl_Entry_name_list.append(aux_en)
        sprot_trembl_ProteinName_list.append(aux_pn)
        sprot_trembl_GeneName_list.append(aux_gn)
        sprot_trembl_Organism.append(aux_o)
        sprot_trembl_KEGG_list.append(aux_k)
        sprot_trembl_FlyBase_ID_list.append(aux_fb)
        sprot_trembl_GO_BioProcess_list.append(aux_goBP)
        sprot_trembl_GO_MolFunction_list.append(aux_goMF)
        sprot_trembl_GO_CelComponent_list.append(aux_goCC)
        sprot_trembl_Pfam_list.append(aux_pf)
        #print(aux_en, code,aux_pn,aux_gn,aux_o,aux_k,aux_fb,aux_goBP,aux_pf)






### extract Trinity names ; extract associated gene names found by Trinotate; store in dictionary

Trinity_DIC = dict()             # dictionary: associates each Trinity gene to a Drosophila or Lep gene (or another, if no Dro/Lep gene found) 
Trinity_DIC_shadow = dict()      # keeps track of gene hits not associated to GO-BioProc terms

trinity_list = list()            #list of all individual Trinity genes

#org_List = ("inv","inv","Drosophila melanogaster", "Bombyx mori", "Danaus plexippus", "Papilio machaon","Drosophila melanogaster", "Bombyx mori", "Danaus plexippus", "Papilio machaon")                #list of species used during BLASTX and BLASTP

flag1 = 0                         #flag used to skip header
 

for line in trinotate_list:
    #
    gene_REF = "."
    geneBasket_list = list()         #list of hits (genes) 
    if flag1 == 1 :
        #
        trinity_name = line[0]                                  #Trinity gene name, as listed in Trinotate
        if trinity_name not in trinity_list :
            trinity_list.append(trinity_name)                   #saves Trinity name to list
        #
        sprot_Top_BLASTX_hit = line[2]
        if sprot_Top_BLASTX_hit != "." :
            mark_g = sprot_Top_BLASTX_hit.find('^')
            if mark_g != -1 :
                sprot_Top_BLASTX_hit = sprot_Top_BLASTX_hit[:mark_g]
        geneBasket_list.append(sprot_Top_BLASTX_hit)
        #
        sprot_Top_BLASTP_hit = line[6]
        if sprot_Top_BLASTP_hit != "." :
            mark_g = sprot_Top_BLASTP_hit.find('^')
            if mark_g != -1 :
                sprot_Top_BLASTP_hit = sprot_Top_BLASTP_hit[:mark_g]
        geneBasket_list.append(sprot_Top_BLASTP_hit)
        #
        Drosophila_Uniprot_blastx_hit = line[10]
        Drosophila_Uniprot_blastx_hit = readGene(Drosophila_Uniprot_blastx_hit)
        geneBasket_list.append(Drosophila_Uniprot_blastx_hit)
        #
        Bombyx_blastx_hit = line[14]
        Bombyx_blastx_hit = readGene(Bombyx_blastx_hit)
        geneBasket_list.append(Bombyx_blastx_hit)
        #
        Danaus_blastx_hit = line[16]
        Danaus_blastx_hit = readGene(Danaus_blastx_hit)
        geneBasket_list.append(Danaus_blastx_hit)
        #
        Papilio_blastx_hit = line[18]
        Papilio_blastx_hit = readGene(Papilio_blastx_hit)
        geneBasket_list.append(Papilio_blastx_hit)
        #
        Drosophila_Uniprot_blastp_hit = line[21]
        Drosophila_Uniprot_blastp_hit = readGene(Drosophila_Uniprot_blastp_hit)
        geneBasket_list.append(Drosophila_Uniprot_blastp_hit)
        #
        Bombyx_blastp_hit = line[25]
        Bombyx_blastp_hit = readGene(Bombyx_blastp_hit)
        geneBasket_list.append(Bombyx_blastp_hit)
        #
        Danaus_blastp_hit = line[27]
        Danaus_blastp_hit = readGene(Danaus_blastp_hit)
        geneBasket_list.append(Danaus_blastp_hit)
        #
        Papilio_blastp_hit = line[29]
        Papilio_blastp_hit = readGene(Papilio_blastp_hit)
        geneBasket_list.append(Papilio_blastp_hit)
        #
        Pfam_Trinotate = line[31]
        if Pfam_Trinotate != "." :
            mark_g = Pfam_Trinotate.find('^')
            if mark_g != -1 :
                Pfam_Trinotate = Pfam_Trinotate[:mark_g]
        #
        Kegg_Trinotate = line[35]
        gene_ontology_blast_Trinotate = line[36]
        gene_ontology_pfam_Trinotate = line[37]
        #
        flag_GO = 0                                                 # indicates whether GO-BioProc terms have been found for the current Trinity gene (0: no; 1: yes)
        #print(trinity_name, geneBasket_list)
        for k in range(len(geneBasket_list)) :
            if k==0 or k==1 : continue                                      # consider only Drosophila and Lep hits
            if geneBasket_list[k] == "." : continue
            if geneBasket_list[k] in sprotKB_Accession_list :               # initial search done on UniProt-KB
                aux_i = sprotKB_Accession_list.index(geneBasket_list[k])
                hit_en = sprotKB_Entry_name_list[aux_i]
                hit_pn = sprotKB_ProteinName_list[aux_i]
                hit_gn = sprotKB_GeneName_list[aux_i]
                hit_on = sprotKB_Organism[aux_i]
                hit_kl = sprotKB_KEGG_list[aux_i]
                hit_fb = sprotKB_FlyBase_ID_list[aux_i]
                hit_goBP = sprotKB_GO_BioProcess_list[aux_i]
                hit_goMF = sprotKB_GO_MolFunction_list[aux_i]                
                hit_goCC = sprotKB_GO_CelComponent_list[aux_i]                
                #hit_pf = sprotKB_Pfam_list[aux_i]
                hit_pf = gene_ontology_pfam_Trinotate                                     #Use Trinotate Pfam hits
                if hit_goBP != '.' and hit_goBP != '' :
                    flag_GO = 1
                    if trinity_name not in Trinity_DIC :
                        Trinity_DIC[trinity_name] = (hit_en, geneBasket_list[k], hit_pn, hit_gn, hit_on, hit_kl, hit_fb, hit_goBP, hit_goMF, hit_goCC, hit_pf, "KB")
                    else :
                        aux_dic_list2 = Trinity_DIC.get(trinity_name ,0)
                        hit_en = updateList(aux_dic_list2[0], hit_en)
                        hit_ac = updateList(aux_dic_list2[1], geneBasket_list[k])                    
                        hit_pn = updateList(aux_dic_list2[2], hit_pn)
                        hit_gn = updateList(aux_dic_list2[3], hit_gn)
                        hit_on = updateList(aux_dic_list2[4], hit_on)
                        hit_kl = updateList(aux_dic_list2[5], hit_kl)
                        hit_fb = updateList(aux_dic_list2[6], hit_fb)
                        hit_goBP = updateList(aux_dic_list2[7], hit_goBP)
                        hit_goMF = updateList(aux_dic_list2[8], hit_goMF)
                        hit_goCC = updateList(aux_dic_list2[9], hit_goCC)
                        hit_pf = updateList(aux_dic_list2[10], hit_pf)
                        Trinity_DIC[trinity_name] = (hit_en, hit_ac, hit_pn, hit_gn, hit_on, hit_kl, hit_fb, hit_goBP, hit_goMF, hit_goCC, hit_pf, "KB")
                else :
                    if flag_GO == 0 :                                      #use shadow dictionary only if no GO terms have been found so far
                        if trinity_name not in Trinity_DIC_shadow :
                            Trinity_DIC_shadow[trinity_name] = (hit_en, geneBasket_list[k], hit_pn, hit_gn, hit_on, hit_kl, hit_fb, hit_goBP, hit_goMF, hit_goCC, hit_pf, "GO-BP_missing")
                        else :
                            aux_dic_list2 = Trinity_DIC_shadow.get(trinity_name ,0)
                            hit_en = updateList(aux_dic_list2[0], hit_en)
                            hit_ac = updateList(aux_dic_list2[1], geneBasket_list[k])                    
                            hit_pn = updateList(aux_dic_list2[2], hit_pn)
                            hit_gn = updateList(aux_dic_list2[3], hit_gn)
                            hit_on = updateList(aux_dic_list2[4], hit_on)
                            hit_kl = updateList(aux_dic_list2[5], hit_kl)
                            hit_fb = updateList(aux_dic_list2[6], hit_fb)
                            hit_goBP = updateList(aux_dic_list2[7], hit_goBP)
                            hit_goMF = updateList(aux_dic_list2[8], hit_goMF)
                            hit_goCC = updateList(aux_dic_list2[9], hit_goCC)
                            hit_pf = updateList(aux_dic_list2[10], hit_pf)
                            Trinity_DIC_shadow[trinity_name] = (hit_en, hit_ac, hit_pn, hit_gn, hit_on, hit_kl, hit_fb, hit_goBP, hit_goMF, hit_goCC, hit_pf, "GO-BP_missing")      
        #
        if flag_GO == 0 :                                # if no GO terms have been found so far (using UniProt-KB), switch to UniProt-trembl
            for k in range(len(geneBasket_list)) :
                if geneBasket_list[k] in sprot_trembl_Accession_list :               
                    aux_i = sprot_trembl_Accession_list.index(geneBasket_list[k])
                    hit_en = sprot_trembl_Entry_name_list[aux_i]
                    hit_pn = sprot_trembl_ProteinName_list[aux_i]
                    hit_gn = sprot_trembl_GeneName_list[aux_i]
                    hit_on = sprot_trembl_Organism[aux_i]
                    hit_kl = sprot_trembl_KEGG_list[aux_i]
                    hit_fb = sprot_trembl_FlyBase_ID_list[aux_i]
                    hit_goBP = sprot_trembl_GO_BioProcess_list[aux_i]
                    hit_goMF = sprot_trembl_GO_MolFunction_list[aux_i]                    
                    hit_goCC = sprot_trembl_GO_CelComponent_list[aux_i]                    
                    #hit_pf = sprotKB_Pfam_list[aux_i]
                    hit_pf = gene_ontology_pfam_Trinotate                                     #Use Trinotate Pfam hits
                    #print(trinity_name,geneBasket_list[k],hit_goBP)
                    if hit_goBP != '.' and hit_goBP != '' :
                        flag_GO = 1
                        if trinity_name not in Trinity_DIC :
                            Trinity_DIC[trinity_name] = (hit_en, geneBasket_list[k], hit_pn, hit_gn, hit_on, hit_kl, hit_fb, hit_goBP, hit_goMF, hit_goCC, hit_pf, "trembl")
                        else :
                            aux_dic_list2 = Trinity_DIC.get(trinity_name ,0)
                            hit_en = updateList(aux_dic_list2[0], hit_en)
                            hit_ac = updateList(aux_dic_list2[1], geneBasket_list[k])                    
                            hit_pn = updateList(aux_dic_list2[2], hit_pn)
                            hit_gn = updateList(aux_dic_list2[3], hit_gn)
                            hit_on = updateList(aux_dic_list2[4], hit_on)
                            hit_kl = updateList(aux_dic_list2[5], hit_kl)
                            hit_fb = updateList(aux_dic_list2[6], hit_fb)
                            hit_goBP = updateList(aux_dic_list2[7], hit_goBP)
                            hit_goMF = updateList(aux_dic_list2[8], hit_goMF)
                            hit_goCC = updateList(aux_dic_list2[9], hit_goCC)
                            hit_pf = updateList(aux_dic_list2[10], hit_pf)
                            Trinity_DIC[trinity_name] = (hit_en, hit_ac, hit_pn, hit_gn, hit_on, hit_kl, hit_fb, hit_goBP, hit_goMF, hit_goCC, hit_pf, "trembl")
                    else :
                        if flag_GO == 0 :                                      #use shadow dictionary only if no GO terms have been found so far
                            if trinity_name not in Trinity_DIC_shadow :
                                Trinity_DIC_shadow[trinity_name] = (hit_en, geneBasket_list[k], hit_pn, hit_gn, hit_on, hit_kl, hit_fb, hit_goBP, hit_goMF, hit_goCC, hit_pf, "trembl_GO-BP_missing")
                            else :
                                aux_dic_list2 = Trinity_DIC_shadow.get(trinity_name ,0)
                                hit_en = updateList(aux_dic_list2[0], hit_en)
                                hit_ac = updateList(aux_dic_list2[1], geneBasket_list[k])                    
                                hit_pn = updateList(aux_dic_list2[2], hit_pn)
                                hit_gn = updateList(aux_dic_list2[3], hit_gn)
                                hit_on = updateList(aux_dic_list2[4], hit_on)
                                hit_kl = updateList(aux_dic_list2[5], hit_kl)
                                hit_fb = updateList(aux_dic_list2[6], hit_fb)
                                hit_goBP = updateList(aux_dic_list2[7], hit_goBP)
                                hit_goMF = updateList(aux_dic_list2[8], hit_goMF)
                                hit_goCC = updateList(aux_dic_list2[9], hit_goCC)
                                hit_pf = updateList(aux_dic_list2[10], hit_pf)
                                Trinity_DIC_shadow[trinity_name] = (hit_en, hit_ac, hit_pn, hit_gn, hit_on, hit_kl, hit_fb, hit_goBP, hit_goMF, hit_goCC, hit_pf, "trembl_GO-BP_missing")       
    else :
        flag1 = 1   
                      




##### Consolidate dictionaries
for gene in trinity_list:
    if gene not in Trinity_DIC :
        Trinity_DIC[gene] = Trinity_DIC_shadow.get(gene ,0)




### remove duplicate GO entries from Trinity_DIC dictionary
for gene in trinity_list :
    entry_aux_list = Trinity_DIC.get(gene ,0)      
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
        entry_aux_list.append(new_goBP)
        entry_aux_list.append(new_goMF)
        entry_aux_list.append(new_goCC)
        entry_aux_list.append(new_gopfam)
        entry_aux_list.append(new_allGo)
        #
        Trinity_DIC[gene] = entry_aux_list         # new columns containing clean GO terms add to each entry
    


#print dictionaries
#for gene in trinity_list:
#    aux_dic_list2 = Trinity_DIC.get(gene ,0)
#    print()
#    print(gene, aux_dic_list2)







######################################################## SAVE RESULTS TO FILE


outFileName1 = outFileName_user
outfile1 = open(outFileName1, 'w') 


aux_save = '"Trinity_gene_ID"' + '\t' + '"Entry_name"' + '\t' + '"Accession_code"' + '\t' + '"Protein_name"' + '\t' + '"Gene_name"' + '\t' + '"Organism_name"' + '\t' + '"KEGG"' + '\t' + '"FlyBase"' + '\t' + '"GO_BioProcess"' + '\t' + '"GO_MolFunction"' + '\t' + '"GO_CelComponent"' + '\t' + '"GO_Pfam"' + '\t' + '"KB_or_trembl"' + '\t' + '"GO_BioProcess_clean"' + '\t' + '"GO_MolFunction_clean"' + '\t' + '"GO_CelComponent_clean"' + '\t' + '"GO_Pfam_clean"' + '\t' + '"All_GO_terms_clean"' + '\n'          #write header

outfile1.write(aux_save)

for gene in trinity_list :
    entry_aux_list = Trinity_DIC.get(gene ,0)
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


outfile1.close()



print()
print('\n Done!')



