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
# Script used to glean data from UniProt database
# (UniProt file downloaded from ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/)  
#
##################################################################################################################################################################### 




print('\n Gleaning data ...')

######################################################### USAGE

error_message1 = ' \n \
 USAGE: $python3 glean_uniprot_spot_DB.py <uniprot_sprot_file.dat> <organism type>'


# Example: $python3 glean_uniprot_spot_DB.py uniprot_sprot_file.dat Insecta
# In this example, the script gleans records only for organisms within Insecta.


######################################################### LOAD STANDARD MODULES

import sys
import re                                                 # 're' module provides regular expression matching operations
import ast                                                # required to convert reference tree from string to list
import os                                                 # required to create new folders
import time



######################################################## OPEN INPUT FILES

try:
    inputfile1 = open(sys.argv[1], 'r')                  # open uniprot file
except:
    print('\n Error: input files missing.')
    print(error_message1) 
    exit() 

try:
    inputOrg = sys.argv[2]                             # read organism type
    if inputOrg != '' :
        flagAll = 0
    print()
    print('Organisms included: ', inputOrg)
except:
    flagAll = 1
    print('\n All species included.')
    print() 
 



    
######################################################## GLEAN DATA; SAVE DATA TO FILE 





### output file
outFileName = 'sprot_OUT.txt'
outfile1 = open(outFileName, 'w') 
headers = 'ID_Entry_name' + '\t' + 'AC_Accession' + '\t' + 'DE_Protein' + '\t' + 'DE_Protein_alternative_name' + '\t' + 'GN_Gene' + '\t' + 'OS_Organism' + '\t' + 'DR_Ensembl' + '\t' + 'KEGG' + '\t' + 'FlyBase_ID' + '\t' + 'FlyBase_Symbol' + '\t' + 'OrthoDB' + '\t' + 'GO_P_Biological_process' + '\t' + 'GO_F_Molecular_Function' + '\t' + 'GO_C_Cellular_Component' + '\t' + 'Pfam' + '\n'
outfile1.write(headers)


### glean data, one entry at a time

flag_start = 1
counter_main = 0

for line in inputfile1 :
    #print(flag_start)
    counter_main += 1
    line = line.strip()
    line_aux = line[0:2]
    if 	line_aux != '//' :	
        if flag_start == 1 :						
            Entry_name = ''
            Accession = '' 
            Protein = ''
            Protein_alt = ''
            Gene = ''
            Organism = ''
            Ensembl = ''
            KEGG = ''
            FlyBase_ID = ''
            FlyBase_Symbol = ''
            OrthoDB = ''
            GO_P = ''
            GO_F = ''
            GO_C = ''
            Pfam = ''
            flag_start = 0
            FlagSave = 0
            try :
                if line_aux == 'ID' :
                    ename_aux = line[5:]
                    ename_mark =  ename_aux.find(' ')
                    Entry_name = ename_aux[:ename_mark]
            except:
                print('\n Error: Cannot find entry name.')
                print('\n', line, '\n (line ',  counter_main,') \n') 
                exit() 
        else :
            string_aux = 'AC   '
            line_aux2 = line[0:len(string_aux)]
            if line_aux2 == string_aux :
                if Accession == '' : 
                    Accession = line[5:(len(line)-1)]
                else :
                    Accession = Accession + '; ' + line[5:(len(line)-1)]
            #
            string_aux = 'DE   RecName: Full='                            #UniProt-KB
            line_aux2 = line[0:len(string_aux)]
            if line_aux2 == string_aux :
                Protein = line[len(string_aux):(len(line)-1)]
            #
            string_aux = 'DE   AltName: Full='                            #UniProt-KB
            line_aux2 = line[0:len(string_aux)]
            if line_aux2 == string_aux :
                Protein_alt = line[len(string_aux):(len(line)-1)]
            #
            if Protein == '' :                                            #UniProt trembl
                string_aux = 'DE   SubName: Full='
                line_aux2 = line[0:len(string_aux)]
                if line_aux2 == string_aux :
                    Protein = line[len(string_aux):(len(line)-1)]
            #
            string_aux = 'GN   '
            line_aux2 = line[0:len(string_aux)]
            if line_aux2 == string_aux :
                if Gene == '' :
                    Gene = line[len(string_aux):(len(line)-1)]
                else :
                    Gene = Gene + '; ' + line[len(string_aux):(len(line)-1)]
            #
            string_aux = 'OS   '
            line_aux2 = line[0:len(string_aux)]
            if line_aux2 == string_aux :
                Organism = line[len(string_aux):(len(line)-1)]
            #
            string_aux = 'OC   '
            line_aux2 = line[0:len(string_aux)]
            if line_aux2 == string_aux and flagAll == 0 :
                if inputOrg in line :
                    FlagSave = 1
            #
            string_aux = 'DR   Ensembl'
            line_aux2 = line[0:len(string_aux)]
            if line_aux2 == string_aux :
                line_mark =  line.find(';')
                if Ensembl == '' :
                    Ensembl = line[(line_mark+2):(len(line)-1)]
                else :
                    Ensembl = Ensembl + '; ' + line[(line_mark+2):(len(line)-1)]
            #
            string_aux = 'DR   KEGG;'
            line_aux2 = line[0:len(string_aux)]
            if line_aux2 == string_aux :
                KEGG = line[(len(string_aux)+1):(len(line)-1)]
            #
            string_aux = 'DR   FlyBase; '
            line_aux2 = line[0:len(string_aux)]
            if line_aux2 == string_aux :
                line_aux3 = line[len(string_aux):]
                line_mark =  line_aux3.find(';')
                FlyBase_ID = line_aux3[:line_mark]
                FlyBase_Symbol = line_aux3[(line_mark+2):(len(line_aux3)-1)]
            #
            string_aux = 'DR   OrthoDB;'
            line_aux2 = line[0:len(string_aux)]
            if line_aux2 == string_aux :
                OrthoDB = line[(len(string_aux)+1):(len(line)-1)]
             #
            string_aux = 'DR   GO;'
            line_aux2 = line[0:len(string_aux)]
            if line_aux2 == string_aux :
                line_aux3 = line[len(string_aux):]
                line_mark = line_aux3.find(';')
                line_aux4 = line_aux3[(line_mark+2):(line_mark+4)]
                if line_aux4 == 'P:' :
                    if GO_P == '' :
                        #GO_P = line[(len(string_aux)+1):(len(string_aux)+11)]
                        GO_P = line[(len(string_aux)+1):]
                        #print(Entry_name,GO_P)
                    else :
                        #GO_P = GO_P + "^NA`" + line[(len(string_aux)+1):(len(string_aux)+11)]
                        GO_P = GO_P + "`" + line[(len(string_aux)+1):]
                if line_aux4 == 'F:' :
                    if GO_F == '' :
                        #GO_F = line[(len(string_aux)+1):(len(string_aux)+11)]
                        GO_F = line[(len(string_aux)+1):]
                    else :
                        #GO_F = GO_F + "^NA`" + line[(len(string_aux)+1):(len(string_aux)+11)]
                        GO_F = GO_F + "`" + line[(len(string_aux)+1):]
                if line_aux4 == 'C:' :
                    if GO_C == '' :
                        #GO_C = line[(len(string_aux)+1):(len(string_aux)+11)]
                        GO_C = line[(len(string_aux)+1):]
                    else :
                        #GO_C = GO_C + "^NA`" + line[(len(string_aux)+1):(len(string_aux)+11)]
                        GO_C = GO_C + "`" + line[(len(string_aux)+1):]
            #
            string_aux = 'DR   Pfam;'
            line_aux2 = line[0:len(string_aux)]
            if line_aux2 == string_aux :
                line_aux3 = line[(len(string_aux)+1):]
                line_mark = line_aux3.find(';')
                if Pfam == '' :
                    Pfam = line_aux3[:line_mark]  
                else :
                    Pfam = Pfam + ';' + line_aux3[:line_mark]          
    else :
        aux_summary = Entry_name + '\t' + Accession + '\t' + Protein + '\t' + Protein_alt + '\t' + Gene + '\t' + Organism + '\t' + Ensembl + '\t' + KEGG + '\t' + FlyBase_ID + '\t' + FlyBase_Symbol + '\t' + OrthoDB + '\t' + GO_P + '\t' + GO_F + '\t' + GO_C + '\t' + Pfam + '\n'
        if flagAll == 1 or (flagAll == 0 and FlagSave == 1) :
            outfile1.write(aux_summary)
        flag_start = 1
    
print('\n Done! \n')    
outfile1.close()									#close output file    



