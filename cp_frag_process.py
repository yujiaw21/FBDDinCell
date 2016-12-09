#! /usr/bin/python 

import sys,os,copy,math

input_1 = open(sys.argv[1])

dic = {}

for line in input_1.readlines():
        if line[0] =="-":
                compound = line[1:].strip()
        if line != "\n" and line[0]!="-":
                experiment_path= line.strip()

                if compound not in dic.keys():
                        dic[compound] = [experiment_path]

                else:
                        dic[compound].append(experiment_path)


reactivity_results = {}
tryp_clean = []

#prot_pep_dic = {} # this needs to be revised 

prot_name = {}
#print dic 

master_dic={}
exp_count = 0
compound_counter = 0 
header_list=[]

for compound in dic.keys():
	header_list.append(compound)
	compound_counter +=1 	
        #out_notryp = open("no_tryp.txt","w") #fis this to change with the name of the file open
        prot_pep_dic = {}
	for experiment in dic[compound]:
                exp_count +=1
                #print experiment 
                file_imp = open(experiment,"r")
                element_holder = []
                tryp_clean = []

                #print file_imp
                for line in file_imp.readlines():
                        if line[0] != "i" and line[0] != " ":
                                #print "line", line.split()
                                data = line.split()
                                #print element_holder
                                element_holder_copy = copy.deepcopy(element_holder)
                                if len(element_holder) > 0:
                                        #check and remove tryptics

                                        for peptide in element_holder_copy:
                                                isotope = 0
                                                peptide_split = peptide.split("\t")
                                                peptide_seq = peptide_split[4]
                                                #check half tryptic
                                                #print peptide_seq

                                                if peptide_seq[0] not in ["K","R","-"]: # nterminus
                                                        element_holder.remove(peptide)
                                                        #print "half",peptide_seq
                                                if  peptide_seq[-3] not in ["K","R"] and peptide_seq[-1] != "-": # cterminus    
                                                        #print "half",peptide_seq
                                                        if peptide in element_holder:
                                                                element_holder.remove(peptide)

                                                #remove peptides that should not have a silac ratio attributed since they don't have K or R
                                                peptide_core = peptide_seq.split('.')[1]
                                                for j in ["K","R"]:

                                                        if j in peptide_core:
                                                                isotope = 1
                                                        #print j, peptide_core, str(isotope)    


                                                if isotope == 0 :
                                                        #print "Remove",peptide
                                                        if peptide in element_holder:
                                                                #print "Remove_2",peptide
                                                                element_holder.remove(peptide)


                                                #if  peptide_seq[-3] == "*" and peptide_seq[-1] != "-": # cterminus
                                                #       element_holder.remove(peptide)
                                                #       print "half", peptide_seq 


                                        if len(element_holder)>0:
                                                tryp_clean.append(master_peptide_info)
                                                master_info_split = master_peptide_info.split("\t")
                                                up = master_info_split[1]
                                                name = master_info_split[3]
                                                description = master_info_split[2]
					
                                                prot_name[up] = [name,description]
                                                if up not in prot_pep_dic.keys():
                                                        prot_pep_dic[up] = []


                                                for elem in element_holder:
                                                        elem_split = elem.split("\t")
                                                        ratio = elem_split[6]
                                                        seq = elem_split[4]
                                                        step = elem_split[len(elem_split)-2]
                                                        charge= elem_split[len(elem_split)-3]

                                                        #add step and charge 
                                                        #print elem_split
                                                        tryp_clean.append(elem)
                                                        if float(ratio) > 0:
                                                                prot_pep_dic[up].append((seq,ratio,str(exp_count),charge,step))

                                element_holder = []
                                master_peptide_info = line
                        #else:
                        if line[0] != "i" and line[0] == " ":
                                element_holder.append(line)



	#print prot_pep_dic
	cleaned_prot_pep_dic = {}
	for entry in prot_pep_dic.keys():
        	peptide_pool = []
        	if len(prot_pep_dic[entry]) > 2:
                	for peptide in prot_pep_dic[entry]:
				if type(peptide) is tuple:
					#print peptide,entry,compound,prot_pep_dic[entry]
                       			peptide_seq = peptide[0]
                        		peptide_seq_clean = peptide_seq.split('*')
                        		peptide_seq_clean_b = "".join(peptide_seq_clean)
                        		peptide_seq_clean_c = peptide_seq_clean_b.split("+")
                        		peptide_seq_clean_d ="".join(peptide_seq_clean_c)

                        		if peptide_seq_clean_d not in peptide_pool:
                                		peptide_pool.append(peptide_seq_clean_d)

                	if len(peptide_pool) >2:
                        	cleaned_prot_pep_dic[entry] = prot_pep_dic[entry]
                        	#cleaned_prot_pep_dic[entry].append(len(peptide_pool))

	exp_count = 0
	for experiment in dic[compound]:
		exp_count += 1 
		file_imp = open(experiment,"r")
		
		for line in file_imp.readlines():
                        if line[0] != "i" and line[0] != " ":
                                line_split= line.split("\t")
                                #print line_split
                                ratio = line_split[6]
                                sd = line_split[7]
                                prot_up_b = line_split[1]
				if prot_up_b in cleaned_prot_pep_dic.keys():
					if isinstance(float(sd),float):
                                        #print "I'm in "        
                                                if float(sd) < 55: #number bigger than 0
                                                        if prot_up_b in master_dic.keys(): 
								if len(master_dic[prot_up_b]) == compound_counter-1:
									master_dic[prot_up_b].append([ratio])
								if len(master_dic[prot_up_b]) == compound_counter and exp_count>1:
									#print compound_counter
									#print master_dic[prot_up_b]
									master_dic[prot_up_b][compound_counter-1].append(ratio)
							
								
								if len(master_dic[prot_up_b]) < compound_counter-1:
									missing_fields = (compound_counter-1) -len(master_dic[prot_up_b])
									
									for k in range(missing_fields):
										master_dic[prot_up_b].append(["-"])
									master_dic[prot_up_b].append([ratio])		
 							
							if prot_up_b not in master_dic.keys():
								master_dic[prot_up_b] = []

								if len(master_dic[prot_up_b]) == compound_counter-1:
                                                                        master_dic[prot_up_b].append([ratio]) 	
								
								if len(master_dic[prot_up_b]) < compound_counter-1: 
                                                                        missing_fields = (compound_counter-1) -len(master_dic[prot_up_b])                                                                                                  
									for k in range(missing_fields):
                                                                                master_dic[prot_up_b].append(["-"])             
                                                                        master_dic[prot_up_b].append([ratio])


#print master_dic 

for key in master_dic.keys():
	missing_length = compound_counter - len(master_dic[key])

	for j in range(missing_length):
		master_dic[key].append(["-"]) 	


final_header= "Identifier\tprot name\t"+"\t".join(header_list)
print final_header

text_table = open("entry_table.txt","w")

text_table.write(final_header+"\n")
for key in master_dic.keys():
	entry_str= ""
	
	for k in master_dic[key]:
		entry_str+= " ".join(k)+"\t"
	final_str = key+"\t"+prot_name[key][1]+"\t"+entry_str[:-1]+"\n"
	text_table.write(final_str)
