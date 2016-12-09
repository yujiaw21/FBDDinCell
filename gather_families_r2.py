#! /usr/bin/python

import sys

input_1 = open(sys.argv[1]) #uniprot database

input_2 = open(sys.argv[2]) #data 

kw_annot ={}
enz_annot = {}
sequence_dict = {}
go_annot_dict = {}

KW_terms = ["DNA-binding","Transcription","Receptor","Porin","Ion transport","Ion channel","Transport","Receptor","SH2 domain","SH3 domain", "Nucleotide-binding", "RNA-binding", "G-protein coupled receptor", "Transcription regulation","Chaperone","Receptor","Repressor","Lipid-binding","Heme","Fatty acid metabolism","Lipid metabolism","Cytoskeleton","Microtubule"]

GO_terms = ["protein binding","enzyme binding","enzyme regulator activity","signaling pathway","lipid binding", "lipid transport","RNA binding","DNA binding","carbohydrate binding", "heme biosynthetic process", "small molecule biosynthetic process", "cytochrome c oxidase activity", "signal transduction", "signal transducer activity","structural molecule activity", "regulation of catalytic activity", "protein targeting", "protein-plasma membrane targeting", "sphingosine N-acyltransferase activity","cell adhesion" ]


ac_line=0
line_counter = 0
seq_read = 0
for line in input_1.readlines():
        line_counter += 1

        if line[0:2] == "AC":
                #print line.strip() 
                ac_split = line[4:].split(";")

                if line_counter == ac_line +1:
                        #print "ac_secon", ac_split
                        ac_split_clean = ac_split[:-1]
                        if ac in kw_annot.keys():
                                for code in  ac_split_clean:
                                        kw_annot[ac][0].append(code.strip())
                                        enz_annot[ac][0].append(code.strip())
                                        sequence_dict[ac][0].append(code.strip())
       					go_annot_dict[ac][0].append(code.strip())
				ac_line = line_counter
		  	else:
                           	print "MAYDAY MAYDAY WERE GOING To Crash"


                else:# its the first line of entries
                        ac = ac_split[0].strip()
                        ac_list = ac_split[1:-1]
                        ac_list_clean = [k.strip() for k in ac_list ]

                        #print "ac pimary",ac ,ac_list_clean 
                        ac_line = line_counter

                        if len(ac_list_clean) == 0 and  ac not in kw_annot.keys() :
                                kw_annot[ac] = [[]]
                                enz_annot[ac] = [[]]
                                sequence_dict[ac] = [[]]
				go_annot_dict[ac] =[[]]				

                        if len(ac_list_clean) > 0 and  ac not in kw_annot.keys() :
                                kw_annot[ac] = [ac_list_clean]
                                enz_annot[ac] = [ac_list_clean]
                                sequence_dict[ac] = [ac_list_clean]
				go_annot_dict[ac] = [ac_list_clean]

                        #if len(ac_list_clean)> 0 and ac in up_annot.keys():
                        #       for k in ac_list_clean:
                        #               up_annot[ac][0].append(k)
				

        if line [0:2] == "DE":
		line_split = line.split("EC=")
		if len(line_split) >1:
 			ec_number = line_split[1].strip()[:-1]
                        enz_annot[ac].append(ec_number)

		


#print enz_annot

        if line [0:2] == "KW":
                line_split = line[5:-1].split(";")
		for k in line_split:
			k_clean = k.strip().strip(".")
               		if k_clean in KW_terms:
				kw_annot[ac].append(k_clean)
                #if line_split[3] == "phenotype." or line_split[3] == 'gene+phenotype.' :
                #       omim_annot[ac].append(line_split[2][:-1])

	if line [5:7] =="GO" and line[0:2] == "DR":
		line_split = line.split(";")
		#print line_split
		expression = line_split[2].split(":")
		
		if  expression[1] in GO_terms:
			go_annot_dict[ac].append(expression[1]) 		

		

print  go_annot_dict

no_category = open("no_category.txt","w")
enzyme = open("KW_enzymes.txt","w")
tf = open("KW_tf_and_regulators.txt","w") #one can create a separat category for the regulators
channel = open("KW_channels.txt","w")
transporter = open("KW_transporter.txt","w")
receptors = open("KW_receptors.txt","w")
adaptors = open("KW_adatptors.txt","w")
chaperones = open("KW_chaperones.txt","w")
nucbinding = open("KW_nucbinding.txt","w")
dna_rna_bind = open("KW_dna_rna_bind.txt","w")
gpcrs = open("KW_gpcrs.txt","w")
repressors = open("KW_repressors.txt","w")
lipid_bind = open("KW_lipid_bind.txt","w")
heme_related = open ("KW_heme_related.txt","w")
fatty_acid_lip_met = open ("KW_fatty_acid_lip_meta.txt","w")
structural_proteins = open("KW_structural_proteins.txt","w")


GO_prot_enz_bind = open ("GO_prot_enz_bind.txt","w")
GO_enzyme_regulators = open("GO_enzyme_regulators.txt","w")
GO_signaling_path = open("GO_signaling.txt","w")
GO_lipid_bind_transport = open("GO_lipid_bind_transport.txt","w")
GO_dna_rna_bind = open("GO_dna_rna_bind.txt","w")
GO_carbohydratebind = open("GO_carbohydrate_bind.txt","w")
GO_heme_bio_process = open("GO_heme_bio_process.txt","w")
GO_small_molecule_bio_pro = open("GO_small_molecule_bio_process.txt","w")
GO_cytochrome_c_oxidase_act = open ("GO_cyto_c_ox_act.txt","w")
GO_structural_molecule = open("GO_struct_mol.txt","w")
GO_regulation_cat_act = open("GO_reg_cat_act.txt","w")
GO_protein_targ = open("GO_prot_targ.txt","w")
GO_sphingo_nacyl_act = open("GO_sphingo_acyl_activity.txt","w")
GO_cell_adhesion = open("GO_cell_adhesion.txt","w")



# "DNA-binding","Transcription","Receptor","Porin","Ion transport","Transport","Receptor","SH2 domain","SH3 domain"

input_dic = {}
for line in input_2.readlines():
	line_split = line.split("\t")
	up = line_split[0]
	name = line_split[1].split(" ")[0]
	if up not in input_dic.keys():
		input_dic[up]=name

for entry in input_dic.keys():
	enz = 0
	other = 0
	if entry in enz_annot.keys():
		if len(enz_annot[entry]) > 1:
			enzyme.write(entry+"\t"+input_dic[entry]+"\n")
			enz = 1 	
	
	if entry in kw_annot.keys():
		terms = kw_annot[entry][1:]
		if "DNA-binding" in terms and "Transcription" in terms:
			if "Receptor" not in terms:
				tf.write(entry+"\t"+input_dic[entry]+"\n")
				other =1 
		if "Receptor" in terms:
			receptors.write(entry+"\t"+input_dic[entry]+"\n")		
			other =1

		if "G-protein coupled receptor" in terms:
			gpcrs.write(entry+"\t"+input_dic[entry]+"\n")
 			other = 1 
		
		if "DNA-binding" not in terms and "Transcription" in terms:
			if "Receptor" not in terms and "Transcription regulation" in terms :
				tf.write(entry+"\t"+input_dic[entry]+"\n")
                                other =1


		if "Porin" in terms or "Ion transport" in terms or "Ion channel" in terms:
			channel.write(entry+"\t"+input_dic[entry]+"\n") 	 
			other =1 

		if "Transport" in terms:
			transporter.write(entry+"\t"+input_dic[entry]+"\n")
			other =1 			
			
		if "SH2 domain" in terms or "SH3 domain" in terms:
			if enz ==0:
				adaptors.write(entry+"\t"+input_dic[entry]+"\n") 
				other = 1 
		if "Chaperone" in terms :
			chaperones.write(entry+"\t"+input_dic[entry]+"\n")
			other = 1	
		
		if "Nucleotide-binding" in terms :
			nucbinding.write(entry+"\t"+input_dic[entry]+"\n")
			other = 1
		
		if "RNA-binding" in terms or "DNA-binding" in terms:
			dna_rna_bind.write(entry+"\t"+input_dic[entry]+"\n") 
		 	other = 1

		if "Repressor" in terms:
			repressors.write(entry+"\t"+input_dic[entry]+"\n")
			other =1 
		
		if "Enzyme regulator activity" in terms:
			enzyme_regulators.write(entry+"\t"+input_dic[entry]+"\n")
			other = 1 
		
		if "Lipid-binding" in terms:  # keep going form here 

			lipid_bind.write(entry+"\t"+input_dic[entry]+"\n")
			other =1 
		if "Heme" in terms: 
			heme_related.write(entry+"\t"+input_dic[entry]+"\n")
			other =1 

		if "Fatty acid metabolism" in terms or "Lipid metabolism" in terms:
			fatty_acid_lip_met.write(entry+"\t"+input_dic[entry]+"\n")
			other = 1 

		if "Cytoskeleton" in terms or "Microtubule" in terms:
			structural_proteins.write(entry+"\t"+input_dic[entry]+"\n")
			other=1 

	if entry in go_annot_dict.keys():
                terms = go_annot_dict[entry][1:]		
		
		if "enzyme regulator activity" in terms:
			GO_enzyme_regulators.write(entry+"\t"+input_dic[entry]+"\n")
			other=1 
		if "signaling pathway"  in terms or "signal transduction" in terms or "signal transducer activity" in terms:
			GO_signaling_path.write(entry+"\t"+input_dic[entry]+"\n")
			other=1

		if "protein binding" in terms or "enzyme binding" in terms:
                        GO_prot_enz_bind.write(entry+"\t"+input_dic[entry]+"\n")
                        other = 1		

		if "lipid binding" in terms or "lipid transport" in terms:
			GO_lipid_bind_transport.write(entry+"\t"+input_dic[entry]+"\n")
			other=1 
		if "RNA binding" in terms or "DNA binding" in terms:
			GO_dna_rna_bind.write(entry+"\t"+input_dic[entry]+"\n")
			other =1 
		if "carbohydrate binding" in terms:
			GO_carbohydratebind.write(entry+"\t"+input_dic[entry]+"\n")
			other =1 

		if "heme biosynthetic process" in terms:
			GO_heme_bio_process.write(entry+"\t"+input_dic[entry]+"\n")
			other =1 
		if "small molecule biosynthetic process" in terms: 
			GO_small_molecule_bio_pro.write(entry+"\t"+input_dic[entry]+"\n")
			other =1 
		if "cytochrome c oxidase activity" in terms:
			GO_cytochrome_c_oxidase_act.write(entry+"\t"+input_dic[entry]+"\n")
			other =1 
		if "structural molecule activity" in terms:
			GO_structural_molecule.write(entry+"\t"+input_dic[entry]+"\n")
			other=1 

		if "regulation of catalytic activity" in terms:
			GO_regulation_cat_act.write(entry+"\t"+input_dic[entry]+"\n")
			other=1 
		if "protein targeting" in terms or "protein-plasma membrane targeting" in terms :
			GO_protein_targ.write(entry+"\t"+input_dic[entry]+"\n")
			other=1
		if "sphingosine N-acyltransferase activity" in terms:
			GO_sphingo_nacyl_act.write(entry+"\t"+input_dic[entry]+"\n")
			other =1 
		
		if "cell adhesion" in terms:
			GO_cell_adhesion.write(entry+"\t"+input_dic[entry]+"\n")
			other =1 

			

	if enz == 0 and other ==0:
		no_category.write(entry+"\t"+input_dic[entry]+"\n")


"""
        if line[0:2] == "SQ":
                seq_read =1
                seq=""

        if seq_read == 1 and line[0:2] not in ["//","SQ"]:
                seq += line.strip().replace(" ", "")

        if line[0:2]=="//":
                seq_read = 0
                sequence_dict[ac].append(seq)


"""	 
