#! /usr/bin/python

import sys 

input_comp = open(sys.argv[1])

input_table = open(sys.argv[2])


competition_dic = {}
 
compound_dic = {}
for line in input_comp.readlines():
	line_split = line.split("\t")
	
	if line_split[0] == "Identifier":
		for k in range(2,len(line_split)-1):
			name = line_split[k].strip()
			compound_dic[k] = name	
	else:
		up = line_split[0]
		for k in range(2,len(line_split)-1):	
			value = line_split[k].strip()
			print value
			if value != "-":
				if float(value) >= 2.5:
					if up in competition_dic.keys():
						competition_dic[up].append(compound_dic[k]+"_"+str(value))
					else:
						competition_dic[up]=[compound_dic[k]+"_"+str(value)] 
						
out = open("table_competitior.txt","w")


for line in input_table.readlines():
        line_split = line.split("\t")
        
        if line_split[0] == "Identifier":

		out.write(line[:-1]+"\t"+"Competitor"+"\n")
        else:
		up_line = line_split[0]
		
		if up_line in competition_dic.keys():
			out.write(line[:-1]+"\t"+" ".join(competition_dic[up_line])+"\n")
		else:
			out.write(line[:-1]+"\t-\n")
		
