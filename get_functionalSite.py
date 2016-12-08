import pandas as pd
from Bio import ExPASy
from Bio import SwissProt
import os

sol = pd.read_csv("SOL_original.csv")

feature_of_interest = ["ACT_SITE", "METAL", "BINDING", "SITE", "NP_BIND"]


uniprot_id = sol.Uniprot
res = dict()

with open("uniprot_sprot_human.dat") as handle:
    for record in SwissProt.parse(handle):
        for id in sol.Uniprot:
            if id in record.accessions:
                tmp = []
                for feature in record.features:
                    if feature[0] in feature_of_interest:
                        tmp.append(
                            feature[0] + " " + str(feature[1]) + " " + str(feature[2]))

                res[id] = [record.sequence, tmp]


with open("site.txt", "w") as file:

    for key, value in res.items():

        file.write("{}\t{}\t{}\n".format(key, value[0], ", ".join(value[1])))

