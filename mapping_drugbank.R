all <-
  read.csv(
    "./drugbank_all_target_polypeptide_ids.fasta/all.fasta",
    header = TRUE,
    stringsAsFactors = FALSE
  )
nutraceutical <-
  read.csv(
    "./drugbank_nutraceutical_target_polypeptide_ids.fasta/all.fasta",
    header = TRUE,
    stringsAsFactors = FALSE
  )

target.list <-
  read.csv("./target.list.csv",
           header = TRUE,
           stringsAsFactors = FALSE)

#exclude drugbank targets with only endogenous ligands
all$No.Cofactor <- FALSE

in_nutraceutical <- all$UniProt.ID %in% nutraceutical$UniProt.ID

all.drug.id <- strsplit(all$Drug.IDs, ";")

for (i in 1:dim(all)[1]) {
  if (!in_nutraceutical[i]) {
    all$No.Cofactor[i] <- TRUE
  }
  else if (length(all.drug.id[[i]]) > length(unlist(nutraceutical$Drug.ID[which(nutraceutical$UniProt.ID == all$UniProt.ID[i])]))) {
    all$No.Cofactor[i] <- TRUE
  }
}


#exclude targets with only endogenous ligands
target.list$DrugBank <-
  (target.list$Identifier %in% all[all$No.Cofactor, ]$UniProt.ID)

#not exclude above
#target.list$DrugBank <-target.list$Identifier %in% all$UniProt.ID

#write output data
write.csv(target.list, file = "target.list.w.drugbank.csv", row.names = FALSE)
