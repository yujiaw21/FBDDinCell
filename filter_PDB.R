require(bio3d)

sol <-
  read.csv("./SOL_for_fpocket.csv",
           stringsAsFactors = FALSE,
           header = TRUE)

sol$PDB_clean <- ""
i <- 1
while (i <= dim(sol)[1]) {
  #those have multiple PDB
  if (nchar(sol$PDB[i]) > 4) {
    #get isoTOP peptides
    pep_seq <- sol[sol$Uniprot == sol$Uniprot[i], ]$PEPTIDE.SEQUENCE
    #get PDBs associated with this Uniprot ID
    pdbs_to_check <- strsplit(sol$PDB[i], ", ")[[1]]
    number_of_peptide_in_PDB <- numeric()
    
    for (j in 1:length(pdbs_to_check)) {
      #find the PDB that contains most of isoTOP peptide
      number_of_peptide_in_PDB[j] <- tryCatch({
        #read PDB file, get its chains and sequences
        pdb <-
          read.pdb(file.path(getwd(), "pdb", paste(
            tolower(pdbs_to_check[j]), ".pdb", sep = ""
          )))
        pdb_info <- list(chain = character(),
                         seq = character())
        pdb_info$chain <-
          unique(pdb$atom$chain[pdb$atom$type == "ATOM"])
        for (k in 1:length(pdb_info$chain)) {
          pdb_info$seq[k] <-
            paste(pdbseq(pdb, inds = atom.select(pdb, "calpha", chain = pdb_info$chain[k])), collapse = "")
        }
        
        #determine whether or not each isoTOP peptide is contained in this PDB
        in_pdb <- logical()
        for (n in 1:length(pep_seq)) {
          in_pdb[n] <-
            ifelse(length(grep(pep_seq[n], pdb_info$seq, fixed = TRUE)) > 0, TRUE, FALSE)
        }
        sum(in_pdb)
      }, error = function(e) {
        0
      })
    }
    
    sol[sol$Uniprot == sol$Uniprot[i], ]$PDB_clean <-
      paste(pdbs_to_check[number_of_peptide_in_PDB == max(number_of_peptide_in_PDB)], collapse = " ")
    
    i <- i + length(pep_seq)
  }
  else{
    sol$PDB_clean[i] <- sol$PDB[i]
    i <- i + 1
  }
}

write.csv(sol, file = "SOL_w_PDB_filter.csv", row.names = FALSE)
