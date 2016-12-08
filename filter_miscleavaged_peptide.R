require(UniProt.ws)
up <- UniProt.ws(taxId = 9606)

sol_original$valid <- "-"
res <-
  select(
    up,
    keys = sol_original$Uniprot,
    columns = c("LENGTH", "SEQUENCE"),
    keytype = "UNIPROTKB"
  )

for (i in 1:dim(sol_original)[1]) {
  seq <- strsplit(sol_original$PEPTIDE.SEQUENCE[i], "")[[1]]
  
  #End with K or R but not C-term peptide
  if (seq[length(seq)] != "K" & seq[length(seq)] != "R") {
    try(if (res[i, 2] != as.numeric(strsplit(sol_original$TARG_PEPRANGE[i], "_")[[1]][3]) -
            1) {
      sol_original$valid[i] <- FALSE
    })
  }
  #contain miscleavaged K or R
  if (("R" %in% seq[1:(length(seq) - 1)]) |
      ("K" %in% seq[1:(length(seq) - 1)])) {
    sol_original$valid[i] <- FALSE
  }
}

for (i in 1:dim(sol_original)[1]) {
  if (sol_original$valid[i] == "-") {
    if (strsplit(sol_original$TARG_PEPRANGE[i], "_")[[1]][2] > 1) {
      pos <-
        as.numeric(strsplit(sol_original$TARG_PEPRANGE[i], "_")[[1]][2]) - 1
      try(if (strsplit(res[i, 3], "")[[1]][pos] != "K" &
              strsplit(res[i, 3], "")[[1]][pos] != "R") {
        sol_original$valid[i] <- FALSE
      })
    }
  }
}
