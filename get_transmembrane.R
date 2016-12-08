require(UniProt.ws)
sol_byProtein$transmembrane <- ""

up <- UniProt.ws(taxId = 9606)
uniprot.id <- sol_byProtein$V1
res <-
  select(up,
         keys = uniprot.id,
         columns = c("FEATURES"),
         keytype = "UNIPROTKB")


for (i in 1:dim(sol_byProtein)[1]) {
  if (length(grep("transmembrane", tolower(res$FEATURES[i]))) > 0) {
    sol_byProtein$transmembrane[i] <- TRUE
  } else{
    sol_byProtein$transmembrane[i] <- FALSE
  }
}
