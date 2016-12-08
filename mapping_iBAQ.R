require(parallel)

target.list <-
  read.csv("target.list.csv",
           header = TRUE,
           stringsAsFactors = FALSE)

iBAQ <-
  read.csv("./IBAQ.csv", header = TRUE, stringsAsFactors = FALSE)

#clean the iBAQ data
iBAQ[iBAQ$median == "#NUM!", "median"] <- "-"

#find position of uniprot ID from target list in iBAQ data
iBAQ.Uniprot <- strsplit(iBAQ$Uniprot, ";")
pos <-
  mclapply(
    target.list$Identifier,
    grep,
    iBAQ.Uniprot,
    fixed  = TRUE,
    mc.cores = detectCores()
  )

#assign iBAQ values to target list
for (i in 1:length(pos)) {
  #no match
  if (length(pos[[i]]) == 0) {
    target.list$iBAQ[i] <- "-"
  }
  #one match
  else if (length(pos[[i]]) == 1) {
    target.list$iBAQ[i] <- iBAQ$median[pos[[i]]]
  }
  #more than one match, use the median value
  else if (length(pos[[i]]) > 1) {
    len <- length(pos[[i]])
    value <- character()
    for (j in 1:len) {
      value <- c(value, iBAQ$median[pos[[i]][j]])
    }
    target.list$iBAQ[i] <- median(as.numeric(value), na.rm = TRUE)
  }
}


write.csv(target.list, file = "target.list.w.iBAQ.csv", row.names = FALSE)