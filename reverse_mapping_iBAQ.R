require(parallel)


target.list.20uM <-
  read.csv(
    "./20uMtargets_found_200uMtargets_w_full_iBAQ_DB.csv",
    header = TRUE,
    stringsAsFactors = FALSE
  )
target.list <-
  read.csv("./target_w_full_iBAQ_DB.csv",
           header = TRUE,
           stringsAsFactors = FALSE)
competition.20uM <-
  read.csv(
    "./20uM_competition_w_full_iBAQ_DB.csv",
    header = TRUE,
    stringsAsFactors = FALSE
  )
competition.200uM <-
  read.csv(
    "./200uM_competition_w_full_iBAQ_DB.csv",
    header = TRUE,
    stringsAsFactors = FALSE
  )
probe.vs.probe <-
  read.csv(
    "./probe_v_probe_w_full_iBAQ_DB.csv",
    header = TRUE,
    stringsAsFactors = FALSE
  )
iBAQ <-
  read.csv("./iBAQ.csv", header = TRUE, stringsAsFactors = FALSE)

#find position of uniprot ID from target list in iBAQ data
iBAQ.Uniprot <- strsplit(iBAQ$Uniprot, ";")
pos.target.list.20uM <-
  mclapply(
    target.list.20uM$Identifier,
    grep,
    iBAQ.Uniprot,
    fixed  = TRUE,
    mc.cores = detectCores()
  )
pos.target.list <-
  mclapply(
    target.list$Identifier,
    grep,
    iBAQ.Uniprot,
    fixed  = TRUE,
    mc.cores = detectCores()
  )
pos.20uM <-
  mclapply(
    competition.20uM$Identifier,
    grep,
    iBAQ.Uniprot,
    fixed  = TRUE,
    mc.cores = detectCores()
  )
pos.200uM <-
  mclapply(
    competition.200uM$Identifier,
    grep,
    iBAQ.Uniprot,
    fixed  = TRUE,
    mc.cores = detectCores()
  )
pos.probe.vs.probe <-
  mclapply(
    probe.vs.probe$Identifier,
    grep,
    iBAQ.Uniprot,
    fixed  = TRUE,
    mc.cores = detectCores()
  )

iBAQ$in_20uM_targets_list <- FALSE
iBAQ$in_targets_list <- FALSE
iBAQ$in_20uM_competition <- FALSE
iBAQ$in_200uM_competition <- FALSE
iBAQ$in_probe_vs_probe <- FALSE

for (i in 1:length(pos.target.list.20uM)) {
  #at least one match
  if (length(pos.target.list.20uM[[i]]) > 0) {
    iBAQ$in_20uM_targets_list[pos.target.list.20uM[[i]]] <- TRUE
  }
}

for (i in 1:length(pos.target.list)) {
  #at least one match
  if (length(pos.target.list[[i]]) > 0) {
    iBAQ$in_targets_list[pos.target.list[[i]]] <- TRUE
  }
}

for (i in 1:length(pos.20uM)) {
  #at least one match
  if (length(pos.20uM[[i]]) > 0) {
    iBAQ$in_20uM_competition[pos.20uM[[i]]] <- TRUE
  }
}

for (i in 1:length(pos.200uM)) {
  #at least one match
  if (length(pos.200uM[[i]]) > 0) {
    iBAQ$in_200uM_competition[pos.200uM[[i]]] <- TRUE
  }
}

for (i in 1:length(pos.probe.vs.probe)) {
  #at least one match
  if (length(pos.probe.vs.probe[[i]]) > 0) {
    iBAQ$in_probe_vs_probe[pos.probe.vs.probe[[i]]] <- TRUE
  }
}

write.csv(iBAQ, file = "reverse_mapped_iBAQ.csv", row.names = FALSE)