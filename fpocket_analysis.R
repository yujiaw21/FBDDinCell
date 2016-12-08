require(bio3d)

fpocket$resid_by_fpocket <- ""
fpocket$shifted_isotop <- ""
fpocket$overlapped <- ""
fpocket$isotop.covered <- ""
fpocket$no.filtered.pockets <- ""
fpocket$number.of.pockets <- "-"
fpocket$pocket.overlapped <- "-"

vol_cutoff <- 500

#get the residues around the binding pockets
get_resid <- function(path,
                      pdb_id,
                      pep_seq,
                      isotop_start,
                      isotop_end) {
  #read the pdb file and obtain information of every chain and sequence
  pdb <-
    read.pdb(file.path(dirname(path), tolower(paste(
      pdb_id, ".pdb", sep = ""
    ))), maxlines = -1)
  pdb_info <- list(
    chain = character(),
    seq = character(),
    shift = numeric(),
    target = c()
  )
  pdb_info$chain <- unique(pdb$atom$chain[pdb$atom$type == "ATOM"])
  for (i in 1:length(pdb_info$chain)) {
    pdb_info$seq[i] <-
      paste(pdbseq(pdb, inds = atom.select(pdb, "calpha", chain = pdb_info$chain[i])), collapse = "")
    pdb_info$target[i] <-
      ifelse(length(grep(pep_seq, pdb_info$seq[i], fixed = TRUE)) != 0, TRUE, FALSE)
    if (pdb_info$target[i] == TRUE) {
      pdb_info$shift[i] <-
        as.numeric(names(pdbseq(
          pdb, inds = atom.select(pdb, "calpha", chain = pdb_info$chain[i])
        )[regexpr(pep_seq, pdb_info$seq[i])])) - isotop_start
    } else{
      pdb_info$shift[i] <- NA
    }
  }
  
  if (any(pdb_info$target) == FALSE) {
    return(paste("Isotop Peptide Not Found in ", pdb_id, sep = ""))
  }
  
  #only consider those chains that contain isotop peptide sequence
  chain_of_target <- pdb_info$chain[pdb_info$target]
  number_of_pockets <- length(list.files(file.path(path, "pockets"))) /
    2
  
  if (number_of_pockets > 0) {
    resid_by_fpocket <- replicate(length(chain_of_target), character())
    names(resid_by_fpocket) <- chain_of_target
    
    #filter the predicted pockets based on volume
    pocket.info <-
      read.delim(
        file.path(path, tolower(paste(
          pdb_id, "info.txt", sep = "_"
        ))),
        header = FALSE,
        stringsAsFactors = FALSE,
        strip.white = TRUE
      )
    pocket_vol <- numeric()
    for (i in 1:number_of_pockets) {
      pocket_vol <- c(pocket_vol, pocket.info$V5[20 * i - 12])
    }
    pockets_ind <- (1:number_of_pockets)[pocket_vol >= vol_cutoff]
    
    if (length(pockets_ind) > 0) {
      for (i in pockets_ind) {
        pdb <-
          read.pdb(file.path(
            path,
            "pockets",
            paste("pocket", i - 1, "_atm.pdb", sep = "")
          ))
        chain_of_pocket <-
          unique(pdb$atom$chain[pdb$atom$type == "ATOM"])
        if ("FALSE" %in% as.character(chain_of_pocket)) {
          chain_of_pocket[as.character(chain_of_pocket) == "FALSE"] <- "F"
          chain_to_go <- intersect(chain_of_target, chain_of_pocket)
          if (length(chain_to_go) != 0) {
            for (j in 1:length(chain_to_go)) {
              if (chain_to_go[j] == "F") {
                res <- pdbseq(pdb, inds = atom.select(pdb, chain = "FALSE", type = "ATOM"))
                res <- res[!duplicated(res)]
                resid_by_fpocket[[chain_to_go[j]]] <-
                  union(resid_by_fpocket[[chain_to_go[j]]],
                        paste(names(res), chain_to_go[j], sep = "."))
              } else{
                res <-
                  pdbseq(pdb,
                         inds = atom.select(pdb, chain = chain_to_go[j], type = "ATOM"))
                res <- res[!duplicated(res)]
                resid_by_fpocket[[chain_to_go[j]]] <-
                  union(resid_by_fpocket[[chain_to_go[j]]],
                        paste(names(res), chain_to_go[j], sep = "."))
              }
            }
          }
        } else{
          chain_to_go <-
            as.character(intersect(chain_of_target, chain_of_pocket))
          if (length(chain_to_go) != 0) {
            for (j in 1:length(chain_to_go)) {
              res <-
                pdbseq(pdb,
                       inds = atom.select(pdb, chain = chain_to_go[j], type = "ATOM"))
              res <- res[!duplicated(res)]
              resid_by_fpocket[[chain_to_go[j]]] <-
                union(resid_by_fpocket[[chain_to_go[j]]],
                      paste(names(res), chain_to_go[j], sep = "."))
            }
          }
        }
      }
    } else{
      return("No Pockets Pass Filter")
    }
  } else{
    return("No Predicted Pockets")
  }
  
  overlapped <- replicate(length(resid_by_fpocket), character())
  names(overlapped) <- names(resid_by_fpocket)
  shifted_isotop <- replicate(length(resid_by_fpocket), character())
  names(shifted_isotop) <- names(resid_by_fpocket)
  
  for (i in 1:length(resid_by_fpocket)) {
    shifted_isotop[[names(resid_by_fpocket[i])]] <-
      paste((isotop_start:(isotop_end)) + pdb_info$shift[pdb_info$chain == names(resid_by_fpocket[i])],
            names(resid_by_fpocket[i]),
            sep = ".")
    overlapped[[names(resid_by_fpocket[i])]] <-
      intersect(shifted_isotop[[names(resid_by_fpocket[i])]], resid_by_fpocket[[i]])
  }
  
  resid_by_fpocket <-
    lapply(Filter(length, resid_by_fpocket), paste, collapse = ",")
  overlapped <- lapply(Filter(length, overlapped), paste, collapse = ",")
  shifted_isotop <-
    lapply(Filter(length, shifted_isotop), paste, collapse = ",")
  
  if (length(overlapped) != 0) {
    return(c(
      paste(paste(resid_by_fpocket, collapse = ","), sep = " "),
      paste(paste(overlapped, collapse = ","), sep = " "),
      paste(paste(shifted_isotop, collapse = ","), sep = " "),
      length(pockets_ind)
    ))
  } else{
    return(c(
      paste(paste(resid_by_fpocket, collapse = ","), sep = " "),
      "No Overlap",
      paste(paste(shifted_isotop, collapse = ","), sep = " "),
      length(pockets_ind)
    ))
  }
}

get_pocket <- function(path, resid.overlapped) {
  #determine the number of pockets of the pdb file
  number_of_pockets <- length(list.files(file.path(path, "pockets"))) /
    2
  
  pocket.info <-
    read.delim(
      file.path(path, tolower(paste(
        pdb_id, "info.txt", sep = "_"
      ))),
      header = FALSE,
      stringsAsFactors = FALSE,
      strip.white = TRUE
    )
  pocket_vol <- numeric()
  for (i in 1:number_of_pockets) {
    pocket_vol <- c(pocket_vol, pocket.info$V5[20 * i - 12])
  }
  pockets_ind <- (1:number_of_pockets)[pocket_vol >= vol_cutoff]
  
  resid.overlapped <- strsplit(resid.overlapped, ",")[[1]]
  overlap_w_pocket <- c()
  for (i in pockets_ind) {
    pdb <-
      read.pdb(file.path(path, "pockets", paste("pocket", i - 1, "_atm.pdb", sep = "")))
    chain_of_pocket <- unique(pdb$atom$chain[pdb$atom$type == "ATOM"])
    resid_by_pocket <- c()
    for (j in 1:length(chain_of_pocket)) {
      res <-
        pdbseq(pdb, inds = atom.select(
          pdb,
          chain = as.character(chain_of_pocket[j]),
          type = "ATOM"
        ))
      res <- res[!duplicated(res)]
      if (as.character(chain_of_pocket[j]) == "FALSE") {
        resid_by_pocket <- c(resid_by_pocket, paste(names(res), "F", sep = "."))
      } else{
        resid_by_pocket <-
          c(resid_by_pocket,
            paste(names(res), chain_of_pocket[j], sep = "."))
      }
    }
    if (length(intersect(resid.overlapped, resid_by_pocket)) > 0) {
      overlap_w_pocket <- c(overlap_w_pocket, i)
    }
  }
  return(overlap_w_pocket)
}

for (i in 1:dim(fpocket)[1]) {
  path <-
    file.path("./Site_of_Labeling/version_5", "pdb", tolower(paste(fpocket$PDB[i], "out", sep = "_")))
  isotop_start <-
    as.numeric(strsplit(fpocket$Labeled.Peptide[i], "-")[[1]][1])
  isotop_end <-
    as.numeric(strsplit(fpocket$Labeled.Peptide[i], "-")[[1]][2])
  pdb_id <- fpocket$PDB[i]
  pep_seq <- fpocket$Peptide.Sequence[i]
  
  result <- get_resid(path, pdb_id, pep_seq, isotop_start, isotop_end)
  fpocket$resid_by_fpocket[i] <- result[1]
  fpocket$overlapped[i] <- result[2]
  fpocket$shifted_isotop[i] <- result[3]
  fpocket$no.filtered.pockets[i] <- result[4]
  fpocket$overlapped[is.na(fpocket$overlapped)] <- "No Overlap"
  
  if (fpocket$overlapped[i] != "No Overlap") {
    overlapped_resid <- strsplit(fpocket$overlapped[i], ",")[[1]]
    
    mean_overlapped <-
      mean(table(gsub("[0-9]+.", "", overlapped_resid)))
    isotop_peptide_length <- nchar(fpocket$Peptide.Sequence[i])
    
    fpocket$isotop.covered[i] <-
      round(mean_overlapped / isotop_peptide_length, digits = 3)
  } else{
    fpocket$isotop.covered[i] <- 0
  }
  
  fpocket$number.of.pockets[i] <-
    length(list.files(file.path(path, "pockets"))) / 2
  
  if (fpocket$overlapped[i] != "No Overlap") {
    fpocket$pocket.overlapped[i] <-
      toString(get_pocket(path, fpocket$overlapped[i]))
  }
}

hist(as.numeric(fpocket$isotop.covered), breaks = 20)

write.table(
  fpocket,
  file = "fpocket.vol_500.txt",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
