require(bio3d)
sol$spatial_distance <- ""
sol$shift <- ""
sol$site_total <- ""
sol$site_in_pdb <- ""
sol$site_near_isotop <- ""

for (i in 1:dim(sol)[1]) {
  max_shift <- 0
  try(if (sol$active_site[i] != "") {
    pdb <-
      read.pdb(file.path(dirname(getwd()), "pdb", paste(tolower(sol$PDB[i]), ".pdb", sep = "")))
    isotop <- sol$PEPTIDE.SEQUENCE[i]
    isotop_start <-
      as.numeric(strsplit(sol$TARG_PEPRANGE[i], "_")[[1]][2])
    ref_seq <- strsplit(sol$sequence[i], "")[[1]]
    
    site_list <-
      strsplit(unlist(strsplit(sol$active_site[i], ", ")), " ")
    site <- character()
    for (j in 1:length(site_list)) {
      site <-
        c(site, as.numeric(site_list[[j]][2]):as.numeric(site_list[[j]][3]))
    }
    site <- as.numeric(site)
    
    pdb_chain <- unique(pdb$atom$chain[pdb$atom$type == "ATOM"])
    min_distance <- character()
    site_distance <- vector(mode = "list", length = length(site))
    names(site_distance) <- site
    
    
    for (j in 1:length(pdb_chain)) {
      chain_seq <-
        pdbseq(pdb, inds = atom.select(pdb, "calpha", chain = pdb_chain[j]))
      pos <- regexpr(isotop, paste(chain_seq, collapse = ""))
      if (pos != -1)
      {
        seq_pos <- chain_seq[pos:(pos + nchar(isotop) - 1)]
        isotop_atoms <-
          atom.select(pdb,
                      "protein",
                      chain = pdb_chain[j],
                      resno = as.numeric(names(seq_pos)))
        shift <- as.numeric(names(seq_pos)[1]) - isotop_start
        
        max_shift <-
          ifelse(abs(shift) >= abs(max_shift),
                 abs(shift),
                 abs(max_shift))
        
        mis_alignment <-
          tryCatch(
            sum(!chain_seq == ref_seq[as.numeric(names(chain_seq)) - shift]),
            error = function(e) {
              9999
            }
          )
        #show the mismathch residue
        #which(!chain_seq==ref_seq[as.numeric(names(chain_seq))-shift])
        
        if (ifelse(mis_alignment >= 4, TRUE, FALSE)) {
          min_distance[j] <-
            paste("Mismatch", sum(!chain_seq == ref_seq[as.numeric(names(chain_seq)) -
                                                          shift]), sep = ": ")
        }
        else{
          site_atoms <-
            atom.select(pdb,
                        "protein",
                        chain = pdb_chain[j],
                        resno = site + shift)
          if (length(site_atoms$atom) == 0) {
            min_distance[j] <- "Active Site Not Found"
          }
          else{
            min_distance[j] <-
              tryCatch(
                min(dist.xyz(
                  pdb$xyz[isotop_atoms$xyz], pdb$xyz[site_atoms$xyz], all.pairs = TRUE
                )),
                error = function(e) {
                  "ERROR"
                }
              )
            for (k in 1:length(site)) {
              site_sle <-
                atom.select(pdb,
                            "protein",
                            chain = pdb_chain[j],
                            resno = site[k] + shift)
              site_distance[[k]][j] <-
                tryCatch(
                  min(dist.xyz(
                    pdb$xyz[isotop_atoms$xyz], pdb$xyz[site_sle$xyz], all.pairs = TRUE
                  )),
                  error = function(e) {
                    9999
                  }
                )
            }
          }
        }
      }
      else{
        min_distance[j] <- "isoTOP Peptide Not Found"
      }
    }
    
    site_distance <- lapply(site_distance, min)
    
    if (any(min_distance == "ERROR")) {
      sol$spatial_distance[i] <- "Something's Wrong!"
    } else if (any(!is.na(as.numeric(min_distance)))) {
      sol$spatial_distance[i] <-
        round(min(as.numeric(min_distance), na.rm = TRUE), digits = 3)
      sol$shift[i] <- max_shift
      sol$site_total[i] <- length(site)
      sol$site_in_pdb[i] <- sum(site_distance != 9999)
      sol$site_near_isotop[i] <- sum(site_distance <= 6)
    } else{
      sol$spatial_distance[i] <- paste(min_distance, collapse = ";")
      sol$shift[i] <- max_shift
    }
  })
}

#write.csv(sol,file = "SOL_spatial_distance.csv",row.names = FALSE)