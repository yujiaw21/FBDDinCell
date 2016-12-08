target.list <-
  read.csv(
    "./processed.200uM.competition.csv",
    header = TRUE,
    stringsAsFactors = FALSE
  )

target.list$competited <- FALSE

start <- grep("_w_", colnames(target.list))[1]
end <-
  grep("_w_", colnames(target.list))[length(grep("_w_", colnames(target.list)))]
cutoff_ratio <- 3

for (i in 1:dim(target.list)[1]) {
  if (sum(as.numeric(target.list[i, start:end]) > cutoff_ratio, na.rm = TRUE) >=
      1) {
    target.list$competited[i] <- TRUE
  }
}

write.csv(target.list, file = "processed.200uM.competition.filtered.csv", row.names = FALSE)
