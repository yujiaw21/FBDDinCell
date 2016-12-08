data <-
  read.delim(
    "./crossed_table_HEKS_probemethyl_probeprobe_r4.txt",
    header = FALSE,
    stringsAsFactors = FALSE
  )
colnames(data) <- data[1, ]
data <- data[-1, ]

probe.pair <- strsplit(gsub("200uM", "", colnames(data)[14:22]), "_")

#filter probe targets
for (i in 1:length(probe.pair))
{
  probe.1 <- grep(probe.pair[[i]][1], colnames(data)[3:13], value = TRUE)
  probe.2 <-
    grep(probe.pair[[i]][2], colnames(data)[3:13], value = TRUE)
  data[data[, probe.1] == 0 & data[, probe.2] == 0, i + 13] <- "-"
}

#check redundancy
for (i in 1:dim(data)[1]) {
  if (sum(data[i, 14:22] == "-") == 9) {
    print(data$Identifier[i])
  }
}

write.csv(data, file = "processed.probe.vs.probe.csv", row.names = FALSE)