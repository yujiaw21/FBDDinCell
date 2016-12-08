data <-
  read.csv("./20uMtargets_found_200uMtargets.csv", header = TRUE)

ratio_cutoff <- 5

data[data$coumarin_20uM_293T < ratio_cutoff, grep("coumarin", colnames(data), value = TRUE)[-1]] <-
  "-"
data[data$pipphen_20uM_293T < ratio_cutoff, grep("pipphen", colnames(data), value = TRUE)[-1]] <-
  "-"
data[data$hydrooxoquin_20uM_293T < ratio_cutoff, grep("hydrooxoquin", colnames(data), value = TRUE)[-1]] <-
  "-"

data <-
  data[!(
    data$coumarin_20uM_293T < ratio_cutoff &
      data$pipphen_20uM_293T < ratio_cutoff &
      data$hydrooxoquin_20uM_293T < ratio_cutoff &
      data$X20phenpip_w_B00174059 == 0 & data$X20phenpip_w_B00174040 == 0
  ), ]

write.csv(data, file = "20uMtargets_found_200uMtargets_filtered.csv", row.names = FALSE)
