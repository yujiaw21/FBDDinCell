data <- read.csv("./crossed_table_HEKS_probemethyl_comp200_r4.csv")

ratio_cutoff <- 5

data[data$coumarin_200uM_293T < ratio_cutoff, grep("coumarin", colnames(data), value = TRUE)[-1]] <-
  "-"
data[data$pipphen_200uM_293T < ratio_cutoff, grep("pipphen", colnames(data), value = TRUE)[-1]] <-
  "-"
data[data$hydrooxoquin_200uM_293T < ratio_cutoff, grep("hydrooxoquin", colnames(data), value = TRUE)[-1]] <-
  "-"

data <-
  data[!(
    data$coumarin_200uM_293T < ratio_cutoff &
      data$pipphen_200uM_293T < ratio_cutoff &
      data$hydrooxoquin_200uM_293T < ratio_cutoff
  ), ]

write.csv(data, file = "crossed_table_HEKS_probemethyl_comp200_r4.filtered.csv")
