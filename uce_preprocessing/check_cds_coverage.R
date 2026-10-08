library(tidyverse)

cov10 <- read.csv("uce_preprocessing/windows_50_dp10_samples.csv") #windows with dp >=10


pops <- cov10$samples |>
  str_extract_all("[A-Z]+\\d{2}") |>
  unlist() |>
  unique()

grep("NTKG", cov10$samples)

#pop groups
outgroup <- "NTKG00"
LH_s1 <- c("LH00", "TWC99")
LH_s1 <- c("LH23", "TWC22")
LSC_changed_s1 <- c("LA98", "LJC00", "LJM98", "SB98", "SJ98", "PA99")
LSC_changed_s2 <- c("LA22", "LJC22", "LJM22", "SB22", "SJ22", "PA24")
LSC_unchanged_s1 <- c(
  "KC98",
  "MP98",
  "SAU98",
  "SAL98",
  "SM98",
  "LP04",
  "MT98",
  "WL03",
  "MB98",
  "UK98",
  "SLO98"
)
lsc_unchanged_s2 <- c(
  "KC23",
  "MP23",
  "SAU22",
  "SAL22",
  "SM22",
  "LP23",
  "MT23",
  "WL22",
  "MB23",
  "UK22",
  "SLO22"
)
sk_s1 <- "LS97"
sk_s2 <- "LS23"


#add a list column with pop IDs
cov10_popcount <- cov10 |>
  mutate(pops = str_extract_all(samples, "[A-Z]+\\d{2}")) |>
  mutate(n_pops = unlist(lapply(pops, n_distinct)))
cov10_start <- cov10_popcount |>
  filter(start == 1)


ggplot(cov10_start, aes(x = contig, y = n_samples)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x = "Contig", y = "Number of Samples with Coverage >= 10")
