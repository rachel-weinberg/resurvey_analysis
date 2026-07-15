pacman::p_load(
  tidyverse,
  ggplot2,
  rstatix
)

###ADMIXTURE PLOTS####
#Plot cross-validation results
xval_results <- "path/to/admixture_results/CV_error_K1-5.txt" #Replace with filepath for cross-validation results

xval <- read_table(
  xval_results,
  col_names = FALSE,
  col_types = cols(X1 = col_skip(), X2 = col_skip())
)
xval$X3 <- c(1:5)
colnames(xval) <- c("K", "CV_error")

ggplot(xval, aes(x = K, y = CV_error)) +
  geom_point() +
  geom_line() +
  labs(title = "Cross-validation error for k=1-5")


#Read in results for k=3
Q <- "path/to/admixture_results/all_samples_100p.3.Q" #Replace with filepath for admixture results
admix_results_3 <- read_table(Q, col_names = FALSE)


#Set sample order and subsample only two individuals per population
#subset LSC and LH for s1 and s2
LSC_s1 <- sort(c(
  "UK98",
  "KC98",
  "LP04",
  "MB98",
  "MP07",
  "MT98",
  "LJM98",
  "SB98",
  "SM98",
  "SJ98",
  "PA99",
  "SAL98",
  "SAU98",
  "SLO98",
  "LJC00",
  "LA98",
  "WL03"
))
LSC_s2 <- sort(c(
  "UK22",
  "KC23",
  "LP23",
  "MB23",
  "MP23",
  "MT23",
  "SM22",
  "SAL22",
  "SAU22",
  "SLO22",
  "WL22"
))
LH_s1 <- sort(c("LH00", "TWC99"))
LH_s2 <- sort(c(
  "LH23",
  "TWC22",
  "SB22",
  "SJ22",
  "LA22",
  "LJM22",
  "LJC22",
  "PA24"
))
SK_s1 <- sort(c("LS97"))
SK_s2 <- sort(c("LS23"))
historical <- sort(c("MC76", "SBCo"))

LH_both <- c("LH", "TWC")
SK_both <- c("LS")

s2_order <- c(
  "LH23",
  "TWC22",
  "LS23",
  "KC23",
  "LA22",
  "LJC22",
  "LJM22",
  "LP23",
  "MB23",
  "MP23",
  "MT23",
  "PA24",
  "SAL22",
  "SAU22",
  "SB22",
  "SJ22",
  "SLO22",
  "SM22",
  "UK22",
  "WL22"
)

fam_path <- "path/to/admixture_results/all_samples_100p.fam" #Replace with filepath for plink fam file
sample_names <- read_table(
  fam_path,
  col_names = FALSE
)

# Function to find row numbers for partial matches
get_partial_matches <- function(pattern_vector, search_column) {
  # results vector
  row_numbers <- c()

  for (i in seq_along(pattern_vector)) {
    match_row <- which(grepl(pattern_vector[i], search_column, fixed = TRUE))
    row_numbers <- c(row_numbers, match_row)
  }

  return(row_numbers)
}


s1_sample_order <- get_partial_matches(
  c(LH_s1, SK_s1, LSC_s1),
  sample_names$X1
)


s2_sample_order <- get_partial_matches(
  s2_order,
  sample_names$X1
)

historical_sample_order <- get_partial_matches(
  historical,
  sample_names$X1
)

#sample name vectors before removing unequal pops
s1_samples <- sample_names$X1[s1_sample_order]
s2_samples <- sample_names$X1[s2_sample_order]


#Get sample counts per pop per timepoint
s1_samples_df <- data.frame(
  sample = s1_samples,
  pop = str_extract(s1_samples, "^[A-Z]{2,4}")
)
s2_samples_df <- data.frame(
  sample = s2_samples,
  pop = str_extract(s2_samples, "^[A-Z]{2,4}")
)

s1_samples_count <- s1_samples_df |> count(pop)
s2_samples_count <- s2_samples_df |> count(pop)
s1_s2_samples_df <- full_join(
  s1_samples_count,
  s2_samples_count,
  by = "pop",
  suffix = c("_s1", "_s2")
)
unequal_pops <- s1_s2_samples_df |> filter(n_s1 != n_s2)

#Remove 2 samples from LJM s1
#remove 2 samples from MT s2
#Remove 1 sample from PA s1
#Remove 1 sample from SAL s1
#Remove 1 sample form SM s1
rm_samples_s1 <- which(
  s1_samples %in% c("LJM9806", "LJM9802", "PA9904", "SAL9805", "SM9802")
)
rm_samples_s2 <- which(s2_samples %in% c("MT2303", "MT2304"))

s1_samples_2per <- s1_samples_df |>
  group_by(pop) |>
  filter(row_number() %in% c(1:2))
s2_samples_2per <- s2_samples_df |>
  group_by(pop) |>
  filter(row_number() %in% c(1:2))

s1_sample_order_new <- s1_sample_order[which(
  s1_samples %in% s1_samples_2per$sample
)]
s2_sample_order_new <- s2_sample_order[which(
  s2_samples %in% s2_samples_2per$sample
)]


s1_samples <- sample_names$X1[s1_sample_order_new]
s2_samples <- sample_names$X1[s2_sample_order_new]
historical_samples <- sample_names$X1[historical_sample_order]
sample_name_order <- c(
  s1_sample_order_new,
  s2_sample_order_new,
  historical_sample_order
)
sample_names_plot <- c(s1_samples, s2_samples, historical_samples)

#Sort samples
admix_results_3_sorted <- admix_results_3[
  c(s1_sample_order_new, s2_sample_order_new, historical_sample_order),
]

color_beh <- c("#1F78B4", "#FF7F00", "#33A02C", "#6A3D9A", "#B3B3B3")

#Make the barplot
admix_plot3 <- barplot(
  t(as.matrix(admix_results_3_sorted)),
  col = color_beh[c(2, 1, 3)],
  names.arg = sample_names_plot,
  ylab = "Ancestry",
  border = NA,
  beside = F,
  cex.names = 0.3,
  las = 3
)
