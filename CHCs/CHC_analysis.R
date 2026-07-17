pacman::p_load(tidyverse, ggplot2, GCalignR, vegan, patchwork)

#OpenChrom integration and peak calling parameters are described in the Detailed Methods portion of the Supplemental Information
#CSVs with peak retention times and areas can be found in the Dryad Repository:

#data import
folder <- "/Users/rachelweinberg/Documents/UCB/lab_stuff/CHCs/field_sample_analysis/strict_peak_calls" #OpenChrom export folder
pattern <- "_integrated.csv"
files <- list.files(folder, pattern = pattern, full.names = TRUE, recursive = T)

sample_pattern <- "(?<=/\\d{6}_)[^_]+(?=_)"

sample_names <- str_extract(files, pattern = sample_pattern)

CHC_dat <- list()

for (i in 1:length(files)) {
  CHC_dat[[i]] <- read.csv(files[i], header = TRUE, stringsAsFactors = FALSE)
  names(CHC_dat)[i] <- sample_names[i]
  CHC_dat[[i]]$Name <- sample_names[i]
  #Remove periods from column names
  colnames(CHC_dat[[i]]) <- gsub("\\.", "", colnames(CHC_dat[[i]]))
}


#Make RT and Area columns numeric
CHC_dat <- lapply(CHC_dat, function(x) {
  x$RT <- as.numeric(x$RT)
  x$Area <- as.numeric(x$Area)
  #subset to only include RT and Area columns
  x <- x[, c("RT", "Area")]
  return(x)
})

#Examnine the first 5 RT values for each sample to calculate a linear correction to apply so the first RT is the dodecane standard
first_RT_values <- sapply(CHC_dat, function(x) x$RT[1:5])

check_input(CHC_dat, plot = T)

first_RT_values[1, ]


#Check peak interspace to set cutoff for peak2peak alignment
peak_interspace(
  data = CHC_dat,
  rt_col_name = "RT",
  quantile_range = c(0, 0.8),
  quantiles = 0.5
)

#Choose optimal reference
choose_optimal_reference(data = CHC_dat, rt_col_name = "RT")

#linear transform alignment with reference ""231205_RBW140_KingCityLhum""
transformed_peaks <- linear_transformation(
  CHC_dat,
  reference = "NDT932",
  rt_col_name = "RT",
  max_linear_shift = 0.1,
  step_size = 0.01
)
#Put transformed peak data back into the same format as CHC_dat for analysis
transformed_CHCs <- list()
for (i in 1:length(transformed_peaks$chroma_aligned)) {
  transformed_CHCs[[i]] <- transformed_peaks$chroma_aligned[[i]]$shifted
  names(transformed_CHCs)[i] <- names(transformed_peaks$chroma_aligned)[i]
}

#Try aligning from 16-35 min only
aligned_peaks <- align_chromatograms(
  data = transformed_CHCs,
  rt_col_name = "RT",
  rt_cutoff_low = 16,
  rt_cutoff_high = 32,
  reference = "NDT932",
  max_linear_shift = 0.1, # max. shift for linear corrections
  max_diff_peak2mean = 0.2, # max. distance of a peak to the mean across samples
  min_diff_peak2peak = 0.01, # min. expected distance between peaks
  blanks = NULL, # negative control
  delete_single_peak = TRUE, # delete peaks that are present in just one sample
  write_output = NULL
)


gc_heatmap(aligned_peaks)

normed_peaks <- norm_peaks(
  aligned_peaks,
  rt_col_name = "RT",
  conc_col_name = "Area",
  out = "data.frame"
)
#Normalize peaks to C12 standard, the first column in the normed_peaks dataframe (`6.39127272727273)
normed_peaks_t <- t(normed_peaks)
write_csv(
  data.frame(cbind(RT = rownames(normed_peaks_t)), normed_peaks_t),
  file = "transposed_normed_peaks_all_field_samples_250828.csv"
)

#Log-transform normalized peaks
normed_peaks_log <- log(normed_peaks + 1)

heatmap(as.matrix(normed_peaks_log))


#Distance calculations and NMDS plots of normed peaks
#set seed for reproducibility
set.seed(12345)

#normed_peaks_log <- normed_peaks_log[-c(1:4),]

pkdist <- metaMDS(normed_peaks_log, distance = "bray", try = 1000)
pkdist_tidy <- metaMDS(normed_peaks_log, distance = "bray", try = 100, tidy = T)
pkpca <- pca(as.matrix(normed_peaks_log))

dist_pts <- data.frame(pkdist$points)

pca_pts <- data.frame(
  sample = rownames(pkpca[["CA"]][["u"]]),
  PC1 = pkpca[["CA"]][["u"]][, 1],
  PC2 = pkpca[["CA"]][["u"]][, 2]
)
#Add colony (rowname) as factor to dist_pts
dist_pts$colony <- factor(rownames(dist_pts), levels = rownames(dist_pts))
dist_pts$location <- c(
  "KC",
  "MB",
  "PA",
  "WL",
  "LJM",
  "LJC",
  "LA",
  "SB",
  "SJ",
  "CK",
  "AB"
)
pca_pts$colony <- factor(rownames(dist_pts), levels = rownames(dist_pts))
# dist_pts$beh <- c(rep("LSC", 2), "LH", "LSC", rep("LH", 5), "LH", "LSC", "LS")

dist_pts$beh <- c(rep("LSC", 2), "LH", "LSC", rep("LH", 5), "LH", "LSC")
dist_pts$role <- c(rep("sample", nrow(dist_pts) - 2), rep("reference", 2))

pca_pts$beh <- dist_pts$beh
pca_pts$location <- dist_pts$location


#dist_pts$beh <- c(rep("LSC", 2), "PA", "LSC", rep("LJ", 2), rep("LH", 2), "LSC", "SK")
#dist_pts$beh <- c(rep("LSC", 4), "PA", "LSC", rep("LJ", 2), rep("LH", 4), "LSC", "SK")

color_beh <- c("#FF7F00", "#33A02C", "#1F78B4")
color_beh <- c("#FF7F00", "#1F78B4")

CHC_dist_plot <- ggplot(
  dist_pts,
  aes(x = MDS1, y = MDS2, color = beh, shape = role)
) +
  geom_point(size = 3) +
  labs(
    title = "Cuticular Hydrocarbon Bray-Curtis Distance, All Peaks 6-34 Minutes",
    x = "MDS1",
    y = "MDS2"
  ) +
  theme_linedraw() +
  scale_color_manual(values = color_beh) +
  scale_shape_manual(values = c("reference" = 4, "sample" = 16)) +
  theme(panel.grid = element_blank()) +
  geom_label(
    aes(label = location),
    nudge_x = 0.03,
    nudge_y = 0.01,
    size = 3,
    #label.padding = unit(0.1, "lines"),
    label.size = 0.5,
    show.legend = FALSE
  )

ggplot(pca_pts, aes(x = PC1, y = PC2, color = beh)) +
  geom_point(size = 2) +
  labs(title = "PCA of CHC profiles", x = "PC1", y = "PC2") +
  theme_linedraw() +
  scale_color_manual(values = color_beh) +
  geom_label(
    aes(label = location),
    nudge_x = 0.025,
    nudge_y = 0.025,
    size = 3,
    label.size = 0.5,
    show.legend = FALSE
  )
