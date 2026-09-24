# r
# Resurvey_UCE_analysis_selected.R
# Contains code for:
# Analyzing vcftools output: pi, Tajima's D
# Genlight analyses: PCA, isolation by distance, allele sharing distances, heterozygosity, and heirarchical variance components (hierfstat)

# load packages
pacman::p_load(
  BiocManager,
  tidyverse,
  ggplot2,
  vcfR,
  SNPRelate,
  gdsfmt,
  SNPfiltR,
  adegenet,
  spdep,
  ape,
  poppr,
  pegas,
  GenoPop,
  paletteer,
  hierfstat,
  geosphere,
  ggpubr,
  rstatix,
  geodist,
  reshape2,
  viridis,
  gplots
)


set.seed(123)

#plot colors
color_beh <- c("#56B4E9", "#E69F00", "#8F7EE5", "#009E73", "lightgrey")


#Folder with vcftools popgen results
vcftools_results_folder <- "/Users/rachelweinberg/Documents/UCB/lab_stuff/LhumGenomic/UCEanalysis/vcftools_popgen_results_resurvey_noEPOW_MP23ref"
window_pi_file <- file.path(
  vcftools_results_folder,
  "win_pi_1kb_noEPOW_250717.tsv"
)

# Import vcf and behavior/location metadata

behavior_data_path <- "/Users/rachelweinberg/Documents/UCB/lab_stuff/Lhum_aggression/CA_Lhum_graphable2.csv" #Replace with filepath for behavior/location metadata csv
Lhum_behavior <- read_csv(
  behavior_data_path,
  col_names = TRUE
)
colnames(Lhum_behavior) <- c(
  "location",
  "Survey",
  "ID",
  "CHC_GeneticCluster",
  "behavior",
  "lat",
  "lon"
)

Lhum_loc <- data.frame(
  lat = Lhum_behavior$lat,
  lon = Lhum_behavior$lon,
  row.names = Lhum_behavior$ID
)
geo_dist <- geodist(Lhum_loc, measure = "geodesic")
rownames(geo_dist) <- Lhum_behavior$ID
colnames(geo_dist) <- Lhum_behavior$ID

vcf1 <- read.vcfR(
  "/Users/rachelweinberg/Documents/UCB/lab_stuff/LhumGenomic/UCEanalysis/resurvey_MP23_allfilters_no_EPOW_250714.recode.vcf"
)

vcf_hist <- read.vcfR(
  "/Users/rachelweinberg/Documents/UCB/lab_stuff/LhumGenomic/UCEanalysis/vcf_all_samples_ab_filtered_100_0.3_noEPOW.vcf"
)

vcf <- vcf1 #Change vcf depending on whether analyzing resurvey or resurvey + historical samples

# Set population map and behavior assignments (popmap)
popmap <- data.frame(
  id = colnames(vcf@gt)[2:length(colnames(vcf@gt))],
  pop = ""
)

years <- c("s1", "s2")
years_survey1 <- c(97, 98, 99, 00, 03, 04, 07)
for (i in years_survey1) {
  popmap$pop[grep(i, popmap$id)] <- "s1"
}
years_survey2 <- c(22, 23, 24)
for (i in years_survey2) {
  popmap$pop[grep(i, popmap$id)] <- "s2"
}
popmap_yr <- popmap

pops <- c(
  "UK",
  "MT",
  "SB",
  "LS",
  "LJM",
  "SM",
  "SJ",
  "SAL",
  "SAU",
  "SLO",
  "WL",
  "LA",
  "PA",
  "LP",
  "LH",
  "KC",
  "MB",
  "MP",
  "LJC",
  "TWC"
)
for (i in pops) {
  popmap$pop[grep(i, popmap$id)] <- i
}
popmap_loc <- popmap

for (i in pops) {
  for (j in c(years_survey1, years_survey2)) {
    popyr <- paste(i, j, sep = "")
    popmap$pop[grep(popyr, popmap$id)] <- popyr
  }
}
popmap$pop[popmap$id %in% c("LH0001", "LH0003")] <- "LH00"
popmap$pop[popmap$id %in% c("LS9701", "LSS9702", "LS9703", "LS9704")] <- "LS97"
popmap$pop[popmap$id %in% c("MP0703", "MP0704", "MP0701", "MP0702")] <- "MP07"
popmap$pop[popmap$id %in% c("WL0301", "WL0304")] <- "WL03"
popmap$pop[popmap$id %in% c("LJC0008", "LJC0009")] <- "LJC00"
popmap$pop[popmap$id %in% c("LP0401", "LP0402", "LP0406", "LP0407")] <- "LP04"

popmap_yrloc <- popmap

#Assign pops by supercolony
LSC <- c(
  "UK22",
  "UK98",
  "KC23",
  "KC98",
  "LP23",
  "LP04",
  "MB23",
  "MB98",
  "MP23",
  "MP07",
  "MT23",
  "MT98",
  "LJM98",
  "SB98",
  "SM98",
  "SM22",
  "SJ98",
  "PA99",
  "SAL98",
  "SAL22",
  "SAU98",
  "SAU22",
  "SLO98",
  "SLO22",
  "LJC00",
  "LA98",
  "WL03",
  "WL22"
)

LS <- c("LS97", "LS23")

LH <- c(
  "LH00",
  "LH23",
  "TWC99",
  "TWC22",
  "LA22",
  "SB22",
  "SJ22",
  "LJM22",
  "LA22",
  "LJM22",
  "LJC22",
  "PA24"
)
beh_df <- data.frame(
  pop = c(LSC, LS, LH),
  beh = c(rep("LSC", length(LSC)), rep("LS", length(LS)), rep("LH", length(LH)))
)
all_pops <- c(LSC, LS, LH)

#Set up popmap
for (i in 1:nrow(beh_df)) {
  popmap$pop[grep(beh_df$pop[i], popmap$id)] <- beh_df$beh[i]
}
popmap_beh <- popmap
popmap_beh_yr <- left_join(popmap_beh, popmap_yr, by = "id")
popmap_all <- left_join(popmap_beh_yr, popmap_loc, by = "id")
colnames(popmap_all) <- c("id", "beh", "yr", "loc")

#subset pops by timepoint,and stable vs changed status
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
  "EPOW22",
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
LH_s2_changed <- sort(c("SB22", "SJ22", "LA22", "LJM22", "LJC22", "PA24"))

LSC_stable_sites <- sort(c(
  "KC",
  "MP",
  "MB",
  "LP",
  "MT",
  "SAL",
  "SAU",
  "SLO",
  "SM",
  "UK",
  "WL"
))

LSC_unchanged_s1 <- c(
  "KC98",
  "LP04",
  "MB98",
  "MP07",
  "MT98",
  "SAL98",
  "SAU98",
  "SLO98",
  "SM98",
  "UK98",
  "WL03"
)
LSC_changed_s1 <- c("LA98", "LJC00", "LJM98", "PA99", "SB98", "SJ98")


#df to store behavior only data
beh_df <- data.frame(
  pop = c(LSC, LS, LH),
  beh = c(rep("LSC", length(LSC)), rep("LS", length(LS)), rep("LH", length(LH)))
)


# define s1/s2 pop lists used later
s1_pops <- c(LH_s1, LSC_s1, "LS97")
s2_pops <- c(LH_s2, LSC_s2, "LS23")
for (i in 1:nrow(beh_df)) {
  popmap$pop[grep(beh_df$pop[i], popmap$id)] <- beh_df$beh[i]
}
popmap_beh <- popmap

popmap_beh_yr <- left_join(popmap_beh, popmap_yr, by = "id")
popmap_all <- left_join(popmap_beh_yr, popmap_loc, by = "id")
colnames(popmap_all) <- c("id", "beh", "yr", "loc")


### FUNCTIONS ###

#Helper functions to import data
read_tajima_data <- function(tajima_path) {
  tajima_files <- list.files(
    path = tajima_path,
    pattern = "\\.Tajima\\.D$",
    full.names = TRUE
  )
  all_tajima_df <- data.frame()
  for (file in tajima_files) {
    pop <- str_extract(basename(file), "[A-Z]{2,3}[0-9]{2}")
    tajima_data <- read.delim(file, sep = "\t")
    tajima_data$pop <- pop
    all_tajima_df <- rbind(all_tajima_df, tajima_data)
  }
  all_tajima_df <- all_tajima_df |> filter(!is.na(TajimaD))
  return(all_tajima_df)
}


# Function to read and process PI data
# r
read_pi_data <- function(pi_path) {
  pi_files <- list.files(
    path = pi_path,
    pattern = "\\.(sites|windowed)\\.pi$",
    full.names = TRUE
  )

  site_files <- pi_files[grepl("\\.sites\\.pi$", pi_files)]
  win_files <- pi_files[grepl("\\.windowed\\.pi$", pi_files)]

  all_site_pi_df <- data.frame()
  for (file in site_files) {
    pop <- sub("\\.(sites|windowed)\\.pi$", "", basename(file))
    pi_data <- read.delim(file, sep = "\t", stringsAsFactors = FALSE)
    pi_data$pop <- pop
    all_site_pi_df <- rbind(all_site_pi_df, pi_data)
  }

  all_win_pi_df <- data.frame()
  for (file in win_files) {
    pop <- sub("\\.(sites|windowed)\\.pi$", "", basename(file))
    pi_data <- read.delim(file, sep = "\t", stringsAsFactors = FALSE)
    pi_data$pop <- pop
    all_win_pi_df <- rbind(all_win_pi_df, pi_data)
  }

  return(list(site = all_site_pi_df, windowed = all_win_pi_df))
}
#Subsample PCA plots to one sample/pop
subsample_PCA <- function(PCA) {
  PCA <- PCA |> group_by(pop, year) |> filter(row_number() == 2)
}

#Calculate mean pairwise difference by UCE locus
mean_pw_diff_per_locus <- function(distmat) {
  Gst_mean <- distmat |>
    group_by(CHROM) |>
    summarise(across(starts_with("Gst"), ~ mean(.x, na.rm = TRUE)))
  Gst_se <- distmat |>
    group_by(CHROM) |>
    summarise(across(starts_with("Gst"), ~ sd(.x, na.rm = TRUE)))
  n_variants <- distmat |> group_by(CHROM) |> summarise(n_variants = n())

  result1 <- merge(Gst_mean, Gst_se, by = "CHROM", suffixes = c(".mean", ".se"))
  result2 <- merge(result1, n_variants, by = "CHROM")
  return(result2)
}


#standard error function
se <- function(x, na.rm = FALSE) {
  if (na.rm == TRUE) {
    x <- x[is.na(x) == FALSE]
  }
  return(sd(x) / sqrt(length(x)))
}

#pairwise geographic distance
get_pw_geodist <- function(distmat, pop1, pop2) {
  pop1_key <- if (pop1 %in% rownames(distmat)) {
    pop1
  } else {
    str_extract(pop1, "^[A-Z]+")
  }
  pop2_key <- if (pop2 %in% colnames(distmat)) {
    pop2
  } else {
    str_extract(pop2, "^[A-Z]+")
  }

  if (
    is.na(pop1_key) ||
      !pop1_key %in% rownames(distmat) ||
      is.na(pop2_key) ||
      !pop2_key %in% colnames(distmat)
  ) {
    name <- paste(pop1, "_", pop2, sep = "")
    return(c(name, NA_real_))
  }

  dist <- distmat[pop1_key, pop2_key]
  name <- paste(pop1, "_", pop2, sep = "")
  return(c(name, dist))
}


### Convert data formats

#Convert to genlight / genind and create HS / loci objects
genl <- vcfR2genlight(vcf)
ploidy(genl) <- 2
gen <- vcfR2genind(vcf, sep = "[/]", return.alleles = TRUE)
strat <- data.frame(
  loc = popmap_loc$pop,
  yr = popmap_yr$pop,
  beh = popmap_beh$pop
)
strata(gen) <- strat
setPop(gen) <- ~ beh / loc / yr
strata(genl) <- strat
setPop(genl) <- ~ beh / loc / yr

gp <- genind2genpop(gen)
gen_loci <- as.loci(gen, ploidy = 2)
HS_loci <- genind2hierfstat(gen)

gens1_loci <- gen_loci |> filter(grepl("s1", population))
gens2_loci <- gen_loci |> filter(grepl("s2", population))
LSC_loci <- gen_loci |> filter(grepl("LSC", population))
LSC_s1_loci <- LSC_loci |> filter(grepl("s1", population))
LSC_s2_loci <- LSC_loci |> filter(grepl("s2", population))
LH_loci <- gen_loci |> filter(grepl("LH", population))
LS_loci <- gen_loci |> filter(grepl("LS", population))


### Nucleotide diversity analysis

# Read site and window pi outputs from vcftools
pi_path <- paste0(vcftools_results_folder, "/pi_dat")
site_pi <- read_pi_data(pi_path)
site_pi_all <- read_table(
  window_pi_file
) |>
  mutate(
    Timepoint = ifelse(KC23 %in% s1_pops, "S1", "S2"),
    PI = as.numeric(PI)
  ) |>
  filter(!is.na(PI)) |>
  mutate(
    group = case_when(
      KC23 %in% LSC_unchanged_s1 ~ "LSC_unchanged",
      KC23 %in% LSC_changed_s1 ~ "LSC_changed",
      KC23 %in% LSC_s2 ~ "LSC_unchanged",
      KC23 %in% c(LH_s1, "LH23", "TWC22") ~ "LH_unchanged",
      KC23 %in% LH_s2_changed ~ "LH_changed",
      KC23 %in% LS ~ "SK"
    )
  )

colnames(site_pi_all)[1] <- "pop"

site_pi_s1 <- site_pi_all |> filter(pop %in% s1_pops)

ggplot(site_pi_s1, aes(x = group, y = PI)) + geom_boxplot()

pistats_window <- read_delim(
  window_pi_file,
  delim = "\t",
  escape_double = FALSE,
  col_types = cols(
    N_VARIANTS = col_number(),
    PI = col_number(),
    BIN_START = col_number(),
    BIN_END = col_number()
  ),
  trim_ws = TRUE
)

colnames(pistats_window) <- c(
  "POP",
  "CHROM",
  "BIN_START",
  "BIN_END",
  "N_VARIANTS",
  "PI"
)
pistats_window <- pistats_window[is.na(pistats_window$N_VARIANTS) == FALSE, ] |>
  mutate(
    Timepoint = case_when(POP %in% s1_pops ~ "S1", POP %in% s2_pops ~ "S2"),
    pop = case_when(
      POP %in% LH ~ "LH",
      POP %in% LSC ~ "LSC",
      POP %in% LS ~ "LS"
    )
  )

# Summary and plot
pistats_win_summary <- pistats_window |>
  group_by(pop, Timepoint) |>
  summarise(mean_pi = mean(PI), sd_pi = sd(PI), n_windows = n())
ggboxplot(
  data = pistats_window |> group_by(pop, Timepoint),
  x = "pop",
  y = "PI",
  fill = "pop",
  palette = c("#1F78B4", "#FF7F00", "#33A02C"),
  facet.by = "Timepoint",
  ylim = c(0, 0.003)
) +
  labs(
    title = "1kbp window pi across supercolonies",
    subtitle = "vcftools window pi"
  )

## PCA
# Calculate PCA
genlight_PCA <- glPca(genl, nf = 4)
PCA_scores <- data.frame(genlight_PCA$scores[, 1:4])
PCA_scores$pop <- 0
PCA_scores$year <- 0
PCA_scores$supercolony <- ""
for (i in pops) {
  PCA_scores$pop[grep(i, rownames(PCA_scores))] <- i
}
for (i in years_survey1) {
  PCA_scores$year[grep(i, rownames(PCA_scores))] <- "s1"
}
for (i in years_survey2) {
  PCA_scores$year[grep(i, rownames(PCA_scores))] <- "s2"
}
for (i in 1:nrow(beh_df)) {
  PCA_scores$supercolony[grep(
    beh_df$pop[i],
    rownames(PCA_scores)
  )] <- beh_df$beh[i]
}

variance_explained <- (genlight_PCA$eig / length(genlight_PCA$eig)) * 100
colnames(PCA_scores)[1:4] <- c("PC1", "PC2", "PC3", "PC4")

#Subsample one individual per site for PCA visualization
PCA_subsampled <- PCA_scores |>
  group_by(pop, year) |>
  filter(row_number() == 1)
ggplot(
  PCA_subsampled,
  aes(x = PC1, y = PC2, color = supercolony, shape = year)
) +
  geom_point(size = 3) +
  labs(
    title = "PCA PC1 vs PC2",
    x = paste0("PC1 (", round(variance_explained[1], 2), "%)"),
    y = paste0("PC2 (", round(variance_explained[2], 2), "%)")
  ) +
  theme_linedraw()


## Isolation by distance analyses

#Calculate distances (Dch = CSE chord distance)
Dch <- genet.dist(HS_loci, diploid = T, method = "Dch")
Dch_mat <- as.matrix(Dch)
rownames(Dch_mat) <- gsub("^[^_]*_", "", rownames(Dch_mat))
colnames(Dch_mat) <- gsub("^[^_]*_", "", colnames(Dch_mat))

s1_names <- rownames(Dch_mat)[grepl("_s1$", rownames(Dch_mat))]
s2_names <- rownames(Dch_mat)[grepl("_s2$", rownames(Dch_mat))]
s1_base <- gsub("_s1$", "", s1_names)
s1_order <- sort(s1_base)
s1_ordered <- paste0(s1_order, "_s1")
Dch_s1_matrix <- Dch_mat[s1_ordered, s1_ordered]
s2_ordered <- paste0(s1_order, "_s2")
s2_final <- s2_ordered[s2_ordered %in% s2_names]
Dch_s2_matrix <- Dch_mat[s2_final, s2_final]
rownames(Dch_s1_matrix) <- str_remove(rownames(Dch_s1_matrix), "_s1")
colnames(Dch_s1_matrix) <- str_remove(colnames(Dch_s1_matrix), "_s1")
rownames(Dch_s2_matrix) <- str_remove(rownames(Dch_s2_matrix), "_s2")
colnames(Dch_s2_matrix) <- str_remove(colnames(Dch_s2_matrix), "_s2")

# Convert geo_dist to pop-level names
rownames(geo_dist) <- str_extract(rownames(geo_dist), "[A-Z]+")
colnames(geo_dist) <- str_extract(colnames(geo_dist), "[A-Z]+")

# Prepare dataframes
# pairwise geodist LSC stable sites
LSC_stable_sites <- sort(c(
  "KC",
  "MP",
  "MB",
  "LP",
  "MT",
  "SAL",
  "SAU",
  "SLO",
  "SM",
  "UK",
  "WL"
))

pw_geo_df_LSC_stable <- data.frame(pair = "", geo_diff = "")
for (i in LSC_stable_sites) {
  for (j in LSC_stable_sites) {
    gd_pair_stable <- get_pw_geodist(geo_dist, i, j)
    pw_geo_df_LSC_stable <- rbind(pw_geo_df_LSC_stable, gd_pair_stable)
  }
}

pw_geo_df_all <- data.frame(pair = "", geo_diff = "")
for (i in pops) {
  for (j in pops) {
    gd_pair_all <- get_pw_geodist(geo_dist, i, j)
    pw_geo_df_all <- rbind(pw_geo_df_all, gd_pair_all)
  }
}

Dch_df_s1 <- melt(Dch_s1_matrix) |> mutate(timepoint = "S1")
Dch_df_s2 <- melt(Dch_s2_matrix) |> mutate(timepoint = "S2")
Dch_df_all <- rbind(Dch_df_s1, Dch_df_s2)

Dch_1_df_allcols <- melt(Dch_s1_matrix, varnames = c("pop1", "pop2")) |>
  mutate(pair = paste0(pop1, "_", pop2), timepoint = "s1") |>
  inner_join(pw_geo_df_all) |>
  select(pair, timepoint, value, geo_diff)
Dch_2_df_allcols <- melt(Dch_s2_matrix, varnames = c("pop1", "pop2")) |>
  mutate(pair = paste0(pop1, "_", pop2), timepoint = "s2") |>
  inner_join(pw_geo_df_all) |>
  select(pair, timepoint, value, geo_diff)
Dch_gen_geo_dists_all <- rbind(Dch_1_df_allcols, Dch_2_df_allcols) |>
  mutate(geo_diff_km = as.numeric(geo_diff) / 1000) |>
  filter(geo_diff_km > 0)

# Convert chord distances to dataframe
Dch_1_df <- melt(Dch_s1_matrix, varnames = c("pop1", "pop2")) |>
  mutate(pair = paste0(pop1, "_", pop2), timepoint = "s1") |>
  inner_join(pw_geo_df_all) |>
  select(pair, timepoint, value, geo_diff)
Dch_2_df <- melt(Dch_s2_matrix, varnames = c("pop1", "pop2")) |>
  mutate(pair = paste0(pop1, "_", pop2), timepoint = "s2") |>
  inner_join(pw_geo_df_all) |>
  select(pair, timepoint, value, geo_diff)


#All colony matrix 1
Dch_1_df_allcols <- melt(Dch_s1_matrix, varnames = c("pop1", "pop2")) |>
  mutate(pair = paste0(pop1, "_", pop2), timepoint = "s1") |>
  mutate(
    relationship = case_when(
      (pop1 %in%
        c(str_extract(LH_s1, "[A-Z]+"), "LS") |
        pop2 %in% c(str_extract(LH_s1, "[A-Z]+"), "LS")) &
        (pop1 %in%
          str_extract(LSC_s1, "[A-Z]+") |
          pop2 %in% str_extract(LSC_s1, "[A-Z]+")) ~ "Between_supercolonies",
      (pop1 %in%
        str_extract(LH_s1, "[A-Z]+") &
        pop2 %in% str_extract(LH_s1, "[A-Z]+")) ~ "Within_supercolonies_LH",
      .default = "Within_supercolonies_LSC"
    )
  ) |>
  #filter(pop1 != pop2) |>
  inner_join(pw_geo_df_all) |>
  select(pair, timepoint, relationship, value, geo_diff)
#All colony matrix 2
Dch_2_df_allcols <- melt(Dch_s2_matrix, varnames = c("pop1", "pop2")) |>
  mutate(pair = paste0(pop1, "_", pop2), timepoint = "s2") |>
  mutate(
    relationship = case_when(
      ((pop1 %in%
        c(str_extract(LH_s2, "[A-Z]+"), "LS") &
        pop2 %in% str_extract(LSC_s2, "[A-Z]+")) |
        (pop2 %in% c(str_extract(LH_s2, "[A-Z]+"), "LS")) &
          (pop1 %in% str_extract(LSC_s2, "[A-Z]+")) |
        pop1 %in% c("LS") |
        pop2 %in% c("LS")) ~ "Between_supercolonies",
      (pop1 %in%
        str_extract(LH_s2, "[A-Z]+") &
        pop2 %in% str_extract(LH_s2, "[A-Z]+")) ~ "Within_supercolonies_LH",
      .default = "Within_supercolonies_LSC"
    )
  ) |>
  #filter(pop1 != pop2) |>
  inner_join(pw_geo_df_all) |>
  select(pair, timepoint, relationship, value, geo_diff)
#Join both chord distnace timepoints to get df for plotting

Dch_gen_geo_dists_all <- rbind(Dch_1_df_allcols, Dch_2_df_allcols) |>
  mutate(geo_diff_km = as.numeric(geo_diff) / 1000) |>
  filter(geo_diff_km > 0)

#Join together, convert m to km, and remove identical pairs and samples that changed colony identity from S1 to S2
changed_pops <- c(
  "LA_s1",
  "EPOW_s2",
  "SB_s1",
  "SJ_s1",
  "LJC_s1",
  "LJM_s1",
  "PA_s1"
)


# Mantel test for LSC samples with 9999 permutations
Dch_s1_matrix_LSC <- Dch_s1_matrix[LSC_stable_sites, LSC_stable_sites]
mantel_s1_LSC <- mantel.test(
  Dch_s1_matrix_LSC,
  geo_dist[LSC_stable_sites, LSC_stable_sites],
  nperm = 9999,
  graph = T
)

mantel_s2_LSC <- mantel.test(
  Dch_s2_matrix[LSC_stable_sites, LSC_stable_sites],
  geo_dist[LSC_stable_sites, LSC_stable_sites],
  nperm = 9999,
  graph = T
)

#Mantel tests for LH samples
LH_s2_pops <- str_extract(LH_s2, "[A-Z]+")
Dch_s2_matrix_LH <- Dch_s2_matrix[LH_s2_pops, LH_s2_pops]
geo_matrix_LH_s2 <- geo_dist[LH_s2_pops, LH_s2_pops]


mantel_LH_s2 <- mantel.test(
  Dch_s2_matrix_LH,
  geo_matrix_LH_s2,
  nperm = 9999,
  graph = T
)

#Plot genetic vs geographic distance for S1
ggplot(
  Dch_gen_geo_dists_all |> filter(timepoint == "s1"),
  aes(x = geo_diff_km, y = value, color = relationship)
) +
  geom_point() +
  stat_smooth(method = "lm", aes(fill = relationship), alpha = 0.3) +
  labs(
    title = "Genetic (Dch) vs Geographic distance (S1)",
    x = "Distance (km)",
    y = "Cavalli-Sforza & Edwards chord distance"
  ) +
  scale_color_manual(values = c(color_beh[4], color_beh[2], color_beh[1])) +
  scale_fill_manual(values = c(color_beh[4], color_beh[2], color_beh[1])) +
  guides(
    color = guide_legend(
      override.aes = list(fill = c(color_beh[4], color_beh[2], color_beh[1]))
    ),
    fill = "none"
  ) +
  theme_linedraw()

#Plot genetic vs geographic distance for S2
ggplot(
  Dch_gen_geo_dists_all |> filter(timepoint == "s2"),
  aes(x = geo_diff_km, y = value, color = relationship)
) +
  geom_point() +
  stat_smooth(method = "lm", aes(fill = relationship), alpha = 0.3) +
  labs(
    title = "Genetic (Dch) vs Geographic distance (S2)",
    x = "Distance (km)",
    y = "Cavalli-Sforza & Edwards chord distance"
  ) +
  scale_color_manual(values = c(color_beh[4], color_beh[2], color_beh[1])) +
  scale_fill_manual(values = c(color_beh[4], color_beh[2], color_beh[1])) +
  guides(
    color = guide_legend(
      override.aes = list(fill = c(color_beh[4], color_beh[2], color_beh[1]))
    ),
    fill = "none"
  ) +
  theme_linedraw()


# Heterozygosity (using Hierfstat)
stats <- basic.stats(HS_loci)
HS_all <- data.frame(stats$Hs) |>
  pivot_longer(everything(), names_to = "pop", values_to = "Hs") |>
  mutate(Timepoint = ifelse(str_detect(pop, "s1"), "S1", "S2"))
ggboxplot(
  HS_all |> group_by(pop, Timepoint),
  x = "pop",
  y = "Hs",
  fill = "Timepoint"
) +
  labs(title = "Expected heterozygosity (Hs)")

#Calculate Weir and Cockerham's (1984) Fst
#Calculate pairwise Fst for all populations
WC84_Fst_all <- genet.dist(HS_loci, diploid = T, method = "WC84")

#Calculate pairwise Fst for LSC populations
WC84_Fst_LSC_s1 <- genet.dist(
  HS_loci_s1 |> filter(str_detect(pop, "LSC")),
  diploid = T,
  method = "WC84"
)
WC84_Fst_LSC_s2 <- genet.dist(
  HS_loci_s2 |> filter(str_detect(pop, "LSC")),
  diploid = T,
  method = "WC84"
)

#Calculate pairwise Fst for LH populations
WC84_Fst_LH_s1 <- genet.dist(
  HS_loci_s1 |> filter(str_detect(pop, "LH")),
  diploid = T,
  method = "WC84"
)
WC84_Fst_LH_s2 <- genet.dist(
  HS_loci_s2 |> filter(str_detect(pop, "LH")),
  diploid = T,
  method = "WC84"
)

WC84_Fst_all <- genet.dist(HS_loci, diploid = T, method = "WC84")

WC84_Fst_df <- melt(as.matrix(WC84_Fst_all), varnames = c("pop1", "pop2")) |>
  mutate(
    pop1 = str_extract(pop1, pattern = pat1),
    pop2 = str_extract(pop2, pattern = pat1),
    pair = paste0(pop1, "_", pop2)
  )


# Hierarchical variance components and F-statistics (hierfstat)
pop_levels <- data.frame(
  time = factor(popmap_all$yr),
  beh = factor(popmap_beh$pop),
  loc = factor(popmap_all$loc)
)
HS_dat_levels <- data.frame(cbind(pop_levels), HS_loci)
HS_dat_levels_s1 <- HS_dat_levels |> filter(time == "s1")
HS_dat_levels_s2 <- HS_dat_levels |> filter(time == "s2")
HS_dat_levels_only_s1 <- HS_dat_levels_s1[, 2:3]
HS_dat_levels_only_s2 <- HS_dat_levels_s2[, 2:3]

varcomp_s1 <- varcomp.glob(
  levels = HS_dat_levels_only_s1,
  loci = HS_dat_levels_s1[, -c(1:4)]
)
varcomp_s2 <- varcomp.glob(
  levels = HS_dat_levels_only_s2,
  loci = HS_dat_levels_s2[, -c(1:4)]
)
varcomp_perc_s1_beh <- varcomp_s1$overall["beh"] / sum(varcomp_s1$overall) * 100
varcomp_perc_s2_beh <- varcomp_s2$overall["beh"] / sum(varcomp_s2$overall) * 100

# Permutation tests
loc_signif_within_beh_s1 <- test.within(
  data = HS_dat_levels_s1[, -c(1:4)],
  within = HS_dat_levels_only_s1$beh,
  test.lev = HS_dat_levels_only_s1$loc,
  nperm = 1000
)
loc_signif_within_beh_s2 <- test.within(
  data = HS_dat_levels_s2[, -c(1:4)],
  within = HS_dat_levels_only_s2$beh,
  test.lev = HS_dat_levels_only_s2$loc,
  nperm = 1000
)
beh_signif_s1 <- test.between(
  data = HS_dat_levels_s1[, -c(1:4)],
  rand.unit = HS_dat_levels_only_s1$loc,
  test.lev = HS_dat_levels_only_s1$beh,
  nperm = 1000
)
beh_signif_s2 <- test.between(
  data = HS_dat_levels_s2[, -c(1:4)],
  rand.unit = HS_dat_levels_only_s2$loc,
  test.lev = HS_dat_levels_only_s2$beh,
  nperm = 1000
)

# Bootstrap variance components
bootvc_beh_loc_s1 <- boot.vc(
  levels = HS_dat_levels_only_s1,
  loci = HS_dat_levels_s1[, -c(1:4)],
  nboot = 1000
)
bootvc_beh_loc_s2 <- boot.vc(
  levels = HS_dat_levels_only_s2,
  loci = HS_dat_levels_s2[, -c(1:4)],
  nboot = 1000
)


bootvc_beh_loc_df <- data.frame(
  rbind(
    data.frame(bootvc_beh_loc_s1[["ci"]]) |>
      mutate(timepoint = "s1", ci = rownames(bootvc_beh_loc_s1[["ci"]])) |>
      select(ci, everything()) |>
      pivot_longer(
        cols = -c(ci, timepoint),
        names_to = "level",
        values_to = "value"
      ) |>
      pivot_wider(names_from = ci, values_from = value, names_prefix = "ci_") |>
      rename(
        lower = `ci_2.5%`,
        estimate = `ci_50%`,
        upper = `ci_97.5%`
      ),
    data.frame(bootvc_beh_loc_s2[["ci"]]) |>
      mutate(timepoint = "s2", ci = rownames(bootvc_beh_loc_s2[["ci"]])) |>
      select(ci, everything()) |>
      pivot_longer(
        cols = -c(ci, timepoint),
        names_to = "level",
        values_to = "value"
      ) |>
      pivot_wider(names_from = ci, values_from = value, names_prefix = "ci_") |>
      rename(
        lower = `ci_2.5%`,
        estimate = `ci_50%`,
        upper = `ci_97.5%`
      )
  )
)

ggplot(
  data = bootvc_beh_loc_df,
  aes(x = level, y = estimate, color = timepoint)
) +
  geom_point(size = 3, position = position_dodge(width = 0.3)) +
  geom_errorbar(
    aes(ymin = lower, ymax = upper),
    width = 0.2,
    position = position_dodge(width = 0.3)
  ) +
  theme_linedraw() +
  theme(
    axis.text.x = element_text(angle = 90),
    panel.grid.minor = element_blank()
  ) +
  labs(
    title = "Bootstrapped variance components",
    x = "Variance component level",
    y = "Fst"
  ) +
  scale_color_manual(values = col_s1_v_s2)

# ---- Tajima's D ----
Tajima_path <- file.path(vcftools_results_folder, "TajimaD_by_UCE")

Tajima_df <- read_tajima_data(Tajima_path) |>
  mutate(
    Timepoint = case_when(pop %in% s1_pops ~ "S1", pop %in% s2_pops ~ "S2"),
    POP = case_when(
      pop %in% c("TWC99", "TWC22", "LH00", "LH23") ~ "LH",
      pop %in% LSC_changed_s1 ~ "LSC_changed",
      pop %in% LH_s2_changed ~ "LH_changed",
      pop %in% LS ~ "SK",
      .default = "LSC"
    ),
    colony = case_when(
      grepl("LSC", POP) ~ "LSC",
      grepl("LH", POP) ~ "LH",
      grepl("SK", POP) ~ "SK"
    ),
    changed = ifelse(grepl("changed", POP), TRUE, FALSE)
  ) |>
  filter(N_SNPS >= 2)


#Plot Tajima's D for known colonies
ggboxplot(
  data = Tajima_df |>
    mutate(colxtime = paste0(colony, "_", Timepoint)) |>
    group_by(colony, Timepoint),
  x = "colony",
  y = "TajimaD",
  fill = "colxtime",
  palette = c("#FDBF6F", "#FF7F00", "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C"),
  ylim = c(-2, 3)
) +
  labs(title = "Tajima's D in changed and unchanged sites")

# Plot for changed and unchanged sites
ggboxplot(
  data = Tajima_df |>
    group_by(POP, Timepoint),
  x = "POP",
  y = "TajimaD",
  fill = "POP",
  palette = c("#1F78B4", "#FDBF6F", "#A6CEE3", "#FF7F00", "#33A02C"),
  ylim = c(-2, 3)
) +
  labs(title = "Tajima's D in changed and unchanged sites") +
  facet_wrap(~Timepoint, scales = "free_x") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())


#Historgram of Tajima's D values for each colony and timepoint
ggplot(Tajima_df |> filter(N_SNPS >= 2), aes(x = TajimaD, fill = Timepoint)) +
  geom_histogram(alpha = 0.5) +
  facet_wrap(~colony) +
  theme_bw() +
  scale_fill_brewer(palette = "Set2")


Tajima_summary <- Tajima_df |>
  group_by(POP, Timepoint) |>
  summarise(
    mean_TajimaD = mean(TajimaD),
    sd_TajimaD = sd(TajimaD),
    n_windows = n()
  )

Tajima_site_summary <- Tajima_df |>
  group_by(pop, Timepoint, POP) |>
  summarise(
    mean_TajimaD = mean(TajimaD),
    sd_TajimaD = sd(TajimaD),
    n_windows = n()
  )


Tajima_summary_behavior <- Tajima_df |>
  group_by(colony, Timepoint) |>
  summarise(
    mean_TajimaD = mean(TajimaD),
    sd_TajimaD = sd(TajimaD),
    n_windows = n()
  )

write_csv(Tajima_summary, file = "TajimasD_by_col_timepoint.csv")


Tajima_Wilcox_LH_changed_vs_unchanged_LH <- compare_means(
  TajimaD ~ changed,
  data = Tajima_df[Tajima_df$colony == "LH" & Tajima_df$Timepoint == "S2", ],
  method = "wilcox.test",
  paired = F,
  p.adjust.method = "BH"
)

Tajima_Wilcox_LSC_changed_vs_unchanged_LSC <- compare_means(
  TajimaD ~ changed,
  data = Tajima_df[Tajima_df$colony == "LSC" & Tajima_df$Timepoint == "S1", ],
  method = "wilcox.test",
  paired = F,
  p.adjust.method = "BH"
)

Tajima_Wilcox_beh <- compare_means(
  TajimaD ~ colony,
  data = Tajima_df,
  method = "wilcox.test",
  paired = F,
  group.by = "Timepoint",
  p.adjust.method = "BH"
)

write_csv(Tajima_Wilcox_beh, file = "TajimasD_Wilcox_by_beh_250908.csv")

Tajima_Wilcox_time <- compare_means(
  TajimaD ~ Timepoint,
  data = Tajima_df,
  method = "wilcox.test",
  paired = F,
  group.by = "colony",
  p.adjust.method = "BH"
)

# Allelic Richness
AR_all <- allelic.richness(HS_loci)

HS_loci_beh <- HS_loci |> mutate(pop = str_replace(pop, "_[A-Z]+_", "_"))

AR_beh_time <- allelic.richness(HS_loci_beh)


AR_beh_df <- AR_beh_time$Ar |>
  pivot_longer(cols = everything(), names_to = "pop", values_to = "Ar") |>
  mutate(
    timepoint = str_extract(pop, "s[0-9]"),
    behavior = str_extract(pop, "[^_]+(?=_)")
  )


#Convert AR_s1 to long format to plot boxplot with ggpolot2
pop_pattern <- "(?<=_)[^_]+(?=_)"
AR_df <- AR_all$Ar |>
  pivot_longer(cols = everything(), names_to = "pop", values_to = "Ar") |>
  mutate(
    timepoint = str_extract(pop, "s[0-9]"),
    loc = str_extract(pop, pop_pattern),
    behavior = str_extract(pop, "[^_]+(?=_)")
  )


AR_summary_beh <- AR_beh_df |>
  group_by(timepoint, behavior) |>
  summarise(mean_AR = mean(Ar), n_sites = n(), sd_ar = sd(Ar))

AR_test_between_sites <- compare_means(
  Ar ~ loc,
  data = AR_df,
  method = "wilcox.test",
  paired = T,
  group.by = "timepoint",
  p.adjust.method = "BH"
) |>
  mutate(
    significance = case_when(
      p.adj <= 0.0001 ~ "****",
      p.adj <= 0.001 ~ "***",
      p.adj <= 0.01 ~ "**",
      p.adj <= 0.05 ~ "*",
      TRUE ~ "ns"
    )
  )

Ar_test_between_cols <- compare_means(
  Ar ~ behavior,
  data = AR_df,
  method = "wilcox.test",
  paired = F,
  group.by = "timepoint",
  p.adjust.method = "BH"
) |>
  mutate(
    significance = case_when(
      p.adj <= 0.0001 ~ "****",
      p.adj <= 0.001 ~ "***",
      p.adj <= 0.01 ~ "**",
      p.adj <= 0.05 ~ "*",
      TRUE ~ "ns"
    )
  )

AR_test_timepoint <- compare_means(
  Ar ~ timepoint,
  data = AR_beh_df,
  method = "wilcox.test",
  group.by = "behavior",
  p.adjust.method = "BH"
)


# Allele sharing distance / genetic distance via pegas
gen_loci2 <- genind2loci(gen)
gendist <- dist.gene(gen_loci2)
#allele sharing distance
asd <- dist.asd(gen_loci)

# heatmap of genetic distances by time
heatmap.2(
  Dch_s1_matrix,
  col = viridis(100),
  distfun = as.dist, # Use distance matrix directly
  hclustfun = function(x) hclust(x, method = "average"), # UPGMA cluster method
  density.info = "none",
  trace = "none",
  key = TRUE,
  keysize = 1.5,
  key.title = "Distance",
  margins = c(5, 5)
)

heatmap.2(
  Dch_s2_matrix,
  col = viridis(100),
  distfun = as.dist, # Use distance matrix directly
  hclustfun = function(x) hclust(x, method = "average"), # UPGMA cluster method
  density.info = "none",
  trace = "none",
  key = TRUE,
  keysize = 1.5,
  key.title = "Distance",
  margins = c(5, 5)
)

#ASD tree

# Get individual-level ASD
asd <- dist.asd(gen_loci, pairwise.deletion = TRUE)

# Convert to matrix
asd_mat <- as.matrix(asd)

# Aggregate by population/supercolony
# Extract population from individual names
ind_pops <- popmap_yrloc$pop

# Calculate mean ASD between populations
pop_names <- unique(ind_pops)
asd_pop_mat <- matrix(NA, nrow = length(pop_names), ncol = length(pop_names))
rownames(asd_pop_mat) <- pop_names
colnames(asd_pop_mat) <- pop_names

for (i in 1:length(pop_names)) {
  for (j in 1:length(pop_names)) {
    # Get individuals in each population
    inds_i <- which(ind_pops == pop_names[i])
    inds_j <- which(ind_pops == pop_names[j])

    # Calculate mean ASD between populations
    asd_pop_mat[i, j] <- mean(asd_mat[inds_i, inds_j], na.rm = TRUE)
  }
}
# Build and plot tree
asd_pop_dist <- as.dist(asd_pop_mat)
asd_tree <- nj(asd_pop_dist)
#Set palette for all pops
library(colorspace)
lsc_pal <- sequential_hcl(n = 38, palette = "Blues 3")
lh_pal <- c(
  "#88002D",
  "#AE0732",
  "#D4332A",
  "#EC5D2F",
  "#F9834C",
  "#FFA46C",
  "#FFC08E",
  "#E4739D",
  "#F28891",
  "#F7A086"
)
sk_pal <- c("#009E73", "#65ffd6")
asd_pal <- c(
  lh_pal[1:2],
  sk_pal[2],
  lsc_pal[1:17],
  lh_pal[3:8],
  lsc_pal[18:38],
  sk_pal[1],
  lh_pal[9:10]
)


plot(asd_tree, cex = 0.7, tip.color = asd_pal) #Need to manually recolor in illustrator
title("Population phylogeny - Allele Sharing Distance")
