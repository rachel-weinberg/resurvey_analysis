pacman::p_load(
  readr,
  ggplot2,
  ggmap,
  tidyverse,
  sf,
  readxl,
  sf,
  ggspatial,
  rnaturalearth,
  rnaturalearthdata,
  ggrepel,
  cowplot,
  xml2,
  jsonlite,
  lubridate,
  usmap,
  raster,
  rgdal,
  reshape2,
  geosphere
)

#set theme for all maps and plots
theme_set(theme_bw())

###Mapping###

world <- ne_countries(scale = "medium", returnclass = "sf")

#Read in map data coordinates
#Aggression data can be downloaded from supplemental data
Lhum_behavior <- read_xlsx(
  "AggressionAssayResults_allsites_alltimepoints.xlsx",
  sheet = "Sites"
)

#Read in collection coordinates
Lhum_geo <- read_csv("CA_Lhum_graphable2.csv")

#Extract coordinates and points
coords <- Lhum_geo %>% select(CHC_GeneticCluster, Lat, Long, Category)

points <- coords %>% st_as_sf(coords = c("Lat", "Long"))

#Mapping colors
color_beh <- c("#1F78B4", "#FF7F00", "#33A02C", "#6A3D9A", "#B3B3B3")

#Plot with Survey 2 points only
Lhum_points <- Lhum_geo |> filter(Survey == 2)


Lhum_colors <- c(
  LSC_Changed = "white",
  LSC_NoChange = color_beh[1],
  LH_NoChange = color_beh[2],
  LS_NoChange = color_beh[3]
)

#Get stadia basemap
california_stadia <- get_stadiamap(
  c(left = -124.5, bottom = 31.5, right = -114.5, top = 42),
  zoom = 7,
  color = "bw",
  maptype = "stamen_toner_lines"
)
ggmap(california_stadia)


#Plot statewide map Survey 2
ca_s2 <- ggmap(california_stadia) +
  geom_point(
    data = Lhum_behavior[Lhum_behavior$Survey == 2, ],
    mapping = aes(x = lon, y = lat, color = CHC_GeneticCluster),
    size = 1.5
  ) +
  scale_color_manual(
    breaks = c("LH", "LSC", "LS", "other"),
    values = c("#990F0F", "#2C85B2", "#6B990F", "#8F7EE5")
  )
ca_s2

#Zoomed in map for SoCal sites
SoCal_stadia <- get_stadiamap(
  c(left = -121.5, bottom = 32, right = -116, top = 35.75),
  source = "stamen",
  zoom = 7,
  maptype = "stamen_toner_lines"
)
ggmap(SoCal_stadia)

SoCal_behavior_s2 <- ggmap(SoCal_stadia) +
  ggtitle("L. humile colony identity, survey 2") +
  geom_point(
    data = SoCal_only[SoCal_only$Survey == 2, ],
    mapping = aes(x = lon, y = lat, color = behavior),
    size = 2
  ) +
  scale_color_manual(
    breaks = c("LH", "LSC", "LS", "other"),
    values = c("#990F0F", "#2C85B2", "#6B990F", "#8F7EE5")
  )

SoCal_behavior_s2


###Behavioral Data###

#Boxplots of behavioral data for all sites
all_aggression_boxplot <- ggplot(
  aggression |>
    filter(
      ReferenceColony %in% c("LH", "LSC", "LS"),
      FocalColonyGenetic_or_CHCID %in% c("LH", "LSC", "LS")
    ) |>
    mutate(
      ReferenceColony = factor(ReferenceColony, levels = c("LSC", "LH", "LS")),
      ReferenceColony = recode(ReferenceColony, "LS" = "SK"),
      FocalColonyGenetic_or_CHCID = factor(
        FocalColonyGenetic_or_CHCID,
        levels = c("LSC", "LH", "LS")
      ),
      FocalColonyGenetic_or_CHCID = recode(
        FocalColonyGenetic_or_CHCID,
        "LS" = "SK"
      )
    ),
  aes(
    x = ReferenceColony,
    y = Score,
    fill = ReferenceColony,
    color = ReferenceColony
  )
) +
  geom_boxplot() +
  facet_grid(FocalColonyGenetic_or_CHCID ~ YearBin) +
  scale_fill_manual(values = c(color_beh[1], color_beh[2], color_beh[3])) +
  scale_color_manual(values = c(color_beh[1], color_beh[2], color_beh[3])) +
  theme_linedraw() +
  theme(panel.grid.minor = element_blank(), legend.position = "none")

#Boxplots of behavioral data for changed sites only
site_order <- c("PA", "SJ", "SB", "LA", "LJM-E", "LJC") #Order north to south

beh_boxplot <- ggplot(
  data = aggression |>
    filter(
      ReferenceColony %in% c("LH", "LSC"),
      FocalColony %in% changed_sites
    ) |>
    mutate(FocalColony = factor(FocalColony, levels = site_order)),
  aes(x = factor(YearBin), y = Score, color = ReferenceColony)
) +
  geom_boxplot(size = 1) +
  facet_wrap(~FocalColony, ncol = 1, strip.position = "right") +
  scale_color_manual(values = c(color_beh_new[2], color_beh_new[1])) +
  theme_bw()

beh_boxplot

#La Jolla Mesa and La Jolla Coast aggression results by date
LJ_agg <- read_excel(
  "AggressionAssayResults_allsites_alltimepoints.xlsx",
  sheet = 2
) |>
  filter(
    FocalColony %in% c("LJM-E", "LJC"),
    YearBin == "2022-2024",
    Date != as.POSIXct("2023-10-27", tz = "UTC")
  ) |>
  filter(ReferenceColony %in% c("LH", "LSC")) |>
  filter(Score > 0)

LJ_summary <- LJ_agg |>
  group_by(Date, ReferenceColony, FocalColony) |>
  summarise(meanScore = mean(Score), sd = sd(Score))

LJ_summary$Date <- as.Date(LJ_summary$Date)

LJ_agg_plot <- ggplot(
  data = LJ_summary,
  aes(x = Date, y = meanScore, color = ReferenceColony)
) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin = meanScore - sd / 2, ymax = meanScore + sd / 2)) +
  scale_color_manual(values = c(color_beh[2], color_beh[1])) +
  facet_wrap(~FocalColony, scales = "free_x") +
  scale_x_date(
    breaks = unique(LJ_summary$Date),
    date_labels = "%Y-%m-%d",
    guide = guide_axis(check.overlap = TRUE)
  ) +
  labs(
    title = "Aggression scores for La Jolla Mesa and La Jolla Coast 2022-2024"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

LJ_agg_plot
