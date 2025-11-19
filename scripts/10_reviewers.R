# 10 - Additional figures for GEB reviewers
# November 2025

# PACKAGES ----

# Load in the required packages
library(tidyverse)
library(ggrepel)
library(Taxonstand)
library(rworldmap) # For mapping
library(sp) # For mapping
library(raster) # For mapping
# devtools::install_github("eliocamp/ggalt@new-coord-proj")
library(ggalt) # For mapping
library(gridExtra)
# devtools::install_github("MikkoVihtakari/ggOceanMapsData")
# devtools::install_github("MikkoVihtakari/ggOceanMaps")
library(ggOceanMaps)
library(viridis)


# THEME 4: MAPS produced using ggOceanMaps ----

theme_4_map <- function(){
  theme(text = element_text(family = "Helvetica Light"),
        plot.margin = unit(c(1, 1, 1, 1), units = , "cm"),
        plot.title = element_text(size = 16, vjust = 1, hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(size = 14, hjust = 0.5, face = "plain"),
        plot.caption = element_text(size = 12, hjust = 1, face = "plain"),
        legend.text = element_text(size = 12),          
        legend.title = element_text(size = 10, face = "bold", hjust = 0.5),
        legend.position = "right", 
        legend.key = element_blank(),
        legend.background = element_rect(color = "black", 
                                         fill = "transparent", 
                                         size = 2, linetype = "blank"),
        plot.background = element_rect(fill = "white", color = NA))
}


# MEAN AND SD TRAIT TABLE ----

# Load in trait data
funcgroups <- read.csv("data/output_fg_export_traits.csv")

traits <- read.csv("data/output_06_TRY_TTT_clean.csv") %>% 
  left_join(., funcgroups, by = c("NAME" = "NAME")) %>% 
  mutate(FuncGroup = ifelse(str_detect(FuncGroup, pattern = "hrub"), "Shrub", FuncGroup)) %>% 
  filter(!is.na(FuncGroup)) %>% 
  group_by(FuncGroup, TraitName) %>% 
  summarise(Mean = mean(TraitValue, na.rm = TRUE),
            SD = sd(TraitValue, na.rm = TRUE)) %>% 
  ungroup()

mean.traits <- traits %>% 
  dplyr::select(-SD) %>% 
  pivot_wider(names_from = "TraitName", values_from = "Mean")

sd.traits <- traits %>% 
  dplyr::select(-Mean) %>% 
  pivot_wider(names_from = "TraitName", values_from = "SD")

# Export to table
write.csv(sd.traits, file = "data/AAA_sd_table.csv", row.names = FALSE)


# LENGTH OF TIME PLOT ----

# Load final temporal data
temporal <- read.csv("data/AAA_temporal_input.csv")

# Work out legnth of years each sampled for
years <- temporal %>% 
  dplyr::select(SiteSubsitePlot, Duration) %>% 
  unique()

# Create histogram
(duration.plot <- ggplot(years) +
    geom_histogram(aes(x = Duration), stat = "bin", binwidth = 1, fill = "#D02090", colour = "#000000", alpha = 0.65) +
    geom_vline(aes(xintercept = mean(Duration)), colour = "#000000", linetype = "dashed", size = 0.5) +
    labs(y = "Plots \n",
         x = "\nTime Series Duration (Years)") +
    theme_1())

# Export to png
ggsave(duration.plot, filename = "figures/outputs_new/manuscript/duration_plot.png", width = 8, height = 5)


# MAP OF ITEX SITES - SPATIAL ----

# Import the ITEX and FD data for the latest year at each plot ONLY
spatial <- read.csv("data/output_08_fd_output_combined_nonPCA_3_latest_years.csv")

spatial.map <- spatial %>% 
  mutate(SITE = ifelse(SITE %in% c("BARROW"), "Utqiagvik", SITE)) %>% 
  group_by(SITE) %>% 
  mutate(Number_Plots = length(unique(PLOT))) %>% 
  ungroup() %>% 
  dplyr::select(LAT, LONG, SITE, Region, Number_Plots) %>% 
  distinct(SITE, Number_Plots, .keep_all = TRUE) %>% 
  arrange(Number_Plots) %>% 
  mutate(Number_Plots_Binned = case_when(Number_Plots < 10 ~ "1 - 9",
                                         Number_Plots >= 10 & Number_Plots < 50 ~ "10 - 49",
                                         Number_Plots >= 50 & Number_Plots < 100 ~ "50 - 99",
                                         Number_Plots >= 100 & Number_Plots < 200 ~ "100 - 199",
                                         Number_Plots >= 200 ~ "200 +")) %>% 
  mutate(Number_Plots_Binned = factor(.$Number_Plots_Binned,
                                      levels = c("1 - 9", "10 - 49", "50 - 99", "100 - 199", "200 +"))) %>% # Reorder factor
  mutate(SITE = str_to_title(SITE),
         Region = case_when(Region == "North America-East" ~ "North America East",
                            Region == "North America-West" ~ "North America West",
                            Region == "GreenIceland" ~ "Greenland & Iceland",
                            Region == "Eurasia" ~ "Eurasia",
                            NA ~ Region))

# Transform coordinates to UTM
spatial.map.transform <- transform_coord(spatial.map, lon = "LONG", lat = "LAT", bind = TRUE)

# Map based on size of circles
(output.spatial <- basemap(limits = 55, grid.size = 0.05, grid.col = "#949494",
                                land.size = 0.05, land.col = "#dcdcdc", land.border.col = "#000000") +
    geom_point(data = spatial.map.transform,
               aes(x = lon.proj, y = lat.proj, size = Number_Plots_Binned, fill = Region),
               alpha = 0.85, colour = "#000000", shape = 21) +
    scale_size_discrete(range = c(4,12,20,28,36)) +
    scale_fill_manual(values = c("#D02090", "#228B22", "#EE7600", "#1C86EE")) +
    geom_label_repel(data = subset(spatial.map.transform, SITE %in% c("Utqiagvik", "Thufuver", "Alexfiord", "Pyramiden")),
                     aes(lon.proj, lat.proj, label = SITE), color = "black", box.padding = 2,
                     segment.color = "black", segment.size = 0.7, fill = "white", label.size = 0.4,  size = 4.5) +
    labs(title = "ITEX Site Locations (Spatial)",
         subtitle = "",
         size = "Number of Plots") +
    theme_4_map())

# Export map
ggsave(output.spatial, file = "figures/itex_sites_spatial.png", height = 8, width = 8)


# MAP OF ITEX SITES - TEMPORAL ----

# Load final temporal data
temporal <- read.csv("data/AAA_temporal_input.csv")

# Modify data
temporal.map <- temporal %>% 
  mutate(SITE = ifelse(SITE %in% c("BARROW"), "Utqiagvik", SITE)) %>% 
  group_by(SITE) %>% 
  mutate(Number_Plots = length(unique(PLOT))) %>% 
  ungroup() %>% 
  dplyr::select(LAT, LONG, SITE, Region, Number_Plots) %>% 
  distinct(SITE, Number_Plots, .keep_all = TRUE) %>% 
  arrange(Number_Plots) %>% 
  mutate(Number_Plots_Binned = case_when(Number_Plots < 10 ~ "1 - 9",
                                         Number_Plots >= 10 & Number_Plots < 50 ~ "10 - 49",
                                         Number_Plots >= 50 & Number_Plots < 100 ~ "50 - 99",
                                         Number_Plots >= 100 & Number_Plots < 200 ~ "100 - 199",
                                         Number_Plots >= 200 ~ "200 +")) %>% 
  mutate(Number_Plots_Binned = factor(.$Number_Plots_Binned,
                                      levels = c("1 - 9", "10 - 49", "50 - 99", "100 - 199", "200 +"))) %>% # Reorder factor
  mutate(SITE = str_to_title(SITE),
         Region = case_when(Region == "North America-East" ~ "North America East",
                            Region == "North America-West" ~ "North America West",
                            Region == "GreenIceland" ~ "Greenland & Iceland",
                            Region == "Eurasia" ~ "Eurasia",
                            NA ~ Region))

# Transform coordinates to UTM
temporal.map.transform <- transform_coord(temporal.map, lon = "LONG", lat = "LAT", bind = TRUE)

# Map based on size of circles
(output.temporal <- basemap(limits = 55, grid.size = 0.05, grid.col = "#949494",
                           land.size = 0.05, land.col = "#dcdcdc", land.border.col = "#000000") +
    geom_point(data = temporal.map.transform,
               aes(x = lon.proj, y = lat.proj, size = Number_Plots_Binned, fill = Region),
               alpha = 0.85, colour = "#000000", shape = 21) +
    scale_size_discrete(range = c(4,12,20,28,36)) +
    scale_fill_manual(values = c("#D02090", "#228B22", "#EE7600", "#1C86EE")) +
    geom_label_repel(data = subset(temporal.map.transform, SITE %in% c("Utqiagvik", "Thufuver", "Alexfiord", "Abisko")),
                     aes(lon.proj, lat.proj, label = SITE), color = "black", box.padding = 2,
                     segment.color = "black", segment.size = 0.7, fill = "white", label.size = 0.4,  size = 4.5) +
    labs(title = "ITEX Site Locations (Temporal)",
         subtitle = "",
         size = "Number of Plots") +
    theme_4_map())

# Export map
ggsave(output.temporal, file = "figures/itex_sites_temporal.png", height = 8, width = 8)


