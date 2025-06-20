# 03 - Permafrost raster data extraction
# Joseph Everest (with help from M. Garcia Criado and J. Assmann)
# March 2021, adapted October 2021, February 2022


# PACKAGES ----

# Load required packages
library(tidyverse)
library(raster)
library(rgdal)
library(rasterVis)
library(sp)
library(broom)


# DATA IMPORT ----

# Load the ITEX coordinates
itex <- read.csv("data/output_01_itex_coords.csv") %>% 
  dplyr::select(-X)

# Permafrost Map:
  # Permafrost map (100m resolution) (.tif) stored locally on a hard drive
  # Contact Joe Everest (joseph.everest@ed.ac.uk) for more information on access...
  # ... or visit Obu et al. (2019) at https://www.sciencedirect.com/science/article/pii/S0012825218305907

# Defining the filepath to the files
# folderpath.permafrost <- ("/Volumes/Everest_3/Permafrost/") # Hard drive on Mac
folderpath.permafrost <- ("D:/Permafrost/") # Hard drive on Windows
filename.permafrost <- list.files(folderpath.permafrost, pattern = "*.tif")
filepath.permafrost = paste0(folderpath.permafrost, filename.permafrost)


# EXTRACTION ----

# Create SpatialPoints (sp) object of unique coordinates
itex.coord <- SpatialPoints(itex)

# create raster object for permafrost tif
permafrost.raster <- stack(filepath.permafrost)

# Extract variables values for each pair of coordinates
permafrost.extract <- raster::extract(permafrost.raster, itex.coord, df = TRUE)


# COMBINED DATAFRAMES ----

# Convert the SpatialPoints (sp) object into a dataframe 
itex.coord.df <- as.data.frame(itex.coord)

# Reassign the 'ID' to the ITEX coordinates dataframe
itex.coord.df$ID <- row.names(itex.coord.df)
itex.coord.df$ID <- as.numeric(itex.coord.df$ID) # Make numeric

# Merge the two dataframes: extracted permafrost values and ITEX coordinates
coord.permafrost.combo <- left_join(permafrost.extract, itex.coord.df, by = c("ID" = "ID")) %>% 
  rename(permafrost = permafrost_obu.2018) %>% 
  dplyr::select(-ID) %>% 
  mutate(permafrost = case_when(permafrost == 4 ~ "continuous",
                                permafrost == 3 ~ "discontinuous",
                                permafrost == 2 ~ "sporadic",
                                permafrost == 1 ~ "none")) %>% 
  mutate(permafrost = ifelse(is.na(permafrost), "none", permafrost))


# EXPORT TO CSV ----

# Export dataframe to combine with ITEX data
write.csv(coord.permafrost.combo, "data/output_03_permafrost")
