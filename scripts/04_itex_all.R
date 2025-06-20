# 04 - Adding Climate and Permafrost Data to ITEX dataframe
# Joseph Everest
# March 2021, adapted October 2021, February 2022


# PACKAGES ----

# Load required packages
library(tidyverse)


# ITEX INPUT: Base data (itex.1) ----

# Import full, cleaned ITEX dataset
itex.1 <- read.csv("data/output_01_itex_full.csv") %>% 
  dplyr::select(-X)


# Export growth form list
itex.gf <- itex.1 %>% 
  dplyr::select(SPECIES, FuncGroup) %>% 
  distinct()

# Export to .csv
write.csv(itex.gf, file = "data/output_fg_export_itex.csv", row.names = FALSE)


# CHELSA INPUT: Climate data (itex.2) ----

# Read in extracted CHELSA variables for each unique coordinate pair (Script: 02)
chelsa <- read.csv("data/output_02_chelsa.csv") %>% 
  dplyr::select(-X)

# Create vector of column names (without lat and long) for reordering below
colnames_chelsa <- chelsa %>% 
  dplyr::select(-LAT, -LONG) %>% 
  colnames(.)

# Combine with main ITEX database
itex.2 <- left_join(itex.1, chelsa, by = c("LAT", "LONG")) %>% 
  relocate(all_of(colnames_chelsa), .after = MOISTURE)


# PERMAFROST INPUT: Permafrost data (itex.3) ----

# Read in extracted permafrost variable for each unique coordinate pair (Script: 03)
permafrost <- read.csv("data/output_03_permafrost.csv") %>% 
  dplyr::select(-X)

# Combine with main ITEX database
itex.3 <- left_join(itex.2, permafrost, by = c("LAT", "LONG")) %>% 
  relocate(permafrost, .after = SnowAnn)


# WRITE CSV: Output csv for combination with trait data ----

# Write CSV
write.csv(itex.3, "data/output_04_itex_all.csv")


# UNIQUE SPECIES: Output list of unique ITEX species (itex.4) ----

# Get unique species as dataframe
itex.4 <- data.frame(unique(itex.3$SPECIES)) %>% 
  rename(Species = unique.itex.3.SPECIES.) %>% 
  filter(!str_detect(Species, pattern = "XXX")) %>% 
  arrange(Species)

# Write CSV
write.csv(itex.4, "data/output_04_itex_final_species_list.csv")
