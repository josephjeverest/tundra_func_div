# 05 - Cleaning the TRY v5.0 database for use in trait gap-filling
# Joseph Everest
# November 2021, adapted February 2022


# LOAD PACKAGES ----

library(tidyverse)
library(Taxonstand)


# TRY - LOAD DATA: (try.1) ----

  # Load in the original (19 GB) TRY v5 .txt file [HARD DRIVE]
  # try.1 <- fread("E:/TRY/8400.txt", header = T, sep = "\t", dec = ".", quote = "", data.table = T, verbose = TRUE) # Windows
  # try.1 <- fread("/Volumes/Everest_3/TRY/8400.txt", header = T, sep = "\t", dec = ".", quote = "", data.table = T, verbose = TRUE) # Mac


  # Export as RData file for smaller storage [HARD DRIVE]
  # save(try.1, file = "E:/TRY/TRYv5.RData") # Windows
  # save(try.1, file = "/Volumes/Everest_3/TRY/TRYv5.RData") # Mac

# Load the smaller (~2 GB) TRY v5 .RData file [HARD DRIVE]
try.1 <- get(load("E:/TRY/TRYv5.RData")) # Windows
# try.1 <- get(load("/Volumes/Everest_3/TRY/TRYv5.RData")) # Mac



# TRY - LOCATION INFORMATION: (try.2) ----

  # Load in the TRY location data [HARD DRIVE]
  # try.location <- read.csv("E:/TRY/TRY_5_Site_Climate_Soil_2019-03-25.csv") # Windows
  # try.location <- read.csv("/Volumes/Everest_3/TRY/TRY_5_Site_Climate_Soil_2019-03-25.csv") # Mac


  # Export as RData file for smaller storage [TundraDivHub]
  # save(try.location, file = "data/output_05_TRY_locations.RData")

# Load in TRY locations data [TundraDivHub]
try.location <- get(load("data/output_05_TRY_locations.RData"))

# Trim to useful columns 
try.location.cut <- try.location %>% 
  dplyr::select(LAT_site, LON_site, observationId) %>% 
  rename(LAT = LAT_site,
         LONG = LON_site,
         ObservationID = observationId)

# Join location data to TRY database
try.2 <- left_join(try.1, try.location.cut, by = c("ObservationID" = "ObservationID")) %>% 
  relocate(LAT, LONG, .after = ObservationID)


# TRY - TAXONOMY CHECK: (try.3) ----

  # NOTE: as per Anne Bjorkman's earlier work, appears the 'AccSpeciesName' column in TRY is as per the The Plant List already
    # Running taxonomy checker on this column to make sure its up to date with the latest TaxonomyChecker version (2.3)

# Extract vector of unique species names from TRY
try.species <- unique(try.2$AccSpeciesName)

# Run the taxonomy checker on the TRY names
# try.tpl <- TPL(try.species)

  # Export the taxonomy checker results to RData file
  # save(try.tpl, file = "data/output_05_tpl_species_try.RData")

  # Import the taxonomy checker results for TRY from the RData file
  try.tpl <- get(load("data/output_05_tpl_species_try.RData"))

# Subset the TPL output to useful columns
try.tpl.cut <- try.tpl %>% 
  dplyr::select(Taxon, Family, New.Genus, New.Species, New.Taxonomic.status) %>% 
  mutate(NAME = paste(New.Genus, New.Species, sep = " ")) %>% # Generate new binomial name
  rename(GENUS = New.Genus,
         FAMILY = Family) %>% 
  relocate(NAME, GENUS, FAMILY, .before = New.Taxonomic.status) %>% # Relocate column
  dplyr::select(-New.Species)

# Counts of the number of Accepted, Unresolved and NA records there are
try.tpl.count <- try.tpl.cut %>% 
  group_by(New.Taxonomic.status) %>% 
  tally() %>% 
  ungroup()
    # vast majority of species now have accepted records after Taxonomy check

# Join new names to the TRY database
try.3 <- try.tpl.cut %>% 
  dplyr::select(Taxon, NAME, GENUS, FAMILY) %>% 
  left_join(try.2, ., by = c("AccSpeciesName" = "Taxon")) %>% 
  relocate(NAME, GENUS, FAMILY, .after = AccSpeciesName)


# TRY - MODIFY DATAFRAME (try.4) ----

# Add contributor column and remove superfluous columns to shrink dataset size
try.4 <- try.3 %>% 
  mutate(DataContributor = paste(FirstName, LastName, sep = " ")) %>% 
  dplyr::select(DatasetID, Dataset, DataContributor, AccSpeciesID, NAME, GENUS, FAMILY, ObservationID,
                LAT, LONG, TraitID, TraitName, DataID, OrigValueStr, OrigObsDataID, ValueKindName, UnitName, Reference, Comment) %>% 
  rename(TraitValue = OrigValueStr)

# Export as RData file for smaller storage [HARD DRIVE]
save(try.4, file = "E:/TRY/output_05_TRY_trimmed.RData") # Windows
# save(try.4, file = "/Volumes/Everest_3/TRY/output_05_TRY_trimmed.RData") # Mac

