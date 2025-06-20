# 06 - Trait data manipulation and gap-filling
# Joseph Everest
# March 2021, adapted November 2021, February 2022, March 2022

  # NOTE: TRY data is taxon checked and trimmed to the relevant columns in script 05 to make it easier to work with


# LOAD PACKAGES ----

library(tidyverse)
library(ggrepel)
library(Taxonstand)
library(rworldmap) # For mapping
library(sp) # For mapping
library(raster) # For mapping
# devtools::install_github("eliocamp/ggalt@new-coord-proj")
library(ggalt) # For mapping
library(gridExtra)
library(rgeos)
library(factoextra)
library(ggOceanMaps)
library(viridis)

# CUSTOM FUNCTIONS ----

# Not included in...
`%notin%` <- Negate(`%in%`)

# Function to convert pair of coordinates into country names
coords2country = function(r.pts){
  # Takes single argument: data frame with two cols (col 1 = LONG in deg, col 2 = LAT in deg.)
  # Sourced from Phoebe Stewart-Sinclair at https://stackoverflow.com/questions/41105042/using-coords2country-function-in-r-on-exclusive-economic-zones-not-country-bound
  countriesSP <- getMap(resolution='high') # You could use high res map from rworldxtra if you were concerned about detail
  r.pts = SpatialPoints(r.pts, proj4string=CRS(proj4string(countriesSP))) # Setting CRS directly to that from rworldmap
  indices = over(r.pts, countriesSP) # Use 'over' to get indices of the Polygons object containing each point 
  indices$ADMIN    # Returns the ADMIN names of each country
  #indices$ISO3    # Returns the ISO3 code
  #indices$continent   # Returns the continent (6 continent model)
  # indices$REGION   # Returns the continent (7 continent model)
}

# Function to convert pair of coordinates into continent names
coords2continent = function(r.pts){
  # Takes single argument: data frame with two cols (col 1 = LONG in deg, col 2 = LAT in deg.)
  # Sourced from Phoebe Stewart-Sinclair at https://stackoverflow.com/questions/41105042/using-coords2country-function-in-r-on-exclusive-economic-zones-not-country-bound
  countriesSP <- getMap(resolution='high') # You could use high res map from rworldxtra if you were concerned about detail
  r.pts = SpatialPoints(r.pts, proj4string=CRS(proj4string(countriesSP))) # Setting CRS directly to that from rworldmap
  indices = over(r.pts, countriesSP) # Use 'over' to get indices of the Polygons object containing each point 
  # indices$ADMIN    # Returns the ADMIN names of each country
  # indices$ISO3    # Returns the ISO3 code
  # indices$continent   # Returns the continent (6 continent model)
  indices$REGION   # Returns the continent (7 continent model)
}

# THEMES ----

# Load ggplot themes from separate source script
source("scripts/00_ggplot_themes.R")


# TRY - LOAD IN REQUIRED DATA (try.1) ----

# Load the smaller TRY v5 .RData with locations,taxonomy checked and trimmed cols [HARD DRIVE]
try.1 <- get(load("E:/TRY/output_05_TRY_trimmed.RData")) # Windows
# try.1 <- get(load("/Volumes/Everest_3/TRY/output_05_TRY_trimmed.RData")) # Mac

# Remove unwanted second .RData import (try.4)
rm(try.4)


# TRY - FILTER FOR REQUIRED TRAITS, NON-EXPERIMENTAL & MEASUREMENT TYPE: (try.2) ----
  
# Create vectors of TraitIDs for each of the seven required traits using TRY website (https://www.try-db.org/de/TabDetails.php)
SLA <- c(3115) # Removed 125 (fresh mass), removed 3116 and 3117 (do/may include petioles)
Height <- c(3106, 3107)
SDM <- c(26)
LDMC <- c(47)
LeafN <- c(14) # Removed 660 (as some organic thing)
LeafP <- c(15)
LeafArea <- c(3108) # Removed 3114, 3112, 3113 (may include petioles), removed 3110, 3111 (petiole included)
  # 3108 - leaf, 3109 - petiole

# Filter for only the required traits and standardise names
try.seven.traits <- try.1 %>% 
  filter(TraitID %in% SLA | TraitID %in% Height | TraitID %in% SDM | # Filter for required traits by TraitID
         TraitID %in% LDMC | TraitID %in% LeafN | TraitID %in% LeafP |
         TraitID %in% LeafArea) %>% 
  mutate(TraitName_NEW = TraitName) %>% 
  mutate(TraitName_NEW = ifelse(TraitID %in% SLA, "SLA", TraitName_NEW), # Standardise names within traits
         TraitName_NEW = ifelse(TraitID %in% Height, "Plant_Height", TraitName_NEW),
         TraitName_NEW = ifelse(TraitID %in% SDM, "SDM", TraitName_NEW),
         TraitName_NEW = ifelse(TraitID %in% LDMC, "LDMC", TraitName_NEW),
         TraitName_NEW = ifelse(TraitID %in% LeafN, "LeafN", TraitName_NEW),
         TraitName_NEW = ifelse(TraitID %in% LeafP, "LeafP", TraitName_NEW),
         TraitName_NEW = ifelse(TraitID %in% LeafArea, "Leaf_Area", TraitName_NEW)) %>% 
  relocate(TraitName_NEW, .after = TraitName)

# Remove traits recorded from plants under experimental conditions, if known (DataID == 327)
try.non.experiment <- try.seven.traits %>% 
  filter(!DataID == 327) # No experimental records remaining

# See what types of measurement are present in data
try.value.types <- try.non.experiment %>% 
  dplyr::select(ValueKindName) %>% 
  unique()

# Create vector of measurement types to retain
value.types <- c("Single", "Mean", "Best estimate", "Median", "Site specific mean", "Species mean")

# Retain only records that are mean or individual (remove max and min etc.)
try.2 <- try.non.experiment %>% 
  filter(ValueKindName %in% value.types)

# # Add temporary export here for corrections
# write.csv(try.2, "data/output_06_TRY_initial_cut.csv")
# 
# # Read in temporary export
# try.2 <- read.csv("data/output_06_TRY_initial_cut.csv") %>% dplyr::select(-X)


# TRY - FILTER RECORDS: (try.3) ----

# Remove NA values
try.na <- try.2 %>% 
  filter(!is.na(LAT),
         !is.na(TraitValue))

# Remove any records not sourced from 60 degN or above
try.3 <- try.na %>% 
  filter(LAT >= 60)


# TRY - REMOVE DUPLICATES: (try.4) ----

  # NOTE: ObservationID is unique to one site, no repetition of ObservationIDs between sites

# First remove any duplicates as identified by TRY (OrigObsDataID)
try.dupl.cut <- try.3 %>% 
  filter(is.na(OrigObsDataID)) # Retain records that aren't identified as having an OrigObsDataID

# Convert TraitValue to numeric
try.numeric <- try.dupl.cut %>% 
  mutate(TraitValue = trimws(TraitValue)) %>% # Removes leading or trailing spaces
  mutate(TraitValue = as.numeric(TraitValue))

# Make columns of lat, long and trait value to 6.d.p. (as per A. Bjorkman) for identifying duplicates
try.round <- try.numeric %>%
  mutate(LAT_r = round(LAT, digits = 6),
         LONG_r = round(LONG, digits = 6),
         TraitValue_r = round(TraitValue, digits = 6))

# First remove duplicates where same trait (new name) measured multiple times on same individual observation (e.g. leaf/leaflet or height veg/reprod.)
try.duplicates.1 <- try.round %>% 
  group_by(ObservationID, NAME, TraitName_NEW, TraitValue_r, LAT_r, LONG_r) %>% # Using new trait name, rather than old as measuring same thing on same species
  mutate(Duplicate = n()>1) %>% 
  ungroup() %>% # Removes 167 clear duplicate records
  distinct(., ObservationID, NAME, TraitName_NEW, TraitValue_r, LAT_r, LONG_r, .keep_all = TRUE) %>% 
  dplyr::select(-Duplicate)

# Make a column for less obvious 'duplicates' where the ObservationID is different
try.duplicates.2 <- try.duplicates.1 %>% 
  group_by(NAME, TraitName_NEW, TraitValue_r, LAT_r, LONG_r) %>%
  mutate(Duplicate = n()>1) %>% 
  ungroup()

# Assign unique identifier to each duplicate (based on name lat, long, trait and trait value)
try.duplicates.3 <- try.duplicates.2 %>% 
  mutate(Duplicate_ID = ifelse(Duplicate == "TRUE", paste(NAME, TraitName_NEW, TraitValue_r, LAT_r, LONG_r,  sep = ":"), NA))

# Create dataframe for the records that are not duplicates (all records to be carried forwards)
try.non.duplicates <- try.duplicates.3 %>% 
  filter(Duplicate %in% c("FALSE"))

# Create dataframe for the records that are identified as duplicates
try.duplicates.all <- try.duplicates.3 %>% 
  filter(Duplicate %in% c("TRUE"))
    # 52 records here, suggesting 26 duplicates

# From duplicates dataframe, create a tally of no. of datasets in which each duplicate occurs
try.duplicates.all.sets <- try.duplicates.all %>% 
  group_by(Duplicate_ID) %>% 
  mutate(Duplicate_SETS = length(unique(Dataset))) %>% 
  ungroup()

# Check the number of datasets in which duplicates occur
unique(try.duplicates.all.sets$Duplicate_SETS)
  # No record is duplicated across three or more datasets, only across one or two datasets

# Filter out any records where only duplicate within one site (875 records)
try.duplicates.sets.one <- try.duplicates.all.sets %>% 
  filter(Duplicate_SETS == 1) %>% 
  group_by(Duplicate_ID) %>% 
  mutate(Occurrences = length(TraitValue_r)) %>% # No. of occurrences of each duplicate
  ungroup() %>% 
  arrange(Duplicate_ID) # Order by Duplicate_ID to group the duplicates together

# Now check that each 'duplicate' record within the dataset is for a different ObservationID
  # These records shouldn't actually be duplicates, just same values for different plants measured at the same location
try.duplicates.sets.one.id <- try.duplicates.sets.one %>% 
  group_by(Duplicate_ID) %>% 
  mutate(Duplicate_ObsID = length(unique(ObservationID))) %>% 
  ungroup() %>% 
  mutate(Match = ifelse(Occurrences == Duplicate_ObsID, "YES", "NO")) %>% 
  filter(Match %in% c("NO"))
    # No records, e.g. every duplicate from within one site ONLY is actually measured on all different individuals (ObservationID)
      # Hence, no actual duplicates so no need to filter anything from try.duplicates.sets.one

# Remove intermediary calculation columns from try.duplicates.sets.one
try.duplicates.sets.one <- dplyr::select(try.duplicates.sets.one, -Duplicate_SETS, -Occurrences)

# Now filter records where duplicates occur across two datasets (likely to be actual duplicates)
try.duplicates.sets.two <- try.duplicates.all.sets %>% 
  filter(Duplicate_SETS == 2) %>% 
  group_by(Duplicate_ID) %>% 
  mutate(Occurrences = length(TraitValue_r)) %>% # No. of occurrences of each duplicate
  ungroup() %>% 
  arrange(Duplicate_ID) # Order by Duplicate_ID to group the duplicates together

# See how many occurrences of each duplicate that occurs in different datasets there are
  # IMPORTANT: need to do this to make sure there aren't any records duplicated within a dataset AND between datasets
  # Then we will lose information we don't want to lose in the filter
unique(try.duplicates.sets.two$Occurrences)
  # No duplicate occurs more than twice so NO records duplicated within a dataset AND between datasets

# Therefore, can filter the real duplicates that occur in different datasets
try.duplicates.sets.two <- try.duplicates.sets.two %>% 
  distinct(., Duplicate_ID, .keep_all = TRUE) %>% # Remove identical records between sites
  dplyr::select(-Duplicate_SETS, -Occurrences) # Remove temporary columns for identifying duplicates

# Now join together all three dataframes:
  # try.non.duplicates - majority of records, never identified as duplicates
  # try.duplicates.sets.one - no records removed, not actually duplicates
  # try.duplicates.sets.two - records removed, actual duplicates in different datasets

try.4 <- rbind(try.non.duplicates, try.duplicates.sets.one, try.duplicates.sets.two) %>% 
  dplyr::select(-LAT_r, -LONG_r, -TraitValue_r, -Duplicate, -Duplicate_ID) %>% # Remove additional intermediary columns
  arrange(NAME)


# TRY - FINAL TAXONOMY CORRECTIONS (try.5) ----

# Generate vector of unique species retained in TRY after filtering and duplicate removal
species.try.4 <- unique(try.4$NAME)

# # Run taxonomy checker on these species
# species.try.checked <- TPL(species.try.4)
# 
# # Write csv of taxonomy checker results
# write.csv(species.try.checked, "data/output_06_tpl_species_try.csv")

# Read csv of taxonomy checker results
species.try.checked <- read.csv("data/output_06_tpl_species_try.csv") %>% 
  dplyr::select(-X)

# Create a trimmed dataframe of the checked species
species.try.checked.2 <- species.try.checked %>%
  dplyr::select(Taxon, New.Genus, New.Species, New.Taxonomic.status, Family) %>%
  mutate(Name_TPL = paste(New.Genus, New.Species, sep = " ")) %>%
  relocate(Name_TPL, New.Genus, Family, .before = New.Taxonomic.status) %>%
  dplyr::select(-New.Species) %>% 
  rename(New.Species = Name_TPL, New.Family = Family) %>% 
  arrange(Taxon)

# Records to amend (Blank and NA status): accept the unresolved records as should be consistent throughout TPL checks
species.try.checked.3 <- species.try.checked.2 %>% 
  filter(New.Taxonomic.status %in% c("", NA))

# Manually fix the broken records in the TPL dataframe
species.try.checked.4 <- species.try.checked.2 %>% 
  mutate(New.Species = ifelse(Taxon %in% c("Arctagrostis NA"), "XXXArctagrostis", New.Species), # No species name
         New.Genus = ifelse(Taxon %in% c("Arctagrostis NA"), "Arctagrostis", New.Genus),
         New.Family = ifelse(Taxon %in% c("Arctagrostis NA"), "Poaceae", New.Family),
         New.Taxonomic.status = ifelse(Taxon %in% c("Arctagrostis NA"), "Fixed", New.Taxonomic.status),
         New.Species = ifelse(Taxon %in% c("Arctostaphylos NA"), "XXXArctostaphylos", New.Species), # No species name
         New.Genus = ifelse(Taxon %in% c("Arctostaphylos NA"), "Arctostaphylos", New.Genus),
         New.Family = ifelse(Taxon %in% c("Arctostaphylos NA"), "Ericaceae", New.Family),
         New.Taxonomic.status = ifelse(Taxon %in% c("Arctostaphylos NA"), "Fixed", New.Taxonomic.status),
         New.Species = ifelse(Taxon %in% c("Betula NA"), "XXXBetula", New.Species), # No species name
         New.Genus = ifelse(Taxon %in% c("Betula NA"), "Betula", New.Genus),
         New.Family = ifelse(Taxon %in% c("Betula NA"), "Betulaceae", New.Family),
         New.Taxonomic.status = ifelse(Taxon %in% c("Betula NA"), "Fixed", New.Taxonomic.status),
         New.Species = ifelse(Taxon %in% c("Carex NA"), "XXXCarex", New.Species), # No species name
         New.Genus = ifelse(Taxon %in% c("Carex NA"), "Carex", New.Genus),
         New.Family = ifelse(Taxon %in% c("Carex NA"), "Cyperaceae", New.Family),
         New.Taxonomic.status = ifelse(Taxon %in% c("Carex NA"), "Fixed", New.Taxonomic.status),
         New.Species = ifelse(Taxon %in% c("Dryas NA"), "XXXDryas", New.Species), # No species name
         New.Genus = ifelse(Taxon %in% c("Dryas NA"), "Dryas", New.Genus),
         New.Family = ifelse(Taxon %in% c("Dryas NA"), "Rosaceae", New.Family),
         New.Taxonomic.status = ifelse(Taxon %in% c("Dryas NA"), "Fixed", New.Taxonomic.status),
         New.Species = ifelse(Taxon %in% c("Eri NA", "Eri sch"), "XXXUnknown", New.Species), # Unknown species/genus/family
         New.Genus = ifelse(Taxon %in% c("Eri NA", "Eri sch"), "Unknown", New.Genus),
         New.Family = ifelse(Taxon %in% c("Eri NA", "Eri sch"), "Unknown", New.Family),
         New.Taxonomic.status = ifelse(Taxon %in% c("Eri NA", "Eri sch"), "Fixed", New.Taxonomic.status),
         New.Species = ifelse(Taxon %in% c("Eria NA"), "XXXEria", New.Species), # No species name
         New.Genus = ifelse(Taxon %in% c("Eria NA"), "Eria", New.Genus),
         New.Family = ifelse(Taxon %in% c("Eria NA"), "Orchidaceae", New.Family),
         New.Taxonomic.status = ifelse(Taxon %in% c("Eria NA"), "Fixed", New.Taxonomic.status),
         New.Species = ifelse(Taxon %in% c("Festuca NA"), "XXXFestuca", New.Species), # No species name
         New.Genus = ifelse(Taxon %in% c("Festuca NA"), "Festuca", New.Genus),
         New.Family = ifelse(Taxon %in% c("Festuca NA"), "Poaceae", New.Family),
         New.Taxonomic.status = ifelse(Taxon %in% c("Festuca NA"), "Fixed", New.Taxonomic.status),
         New.Species = ifelse(Taxon %in% c("Gentiaglauca NA"), "Gentiana glauca", New.Species), # Incorrect species name
         New.Genus = ifelse(Taxon %in% c("Gentiaglauca NA"), "Gentiana", New.Genus),
         New.Family = ifelse(Taxon %in% c("Gentiaglauca NA"), "Gentianaceae", New.Family),
         New.Taxonomic.status = ifelse(Taxon %in% c("Gentiaglauca NA"), "Fixed", New.Taxonomic.status),
         New.Species = ifelse(Taxon %in% c("Kobresia NA"), "XXXKobresia", New.Species), # No species name
         New.Genus = ifelse(Taxon %in% c("Kobresia NA"), "Kobresia", New.Genus),
         New.Family = ifelse(Taxon %in% c("Kobresia NA"), "Cyperaceae", New.Family),
         New.Taxonomic.status = ifelse(Taxon %in% c("Kobresia NA"), "Fixed", New.Taxonomic.status),
         New.Species = ifelse(Taxon %in% c("Loiselura procumbens"), "Loiseleuria procumbens", New.Species), # Incorrect species name
         New.Genus = ifelse(Taxon %in% c("Loiselura procumbens"), "Loiseleuria", New.Genus),
         New.Family = ifelse(Taxon %in% c("Loiselura procumbens"), "Ericaceae", New.Family),
         New.Taxonomic.status = ifelse(Taxon %in% c("Loiselura procumbens"), "Fixed", New.Taxonomic.status),
         New.Species = ifelse(Taxon %in% c("Pedicularis NA"), "XXXPedicularis", New.Species), # No species name
         New.Genus = ifelse(Taxon %in% c("Pedicularis NA"), "Pedicularis", New.Genus),
         New.Family = ifelse(Taxon %in% c("Pedicularis NA"), "Orobanchaceae", New.Family),
         New.Taxonomic.status = ifelse(Taxon %in% c("Pedicularis NA"), "Fixed", New.Taxonomic.status),
         New.Species = ifelse(Taxon %in% c("Picea NA"), "XXXPicea", New.Species), # No species name
         New.Genus = ifelse(Taxon %in% c("Picea NA"), "Picea", New.Genus),
         New.Family = ifelse(Taxon %in% c("Picea NA"), "Pinaceae", New.Family),
         New.Taxonomic.status = ifelse(Taxon %in% c("Picea NA"), "Fixed", New.Taxonomic.status),
         New.Species = ifelse(Taxon %in% c("Poa NA"), "XXXPoa", New.Species), # No species name
         New.Genus = ifelse(Taxon %in% c("Poa NA"), "Poa", New.Genus),
         New.Family = ifelse(Taxon %in% c("Poa NA"), "Poaceae", New.Family),
         New.Taxonomic.status = ifelse(Taxon %in% c("Poa NA"), "Fixed", New.Taxonomic.status),
         New.Species = ifelse(Taxon %in% c("Rhododendron laponica"), "Rhododendron lapponicum", New.Species), # Incorrect species name
         New.Genus = ifelse(Taxon %in% c("Rhododendron laponica"), "Rhododendron", New.Genus),
         New.Family = ifelse(Taxon %in% c("Rhododendron laponica"), "Ericaceae", New.Family),
         New.Taxonomic.status = ifelse(Taxon %in% c("Rhododendron laponica"), "Fixed", New.Taxonomic.status),
         New.Species = ifelse(Taxon %in% c("Ribes NA"), "XXXRibes", New.Species), # No species name
         New.Genus = ifelse(Taxon %in% c("Ribes NA"), "Ribes", New.Genus),
         New.Family = ifelse(Taxon %in% c("Ribes NA"), "Grossulariaceae", New.Family),
         New.Taxonomic.status = ifelse(Taxon %in% c("Ribes NA"), "Fixed", New.Taxonomic.status),
         New.Species = ifelse(Taxon %in% c("Salix NA"), "XXXSalix", New.Species), # No species name
         New.Genus = ifelse(Taxon %in% c("Salix NA"), "Salix", New.Genus),
         New.Family = ifelse(Taxon %in% c("Salix NA"), "Salicaceae", New.Family),
         New.Taxonomic.status = ifelse(Taxon %in% c("Salix NA"), "Fixed", New.Taxonomic.status),
         New.Species = ifelse(Taxon %in% c("Vaccinium NA"), "XXXVaccinium", New.Species), # No species name
         New.Genus = ifelse(Taxon %in% c("Vaccinium NA"), "Vaccinium", New.Genus),
         New.Family = ifelse(Taxon %in% c("Vaccinium NA"), "Ericaceae", New.Family),
         New.Taxonomic.status = ifelse(Taxon %in% c("Vaccinium NA"), "Fixed", New.Taxonomic.status))

# Rejoin new species information to full ITEX dataset
try.5 <- left_join(try.4, species.try.checked.4, by = c("NAME" = "Taxon")) %>% 
  relocate(New.Species, New.Genus, New.Family, .after = FAMILY) %>% 
  dplyr::select(-c(GENUS, FAMILY)) %>%
  rename(ORIGINAL_NAME = NAME,
         NAME = New.Species,
         GENUS = New.Genus,
         FAMILY = New.Family)


# TRY - UNIT EXTRACTION ----

# Check what units each trait has
only_SLA_try <- filter(try.5, TraitName_NEW == "SLA")
units_SLA_try <- unique(only_SLA_try$UnitName)

only_Height_try <- filter(try.5, TraitName_NEW == "Plant_Height")
units_Height_try <- unique(only_Height_try$UnitName)

only_SDM_try <- filter(try.5, TraitName_NEW == "SDM")
units_SDM_try <- unique(only_SDM_try$UnitName)

only_LDMC_try <- filter(try.5, TraitName_NEW == "LDMC")
units_LDMC_try <- unique(only_LDMC_try$UnitName)

only_LeafN_try <- filter(try.5, TraitName_NEW == "LeafN")
units_LeafN_try <- unique(only_LeafN_try$UnitName)

only_LeafP_try <- filter(try.5, TraitName_NEW == "LeafP")
units_LeafP_try <- unique(only_LeafP_try$UnitName)

only_LeafArea_try <- filter(try.5, TraitName_NEW == "Leaf_Area")
units_LeafArea_try <- unique(only_LeafArea_try$UnitName)


# TRY - EXPORT CLEANED FILE ----

# Export .csv of cleaned TRY data
write.csv(try.5, "data/output_06_TRY_clean.csv")

# Remove unwanted intermediary variables
rm(try.dupl.cut, try.duplicates.1, try.duplicates.2, try.duplicates.3, try.duplicates.all,
   try.duplicates.all.sets, try.duplicates.sets.one, try.duplicates.sets.one.id,
   try.duplicates.sets.two, try.na, try.non.duplicates, try.non.experiment, try.numeric,
   try.round, try.seven.traits, try.value.types, species.try.checked, species.try.checked.2,
   species.try.checked.3, species.try.checked.4, only_SLA_try, only_Height_try, only_SDM_try,
   only_LDMC_try, only_LeafN_try, only_LeafP_try, only_LeafArea_try)


# TTT - LOAD DATA: (ttt.1) ----

# NOTE: as per Anne Bjorkman's earlier work, appears the 'AccSpeciesName' column in TTT is as per the The Plant List already
  # Running taxonomy checker later to make sure its up to date with the latest TaxonomyChecker version (2.3)

# Load in cleaned TTT database
ttt.1 <- read.csv("data/input_TTT_full_cleaned.csv", stringsAsFactors = FALSE) %>% 
  dplyr::select(-X)


# TTT - FILTER FOR REQUIRED TRAITS & MEASUREMENT TYPE: (ttt.2) ----

# Remove superfluous columns to shrink dataset size
ttt.cut <- ttt.1 %>% 
  dplyr::select(-DayOfYear) %>% 
  rename(TraitValue = Value)


ttt.seven.traits <- ttt.cut %>% 
  rename(TraitName = Trait) %>% 
  filter(TraitName %in% c("Leaf dry mass per leaf fresh mass (Leaf dry matter content, LDMC)",
                          "Leaf area",
                          "Leaf area per leaf dry mass (specific leaf area, SLA)",
                          "Leaf phosphorus (P) content per leaf dry mass",
                          "Leaf nitrogen (N) content per leaf dry mass",
                          "Plant height, vegetative",
                          "Plant height, reproductive",
                          "Seed dry mass")) %>%
  mutate(TraitName_NEW = TraitName) %>% 
  mutate(TraitName_NEW = ifelse(TraitName %in% c("Leaf dry mass per leaf fresh mass (Leaf dry matter content, LDMC)"), "LDMC", TraitName_NEW),
         TraitName_NEW = ifelse(TraitName %in% c("Leaf area"), "Leaf_Area", TraitName_NEW),
         TraitName_NEW = ifelse(TraitName %in% c("Leaf area per leaf dry mass (specific leaf area, SLA)"), "SLA", TraitName_NEW),
         TraitName_NEW = ifelse(TraitName %in% c("Leaf phosphorus (P) content per leaf dry mass"), "LeafP", TraitName_NEW),
         TraitName_NEW = ifelse(TraitName %in% c("Leaf nitrogen (N) content per leaf dry mass"), "LeafN", TraitName_NEW),
         TraitName_NEW = ifelse(TraitName %in% c("Plant height, vegetative", "Plant height, reproductive"), "Plant_Height", TraitName_NEW),
         TraitName_NEW = ifelse(TraitName %in% c("Seed dry mass"), "SDM", TraitName_NEW)) %>% 
  relocate(TraitName_NEW, .after = TraitName)

# Remove unwanted measurement methods
ttt.2 <- ttt.seven.traits %>% 
  filter(!ValueKindName %in% c("Maximum in plot")) # Only want individuals or means


# TTT - FILTER BY LATITUDE: (ttt.3) ----

# Identify sites with missing coordinates
ttt.coords.NA <- ttt.2 %>% 
  filter(is.na(Latitude)) %>% 
  dplyr::select(SiteName) %>% 
  unique()
    # Four sites have NAs - Eastern Alps, Svalbard, Urals Picos de Urals - but all in specified NH range
      # Want to retain these sites

# Retain all records over 30 degN and the E.Alps/Svalbard/Urals records missing LAT information
ttt.3 <- ttt.2 %>% # Assign missing LAT values for three sites as 999 to retain in > 30 filter [BELOW]
  mutate(Latitude = ifelse(SiteName %in% c("Eastern Alps", "Svalbard", "Urals", "Picos de Europa, Spain") & is.na(Latitude), 999, Latitude)) %>% 
  filter(Latitude > 0) %>% # Filter for records in NH
  mutate(Latitude = ifelse(Latitude == 999, NA, Latitude)) %>% # Change modified latitudes back to NA
  mutate(SiteName = ifelse(is.na(SiteName), "UNKNOWN", SiteName)) # 49 records have no SiteName, rename as UNKNOWN
    # Removes SH (Marion Islands, Australian and NZ) sites


# TTT - REMOVE DUPLICATES: (ttt.4) ----

  # NOTE: IndividualID is unique to individual sites, no repetition of Individual IDs between sites

# Make columns of lat, long and trait value to 6.d.p. (as per A. Bjorkman) for identifying duplicates
ttt.round <- ttt.3 %>%
  mutate(LAT_r = round(Latitude, digits = 6),
         LONG_r = round(Longitude, digits = 6),
         TraitValue_r = round(TraitValue, digits = 6))

# First remove duplicates where same trait (new name) measured multiple times on same individual observation (e.g. leaf/leaflet or height veg/reprod.)
ttt.duplicates.1 <- ttt.round %>% 
  group_by(IndividualID, AccSpeciesName, TraitName_NEW, TraitValue_r, LAT_r, LONG_r, Year) %>% # Using new trait name, rather than old as measuring same thing on same species
  mutate(Duplicate = n()>1) %>%
  ungroup() %>% # [BELOW] Removes one instance of an identical record
  distinct(., IndividualID, AccSpeciesName, TraitName_NEW, TraitValue_r, LAT_r, LONG_r, Year, .keep_all = TRUE) %>%
  dplyr::select(-Duplicate)
    # Removes 188 records (all where both types of plant height measured on the same individual)

# Make a column to check for less obvious potential 'duplicates' where the IndividualID is different
ttt.duplicates.2 <- ttt.duplicates.1 %>% 
  group_by(AccSpeciesName, TraitName_NEW, TraitValue_r, LAT_r, LONG_r, Year) %>%
  mutate(Duplicate = n()>1) %>% 
  ungroup()

# Assign unique identifier to each duplicate (based on name, year, lat, long, trait and trait value)
ttt.duplicates.3 <- ttt.duplicates.2 %>% 
  mutate(Duplicate_ID = ifelse(Duplicate == "TRUE", paste(AccSpeciesName, Year, TraitName_NEW, TraitValue_r, LAT_r, LONG_r,  sep = ":"), NA))

# Create dataframe for the records that are not duplicates (all records to be carried forwards)
ttt.non.duplicates <- ttt.duplicates.3 %>% 
  filter(Duplicate %in% c("FALSE"))

# Create dataframe for the records that are identified as duplicates
ttt.duplicates.all <- ttt.duplicates.3 %>% 
  filter(Duplicate %in% c("TRUE"))
    # ~12,000 records here, suggesting up to ~6,000 duplicates (seems unlikely...)

# From duplicates dataframe, create a tally of no. of datasets in which each duplicate occurs
ttt.duplicates.all.sets <- ttt.duplicates.all %>% 
  group_by(Duplicate_ID) %>% 
  mutate(Duplicate_SETS = length(unique(SiteName))) %>% 
  ungroup()

# Check the number of datasets in which duplicates occur
unique(ttt.duplicates.all.sets$Duplicate_SETS)
  # Records only duplicated within one site, suggesting not duplicates, but actually different individuals (CHECK BELOW)

# Mutate column to check the number of occurrences of each duplicate#
ttt.duplicates.sets.one <- ttt.duplicates.all.sets %>% 
  group_by(Duplicate_ID) %>% 
  mutate(Occurrences = length(TraitValue_r)) %>% # No. of occurrences of each duplicate
  ungroup() %>% 
  arrange(Duplicate_ID) # Order by Duplicate_ID to group the duplicates together

# Now check that each 'duplicate' record within the dataset is for a different ObservationID
  # These records shouldn't actually be duplicates, just same values for different plants measured at the same location
ttt.duplicates.sets.one.id <- ttt.duplicates.sets.one %>% 
  group_by(Duplicate_ID) %>% 
  mutate(Duplicate_ObsID = length(unique(IndividualID))) %>% 
  ungroup() %>% 
  mutate(Match = ifelse(Occurrences == Duplicate_ObsID, "YES", "NO")) %>% 
  filter(Match %in% c("NO"))
    # No records, e.g. every duplicate from within one site ONLY is actually measured on all different individuals (ObservationID)
      # Hence, no actual duplicates so no need to filter anything from try.duplicates.sets.one

# Remove intermediary calculation columns from try.duplicates.sets.one
ttt.duplicates.sets.one <- dplyr::select(ttt.duplicates.sets.one, -Duplicate_SETS, -Occurrences)

# Now join together both dataframes:
# ttt.non.duplicates - majority of records, never identified as duplicates
# try.duplicates.sets.one - no records removed, not actually duplicates

ttt.4 <- rbind(ttt.non.duplicates, ttt.duplicates.sets.one) %>% 
  dplyr::select(-LAT_r, -LONG_r, -TraitValue_r, -Duplicate, -Duplicate_ID) %>% # Remove additional intermediary columns
  arrange(AccSpeciesName)


# TTT - FINAL TAXONOMY CORRECTIONS (ttt.5) ----

# Generate vector of unique species retained in TRY after filtering and duplicate removal
species.ttt.4 <- unique(ttt.4$AccSpeciesName)

# Run taxonomy checker on these species
# species.ttt.checked <- TPL(species.ttt.4)

# Write csv of taxonomy checker results
# write.csv(species.ttt.checked, "data/output_06_tpl_species_ttt.csv")

# Read csv of taxonomy checker results
species.ttt.checked <- read.csv("data/output_06_tpl_species_ttt.csv")

# Create a trimmed dataframe of the checked species
species.ttt.checked.2 <- species.ttt.checked %>%
  dplyr::select(Taxon, New.Genus, New.Species, New.Taxonomic.status, Family) %>%
  mutate(Name_TPL = paste(New.Genus, New.Species, sep = " ")) %>%
  relocate(Name_TPL, New.Genus, Family, .before = New.Taxonomic.status) %>%
  dplyr::select(-New.Species) %>% 
  rename(New.Species = Name_TPL, New.Family = Family) %>% 
  arrange(Taxon)

# Records to amend (Blank and NA status): accept the unresolved records as should be consistent throughout TPL checks
species.ttt.checked.3 <- species.ttt.checked.2 %>% 
  filter(New.Taxonomic.status %in% c("", NA))

# Manually fix the broken records in the TPL dataframe
species.ttt.checked.4 <- species.ttt.checked.2 %>% 
  mutate(New.Species = ifelse(Taxon %in% c("Alchemilla sp", "Alchemilla sp."), "XXXAlchemilla", New.Species), # No species name
         New.Genus = ifelse(Taxon %in% c("Alchemilla sp", "Alchemilla sp."), "Alchemilla", New.Genus),
         New.Family = ifelse(Taxon %in% c("Alchemilla sp", "Alchemilla sp."), "Rosaceae", New.Family),
         New.Taxonomic.status = ifelse(Taxon %in% c("Alchemilla sp", "Alchemilla sp."), "Fixed", New.Taxonomic.status),
         New.Species = ifelse(Taxon %in% c("Anthyllis sp."), "XXXAnthyllis", New.Species), # No species name
         New.Genus = ifelse(Taxon %in% c("Anthyllis sp."), "Anthyllis", New.Genus),
         New.Family = ifelse(Taxon %in% c("Anthyllis sp."), "Fabaceae", New.Family),
         New.Taxonomic.status = ifelse(Taxon %in% c("Anthyllis sp."), "Fixed", New.Taxonomic.status),
         New.Species = ifelse(Taxon %in% c("Carex sp."), "XXXCarex", New.Species), # No species name
         New.Genus = ifelse(Taxon %in% c("Carex sp."), "Carex", New.Genus),
         New.Family = ifelse(Taxon %in% c("Carex sp."), "Cyperaceae", New.Family),
         New.Taxonomic.status = ifelse(Taxon %in% c("Carex sp."), "Fixed", New.Taxonomic.status),
         New.Species = ifelse(Taxon %in% c("Draba sedenense", "Draba sp."), "XXXDraba", New.Species), # No species name
         New.Genus = ifelse(Taxon %in% c("Draba sedenense", "Draba sp."), "Draba", New.Genus),
         New.Family = ifelse(Taxon %in% c("Draba sedenense", "Draba sp."), "Brassicaceae", New.Family),
         New.Taxonomic.status = ifelse(Taxon %in% c("Draba sedenense", "Draba sp."), "Fixed", New.Taxonomic.status),
         New.Species = ifelse(Taxon %in% c("Dupontia sp."), "XXXDupontia", New.Species), # No species name
         New.Genus = ifelse(Taxon %in% c("Dupontia sp."), "Dupontia", New.Genus),
         New.Family = ifelse(Taxon %in% c("Dupontia sp."), "Poaceae", New.Family),
         New.Taxonomic.status = ifelse(Taxon %in% c("Dupontia sp."), "Fixed", New.Taxonomic.status),
         New.Species = ifelse(Taxon %in% c("Festuca sp."), "XXXFestuca", New.Species), # No species name
         New.Genus = ifelse(Taxon %in% c("Festuca sp."), "Festuca", New.Genus),
         New.Family = ifelse(Taxon %in% c("Festuca sp."), "Poaceae", New.Family),
         New.Taxonomic.status = ifelse(Taxon %in% c("Festuca sp."), "Fixed", New.Taxonomic.status),
         New.Species = ifelse(Taxon %in% c("Graminoid sp."), "XXXGraminoid", New.Species), # Growth form
         New.Genus = ifelse(Taxon %in% c("Graminoid sp."), NA, New.Genus),
         New.Family = ifelse(Taxon %in% c("Graminoid sp."), NA, New.Family),
         New.Taxonomic.status = ifelse(Taxon %in% c("Graminoid sp."), "Fixed", New.Taxonomic.status),
         New.Species = ifelse(Taxon %in% c("Pedicularis sp."), "XXXPedicularis", New.Species), # No species name
         New.Genus = ifelse(Taxon %in% c("Pedicularis sp."), "Pedicularis", New.Genus),
         New.Family = ifelse(Taxon %in% c("Pedicularis sp."), "Orobanchaceae", New.Family),
         New.Taxonomic.status = ifelse(Taxon %in% c("Pedicularis sp."), "Fixed", New.Taxonomic.status),
         New.Species = ifelse(Taxon %in% c("Poa sp."), "XXXPoa", New.Species), # No species name
         New.Genus = ifelse(Taxon %in% c("Poa sp."), "Poa", New.Genus),
         New.Family = ifelse(Taxon %in% c("Poa sp."), "Poaceae", New.Family),
         New.Taxonomic.status = ifelse(Taxon %in% c("Poa sp."), "Fixed", New.Taxonomic.status),
         New.Species = ifelse(Taxon %in% c("Puccinelia sp."), "XXXPuccinellia", New.Species), # No species name; misspelt genus
         New.Genus = ifelse(Taxon %in% c("Puccinelia sp."), "Puccinellia", New.Genus),
         New.Family = ifelse(Taxon %in% c("Puccinelia sp."), "Poaceae", New.Family),
         New.Taxonomic.status = ifelse(Taxon %in% c("Puccinelia sp."), "Fixed", New.Taxonomic.status),
         New.Species = ifelse(Taxon %in% c("Salix glauca-niphoclada", "Salix herbacea-polaris", "Salix sp."), "XXXSalix", New.Species), # ...
         New.Genus = ifelse(Taxon %in% c("Salix glauca-niphoclada", "Salix herbacea-polaris", "Salix sp."), "Salix", New.Genus),
         New.Family = ifelse(Taxon %in% c("Salix glauca-niphoclada", "Salix herbacea-polaris", "Salix sp."), "Salicaceae", New.Family),
         New.Taxonomic.status = ifelse(Taxon %in% c("Salix glauca-niphoclada", "Salix herbacea-polaris", "Salix sp."), "Fixed", New.Taxonomic.status),
         New.Species = ifelse(Taxon %in% c("Taraxacum sp", "Taraxacum sp."), "XXXTaraxacum", New.Species), # ...
         New.Genus = ifelse(Taxon %in% c("Taraxacum sp", "Taraxacum sp."), "Taraxacum", New.Genus),
         New.Family = ifelse(Taxon %in% c("Taraxacum sp", "Taraxacum sp."), "Compositae", New.Family),
         New.Taxonomic.status = ifelse(Taxon %in% c("Taraxacum sp", "Taraxacum sp."), "Fixed", New.Taxonomic.status))

# Rejoin new species information to full ITEX dataset
ttt.5 <- left_join(ttt.4, species.ttt.checked.4, by = c("AccSpeciesName" = "Taxon")) %>% 
  relocate(New.Species, New.Genus, New.Family, New.Taxonomic.status, .after = OriginalName) %>% 
  dplyr::select(-c(AccSpeciesName)) %>%
  rename(NAME = New.Species,
         GENUS = New.Genus,
         FAMILY = New.Family)


# TTT - UNIT EXTRACTION ----

# Check what units each trait has
only_SLA_ttt <- filter(ttt.5, TraitName_NEW == "SLA")
units_SLA_ttt <- unique(only_SLA_ttt$Units)

only_Height_ttt <- filter(ttt.5, TraitName_NEW == "Plant_Height")
units_Height_ttt <- unique(only_Height_ttt$Units)

only_SDM_ttt <- filter(ttt.5, TraitName_NEW == "SDM")
units_SDM_ttt <- unique(only_SDM_ttt$Units)

only_LDMC_ttt <- filter(ttt.5, TraitName_NEW == "LDMC")
units_LDMC_ttt <- unique(only_LDMC_ttt$Units)

only_LeafN_ttt <- filter(ttt.5, TraitName_NEW == "LeafN")
units_LeafN_ttt <- unique(only_LeafN_ttt$Units)

only_LeafP_ttt <- filter(ttt.5, TraitName_NEW == "LeafP")
units_LeafP_ttt <- unique(only_LeafP_ttt$Units)

only_LeafArea_ttt <- filter(ttt.5, TraitName_NEW == "Leaf_Area")
units_LeafArea_ttt <- unique(only_LeafArea_ttt$Units)


# TTT - EXPORT CLEANED FILE ----

# Export .csv of cleaned TTT data
write.csv(ttt.5, file = "data/output_06_TTT_clean.csv")

# Remove unwanted intermediary variables
rm(ttt.coords.NA, ttt.cut, ttt.duplicates.1, ttt.duplicates.2, ttt.duplicates.3,
   ttt.duplicates.all, ttt.duplicates.all.sets, ttt.duplicates.sets.one, ttt.duplicates.sets.one.id,
   ttt.non.duplicates, ttt.round, ttt.seven.traits, species.ttt.checked, species.ttt.checked.2,
   species.ttt.checked.3, species.ttt.checked.4, only_SLA_ttt, only_Height_ttt, only_SDM_ttt,
   only_LDMC_ttt, only_LeafN_ttt, only_LeafP_ttt, only_LeafArea_ttt)


# COMBO - COMBINE THE TRY AND TTT DATA (combo.1) ----

    # If starting from here, import cleaned TRY and TTT dataframes
    # try.5 <- read.csv("data/output_06_TRY_clean.csv")
    # ttt.5 <- read.csv("data/output_06_TTT_clean.csv")

# Cut TRY data so standard subset of columns
try.x <- try.5 %>% 
  dplyr::select(NAME, GENUS, FAMILY, ObservationID, Dataset, DataContributor,
                LAT, LONG, TraitName_NEW, TraitName, TraitValue, ValueKindName) %>% 
  mutate(ObservationID = paste("TRY", ObservationID, sep = ":")) %>% # Add TRY to ID incase same number in TTT
  mutate(TRYorTTT = "TRY") %>% # Add column for what dataset
  relocate(TRYorTTT, .before = ObservationID) %>% 
  rename(Dataset_Sitename = Dataset,
         Singlular_ID = ObservationID,
         TraitName_FULL = TraitName,
         TraitName = TraitName_NEW) %>% 
  mutate(Height_Type = NA) %>% # Add in column for height type
  mutate(Height_Type = ifelse(str_detect(TraitName_FULL, pattern = "reproductive"), "Reproductive", Height_Type),
         Height_Type = ifelse(str_detect(TraitName_FULL, pattern = "generative"), "Reproductive", Height_Type),
         Height_Type = ifelse(str_detect(TraitName_FULL, pattern = "vegetative"), "Vegetative", Height_Type))

# Cut TTT data so standard subset of columns
ttt.x <- ttt.5 %>% 
  dplyr::select(NAME, GENUS, FAMILY, IndividualID, SiteName, DataContributor,
                Latitude, Longitude, TraitName_NEW, TraitName, TraitValue, ValueKindName) %>% 
  mutate(IndividualID = paste("TTT", IndividualID, sep = ":")) %>% # Add TTT to ID incase same number in TRY
  mutate(TRYorTTT = "TTT") %>% # Add column for what dataset
  mutate(GENUS = as.character(GENUS),
         FAMILY = as.character(FAMILY)) %>% 
  relocate(TRYorTTT, .before = IndividualID) %>% 
  rename(Dataset_Sitename = SiteName,
         Singlular_ID = IndividualID,
         TraitName_FULL = TraitName,
         TraitName = TraitName_NEW,
         LAT = Latitude,
         LONG = Longitude) %>% 
  mutate(Height_Type = NA) %>% # Add in column for height type
  mutate(Height_Type = ifelse(str_detect(TraitName_FULL, pattern = "reproductive"), "Reproductive", Height_Type),
         Height_Type = ifelse(str_detect(TraitName_FULL, pattern = "generative"), "Reproductive", Height_Type),
         Height_Type = ifelse(str_detect(TraitName_FULL, pattern = "vegetative"), "Vegetative", Height_Type))

# Join the two dataframes
combo.1 <- rbind(try.x, ttt.x)


# COMBO - CHECK FOR DUPLICATES BETWEEN TRY AND TTT (combo.2) ----

# Make columns of lat, long and trait value to 6.d.p. (as per A. Bjorkman) for identifying duplicates
combo.round <- combo.1 %>%
  mutate(LAT_r = round(LAT, digits = 6),
         LONG_r = round(LONG, digits = 6),
         TraitValue_r = round(TraitValue, digits = 6))

# Check to see if any duplicate records
  # This will flag records in same dataset we have already identified as NOT duplicates
  # Nothing in the same dataset is a duplicate
combo.duplicates.1 <- combo.round %>% 
  group_by(NAME, TraitName, TraitValue_r, LAT_r, LONG_r) %>%
  mutate(Duplicate = n()>1) %>% 
  ungroup()

# Assign unique identifier to each duplicate (based on dataset, name, year, lat, long, trait and trait value)
combo.duplicates.2 <- combo.duplicates.1 %>% 
  mutate(Duplicate_ID = ifelse(Duplicate == "TRUE", paste(NAME, TraitName, TraitValue_r, LAT_r, LONG_r,  sep = ":"), NA))

# Create dataframe for the records that are not duplicates (all records to be carried forwards)
combo.non.duplicates <- combo.duplicates.2 %>% 
  filter(Duplicate %in% c("FALSE"))

# Create dataframe for the records that are identified as duplicates
combo.duplicates.all <- combo.duplicates.2 %>% 
  filter(Duplicate %in% c("TRUE"))
    # ~13,000 records here, suggesting up to ~7,500 duplicates (seems unlikely, expecting zero...)

# From duplicates dataframe, create a tally of no. of datasets in which each duplicate occurs
  # If only occurs in one dataset, can discount from being a duplicate as already cleaned these
  # Only looking for those that occur in two datasets (e.g. TRY and TTT)
combo.duplicates.all.sets <- combo.duplicates.all %>% 
  group_by(Duplicate_ID) %>% 
  mutate(Duplicate_SETS = length(unique(TRYorTTT))) %>% 
  ungroup() %>% 
  filter(Duplicate_SETS > 1)
    # No records occur in both datasets, therefore no duplicates between the two

# **As no duplicates identified, can carry forwards combo.1 as combo.2**
combo.2 <- combo.1


# COMBO - FURTHER CLEANING AND MANUAL CHECKING OF TRAIT VALUES (combo.3) ----

 # NOTE: cleaning methods here employed as per A. Bjorkman and M. Garcia Criado
  # Adapted from TraitHub/scripts/TTT_cleaning_script.R and Anne's "Combine_TRY_ITEX_data_UPDATED.R" script

# Remove records identified by Anne as needing removal
combo.3 <- combo.2 %>% 
  # Remove Niwot mean data (because individual values are added from TTT)
  mutate(REMOVE = ifelse(str_detect(Dataset_Sitename, pattern = "Niwot") & str_detect(ValueKindName, pattern = "ean"), "YES", "NO")) %>% 
  # Remove ALL Niwot SLA data because it's sketchy (apparently)
  mutate(REMOVE = ifelse(str_detect(Dataset_Sitename, pattern = "Niwot") & TraitName %in% c("SLA"), "YES", REMOVE)) %>% 
  # Remove Salix arctica SDM data from Rebecca Klady because units uncertain
  mutate(REMOVE = ifelse(str_detect(DataContributor, pattern = "Rebecca Klady") & NAME %in% c("Salix arctica") & TraitName %in% c("SDM"), "YES", REMOVE)) %>% 
  # Remove Papaver radicatum SDM data from Anne Bjorkman because values seem off (potential units issue)
  mutate(REMOVE = ifelse(str_detect(DataContributor, pattern = "Anne Bjorkman") & NAME %in% c("Papaver radicatum") & TraitName %in% c("SDM"), "YES", REMOVE)) %>% 
  filter(REMOVE %in% c("NO")) %>% # Retain only recorded labelled REMOVE == "NO"
  dplyr::select(-REMOVE) # Remove superfluous column


# COMBO - MANUAL CHECKS OF ERRONEOUS TRAIT VALUES (combo.4) ----

# Create function to generate histogram per species per trait to manually eyeball erroneous values
combo.error.function <- function(species, traits){
  
  # Create loop for generating and saving histogram
  for (j in traits){
    
    for (i in species){
      
      # Create dataframe for that species and trait
      combo.error.df <- input.dataframe %>% 
        filter(NAME == i, TraitName == j)
      
      # Determine number of trait values for that species/trait combination
      combo.error.measurements <- length(combo.error.df$TraitValue)
      
      # Create an if statement dependent on whether any measurements for that species/trait combo
      if (combo.error.measurements > 0){
        
        # Plot histogram
        combo.error.histogram <- ggplot(combo.error.df) +
          geom_histogram(aes(x = TraitValue)) +
          labs(title = i,
               subtitle = j) +
          theme_bw() +
          theme(plot.title = element_text(size = 6),
                plot.subtitle = element_text(size = 6),
                axis.title = element_text(size = 6))
        
        # Create filepath for saving histogram
        combo.error.filepath = paste0("figures/", output.folder, "/", j, "_", i, ".png")
        
        
        # Save histogram
        ggsave(combo.error.histogram, file = combo.error.filepath, width = 2, height = 2)
        
      } else {
        
        # Remove variables
        rm(combo.error.species.name, combo.error.trait.name, combo.error.df, combo.error.measurements)
        
      }
    }
  }
  
  # Remove all variables at end of loop
  rm(combo.error.histogram, combo.error.filepath,combo.error.species.name,
     combo.error.trait.name, combo.error.df, combo.error.measurements)
  
}

# Set parameters to run function by (before filtering)
input.dataframe <- combo.3 # Full unfiltered data frame
output.folder <- "trait_error_checking_pre_filtering" # Folder to output histograms to
combo.error.species <- unique(combo.3$NAME) # Vector of unique species in the trait data
combo.error.traits <- unique(combo.3$TraitName) # Vector of unique traits in the trait data

# # Delete files previously in folder
# unlink(paste0("figures/", output.folder, "/*"), recursive = T, force = T)
# 
# # Run function to generate histograms
# combo.error.function(combo.error.species, combo.error.traits)

    # NOTE: See script EX4_trait_error_checking for details of trait investigations

# Manually remove erroneous values from the combined trait data frame
combo.manual <- combo.3 %>%
  mutate(Remove = "No") %>% # Create column from which to remove values
  # LDMC checks
  mutate(Remove = ifelse(TraitName == "LDMC" & TraitValue > 1, "Yes", Remove)) %>% # Values over 1 not possible
  # LeafN checks
  mutate(Remove = ifelse(TraitName == "LeafN" & Dataset_Sitename == "Abisko, Sweden", "Yes", Remove)) %>% # Lots of erroneous values in dataset
  # LeafP checks
  mutate(Remove = ifelse(TraitName == "LeafP" & Dataset_Sitename == "Reich-Oleksyn Global Leaf N, P Database", "Yes", Remove), # Lots of erroneous values in dataset
         Remove = ifelse(TraitName == "LeafP" & Dataset_Sitename == "Abisko, Sweden", "Yes", Remove)) %>% # Lots of erroneous values in dataset
  # Plant height checks
  mutate(Remove = ifelse(TraitName == "Plant_Height" & Dataset_Sitename == "European North Russia", "Yes", Remove), # Lots of erroneous values in dataset
         Remove = ifelse(TraitName == "Plant_Height" & Dataset_Sitename == "The Xylem/Phloem Database", "Yes", Remove), # Lots of erroneous values in dataset
         Remove = ifelse(TraitName == "Plant_Height" & Dataset_Sitename == "The VISTA Plant Trait Database", "Yes", Remove), # Lots of erroneous values in dataset
         Remove = ifelse(TraitName == "Plant_Height" & Dataset_Sitename == "Siberian shrub allometry", "Yes", Remove)) %>% # Lots of erroneous values in dataset
  # SDM checks
  mutate(Remove = ifelse(TraitName == "SDM" & Dataset_Sitename == "Abisko, Sweden", "Yes", Remove)) %>% # Lots of erroneous values in dataset
  # SLA checks
  mutate(Remove = ifelse(TraitName == "SLA" & Dataset_Sitename == "The Global Leaf Traits", "Yes", Remove), # Lots of erroneous values in dataset
         Remove = ifelse(TraitName == "SLA" & Dataset_Sitename == "GLOPNET - Global Plant Trait Network Database", "Yes", Remove)) %>%  # Lots of erroneous values in dataset
  # Leaf area checks
      # SEE K-MEANS CLUSTERING BELOW
  
  # REMOVE ALL VALUES TO BE REMOVED
  filter(Remove == "No") %>% 
  dplyr::select(-Remove) %>% 
  
  # Standard deviations removal
  group_by(NAME, TraitName) %>% 
  mutate(SD = sd(TraitValue), # s.d. per species per trait
         Mean = mean(TraitValue)) %>% # mean per species per trait
  ungroup() %>% 
  mutate(UpperErrorBound = Mean + 2 * SD, # Values within 2 s.d. of the mean for that species / trait combo
         LowerErrorBound = Mean - 2 * SD) %>% 
  mutate(Remove = "No") %>% 
  mutate(Remove = ifelse(TraitValue <= LowerErrorBound | TraitValue >= UpperErrorBound, "Yes", Remove)) %>% 
  filter(Remove == "No") %>% 
  dplyr::select(-c(Remove, SD, Mean, UpperErrorBound, LowerErrorBound))

# Set parameters to run function by (after filtering)
input.dataframe <- combo.manual # Full filtered data frame  
output.folder <- "trait_error_checking_post_filtering" # Folder to output histograms to
combo.error.species <- unique(combo.manual$NAME) # Vector of unique species in the trait data
combo.error.traits <- unique(combo.manual$TraitName) # Vector of unique traits in the trait data

# # Delete files previously in folder
# unlink(paste0("figures/", output.folder, "/*"), recursive = T, force = T)
# 
# # Run function to generate histograms
# combo.error.function(combo.error.species, combo.error.traits)

# Create list of species with leaf area data
leaf.area.species <- sort(unique(filter(combo.manual, TraitName == "Leaf_Area")$NAME))

# Create blank data frame for outputting values to
combo.leaf.area <- data.frame()


# # Loop to run k-means clustering on each of the leaf_area vs species combos
# for (i in leaf.area.species) {
# 
#   # Filter full data frame to only leaf area and species name
#   leaf.area <- combo.manual %>%
#     filter(TraitName == "Leaf_Area", NAME == i) %>%
#     arrange(TraitValue)
# 
#   # Determine number of observations
#   leaf.area.observations <- as.numeric(nrow(leaf.area))
# 
# 
#   # If only two or less observations, can't run clustering at all
#   if (leaf.area.observations > 2) {
# 
#     # If has three or less observations, can't run automated cluster number detection as: max number of clusters can't be more than the number of obs but MUST be at least 4 (2 more than 2)
#     if (leaf.area.observations > 3) {
# 
#       # if (leaf.area.observations <= 15) {
#       #
#       #
#       #   # Run function to decide ideal number of clusters for each species
#       #   leaf.area.clusters <- NbClust::NbClust(leaf.area$TraitValue, distance = "euclidean", min.nc = 2,
#       #                                          max.nc = leaf_area_observations, method = "kmeans", index = "all")
#       #
#       # } else {
#       #
#       # # Run function to decide ideal number of clusters for each species
#       # leaf.area.clusters <- NbClust::NbClust(leaf.area$TraitValue, distance = "euclidean", min.nc = 2,
#       #                                        max.nc = 15, method = "kmeans", index = "all")
#       #
#       # }
#       #
#       # # Extract the correct number of clusters
#       # number.of.clusters <- as.numeric(length(unique(leaf.area.clusters$Best.partition)))
#       #
#       # # Run the clustering - automated selection of number of clusters
#       # k_means <- kmeans(leaf.area$TraitValue, number.of.clusters, nstart = 25)
# 
#       # Run the clustering - with just two clusters
#       k.means <- kmeans(leaf.area$TraitValue, 2, nstart = 25)
# 
#     } else {
# 
#       # Run the clustering - with just two clusters
#       k.means <- kmeans(leaf.area$TraitValue, 2, nstart = 25)
# 
#     }
# 
# 
#     # Extract the clusters
#     k.means.clusters <- k.means$cluster
# 
#     # Join the clusters back to the original data
#     leaf.area.output <- leaf.area %>%
#       mutate(cluster = as.character(k.means.clusters))
# 
#     # Plot histogram to show data values with the clusters
#     leaf.area.histogram <-
#       ggplot(leaf.area.output, aes(x = TraitValue, fill = cluster)) +
#       geom_histogram() +
#       theme_bw() +
#       theme(plot.title = element_text(size = 6),
#             plot.subtitle = element_text(size = 6),
#             axis.title = element_text(size = 6))
# 
#     # Determine filepath
#     leaf.area.filepath <-
#       paste0("figures/trait_error_checking_clustering/",
#              str_replace_all(as.character(i), pattern = " ", replacement = "_"), ".png")
# 
#     # Output histograms
#     ggsave(leaf.area.histogram, file = leaf.area.filepath)
# 
# 
#   } else {
# 
#     # Add NAs to cluster column
#     leaf.area.output <- mutate(leaf.area, cluster = NA)
# 
#   }
# 
#   # Bind outputs back to main data frame
#   combo.leaf.area <- rbind(combo.leaf.area, leaf.area.output)
# 
#   # Remove intermediate variables
#   rm(leaf.area, leaf.area.observations, leaf.area.clusters, leaf.area.histogram,
#      leaf.area.output,  k.means)
# 
# }
# 
# # Save output to prevent having to run loop each time
# write.csv(combo.leaf.area, file = "data/output_06_leaf_area_clusters.csv")

# Import the csv file to prevent having to run loop each time
combo.leaf.area <- read.csv("data/output_06_leaf_area_clusters.csv") %>% dplyr::select(-X) %>% 
  mutate(Height_Type = NA) # Manually add back in Height_Type column


# Define what percentage spread of data you want to encompass in the below cleaning by clusters
combo.leaf.area.percent <- 0.5


# If the spread between the two means is over x% of the range of whole data, considered distinct, retain only cluster of smaller average mean
combo.leaf.area.cut <- combo.leaf.area %>% 
  arrange(NAME, cluster) %>% 
  group_by(NAME) %>% 
  mutate(min = min(TraitValue), # Work out minimum and maximum trait values
         max = max(TraitValue)) %>% 
  ungroup() %>%
  mutate(range = max - min, # Work out spread of trait values
         range_percent = range * combo.leaf.area.percent) %>% # Work out 50/60/70/80% etc. of that spread
  group_by(NAME, cluster) %>% 
  mutate(mean_cluster = mean(TraitValue)) %>% # Work out the mean of each cluster
  ungroup() %>% 
  group_by(NAME) %>% 
  mutate(range_cluster = max(mean_cluster) - min(mean_cluster), # The spread of the clustered data
         larger_mean_cluster = max(mean_cluster)) %>% 
  ungroup() %>% 
  mutate(over_percent = ifelse(range_cluster >= range_percent, "YES", "NO"), # Is the spread of the cluster means bigger than x% of the spread of the overall data
         bigger_cluster = ifelse(larger_mean_cluster == mean_cluster, "YES", "NO")) %>% 
  mutate(REMOVE = ifelse(over_percent == "YES" & bigger_cluster == "YES", "YES", "NO")) %>% 
  filter(REMOVE == "NO") %>% 
  dplyr::select(-c(min, max, range, range_percent, mean_cluster, range_cluster,
                   larger_mean_cluster, over_percent, bigger_cluster))
         

# # Now run another loop to create histograms of the filtered leaf area data
# for (i in leaf.area.species){
# 
#   # Filter to required species (and replace cluster NA as 0)
#   leaf.area.cut <- combo.leaf.area.cut %>%
#     filter(NAME == i) %>%
#     mutate(cluster = ifelse(is.na(cluster), 0, cluster))
# 
#   # Create histograms
#   leaf.area.cut.histogram <- ggplot(leaf.area.cut, aes(x = TraitValue, fill = cluster)) +
#     geom_histogram() +
#     theme_bw() +
#     theme(plot.title = element_text(size = 6),
#           plot.subtitle = element_text(size = 6),
#           axis.title = element_text(size = 6))
# 
#   # Determine filepath
#   leaf.area.cut.filepath <- paste0("figures/trait_error_checking_cut/",
#                                    str_replace_all(as.character(i), pattern = " ", replacement = "_"), "_CUT.png")
# 
#   # Output histograms
#   ggsave(leaf.area.cut.histogram, file = leaf.area.cut.filepath)
# 
#   # Remove intermediate objects
#   rm(leaf.area.cut, leaf.area.cut.histogram)
# 
# }

# Rejoin the leaf area data back to the main data frame
combo.4 <- rbind(filter(combo.manual, !TraitName == "Leaf_Area"),
                 dplyr::select(combo.leaf.area.cut, -c(cluster, REMOVE)))


# COMBO - PREPARE GEOGRAPHICAL INFORMATION FOR GAP-FILLING (combo.5) ----
  
# Remove unwanted columns and NA trait values from trait database
combo.cut <- combo.4 %>% 
  dplyr::select(NAME, GENUS, FAMILY, TraitName, TraitValue, Dataset_Sitename, LAT, LONG, Height_Type) %>% 
  filter(!is.na(TraitValue), !is.na(LAT), !is.na(LONG)) # Remove NAs (6 obs.) from TraitValue 

# Create a dataframe of coordinate pairs with dataset_sitename attached
coords.combo <- combo.cut %>% 
  dplyr::select(LONG, LAT) %>% 
  distinct(LONG, LAT, .keep_all = TRUE)

# Get country names for each pair of points
coords.country <- coords2country(coords.combo)

# Get continent names for each pair of points
coords.continent <- coords2continent(coords.combo)

# Join country and continent names to coordinate pairs dataframe
coords.c.c <- data.frame(coords.combo, coords.country, coords.continent) %>% 
  rename(Country = coords.country, 
         Continent = coords.continent) %>% 
  mutate(Country = as.character(Country),
         Continent = as.character(Continent)) %>% 
  mutate(LAT_r = round(LAT, digits = 1), # For manually fixing NAs
         LONG_r = round(LONG, digits = 1))

# Check NAs for missing country names (when point falls into see just off coast)
coords.c.c.na <- coords.c.c %>% 
  filter(is.na(Country)) %>% 
  distinct(LAT_r, LONG_r, .keep_all = TRUE)

# Fill in the missing country and continent data based on the rounded LAT and LONG coords
  # Some fall in middle of ocean, marked to REMOVE
coords.c.c.fix <- coords.c.c %>% 
  mutate(Country = case_when(LAT_r == 74.8 & LONG_r == 6.5 ~ "REMOVE",
                             LAT_r == 61.8 & LONG_r == 20.7 ~ "REMOVE",
                             LAT_r == 64.9 & LONG_r == 22.8 ~ "REMOVE",
                             LAT_r == 78.5 & LONG_r == 11.8 ~ "Norway",
                             LAT_r == 72.0 & LONG_r == -120.0 ~ "Canada",
                             LAT_r == 80.0 & LONG_r == 53.0 ~ "Russia",
                             LAT_r == 65.0 & LONG_r == 40.0 ~ "Russia",
                             LAT_r == 69.6 & LONG_r == -138.9 ~ "Canada",
                             LAT_r == 56.5 & LONG_r == -76.5 ~ "Canada",
                             LAT_r == 69.6 & LONG_r == -139.0 ~ "Canada",
                             LAT_r == 65.6 & LONG_r == -37.7 ~ "Greenland",
                             LAT_r == 62.1 & LONG_r == -74.7 ~ "Canada",
                             LAT_r == 78.6 & LONG_r == 16.4 ~ "Norway",
                             LAT_r == 56.6 & LONG_r == -76.5 ~ "Canada",
                             LAT_r == 56.5 & LONG_r == -61.7 ~ "Canada",
                             LAT_r == 56.5 & LONG_r == -76.5 ~ "Canada",
                             LAT_r == 74.3 & LONG_r == -19.9 ~ "Greenland",
                             LAT_r == 69.2 & LONG_r == -51.1 ~ "Greenland",
                             LAT_r == 69.5 & LONG_r == -53.8 ~ "Greenland",
                             LAT_r == 69.5 & LONG_r == -54.0 ~ "Greenland",
                             TRUE ~ Country)) %>%
  mutate(Continent = ifelse(Country %in% c("Norway", "Russia"), "Europe", Continent),
         Continent = ifelse(Country == "Canada", "North America", Continent),
         Continent = ifelse(Country == "Greenland", "Greenland", Continent),
         Continent = ifelse(Country == "REMOVE", "REMOVE", Continent)) %>% 
  dplyr::select(-LAT_r, -LONG_r)

# Create spatial points object of the unique coordinate pairs
lgm.coords <- coords.c.c.fix %>% 
  dplyr::select(LONG, LAT) %>% 
  SpatialPoints(.)

# Check the projection of the coordinates layer
projection(lgm.coords)

# Assign # WGS84 projection
projection(lgm.coords) <- CRS("+init=epsg:4326")

# Recheck the projection of the coordinates layer
projection(lgm.coords)

# Load in the shapefile of the last LGM
lgm <- shapefile("data/lgm/lgm.shp")

# Check the projection of the shapefile layer
projection(lgm) # WGS84

# Extract the ice-covered or not extent for each LAT-LONG pair
lgm.extract <- raster::extract(lgm, lgm.coords, df = TRUE) %>% 
  rename(glaciated = ID)

# Check that duplicate records for certain LAT-LONG pairs are due to overlapping polygons and give same result
lgm.dupl.count <- lgm.extract %>% 
  group_by(point.ID) %>% 
  summarise(records = length(unique(glaciated))) %>% 
  ungroup()

# Check that all pairs only have one value for ice covered or not
unique(lgm.dupl.count$records) # CHECK

# All have one record for each point, so can remove duplicate ID records
lgm.extract.clean <- lgm.extract %>% 
  distinct(point.ID, .keep_all = TRUE)

# Convert the SpatialPoints (sp) object into a dataframe 
lgm.coords.df <- as.data.frame(lgm.coords)

# Reassign the 'ID' to the coordinates dataframe
lgm.coords.df$ID <- row.names(lgm.coords.df)
lgm.coords.df$ID <- as.numeric(lgm.coords.df$ID) # Make numeric

# Merge the two dataframes: trait coordinates and whether glaciated in LGM
lgm.combo <- left_join(lgm.coords.df, lgm.extract.clean, by = c("ID" = "point.ID")) %>% 
  dplyr::select(-c(ID, poly.ID, FEATURE, COMMENT)) %>% 
  mutate(glaciated = ifelse(glaciated == 0, "Yes", glaciated)) %>% 
  mutate(glaciated = ifelse(is.na(glaciated), "No", glaciated))

# Join information on whether glaciated or now to full coordinates
coords.c.c.lgm <- left_join(coords.c.c.fix, lgm.combo, by = c("LAT" = "LAT", "LONG" = "LONG"))

# Count how many locations were and weren't glaciated by continent
lgm.count <- coords.c.c.lgm %>% 
  group_by(Continent, glaciated) %>% 
  summarise(count = length(LAT)) %>% 
  ungroup() %>% 
  filter(Continent != "REMOVE") %>%
  mutate(Continent = ifelse(Continent == "Europe", "Eurasia", Continent)) %>% 
  add_row(Continent = "Greenland", glaciated = "No", count = 0) %>% # No non-glaciated records from Greenland
  arrange(Continent)

# Create 'Region' column based on the continent and whether glaciated or not
coords.region <- coords.c.c.lgm %>% 
  mutate(Region = case_when(Continent %in% c("Europe", "Asia") & Country != "Iceland" ~ "Eurasia",
                            Continent == "Greenland" | Country == "Iceland" ~ "GreenIceland",
                            Continent == "North America" & glaciated == "Yes" ~ "North America-East",
                            Continent == "North America" & glaciated == "No" ~ "North America-West",
                            Continent == "REMOVE" ~ "REMOVE")) %>% 
  mutate(Region_2 = Region) %>% # Second tier for gap-filling by North America alone
  mutate(Region_2 = ifelse(str_detect(Region, pattern = "North America"), "North America_2", Region_2),
         Region_2 = ifelse(Region %in% c("Eurasia"), "Eurasia_2", Region_2),
         Region_2 = ifelse(Region %in% c("GreenIceland"), "GreenIceland_2", Region_2)) %>% 
  dplyr::select(- c(Country, Continent, glaciated)) # No longer need country, continent or glaciated records

# Join region information to the overall combined traits dataframe
combo.geographic <- left_join(combo.cut, coords.region, by = c("LAT" = "LAT", "LONG" = "LONG")) %>% 
  filter(!Region %in% c("REMOVE")) # Remove records identified as being located in oceans
  
# Check all records have region information
unique(combo.geographic$Region) # No NAs

# Generate Gridcell, Region and Latitudinal Band information
combo.5 <- combo.geographic %>% 
  mutate(LAT_grid = plyr::round_any(LAT, 0.5, f = floor),
         LONG_grid = ifelse(LONG > 0, plyr::round_any(LONG, 0.5, f = floor), plyr::round_any(LONG, 0.5, f = ceiling)),
         LAT_band = plyr::round_any(LAT, 5, f = floor),
         over60degN = ifelse(LAT >= 60, "above_60", "below_60")) %>% 
  mutate(gridcell = paste0("_", LAT_grid, "_", LONG_grid)) %>% 
  mutate(Region_gridcell = paste0(Region, "_", gridcell),
         Region_LAT_band = paste0(Region, "_", LAT_band)) %>% 
  mutate(RecordID = row_number()) %>% # Temporarily adding in RecordID for visually checking trait gap-filling procedure
  relocate(RecordID, .before = NAME)
  
# Gap-filling levels:
  # 1) Species within gridcell
  # 2) Species within LAT_band of Region
  # 3) Species within Region > 60 degN
  # 4) Species within Region
  # 5) Species within Region_2 (combines N.America)
  # 6) Species OVERALL
  # 7-12) REPEAT for Genus-level
  # 13-18) REPEAT for Family-level


# COMBO - MODIFY PLANT HEIGHT VALUES (combo.6) ----

# Work out how many records each for vegetative and reproductive heights
combo.heights <- combo.5 %>%
  filter(TraitName == "Plant_Height") %>% 
  group_by(NAME) %>% 
  mutate(Height_Types = length(unique(Height_Type))) %>% 
  ungroup() %>% 
  group_by(NAME, Height_Type) %>% 
  mutate(Record_Sums = length(TraitValue)) %>% 
  ungroup() %>% 
  dplyr::select(NAME, Height_Type, Height_Types, Record_Sums) %>% 
  distinct(NAME, Height_Type, .keep_all = TRUE) %>% 
  arrange(NAME) %>% 
  pivot_wider(names_from = Height_Type, values_from = Record_Sums) %>% 
  mutate(Vegetative = ifelse(is.na(Vegetative), 0, Vegetative),
         Reproductive = ifelse(is.na(Reproductive), 0, Reproductive))

# Determine whether to retain reproductive records or not
combo.heights.retain <- combo.heights %>% 
  mutate(Keep_Vegetative = NA, Keep_Reproductive = NA) %>% 
  # If there is four or more vegetative, just keep vegetative
  mutate(Keep_Vegetative = ifelse(Vegetative >= 4, 1, Keep_Vegetative),
         Keep_Reproductive = ifelse(Vegetative >= 4, 0, Keep_Reproductive)) %>% 
  # If there is less than four vegetative but four or more reproductive, keep reproductive
  mutate(Keep_Vegetative = ifelse(Vegetative < 4 & Reproductive >= 4, 0, Keep_Vegetative),
         Keep_Reproductive = ifelse(Vegetative < 4 & Reproductive >= 4, 1, Keep_Reproductive)) %>%
  # If there is less than four vegetative and less than four reproductive, keep neither
  mutate(Keep_Vegetative = ifelse(Vegetative < 4 & Reproductive < 4, 0, Keep_Vegetative),
         Keep_Reproductive = ifelse(Vegetative < 4 & Reproductive < 4, 0, Keep_Reproductive)) %>% 
  # Calculate total to check never more than one record being kept
  mutate(TOTAL = Keep_Vegetative + Keep_Reproductive) %>% 
  # Pivot longer to make key
  dplyr::select(NAME, Keep_Vegetative, Keep_Reproductive) %>% 
  pivot_longer(names_to = "Height_Type", values_to = "RETAIN", cols = c(Keep_Vegetative, Keep_Reproductive)) %>%
  # Create key to join back to combo data for what height records to retain
  filter(RETAIN == 1) %>% 
  mutate(RETAIN = ifelse(RETAIN == 1, "YES", "NO"),
         Height_Type = ifelse(str_detect(Height_Type, pattern = "Vegetative"), "Vegetative", "Reproductive"))

# Rejoin key to combo dataset and retain only required height records
combo.6 <- left_join(combo.5, combo.heights.retain, by = c("NAME" = "NAME", "Height_Type" = "Height_Type")) %>% 
  # Complete the retain column
  mutate(RETAIN = case_when(TraitName != "Plant_Height" ~ "YES",
                            TraitName == "Plant_Height" & is.na(RETAIN) ~ "NO",
                            TRUE ~ RETAIN)) %>% 
  # Filter un-required height records
  filter(RETAIN == "YES") %>% 
  # Retain only columns required in the gap-filled dataset
  dplyr::select(RecordID, NAME, GENUS, FAMILY, Region_gridcell, Region_LAT_band, over60degN, Region, Region_2, TraitName, TraitValue)

# Create dataframe with LAT and LONG for map making
combo.6.map <- combo.5 %>% 
  dplyr::select(LAT, LONG, Dataset_Sitename, TraitName, Region)
    # Ready for gap-filling

# Remove unwanted intermediate combo. variables so far
rm(combo.cut, combo.duplicates.1, combo.duplicates.2,  combo.duplicates.all,
   combo.duplicates.all.sets, combo.geographic, combo.non.duplicates,
   combo.round, coords.c.c, coords.c.c.fix, coords.c.c.lgm, coords.c.c.na,
   coords.combo, coords.region, lgm, lgm.combo, lgm.coords, lgm.coords.df,
   lgm.count, lgm.dupl.count, lgm.extract, lgm.extract.clean, coords.country,
   coords.continent, input.dataframe, combo.manual, combo.leaf.area,
   combo.leaf.area.cut, combo.heights, combo.heights.retain)

# Export cleaned, combined, gap-filling ready trait dataframe with geographic information
write.csv(combo.6, "data/output_06_TRY_TTT_clean.csv")
    
# Import cleaned, combined, gap-filling ready trait dataframe with geographic information
    # combo.6 <- read.csv("data/output_06_TRY_TTT_clean.csv") %>% dplyr::select(-X)


# MAP - MAKE A VISUALISATION OF THE TRAIT DATA LOCATIONS ----

# Produce dataframe for making plot
combo.map.df <- combo.6.map %>% 
  group_by(LAT, LONG, Dataset_Sitename) %>% 
  mutate(Number_Traits = length(unique(TraitName))) %>% 
  ungroup() %>% 
  dplyr::select(-TraitName) %>% 
  distinct(LAT, LONG, Dataset_Sitename, Number_Traits, .keep_all = TRUE) %>% 
  filter(LAT > 55) %>% 
  mutate(Number_Traits = as.factor(Number_Traits)) %>% 
  mutate(Number_Traits = factor(.$Number_Traits,
                                levels = c("1", "2", "3", "4", "5", "6", "7"))) %>% # Reorder factor
  arrange(Number_Traits) %>% 
  mutate(Region = case_when(Region == "North America-East" ~ "North America East",
                            Region == "North America-West" ~ "North America West",
                            Region == "GreenIceland" ~ "Greenland & Iceland",
                            Region == "Eurasia" ~ "Eurasia",
                            NA ~ Region))

# Transform coordinates to UTM
combo.map.transform <- transform_coord(combo.map.df, lon = "LONG", lat = "LAT", bind = TRUE)

# Create dataframe just for labelling the plot
combo.map.labels <- combo.map.transform %>% 
  filter(Number_Traits %in% c("5", "6")) %>% 
  distinct(Dataset_Sitename, Number_Traits, .keep_all = TRUE) %>% 
  mutate(Simplified_Names = case_when(Dataset_Sitename == "Qikiqtaruk-Herschel Island, Canada" ~ "Qikiqtaruk-Herschel Island",
                                      Dataset_Sitename == "Caucasus Mountains" ~ "Caucasus",
                                      Dataset_Sitename == "Kilpisjarvi, Finland" ~ "Kilpisjarvi",
                                      Dataset_Sitename == "Umiujaq, QC, Canada" ~ "Umiujaq",
                                      Dataset_Sitename == "Alexandra Fiord, Ellesmere Island, Canada" ~ "Alexandra Fiord"))

# Map based on size of circles
(map.trait.sites <- basemap(limits = 55, grid.size = 0.05, grid.col = "#949494",
                            land.size = 0.05, land.col = "#dcdcdc", land.border.col = "#000000") +
    geom_point(data = combo.map.transform,
               aes(x = lon.proj, y = lat.proj, size = Number_Traits, fill = Region),
               alpha = 0.95, colour = "#000000", shape = 21) +
    scale_size_discrete(range = c(6,12,18,24,30,36)) +
    # scale_fill_viridis(discrete = TRUE, option = "A", begin = 0.1, end = 0.95, direction = -1) +
    geom_label_repel(data = subset(combo.map.labels, Simplified_Names %in% c("Qikiqtaruk-Herschel Island",
                                                                             "Caucasus", "Kilpisjarvi", "Umiujaq",
                                                                             "Alexandra Fiord")),
                     aes(lon.proj, lat.proj, label = Simplified_Names), color = "black", box.padding = 2,
                     segment.color = "black", segment.size = 0.7, fill = "white", label.size = 0.4,  size = 4.5) +
    labs(title = "TTT and TRY Site Locations",
         # subtitle = "Glaciated during the LGM (~18,000 B.P.) (Ehlers and Gibbard, 2004)",
         subtitle = "",
         fill = "Number of Traits per Site") +
    guides(fill = guide_legend(override.aes = list(shape = 21))) + # Fixes bug in the legend
    theme_4_map()
)

# Export map
ggsave(map.trait.sites, file = "figures/traits_sites_polarview.png", height = 8, width = 8)

# Remove intermediate objects
rm(combo.6.map, combo.map.df, combo.map.labels, combo.map.transform, map.trait.sites)

# ITEX - LOAD IN ITEX DATA (itex.1) ----

# Load in manipulated ITEX dataframe
itex.1 <- read.csv("data/output_04_itex_all.csv") %>%
  dplyr::select(-X)


# ITEX - PREPARE FOR GAP-FILLING (itex.2) ----

# Create columns useful for gap-filling
itex.gapfill <- itex.1 %>%
  mutate(gridcell = paste0("_", lat_grid, "_", lon_grid),
         LAT_band = plyr::round_any(LAT, 5, f = floor),
         over60degN = ifelse(LAT >= 60, "above_60", "below_60")) %>% 
  mutate(Region_gridcell = paste0(Region, "_", gridcell),
         Region_LAT_band = paste0(Region, "_", LAT_band)) %>% 
  mutate(Region_2 = Region) %>% # Second tier for gap-filling by North America alone
  mutate(Region_2 = ifelse(str_detect(Region, pattern = "North America"), "North America_2", Region_2),
         Region_2 = ifelse(Region %in% c("Eurasia"), "Eurasia_2", Region_2),
         Region_2 = ifelse(Region %in% c("GreenIceland"), "GreenIceland_2", Region_2))

# Cut to useful columns for gap-filling
itex.2 <- itex.gapfill %>% 
  dplyr::select(ID, SPECIES, GENUS, FAMILY, FuncGroup, Region_gridcell, Region_LAT_band, over60degN, Region, Region_2)

# Remove unwanted intermediary object
rm(itex.gapfill)


# COMBO - ADD FUNCGROUP INFORMATION (combo.7) ----

# Load in dataframe of functional groups generated in script EX1 from Sarah Elmendorf's ITEX spp info.
funcgroup.spp <- read.csv("data/output_EX1_funcgroup.csv") %>% 
  dplyr::select(-c(X, GFNARROWarft, GFNARROWwalker)) %>% 
  mutate(source = "SPP")

# Create dataframe of FuncGroups by SPECIES in ITEX
funcgroup.itex <- itex.1 %>% 
  dplyr::select(SPECIES, FuncGroup) %>% 
  distinct(SPECIES, FuncGroup, .keep_all = TRUE) %>% 
  filter(!str_detect(SPECIES, pattern = "XXX")) %>% 
  mutate(source = "ITEX")

# Combine to form one growth form checklist
funcgroup.all <- rbind(funcgroup.spp, funcgroup.itex) %>% 
  mutate(FuncGroup = ifelse(SPECIES == "Linnaea borealis", "Evergreen_Shrub", FuncGroup)) %>% 
  distinct(SPECIES, FuncGroup, .keep_all = TRUE)

# Check duplicates that occur in both
funcgroup.dupl <- funcgroup.all %>% 
  group_by(SPECIES) %>% 
  mutate(records = length(FuncGroup)) %>% 
  ungroup() %>% 
  arrange(SPECIES) %>%
  filter(records > 1) # No duplicates

# Create dataframe of genuses that are all graminoids (no overlap like forb/shrub)
funcgroup.graminoids <- funcgroup.all %>% 
  separate(SPECIES, into = c("GENUS", "SPECIES"), sep = " ") %>% 
  dplyr::select(-SPECIES, -source) %>% 
  filter(FuncGroup == "Graminoid") %>% 
  distinct(GENUS, .keep_all = TRUE) %>% 
  rename(Graminoid = FuncGroup)

# Check that all graminoid genuses are exclusively allocated to 'Graminoid' FuncGroup
funcgroup.graminoids.checks <- funcgroup.all %>% 
  separate(SPECIES, into = c("GENUS", "SPECIES"), sep = " ") %>% 
  dplyr::select(-SPECIES, -source) %>%
  mutate(GRAMINOID = ifelse(FuncGroup == "Graminoid", "YES", "NO")) %>% 
  group_by(GENUS) %>% 
  mutate(FuncGroups = length(unique(FuncGroup))) %>% 
  ungroup() %>% 
  filter(GRAMINOID == "YES") %>%
  filter(FuncGroups > 1) %>% 
  unique()

# Join FuncGroup information to trait dataset
combo.funcgroup <- left_join(combo.6, funcgroup.all, by = c("NAME" = "SPECIES")) %>% 
  left_join(., funcgroup.graminoids, by = c("GENUS" = "GENUS")) %>% 
  mutate(FuncGroup = ifelse(is.na(FuncGroup) & Graminoid == "Graminoid", "Graminoid", FuncGroup)) %>%
  mutate(FuncGroup = ifelse(is.na(FuncGroup) & FAMILY %in% c("Poaceae", "Cyperaceae", "Juncaceae"), "Graminoid", FuncGroup)) %>% 
  dplyr::select(-c(source, Graminoid)) %>% 
  relocate(FuncGroup, .after = FAMILY)

# Remove unneeded intermediate variables
rm(funcgroup.spp, funcgroup.itex, funcgroup.all, funcgroup.dupl, funcgroup.graminoids, funcgroup.graminoids.checks)

# Check if any species missing FuncGroup information (434!)
combo.funcgroup.missing <- combo.funcgroup %>% 
  dplyr::select(NAME, GENUS, FAMILY, FuncGroup) %>% 
  filter(is.na(FuncGroup)) %>% 
  distinct(NAME, GENUS, .keep_all = TRUE) %>% 
  arrange(NAME)

# Create vectors: genuses that can all be assigned to a certain FuncGroup
missing.g.eshrub <- c("Acinos", "Dryas", "Erica", "Picea")
missing.g.dshrub <- c("Betula", "Larix", "Salix", "Spiraea")
missing.g.forb <- c("Achillea", "Aconitum", "Ajuga", "Alchemilla", "Androsace", "Anemone", "Angelica", "Anthemis",
                    "Anthriscus", "Arnica", "Astragalus", "Athyrium", "Barbarea", "Bupleurum", "Butomus", "Cerastium",
                    "Chaerophyllum", "Cirsium", "Colchicum", "Corydalis", "Cruciata", "Cypripedium", "Dactylorhiza",
                    "Delphinium", "Dianthus", "Draba", "Eremogone", "Eritrichium", "Eryngium", "Erythronium", "Euphrasia",
                    "Gagea", "Gentiana", "Gentianella", "Heracleum", "Herniaria", "Lactuca", "Lagotis", "Lathyrus",
                    "Ligularia", "Ligusticum", "Maianthemum", "Matricaria", "Melampyrum", "Myosotis", "Paeonia",
                    "Papaver", "Parasenecio", "Pedicularis", "Pinguicula", "Platanthera", "Polystichum", "Prunella", "Pulmonaria",
                    "Pulsatilla", "Ranunculus", "Rhinanthus", "Rhodiola", "Rorippa", "Sanguisorba", "Saxifraga",
                    "Seseli", "Stellaria", "Swertia", "Tanacetum", "Taraxacum", "Thymelaea", "Tripleurospermum",
                    "Urtica", "Veratrum", "Anthyllis", "Eria")

# Create vectors: species that can be assigned to a certain FuncGroup
missing.s.graminoid <- c("XXXGraminoid")
missing.s.eshrub <- c("Artemisia alaskana", "Artemisia gmelinii", "Cassiope ericoides", "Chamaedaphne calyculata", 
                      "Daphne laureola", "Dryas grandis", "Eriogonum arcuatum", "Euphorbia flavicoma", "Euphorbia pyrenaica",
                      "Glandora diffusa", "Globularia nudicaulis", "Globularia repens", "Helianthemum apenninum",
                      "Helianthemum oelandicum", "Juniperus sabina", "Linum suffruticosum", "Luetkea pectinata",
                      "Myrtus communis", "Pinus sylvestris", "Rhododendron adansonii", "Rhododendron aureum",
                      "Rhododendron caucasicum", "Rhododendron parviflorum", "Teucrium chamaedrys", "Thymus pulegioides",
                      "Vaccinium oxycoccos")
missing.s.dshrub <- c("Alnus alnobetula", "Alnus hirsuta", "Alnus incana", "Berberis vulgaris", "Caragana arborescens",
                      "Caragana ussuriensis", "Cornus alba", "Cornus canadensis", "Cornus mas", "Cornus suecica",
                      "Corylus avellana", "Genista balansae", "Genista hispanica", "Genista hystrix", "Lonicera caerulea",
                      "Packera cana", "Populus tremula", "Populus tremuloides", "Prunus padus", "Ribes aciculare",
                      "Ribes dicuscha", "Ribes hudsonianum", "Ribes montigenum", "Ribes nigrum", "Ribes procumbens",
                      "Ribes rubrum", "Ribes spicatum", "Ribes triste", "Rubus idaeus", "Rubus saxatilis",
                      "Sambucus sibirica", "Senecio triangularis", "Sorbus aucuparia", "Sorbus intermedia",
                      "Sorbus sambucifolia", "Thymus nummularius")
missing.s.forb <- c("Alisma plantago-aquatica", "Allium schoenoprasum", "Agoseris aurantiaca", "Alyssum montanum",
                    "Anaphalis margaritacea", "Anticlea elegans", "Aquilegia caerulea", "Arctanthemum arcticum",
                    "Arctium lappa", "Arenaria grandiflora", "Arenaria purpurascens", "Arenaria serpyllifolia",
                    "Artemisia vulgaris", "Lycopodium clavatum", "Bellidastrum michelii", "Dolichorrhiza caucasica",
                    "Iranecio taraxacifolius", "Senecio kolenatianus",
                    "Beckwithia glacialis", "Betonica macrantha", "Campanula parryi", "Capsella bursa-pastoris",
                    "Carthamus mitissimus", "Chamaesciadium acaule", "Chenopodium album", "Chrysosplenium alternifolium",
                    "Corallorhiza trifida", "Crepis albida", "Crepis caucasica", "Crepis leontodontoides",
                    "Crepis paludosa", "Crepis sibirica", "Cryptantha cana", "Cyanus cheiranthifolius",
                    "Cymopterus alpinus", "Dethawia splendens", "Diphasiastrum complanatum", "Dicentra peregrina",
                    "Drosera anglica", "Drosera rotundifolia", "Dryopteris filix-mas", "Epipactis atrorubens",
                    "Equisetum sylvaticum", "Erysimum capitatum", "Erysimum cheiranthoides", "Fragaria vesca",
                    "Galeopsis speciosa", "Geranium albiflorum", "Geranium erianthum", "Geum aleppicum", "Geum rossii",
                    "Gnaphalium norvegicum", "Gymnocarpium dryopteris", "Hackelia deflexa", "Hedysarum boreale",
                    "Hepatica nobilis", "Hieracium froelichianum","Hieracium hypoglaucum", "Hieracium laevigatum",
                    "Hieracium mixtum", "Hieracium onosmoides", "Hieracium prenanthoides", "Hieracium umbellatum",
                    "Hypericum linarioides", "Jurinea humilis", "Kemulariella caucasica", "Koenigia hadacii",
                    "Laserpitium halleri", "Linaria vulgaris", "Linum narbonense", "Lomelosia caucasica",
                    "Lotus corniculatus", "Lupinus polyphyllus", "Lysimachia europaea", "Lysimachia vulgaris",
                    "Matteuccia struthiopteris", "Medicago minima", "Melilotus albus", "Mertensia ciliata",
                    "Minuartia recurva", "Moneses uniflora", "Neottia cordata", "Oxalis acetosella", 
                    "Packera crocata", "Paris quadrifolia", "Pedicularis compacta", "Pedicularis condensata",
                    "Pedicularis palustris", "Pedicularis pyrenaica", "Pedicularis resupinata",
                    "Pedicularis sceptrum-carolinum", "Phegopteris connectilis", "Phyteuma spicatum",
                    "Pilosella lactucella", "Pilosella officinarum", "Pimpinella saxifraga", "Plantago lanceolata",
                    "Plantago major", "Plantago maritima", "Plantago media", "Podospermum canum", "Polemonium caeruleum",
                    "Polemonium pulcherrimum", "Polygala calcarea", "Polygonum aviculare", "Polygonum ellipticum",
                    "Potamogeton alpinus", "Potentilla atrosanguinea", "Potentilla concinna", "Potentilla fragiformis",
                    "Potentilla norvegica", "Potentilla rupifraga", "Primula auriculata", "Primula glutinosa",
                    "Primula matthioli", "Primula mazurenkoae", "Pritzelago alpina", "Psephellus caucasicus",
                    "Pseudocymopterus montanus", "Pyrola media", "Rumex acetosella", "Rumex alpinus", "Rumex confertus",
                    "Rumex obtusifolius", "Saussurea mae", "Scheuchzeria palustris", "Scorzoneroides helvetica",
                    "Sedum acre", "Sedum annuum", "Sedum cyaneum", "Sedum kamtschaticum", "Sedum roseum",
                    "Sedum telephium", "Senecio aurantiacus", "Senecio nemorensis", "Serratula tinctoria",
                    "Silene arvatica", "Silene ciliata", "Silene jeniseensis", "Silene nutans", "Silene samojedorum",
                    "Silene stenophylla", "Silene uniflora", "Solidago spathulata", "Sparganium angustifolium",
                    "Tephroseris palustris", "Tetraneuris acaulis", "Teucrium pyrenaicum", "Teucrium siculum",
                    "Thalictrum minus", "Traunsteinera globosa", "Trifolium medium", "Tussilago farfara",
                    "Valeriana alpestris", "Valeriana officinalis", "Valeriana wolgensis", "Veronica longifolia",
                    "Vicia sepium", "Vicia sylvatica", "Viola canina", "Viola epipsiloides")

# Create vectors: species names that weren't able to confidently assign to a certain FuncGroup
missing.uncorrected <- c("Artemisia lagocephala", "Artemisia lagopus", "Artemisia leucophylla",
                         "Campanula gieseckeana", "Catolobus pendulus", "Lycopodium clavatum", 
                         "Lycopodium complanatum", "Minuartia tricostata", "Polygala edmundii",
                         "Potentilla stolonifera", "Ptarmica alpina", "Ptarmica camtschatica",
                         "Sieversia pusilla", "Trifolium brandegei", "XXXArctostaphylos",
                         "XXXRibes", "XXXVaccinium")

# Create vectors: species names that can be removed (moss and unknown)
missing.remove <- c("Polytrichum commune", "XXXUnknown")

# NOTE: certain species of 'shrub' are actually trees (e.g.
# "Populus tremula", "Populus tremuloides", "Pinus sylvestris", "Prunus padus", "Sorbus aucuparia", "Sorbus intermedia"
# Why we put the conditionality on creating gap-filled means

# NOTE: The vast majority of growth forms were collected from the following two resources:
# https://www.usda.gov/
# https://eol.org/

# Manually correct the missing FuncGroups by SPECIES AND GENUS
combo.7 <- combo.funcgroup %>% # CHANGE!
  mutate(FuncGroup = case_when(is.na(FuncGroup) & GENUS %in% missing.g.eshrub ~ "Evergreen_Shrub", # Correct by genus
                               is.na(FuncGroup) & GENUS %in% missing.g.dshrub ~ "Deciduous_Shrub",
                               is.na(FuncGroup) & GENUS %in% missing.g.forb ~ "Forb",
                               TRUE ~ FuncGroup)) %>% 
  mutate(FuncGroup = case_when(is.na(FuncGroup) & NAME %in% missing.s.eshrub ~ "Evergreen_Shrub", # Correct by species
                               is.na(FuncGroup) & NAME %in% missing.s.dshrub ~ "Deciduous_Shrub",
                               is.na(FuncGroup) & NAME %in% missing.s.forb ~ "Forb",
                               is.na(FuncGroup) & NAME %in% missing.s.graminoid ~ "Graminoid",
                               TRUE ~ FuncGroup)) %>% 
  mutate(FuncGroup = ifelse(NAME %in% missing.uncorrected, NA, FuncGroup)) %>% # List as NA if FG not found
  filter(NAME %notin% missing.remove) # Remove species not needed

# Remove unwanted objects
rm(missing.s.eshrub, missing.s.dshrub, missing.s.forb, missing.s.graminoid, missing.g.dshrub, missing.g.eshrub,
   combo.funcgroup, combo.funcgroup.missing, missing.g.forb, missing.remove, missing.uncorrected)

# Create dataframe of growth forms to export
combo.gf.export <- combo.7 %>% 
  dplyr::select(NAME, FuncGroup) %>% 
  mutate(FuncGroup = ifelse(NAME == "Linnaea borealis", "Forb", FuncGroup)) %>% 
  distinct() %>% 
  arrange(NAME)

# Write to csv
write.csv(combo.gf.export, file = "data/output_fg_export_traits.csv", row.names = FALSE)


# ITEX/COMBO/TRAIT - PREPARING TRAIT DATA FRAMES  (itex.3/combo.8/trait....) ----

# Generate unique identifiers for all the levels we want to gap-fill by for ITEX
itex.3 <- itex.2 %>% 
  mutate(sp_gc = paste0(SPECIES, "_", Region_gridcell), # 1) Species within gridcell
         sp_r_b = paste0(SPECIES, "_", Region_LAT_band), # 2) Species within LAT_band of Region
         sp_r_60 = paste0(SPECIES, "_", Region, "_", over60degN), # 3) Species within Region > 60 degN
         sp_r = paste0(SPECIES, "_", Region), # 4) Species within Region
         sp_r_2 = paste0(SPECIES, "_", Region_2), # 5) Species within Region_2 (combines N.America)
         sp = SPECIES, # 6) Species OVERALL
         g_gc = paste0(GENUS, "_", Region_gridcell), # As above, with GENUS
         g_r_b = paste0(GENUS, "_", Region_LAT_band),
         g_r_60 = paste0(GENUS, "_", Region, "_", over60degN),
         g_r = paste0(GENUS, "_", Region),
         g_r_2 = paste0(GENUS, "_", Region_2),
         g = GENUS,
         f_gc = paste0(FAMILY, "_", Region_gridcell), # As above, with FAMILY
         f_r_b = paste0(FAMILY, "_", Region_LAT_band),
         f_r_60 = paste0(FAMILY, "_", Region, "_", over60degN),
         f_r = paste0(FAMILY, "_", Region),
         f_r_2 = paste0(FAMILY, "_", Region_2),
         f = FAMILY,
         fg_gc = paste0(FuncGroup, "_", Region_gridcell), # As above, with FuncGroup
         fg_r_b = paste0(FuncGroup, "_", Region_LAT_band),
         fg_r_60 = paste0(FuncGroup, "_", Region, "_", over60degN),
         fg_r = paste0(FuncGroup, "_", Region),
         fg_r_2 = paste0(FuncGroup, "_", Region_2),
         fg = FuncGroup) %>% 
  dplyr::select(-c(Region_gridcell, Region_LAT_band, over60degN, Region, Region_2)) %>% 
  pivot_longer(., 6:ncol(.), names_to = "UnitType", values_to = "Unit") %>% 
  mutate(Unit = ifelse(str_detect(SPECIES, pattern = "XXX") & str_detect(UnitType, pattern = "sp"), NA, Unit)) %>% 
  dplyr::select(-c(GENUS, FAMILY, FuncGroup)) # ^^ If genus-level morphospecies, removing species-level Unit

# Testing new method
combo.8 <- combo.7 %>% 
  mutate(sp_gc = paste0(NAME, "_", Region_gridcell), # 1) Species within gridcell
         sp_r_b = paste0(NAME, "_", Region_LAT_band), # 2) Species within LAT_band of Region
         sp_r_60 = paste0(NAME, "_", Region, "_", over60degN), # 3) Species within Region > 60 degN
         sp_r = paste0(NAME, "_", Region), # 4) Species within Region
         sp_r_2 = paste0(NAME, "_", Region_2), # 5) Species within Region_2 (combines N.America)
         sp = NAME, # 6) Species OVERALL
         g_gc = paste0(GENUS, "_", Region_gridcell), # As above, with GENUS
         g_r_b = paste0(GENUS, "_", Region_LAT_band),
         g_r_60 = paste0(GENUS, "_", Region, "_", over60degN),
         g_r = paste0(GENUS, "_", Region),
         g_r_2 = paste0(GENUS, "_", Region_2),
         g = GENUS,
         f_gc = paste0(FAMILY, "_", Region_gridcell), # As above, with FAMILY
         f_r_b = paste0(FAMILY, "_", Region_LAT_band),
         f_r_60 = paste0(FAMILY, "_", Region, "_", over60degN),
         f_r = paste0(FAMILY, "_", Region),
         f_r_2 = paste0(FAMILY, "_", Region_2),
         f = FAMILY,
         fg_gc = paste0(FuncGroup, "_", Region_gridcell), # As above, with FuncGroup
         fg_r_b = paste0(FuncGroup, "_", Region_LAT_band),
         fg_r_60 = paste0(FuncGroup, "_", Region, "_", over60degN),
         fg_r = paste0(FuncGroup, "_", Region),
         fg_r_2 = paste0(FuncGroup, "_", Region_2),
         fg = FuncGroup) %>%
  dplyr::select(-c(Region_gridcell, Region_LAT_band, over60degN, Region, Region_2)) %>% 
  pivot_longer(., 8:ncol(.), names_to = "UnitType", values_to = "Unit") %>% 
  mutate(Unit = ifelse(is.na(GENUS) | is.na(FAMILY) | is.na(FuncGroup), NA, Unit)) %>% # If no genus/family/fg, assigning unit as NA
  group_by(Unit, TraitName) %>% # Add columns counting number of records for that value type
  mutate(Records = length(TraitValue)) %>% 
  ungroup() %>% 
  mutate(Records = ifelse(is.na(Unit), NA, Records)) %>% # if unit is NA, assigning Records as NA
  mutate(TraitValue = ifelse(str_detect(Unit, pattern = "below_60"), NA, TraitValue)) # Setting trait values to NA when below 60degN so not included in the 'above60' gap-fill level
  
# Generate error risks per Unit for each trait
combo.error.risk <- combo.8 %>% 
  dplyr::select(RecordID, Unit, UnitType, TraitName, TraitValue, Records) %>% 
  group_by(Unit, TraitName) %>% 
  mutate(sd_unit = sd(TraitValue),
         mean_unit = mean(TraitValue)) %>% 
  ungroup() %>% 
  mutate(Error_Risk = (abs(TraitValue - mean_unit) / sd_unit)) %>% # abs() is absolute value, makes all differences positive
  dplyr::select(-c(sd_unit, mean_unit))

# Generate average trait values per Unit for each trait
combo.medians <- combo.error.risk %>% 
  filter(Records >= 4) %>% # Don't use a unit if 1-3 records
  mutate(Retain = "Remove") %>% # Set default for column to remove
  mutate(Retain = case_when(Records >= 4 & Records < 10 & Error_Risk < 2.25 ~ "Keep", # 4-9 records
                            Records >= 10 & Records < 20 & Error_Risk < 2.75 ~ "Keep", # 10-19 records
                            Records >= 20 & Records < 30 & Error_Risk < 3.25 ~ "Keep", # 20-29 records
                            Records >= 30 & Error_Risk < 4 ~ "Keep", # 30+ records
                            TRUE ~ Retain)) %>%
  filter(Retain == "Keep") %>% # Only keep records that fall within Error Risk stipulations
  dplyr::select(-c(Retain, Records)) %>%
  group_by(Unit, TraitName) %>% # Create new column counting number of records for that Unit type now that values outside error risk band have been removed
  mutate(Records = length(TraitValue)) %>%
  ungroup() %>%
  filter(Records >= 4) %>% # Ensure that Unit still has minimum four units for calculating average trait value
  group_by(Unit, TraitName) %>% # Calculating averages by unit (e.g species/genus etc. in certain region etc.)
  summarise(MedianTraitValue = median(TraitValue, na.rm = TRUE)) %>% 
  ungroup()


# COMBO - GAP-FILLING (combo.9) ----

# Create gap-filling function
gap_fill <- function(df){
  mutate(df, accepted_value = ifelse(!is.na(sp_gc), sp_gc,
                                     ifelse(!is.na(sp_r_b), sp_r_b,
                                            ifelse(!is.na(sp_r_60), sp_r_60,
                                                   ifelse(!is.na(sp_r), sp_r,
                                                          ifelse(!is.na(sp_r_2), sp_r_2,
                                                                 ifelse(!is.na(sp), sp,
                                                                        ifelse(!is.na(g_gc), g_gc,
                                                                               ifelse(!is.na(g_r_b), g_r_b,
                                                                                      ifelse(!is.na(g_r_60), g_r_60,
                                                                                             ifelse(!is.na(g_r), g_r,
                                                                                                    ifelse(!is.na(g_r_2), g_r_2,
                                                                                                           ifelse(!is.na(g), g,
                                                                                                                  ifelse(!is.na(f_gc), f_gc,
                                                                                                                         ifelse(!is.na(f_r_b), f_r_b,
                                                                                                                                ifelse(!is.na(f_r_60), f_r_60,
                                                                                                                                       ifelse(!is.na(f_r), f_r,
                                                                                                                                              ifelse(!is.na(f_r_2), f_r_2,
                                                                                                                                                     ifelse(!is.na(f), f,
                                                                                                                                                            ifelse(!is.na(fg_gc), fg_gc,
                                                                                                                                                                   ifelse(!is.na(fg_r_b), fg_r_b,
                                                                                                                                                                          ifelse(!is.na(fg_r_60), fg_r_60,
                                                                                                                                                                                 ifelse(!is.na(fg_r), fg_r,
                                                                                                                                                                                        ifelse(!is.na(fg_r_2), fg_r_2,
                                                                                                                                                                                               ifelse(!is.na(fg), fg, NA)))))))))))))))))))))))))
}

# # Now create individual ITEX dataframes for each trait
itex.LDMC <- left_join(itex.3, filter(combo.medians, TraitName == "LDMC"), by = c("Unit" = "Unit")) %>% 
  dplyr::select(-c(Unit, TraitName)) %>% 
  pivot_wider(names_from = "UnitType", values_from = "MedianTraitValue") %>% # Trait values vanishing here (using Achillea millefolium as an example)
  gap_fill(.) %>% 
  rename(LDMC = accepted_value) %>% 
  dplyr::select(ID, LDMC)
  
itex.Leaf_Area <- left_join(itex.3, filter(combo.medians, TraitName == "Leaf_Area"), by = c("Unit" = "Unit")) %>% 
  dplyr::select(-c(Unit, TraitName)) %>% 
  pivot_wider(names_from = "UnitType", values_from = "MedianTraitValue") %>% # Trait values vanishing here (using Achillea millefolium as an example)
  gap_fill(.) %>% 
  rename(Leaf_Area = accepted_value) %>% 
  dplyr::select(ID, Leaf_Area)

itex.SLA <- left_join(itex.3, filter(combo.medians, TraitName == "SLA"), by = c("Unit" = "Unit")) %>% 
  dplyr::select(-c(Unit, TraitName)) %>% 
  pivot_wider(names_from = "UnitType", values_from = "MedianTraitValue") %>% # Trait values vanishing here (using Achillea millefolium as an example)
  gap_fill(.) %>% 
  rename(SLA = accepted_value) %>% 
  dplyr::select(ID, SLA)

itex.LeafN <- left_join(itex.3, filter(combo.medians, TraitName == "LeafN"), by = c("Unit" = "Unit")) %>% 
  dplyr::select(-c(Unit, TraitName)) %>% 
  pivot_wider(names_from = "UnitType", values_from = "MedianTraitValue") %>% # Trait values vanishing here (using Achillea millefolium as an example)
  gap_fill(.) %>% 
  rename(LeafN = accepted_value) %>% 
  dplyr::select(ID, LeafN)

itex.LeafP <- left_join(itex.3, filter(combo.medians, TraitName == "LeafP"), by = c("Unit" = "Unit")) %>% 
  dplyr::select(-c(Unit, TraitName)) %>% 
  pivot_wider(names_from = "UnitType", values_from = "MedianTraitValue") %>% # Trait values vanishing here (using Achillea millefolium as an example)
  gap_fill(.) %>% 
  rename(LeafP = accepted_value) %>% 
  dplyr::select(ID, LeafP)

itex.Plant_Height <- left_join(itex.3, filter(combo.medians, TraitName == "Plant_Height"), by = c("Unit" = "Unit")) %>% 
  dplyr::select(-c(Unit, TraitName)) %>% 
  pivot_wider(names_from = "UnitType", values_from = "MedianTraitValue") %>% # Trait values vanishing here (using Achillea millefolium as an example)
  gap_fill(.) %>% 
  rename(Plant_Height = accepted_value) %>% 
  dplyr::select(ID, Plant_Height)

itex.SDM <- left_join(itex.3, filter(combo.medians, TraitName == "SDM"), by = c("Unit" = "Unit")) %>% 
  dplyr::select(-c(Unit, TraitName)) %>% 
  pivot_wider(names_from = "UnitType", values_from = "MedianTraitValue") %>% # Trait values vanishing here (using Achillea millefolium as an example)
  gap_fill(.) %>% 
  rename(SDM = accepted_value) %>% 
  dplyr::select(ID, SDM)

# Join the individual trait values together to generate the complete gap-filled trait dataset
combo.9 <- left_join(itex.LDMC, itex.Leaf_Area, by = c("ID" = "ID")) %>% 
  left_join(., itex.SLA, by = c("ID" = "ID")) %>% 
  left_join(., itex.LeafN, by = c("ID" = "ID")) %>% 
  left_join(., itex.LeafP, by = c("ID" = "ID")) %>% 
  left_join(., itex.Plant_Height, by = c("ID" = "ID")) %>% 
  left_join(., itex.SDM, by = c("ID" = "ID"))

# Remove unwanted intermediary variables
rm(combo.error.risk, itex.LDMC, itex.Leaf_Area, itex.SLA, itex.LeafN, itex.LeafP, itex.Plant_Height, itex.SDM)


# COMBO - ADD TRAIT: Woodiness (combo.10) ----

# Regenerate list of functional groups by ID from itex.1
funcgroup.woodiness <- itex.1 %>% 
  dplyr::select(ID, FuncGroup)

# Create vectors of woody and not woody growth forms
woody <- c("Evergreen_Shrub", "Deciduous_Shrub")
not.woody <- c("Forb", "Graminoid")

# Create column for woodiness
# 0 = NOT woody, or herbaceous (FG: FORBS & GRAMINOIDS)
# 1 = woody (FG: D.SHRUBS & E.SHRUBS)
combo.10 <- left_join(combo.9, funcgroup.woodiness, by = c("ID" = "ID")) %>% 
  relocate(FuncGroup, .after = ID) %>% 
  mutate(Woodiness = case_when(FuncGroup %in% not.woody ~ 0,
                               FuncGroup %in% woody ~ 1))

# Remove unwanted variables
rm(funcgroup.woodiness, woody, not.woody)


# COMBO - ADD TRAIT: Evergreenness (combo.11) ----

# Create vectors of evergreen and not evergreen growth forms
evergreen <- c("Evergreen_Shrub")
not.evergreen <- c("Deciduous_Shrub", "Forb", "Graminoid")

# Create column for evergreenness (evergreen woody species)
# 0 = NOT woody and evergreen (FG: D.SHRUBS, FORBS & GRAMINOIDS)
# 1 = woody and evergreen (FG: E.SHRUBS)
combo.11 <- combo.10 %>% 
  mutate(Evergreenness = case_when(FuncGroup %in% not.evergreen ~ 0,
                                   FuncGroup %in% evergreen ~1)) %>% 
  dplyr::select(-FuncGroup)

# Remove unwanted variables
rm(evergreen, not.evergreen)


# COMBO - EXPORT ----

# Export final trait dataframe
write.csv(combo.11, "data/output_06_traits_final.csv")


# PROP - IDENTIFYING LEVEL OF TRAIT MEDIANS BY TAXONOMIC LEVEL (prop.1) ----

# Create gap-filling function
gapfill_level <- function(df){
  mutate(df, gapfill_level = ifelse(!is.na(sp_gc), "sp_gc",
                                    ifelse(!is.na(sp_r_b), "sp_r_b",
                                           ifelse(!is.na(sp_r_60), "sp_r_60",
                                                  ifelse(!is.na(sp_r), "sp_r",
                                                         ifelse(!is.na(sp_r_2), "sp_r_2",
                                                                ifelse(!is.na(sp), "sp",
                                                                       ifelse(!is.na(g_gc), "g_gc",
                                                                              ifelse(!is.na(g_r_b), "g_r_b",
                                                                                     ifelse(!is.na(g_r_60), "g_r_60",
                                                                                            ifelse(!is.na(g_r), "g_r",
                                                                                                   ifelse(!is.na(g_r_2), "g_r_2",
                                                                                                          ifelse(!is.na(g), "g",
                                                                                                                 ifelse(!is.na(f_gc), "f_gc",
                                                                                                                        ifelse(!is.na(f_r_b), "f_r_b",
                                                                                                                               ifelse(!is.na(f_r_60), "f_r_60",
                                                                                                                                      ifelse(!is.na(f_r), "f_r",
                                                                                                                                             ifelse(!is.na(f_r_2), "f_r_2",
                                                                                                                                                    ifelse(!is.na(f), "f",
                                                                                                                                                           ifelse(!is.na(fg_gc), "fg_gc",
                                                                                                                                                                  ifelse(!is.na(fg_r_b), "fg_r_b",
                                                                                                                                                                         ifelse(!is.na(fg_r_60), "fg_r_60",
                                                                                                                                                                                ifelse(!is.na(fg_r), "fg_r",
                                                                                                                                                                                       ifelse(!is.na(fg_r_2), "fg_r_2",
                                                                                                                                                                                              ifelse(!is.na(fg), "fg", NA)))))))))))))))))))))))))
}

# Determine where the median trait value was sourced from by each taxonomic level
prop.LDMC <- left_join(itex.3, filter(combo.medians, TraitName == "LDMC"), by = c("Unit" = "Unit")) %>% 
  dplyr::select(-Unit, -TraitName) %>% 
  pivot_wider(names_from = "UnitType", values_from = "MedianTraitValue") %>% 
  gapfill_level(.) %>% 
  rename(LDMC = gapfill_level) %>% 
  dplyr::select(ID, LDMC)

prop.Leaf_Area <- left_join(itex.3, filter(combo.medians, TraitName == "Leaf_Area"), by = c("Unit" = "Unit")) %>% 
  dplyr::select(-Unit, -TraitName) %>% 
  pivot_wider(names_from = "UnitType", values_from = "MedianTraitValue") %>% 
  gapfill_level(.) %>% 
  rename(Leaf_Area = gapfill_level) %>% 
  dplyr::select(ID, Leaf_Area)

prop.SLA <- left_join(itex.3, filter(combo.medians, TraitName == "SLA"), by = c("Unit" = "Unit")) %>% 
  dplyr::select(-Unit, -TraitName) %>% 
  pivot_wider(names_from = "UnitType", values_from = "MedianTraitValue") %>% 
  gapfill_level(.) %>% 
  rename(SLA = gapfill_level) %>% 
  dplyr::select(ID, SLA)

prop.LeafN <- left_join(itex.3, filter(combo.medians, TraitName == "LeafN"), by = c("Unit" = "Unit")) %>% 
  dplyr::select(-Unit, -TraitName) %>% 
  pivot_wider(names_from = "UnitType", values_from = "MedianTraitValue") %>% 
  gapfill_level(.) %>% 
  rename(LeafN = gapfill_level) %>% 
  dplyr::select(ID, LeafN)

prop.LeafP <- left_join(itex.3, filter(combo.medians, TraitName == "LeafP"), by = c("Unit" = "Unit")) %>% 
  dplyr::select(-Unit, -TraitName) %>% 
  pivot_wider(names_from = "UnitType", values_from = "MedianTraitValue") %>% 
  gapfill_level(.) %>% 
  rename(LeafP = gapfill_level) %>% 
  dplyr::select(ID, LeafP)

prop.Plant_Height <- left_join(itex.3, filter(combo.medians, TraitName == "Plant_Height"), by = c("Unit" = "Unit")) %>% 
  dplyr::select(-Unit, -TraitName) %>% 
  pivot_wider(names_from = "UnitType", values_from = "MedianTraitValue") %>% 
  gapfill_level(.) %>% 
  rename(Plant_Height = gapfill_level) %>% 
  dplyr::select(ID, Plant_Height)

prop.SDM <- left_join(itex.3, filter(combo.medians, TraitName == "SDM"), by = c("Unit" = "Unit")) %>% 
  dplyr::select(-Unit, -TraitName) %>% 
  pivot_wider(names_from = "UnitType", values_from = "MedianTraitValue") %>% 
  gapfill_level(.) %>% 
  rename(SDM = gapfill_level) %>% 
  dplyr::select(ID, SDM)

# Join the individual types of median trait value together to generate the complete type of trait dataset
prop.1 <- left_join(prop.LDMC, prop.Leaf_Area, by = c("ID" = "ID")) %>% 
  left_join(., prop.SLA, by = c("ID" = "ID")) %>% 
  left_join(., prop.LeafN, by = c("ID" = "ID")) %>% 
  left_join(., prop.LeafP, by = c("ID" = "ID")) %>% 
  left_join(., prop.Plant_Height, by = c("ID" = "ID")) %>%
  left_join(., prop.SDM, by = c("ID" = "ID"))

# Remove unwanted intermediate variables
rm(combo.medians, prop.LDMC, prop.Leaf_Area,  prop.SLA, prop.LeafN, prop.LeafP, prop.Plant_Height, prop.SDM)


# PROP - WORK OUT PERCENTAGE OF EACH TRAIT COMPRISED BY EACH LEVEL (prop.2 & prop.3) ----

# Generate percentage of each trait comprised by each taxonomic level and geographic division
prop.2 <- prop.1 %>% 
  pivot_longer(2:ncol(.), names_to = "Trait", values_to = "Data_Type") %>%
  group_by(Trait) %>% 
  mutate(Total_Plots = length(unique(ID))) %>% # Total number of plot records (32,272 each time)
  ungroup() %>% 
  group_by(Trait, Data_Type, Total_Plots) %>% 
  tally() %>%
  ungroup() %>% 
  mutate(Percentage = (n / Total_Plots) * 100) %>% # Percentage of each trait comprised by each level
  dplyr::select(-n, -Total_Plots)

# Generate rows for mean proportion of trait values comprised by each level across all six traits
prop.3 <- prop.2 %>% 
  group_by(Data_Type) %>% 
  summarise(Percentage = mean(Percentage)) %>% # Calculate mean percentage derived from each level
  ungroup() %>% 
  mutate(Trait = "Mean") %>% # Generate column saying that the trait is in fact 'mean'
  rbind(prop.2, .) %>% # Bind back together with original dataset
  mutate(Trait = case_when(Trait %in% c("LDMC") ~ "LDMC",
                           Trait %in% c("Leaf_Area") ~ "Leaf Area",
                           Trait %in% c("SLA") ~ "SLA",
                           Trait %in% c("LeafN") ~ "Leaf N",
                           Trait %in% c("LeafP") ~ "Leaf P",
                           Trait %in% c("Plant_Height") ~ "Plant Height",
                           Trait %in% c("SDM") ~ "SDM",
                           Trait %in% c("Mean") ~ "Mean"))

# Work out proportion of values at continent of global level for chapter 2
proportion.continent.global <- prop.3 %>% 
  filter(Trait == "Mean") %>% 
  mutate(Retain = TRUE) %>% 
  mutate(Retain = ifelse(str_detect(Data_Type, pattern = "gc"), FALSE, Retain),
         Retain = ifelse(str_detect(Data_Type, pattern = "_r"), FALSE, Retain),
         Retain = ifelse(str_detect(Data_Type, pattern = "_r_2"), TRUE, Retain)) %>% 
  filter(Retain == TRUE)

# Determine total
total.continent.global <- sum(proportion.continent.global$Percentage)


# PROP - STACKED BAR CHART (prop.4) ----

# Create table calculating the total proportion of data drawn from the species-level by trait
prop.table <- prop.3 %>% 
  pivot_wider(names_from = Data_Type, values_from = Percentage) %>% 
  mutate_all(~replace(., is.na(.), 0)) %>% 
  group_by(Trait) %>% 
  mutate(SpeciesLevel = sum(sp_gc, sp_r_b, sp_r, sp_r_2, sp_r_60, sp)) %>% 
  ungroup() %>% 
  dplyr::select(Trait, SpeciesLevel) %>% 
  rename("Proportion of Trait Data calculated at the Species Level" = SpeciesLevel)

# Add in additional row for datatypes not gap-filled at at all (so appear in plot legend)
prop.4 <- prop.3 %>% 
  add_row(Trait = "Mean", Data_Type = "fg_r", Percentage = 0)
  
# Rename and reorder factor levels for the Taxonomic hierarchical column (biggest to smallest)
prop.4$Data_Type <- factor(prop.4$Data_Type,
                           levels = c("fg", "fg_r_2", "fg_r", "fg_r_60", "fg_r_b", "fg_gc",
                                      "f", "f_r_2", "f_r", "f_r_60", "f_r_b", "f_gc",
                                      "g", "g_r_2", "g_r", "g_r_60", "g_r_b", "g_gc",
                                      "sp", "sp_r_2", "sp_r", "sp_r_60", "sp_r_b", "sp_gc"),
                           labels = c("Functional Group",
                                      "Functional Group by Continent",
                                      "Functional Group by Region",
                                      "Functional Group by Region (> 60 degN)",
                                      "Functional Group by Region (in 5 deg Lat. Band)",
                                      "Functional Group by Region (in 0.5 deg Gridcell)",
                                      "Family",
                                      "Family by Continent",
                                      "Family by Region",
                                      "Family by Region (> 60 degN)",
                                      "Family by Region (in 5 deg Lat. Band)",
                                      "Family by Region (in 0.5 deg Gridcell)",
                                      "Genus",
                                      "Genus by Continent",
                                      "Genus by Region",
                                      "Genus by Region (> 60 degN)",
                                      "Genus by Region (in 5 deg Lat. Band)",
                                      "Genus by Region (in 0.5 deg Gridcell)",
                                      "Species",
                                      "Species by Continent",
                                      "Species by Region",
                                      "Species by Region (> 60 degN)",
                                      "Species by Region (in 5 deg Lat. Band)",
                                      "Species by Region (in 0.5 deg Gridcell)"))

# Rename and reorder factor levels for the Traits
prop.4$Trait <- factor(prop.4$Trait,
                       levels = c("Plant Height", "Leaf N", "SLA", "Leaf Area", "LDMC", "SDM", "Leaf P", "Mean"))

# Produce a stacked bar plot of the proportion of each trait comprised by each level
(prop.plot <- ggplot(prop.4, aes(x = Trait, y = Percentage, fill = Data_Type), colour = "black") +
    geom_histogram(stat = "identity", colour = "black") +
    scale_fill_manual(values = c('#f7f7f7','#d9d9d9','#bdbdbd','#969696','#636363','#252525',
                                 '#edf8e9','#c7e9c0','#a1d99b','#74c476','#31a354','#006d2c',
                                 '#eff3ff','#c6dbef','#9ecae1','#6baed6','#3182bd','#08519c',
                                 '#fee5d9','#fcbba1','#fc9272','#fb6a4a','#de2d26','#a50f15')) +
    labs(title = "Proportion of Gap-filled Trait Data",
         subtitle = "Comprised by each Taxonomic & Geographic Level",
         x = "\n Trait",
         y = "Proportion of Median Traits (%) \n",
         fill = "Gap-fill Level") +
    theme_2() +
    theme(legend.position = "right",
          axis.text.x = element_text(face = c("plain", "plain", "plain", "plain", "plain", "plain", "plain", "bold"))))

# Create simplified version just by taxonomic level (with no spatial unit)
prop.5 <- prop.4 %>% 
  filter(!str_detect(Trait, pattern = "Mean")) %>% 
  mutate(Taxonomic_Group = case_when(str_detect(Data_Type, pattern = "Functional Group") ~ "Functional Group",
                                     str_detect(Data_Type, pattern = "Family") ~ "Family",
                                     str_detect(Data_Type, pattern = "Genus") ~ "Genus",
                                     str_detect(Data_Type, pattern = "Species") ~ "Species")) %>% 
  group_by(Taxonomic_Group, Trait) %>% 
  summarise(Percentage = sum(Percentage)) %>% 
  ungroup() %>% 
  rename(Data_Type = Taxonomic_Group)

# Generate means for simplified version
prop.6 <- prop.5 %>% 
  group_by(Data_Type) %>% 
  summarise(Percentage = mean(Percentage)) %>% # Calculate mean percentage derived from each level
  ungroup() %>% 
  mutate(Trait = "Mean") %>% # Generate column saying that the trait is in fact 'mean'
  rbind(prop.5, .) # Bind back together with original dataset

# Rename and reorder factor levels for the Taxonomic hierarchical column (biggest to smallest)
prop.6$Data_Type <- factor(prop.6$Data_Type,
                           levels = c("Functional Group", "Family", "Genus", "Species"))

# Rename and reorder factor levels for the Traits
prop.6$Trait <- factor(prop.6$Trait,
                       levels = c("Plant Height", "SLA", "Leaf N", "Leaf Area", "LDMC", "Leaf P", "SDM", "Mean"))

# Produce a stacked bar plot of the proportion of each trait comprised by each summed Taxonomic level
(prop.plot.simplified <- ggplot(prop.6, aes(x = Trait, y = Percentage, fill = Data_Type), colour = "black") +
    geom_histogram(stat = "identity", colour = "black") +
    scale_fill_manual(values = c('#969696', '#74c476', '#6baed6', '#fb6a4a')) +
    labs(title = "Proportion of Gap-filled Trait Data",
         subtitle = "Comprised by each Taxonomic Level",
         x = "\n Trait",
         y = "Proportion of Median Traits (%) \n",
         fill = "Gap-fill Level") +
    theme_2() +
    theme(legend.position = "right",
          axis.text.x = element_text(face = c("plain", "plain", "plain", "plain", "plain", "plain", "plain", "bold"))))

# Create panel of simplified and non-simplified plots, and table
prop.panel <- grid.arrange(prop.plot, grid.arrange(prop.plot.simplified, tableGrob(prop.table, theme = ttheme_minimal()), widths = c(3,1)), ncol = 1)

# Export panel
ggsave(prop.panel, file = "figures/trait_proportion.png", width = 20, height = 16)


# UNIQUENESS MATRIX: Working out duplicate trait values within cells ----

    # # If starting from here
    # itex.1 <- read.csv("data/output_04_itex_all.csv") %>% dplyr::select(-X)
    # combo.10 <- read.csv("data/output_06_traits_final.csv") %>% dplyr::select(-X)

# Create vectors of columns/values we want to retain for matrix calculation
cols.itex <- c("SiteSubsitePlotYear", "SPECIES")
cols.traits <- c("LDMC", "Leaf_Area", "SLA", "LeafN", "LeafP", "Plant_Height", "SDM")
num.traits <- length(cols.traits)

# Join together toe ITEX and trait dateframes and retain selected columns
uniqueness.matrix <- left_join(itex.1, combo.9, by = c("ID" = "ID")) %>%
  dplyr::select(cols.itex, cols.traits)

# Calculate uniqueness of trait values in each plot by trait
uniqueness.long <- uniqueness.matrix %>%
  group_by(SiteSubsitePlotYear) %>%
  mutate(SR = length(unique(SPECIES)), # Calculate number of species in plot and dif trait values
         LDMC_unique = length(unique(LDMC)),
         Leaf_Area_unique = length(unique(Leaf_Area)),
         SLA_unique = length(unique(SLA)),
         LeafN_unique = length(unique(LeafN)),
         LeafP_unique = length(unique(LeafP)),
         Plant_Height_unique = length(unique(Plant_Height)),
         SDM_unique = length(unique(SDM))) %>%
  ungroup() %>%
  dplyr::select(-cols.traits) %>%
  mutate(LDMC_prop_unique = (LDMC_unique / SR) * 100, # Calculate proportion of trait values in that PlotYear combo unique
         Leaf_Area_prop_unique = (Leaf_Area_unique / SR) * 100,
         SLA_prop_unique = (SLA_unique / SR) * 100,
         LeafN_prop_unique = (LeafN_unique / SR) * 100,
         LeafP_prop_unique = (LeafP_unique / SR) * 100,
         Plant_Height_prop_unique = (Plant_Height_unique / SR) * 100,
         SDM_prop_unique = (SDM_unique / SR) * 100) %>%
  mutate(LDMC_all_unique = ifelse(SR - LDMC_unique == 0, 1, 0), # If all values for that trait are unique, assign them as 1, if not, 0
         Leaf_Area_all_unique = ifelse(SR - Leaf_Area_unique == 0, 1, 0),
         SLA_all_unique = ifelse(SR - SLA_unique == 0, 1, 0),
         LeafN_all_unique = ifelse(SR - LeafN_unique == 0, 1, 0),
         LeafP_all_unique = ifelse(SR - LeafP_unique == 0, 1, 0),
         Plant_Height_all_unique = ifelse(SR - Plant_Height_unique == 0, 1, 0),
         SDM_all_unique = ifelse(SR - SDM_unique == 0, 1, 0)) %>%
  dplyr::select(-SR, -LDMC_unique, -Leaf_Area_unique, -SLA_unique, # Remove superfluous rows
                -LeafN_unique, -LeafP_unique, -Plant_Height_unique, -SDM_unique) %>%
  mutate(TRAITS_num_unique = rowSums(dplyr::select(., LDMC_all_unique, Leaf_Area_all_unique, SLA_all_unique, # Sum number of the six traits unique in each plot
                                                   LeafN_all_unique, LeafP_all_unique, Plant_Height_all_unique, SDM_all_unique,))) %>%
  mutate(LDMC_all_unique = ifelse(LDMC_all_unique == 1, "YES", "NO"), # Convert the 1 to a YES or NO
         Leaf_Area_all_unique = ifelse(Leaf_Area_all_unique == 1, "YES", "NO"),
         SLA_all_unique = ifelse(SLA_all_unique == 1, "YES", "NO"),
         LeafN_all_unique = ifelse(LeafN_all_unique == 1, "YES", "NO"),
         LeafP_all_unique = ifelse(LeafP_all_unique == 1, "YES", "NO"),
         Plant_Height_all_unique = ifelse(Plant_Height_all_unique == 1, "YES", "NO"),
         SDM_all_unique = ifelse(SDM_all_unique == 1, "YES", "NO")) %>%
  distinct(., SiteSubsitePlotYear, .keep_all = TRUE) %>% # Retain one row per PlotYear combo
  dplyr::select(-SPECIES) %>% # No longer need species name as retaining one row per plot
  mutate(TRAITS_all_unique = ifelse(TRAITS_num_unique == num.traits, "YES", "NO")) %>% # Are all trait values in that plot unique
  relocate(., TRAITS_all_unique, TRAITS_num_unique, # Reorder cols into logical order
           LDMC_all_unique, LDMC_prop_unique,
           Leaf_Area_all_unique, Leaf_Area_prop_unique,
           SLA_all_unique, SLA_prop_unique,
           LeafN_all_unique, LeafN_prop_unique,
           LeafP_all_unique, LeafP_prop_unique,
           Plant_Height_all_unique, Plant_Height_prop_unique,
           SDM_all_unique, SDM_prop_unique, .after = SiteSubsitePlotYear)

# Calculate uniqueness summary table
uniqueness.summary <- uniqueness.long %>%
  mutate(LDMC_all_unique = ifelse(LDMC_all_unique == "YES", 1, 0), # Convert the YES/NO back to 1 or 0 for calculating summary stats
         Leaf_Area_all_unique = ifelse(Leaf_Area_all_unique == "YES", 1, 0),
         SLA_all_unique = ifelse(SLA_all_unique == "YES", 1, 0),
         LeafN_all_unique = ifelse(LeafN_all_unique == "YES", 1, 0),
         LeafP_all_unique = ifelse(LeafP_all_unique == "YES", 1, 0),
         Plant_Height_all_unique = ifelse(Plant_Height_all_unique == "YES", 1, 0),
         SDM_all_unique = ifelse(SDM_all_unique == "YES", 1, 0),
         TRAITS_all_unique = ifelse(TRAITS_all_unique == "YES", 1, 0)) %>%
  group_by() %>% # Deliberately not grouping by anything to get overall summary on rows
  summarise(Total_plots = length(unique(SiteSubsitePlotYear)),
            TRAITS_unique = sum(TRAITS_all_unique),
            TRAITS_num_unique = mean(TRAITS_num_unique),
            LDMC_unique = sum(LDMC_all_unique),
            LDMC_prop_unique = mean(LDMC_prop_unique),
            Leaf_Area_unique = sum(Leaf_Area_all_unique),
            Leaf_Area_prop_unique = mean(Leaf_Area_prop_unique),
            SLA_unique = sum(SLA_all_unique),
            SLA_prop_unique = mean(SLA_prop_unique),
            LeafN_unique = sum(LeafN_all_unique),
            LeafN_prop_unique = mean(LeafN_prop_unique),
            LeafP_unique = sum(LeafP_all_unique),
            LeafP_prop_unique = mean(LeafP_prop_unique),
            Plant_Height_unique = sum(Plant_Height_all_unique),
            Plant_Height_prop_unique = mean(Plant_Height_prop_unique),
            SDM_unique = sum(SDM_all_unique),
            SDM_prop_unique = mean(SDM_prop_unique)) %>%
  ungroup() %>%
  pivot_longer(., cols = 1:ncol(.), names_to = "Summary_Stat", values_to = "Value") %>%
  mutate(Value = round(Value, digits = 2))

# Remove unwanted variables
rm(uniqueness.matrix, uniqueness.long, cols.itex, cols.traits)

# Total number of plots
num.plots <- length(unique(itex.1$SiteSubsitePlotYear))


# Create a dataframe for just the overall summary stats
uniqueness.overall <- uniqueness.summary %>%
  filter(Summary_Stat %in% c("TRAITS_unique", "TRAITS_num_unique")) %>%
  mutate(Summary_Stat = case_when(Summary_Stat == "TRAITS_unique" ~ paste("Number of Plots with NO repeated trait values (TOTAL: ", num.plots,")"),
                                  Summary_Stat == "TRAITS_num_unique" ~ paste("Mean Number of Traits with NO repeated values in a plot (TOTAL: ", num.traits, ")"))) %>%
  rename("Overall Statistic" = Summary_Stat, "Plot Value" = Value)

# Create a table for the overall summary stats
uniqueness.table.overall <- tableGrob(uniqueness.overall, theme = ttheme_minimal())

# Remove unneeded intermediate variable
rm(uniqueness.overall, num.traits)

# Create a dataframe looking at the number of plots for each trait with 100% different trait values
uniqueness.traits <- uniqueness.summary %>%
  filter(Summary_Stat %in% c("Total_plots", "LDMC_unique", "Leaf_Area_unique", "SLA_unique",
                             "LeafN_unique", "LeafP_unique", "Plant_Height_unique", "SDM_unique")) %>%
  rename(Trait = Summary_Stat, Number_Plots = Value) %>%
  mutate(Trait = case_when(Trait == "Total_plots" ~ "Total Plots",
                           Trait == "LDMC_unique" ~ "LDMC",
                           Trait == "Leaf_Area_unique" ~ "Leaf Area",
                           Trait == "SLA_unique" ~ "SLA",
                           Trait == "LeafN_unique" ~ "Leaf N",
                           Trait == "LeafP_unique" ~ "Leaf P",
                           Trait == "Plant_Height_unique" ~ "Plant Height",
                           Trait == "SDM_unique" ~ "SDM")) %>%
  mutate(Proportion = round((Number_Plots / num.plots) * 100, digits = 2)) %>%
  mutate(Proportion = ifelse(Trait == "Total Plots", NA, Proportion))

# Modify variable types
uniqueness.traits$Trait <- as.factor(uniqueness.traits$Trait)
uniqueness.traits$Number_Plots <- as.numeric(uniqueness.traits$Number_Plots)

# Reorder bars
uniqueness.traits$Trait <- factor(uniqueness.traits$Trait,
                                  levels = c("Total Plots", "LDMC", "Leaf Area", "SLA",
                                             "Leaf N", "Leaf P", "Plant Height", "SDM"))

# Determine useful plotting stats
mean.num.plots <- mean(uniqueness.traits$Number_Plots)
label.num.plots <- mean(uniqueness.traits$Number_Plots) + (0.05 * (max(uniqueness.traits$Number_Plots) -
                                                                     min(uniqueness.traits$Number_Plots)))

# Create a barplot looking at the number of plots for each trait with 100% different trait values
(uniqueness.plot.traits <- ggplot(uniqueness.traits) +
    geom_bar(aes(x = Trait, y = Number_Plots, fill = Trait, group = Trait),
             stat = "identity", colour = "#000000") +
    scale_y_continuous(limits = c(plyr::round_any(min(uniqueness.traits$Number_Plots), 500, f = floor),
                                  plyr::round_any(max(uniqueness.traits$Number_Plots), 500, f = ceiling)),
                       oob = scales::squish) +
    geom_hline(aes(yintercept = mean.num.plots), colour = "#000000", linetype = "dashed", size = 1) +
    geom_text(aes(x = Trait, y = Proportion, label = paste(Proportion, "%"), group = Trait),
              size = 5, vjust = -16.5) +
    annotate(geom = "text", x = "SDM", y = label.num.plots, label = "Mean") +
    scale_fill_manual(values = c("#3e424b", "#ffffd4", "#fee391", "#fec44f", "#fe9929", "#d95f0e", "#993404", "#800020")) +
    labs(y = "Number of Plots",
         x = "Trait",
         title = "Total number of plots with NO repeated trait values") +
    theme_2() +
    theme(legend.position = "none"))

# Remove unneeded intermediate variable
rm(uniqueness.traits, num.plots)

# Create a dataframe looking at the mean proportion of a plot comprising different trait values
uniqueness.proportion <- uniqueness.summary %>%
  filter(Summary_Stat %in% c("LDMC_prop_unique", "Leaf_Area_prop_unique", "SLA_prop_unique", "LeafN_prop_unique",
                             "LeafP_prop_unique", "Plant_Height_prop_unique", "SDM_prop_unique")) %>%
  rename(Trait = Summary_Stat, Mean_Proportion = Value) %>%
  mutate(Trait = case_when(Trait == "LDMC_prop_unique" ~ "LDMC",
                           Trait == "Leaf_Area_prop_unique" ~ "Leaf Area",
                           Trait == "SLA_prop_unique" ~ "SLA",
                           Trait == "LeafN_prop_unique" ~ "Leaf N",
                           Trait == "LeafP_prop_unique" ~ "Leaf P",
                           Trait == "Plant_Height_prop_unique" ~ "Plant Height",
                           Trait == "SDM_prop_unique" ~ "SDM"))

# Remove unneeded intermediate variable
rm(uniqueness.summary)

# Modify variable types
uniqueness.proportion$Trait <- as.factor(uniqueness.proportion$Trait)
uniqueness.proportion$Mean_Proportion <- as.numeric(uniqueness.proportion$Mean_Proportion)

# Reorder bars
uniqueness.proportion$Trait <- factor(uniqueness.proportion$Trait,
                                      levels = c("LDMC", "Leaf Area", "SLA", "Leaf N",
                                                 "Leaf P", "Plant Height", "SDM"))

# Determine useful plotting stats
mean.prop.plots <- mean(uniqueness.proportion$Mean_Proportion)
label.prop.plots <- mean(uniqueness.proportion$Mean_Proportion) +
  (0.05 * (max(uniqueness.proportion$Mean_Proportion) -
             min(uniqueness.proportion$Mean_Proportion)))

# Create a barplot looking at the mean proportion of a plot comprising different trait values
(uniqueness.plot.proportion <- ggplot(uniqueness.proportion) +
    geom_bar(aes(x = Trait, y = Mean_Proportion, fill = Trait, group = Trait),
             stat = "identity", colour = "#000000") +
    scale_y_continuous(limits = c(plyr::round_any(min(uniqueness.proportion$Mean_Proportion), 1, f = floor),
                                  plyr::round_any(max(uniqueness.proportion$Mean_Proportion), 1, f = ceiling)),
                       oob = scales::squish) +
    geom_hline(aes(yintercept = mean.prop.plots), colour = "#000000", linetype = "dashed", size = 1) +
    geom_text(aes(x = Trait, y = Mean_Proportion, label = paste(Mean_Proportion, "%"), group = Trait),
              size = 5, vjust = -0.5) +
    annotate(geom = "text", x = "SDM", y = label.prop.plots, label = "Mean") +
    scale_fill_manual(values = c("#f1eef6", "#d0d1e6", "#a6bddb", "#74a9cf", "#2b8cbe", "#045a8d", "#001e4e")) +
    labs(y = "Mean Proportion (%)",
         x = "Trait",
         title = "Mean Proportion of unique trait values in plots") +
    theme_2() +
    theme(legend.position = "none"))

# Create a panel of outputs
uniqueness.panel <- grid.arrange(uniqueness.table.overall, uniqueness.plot.traits, uniqueness.plot.proportion, ncol = 3)

# Export panel
ggsave(uniqueness.panel, file = "figures/trait_uniqueness.png", width = 24, height = 5)

# Remove unwanted variables
rm(uniqueness.proportion, mean.num.plots, mean.prop.plots, label.num.plots, label.prop.plots)

