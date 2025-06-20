# 01 - ITEX Dataset Preparation
# October 2021, adapted February 2022, March 2022


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


# CUSTOM FUNCTIONS ----

# Not included in...
`%notin%` <- Negate(`%in%`)

# Convert factor to numeric
as.numeric.factor <- function(x) {as.numeric(levels(x))[x]}

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


# LOAD IN DATA (..._all) ----

# Load in the three required ITEX datasets
load("data/itex_june2022/input_itex_pfxy_all.RData") # Point frame data with xy coordinates
load("data/itex_june2022/input_itex_pfplot_all.RData") # Point frame data with no xy coordinates
load("data/itex_june2022/input_itex_perccov_all.RData") # Percent cover data

# Add in QHI data
qhi <- read.csv("data/itex_june2022/input_itex_qhi_incl2022.csv") %>% 
  dplyr::select(-X.1)

# Bind with the XY data
pfxy_all <- rbind(filter(pfxy_all, SITE != "QHI"), qhi)

# Remove unwanted variable
rm(qhi)


# CLEANING & FILTERING (..._all_2) ----

# We need to clean up separately before converting to cover and then binding because the dataframes have different structures

# Check unique statuses in all three dataframes
unique(perccov_all$STATUS)
unique(pfplot_all$STATUS)
unique(pfxy_all$STATUS)

# Create vector of statuses to retain across all three dataframes - keeps only live plants and NA to correct for later
status <- c("LIVE", NA, "Alive", "Live", "N/A") # Not including UNK, UNKNOWN or OTHER

# # Generating list of sites with NA status data
# status_test_perccov_all <- perccov_all %>% filter(STATUS %in% c(NA, "N/A"))
# status_test_pfplot_all <- pfplot_all %>% filter(STATUS %in% c(NA, "N/A"))
# status_test_pfxy_all <- pfxy_all %>% filter(STATUS %in% c(NA, "N/A"))
# 
# status_sites_perccov_all <- unique(status_test_perccov_all$SITE)
# status_sites_pfplot_all <- unique(status_test_pfplot_all$SITE)
# status_sites_pfxy_all <- unique(status_test_pfxy_all$SITE)
# 
# status_sites <- c(status_sites_perccov_all, status_sites_pfplot_all, status_sites_pfxy_all)
# 
# rm(status_test_perccov_all, status_test_pfplot_all, status_test_pfxy_all, status_sites_perccov_all,
#    status_sites_pfplot_all, status_sites_pfxy_all)

# Check unique treatments in all three dataframes
unique(perccov_all$TREATMENT)
unique(pfplot_all$TREATMENT)
unique(pfxy_all$TREATMENT)

# Create vector of treatments to retain across all three dataframes - keeps only control plots
treatment <- c( "CTL", "CONTROL", "CLT_GRA", NA) # Not including CTLNG and CTL_NGRA as that means grazers excluded, i.e. not control plot

# # Generating list of sites with NA status data or grazing-related status data
# treatment_test_perccov_all <- filter(perccov_all, TREATMENT %in% c("CTLNG", "CLT_GRA", "CTL_NGRA", NA))
# treatment_test_pfplot_all <- filter(pfplot_all, TREATMENT %in% c("CTLNG", "CLT_GRA", "CTL_NGRA", NA))
# treatment_test_pfxy_all <- filter(pfxy_all, TREATMENT %in% c("CTLNG", "CLT_GRA", "CTL_NGRA", NA))
# 
# treatment_sites_perccov_all <- unique(treatment_test_perccov_all$SITE)
# treatment_sites_pfplot_all <- unique(treatment_test_pfplot_all$SITE)
# treatment_sites_pfxy_all <- unique(treatment_test_pfxy_all$SITE)
# 
# treatment_sites <- c(treatment_sites_perccov_all, treatment_sites_pfplot_all, treatment_sites_pfxy_all)
# 
# rm(treatment_test_perccov_all, treatment_test_pfplot_all, treatment_test_pfxy_all, treatment_sites_perccov_all,
#    treatment_sites_pfplot_all, treatment_sites_pfxy_all)

# Check unique functional groups in all three dataframes
unique(perccov_all$GFNARROWwalker)
unique(pfplot_all$GFNARROWwalker)
unique(pfxy_all$GFNARROWwalker)

# # Check which species are unknown functional group
# fg_unknown <- c("UNK", "UNKNOWN")
# 
# fg_unk_perccov_all <- filter(perccov_all, GFNARROWwalker %in% c("UNK", "UNKNOWN"))
# fg_unk_pfplot_all <- filter(pfplot_all, GFNARROWwalker %in% c("UNK", "UNKNOWN"))
# fg_unk_pfxy_all <- filter(pfxy_all, GFNARROWwalker %in% c("UNK", "UNKNOWN"))
# 
# fg_unk_sp_perccov_all <- unique(fg_unk_perccov_all$SPECIES_NAME)
# fg_unk_sp_pfplot_all <- unique(fg_unk_pfplot_all$SPECIES_NAME)
# fg_unk_sp_pfxy_all <- unique(fg_unk_pfxy_all$SPECIES_NAME)
# 
# fg_unk_sp <- c(fg_unk_sp_perccov_all, fg_unk_sp_pfplot_all, fg_unk_sp_pfxy_all)

# Clean up unknown species
perccov_all <- perccov_all %>%
  mutate(GFNARROWwalker = case_when(GFNARROWwalker == "UNK" & SPECIES_NAME == "Cerastium arcticum" ~ "FORB",
                                    TRUE ~ GFNARROWwalker))

pfplot_all <- pfplot_all %>%
  mutate(GFNARROWwalker = case_when(GFNARROWwalker == "UNK" & SPECIES_NAME == "XXXFORB" ~ "FORB",
                                    TRUE ~ GFNARROWwalker))

pfxy_all <- pfxy_all %>%
  mutate(GFNARROWwalker = case_when(GFNARROWwalker == "UNKNOWN" & SPECIES_NAME == "XXXUNKNOWNCOMPOSITE" ~ "FORB",
                                    TRUE ~ GFNARROWwalker))

# Vector for target functional groups - vascular plants, leaves out moss, lichen and liverworts + abiotic stuff
fg <- c("FORB", "SEVER", "SDECI", "SHRUBU", "SHRUB", "GRAMINOIDU", "GRASS", "SEDGE", "RUSH", "GRAMU", "WOODYU", NA) # Not including UNK, UNKNOWN or OTHER

# # Generating list of sites with NA status data or grazing-related status data
# fg_test_perccov_all <- filter(perccov_all, GFNARROWwalker %in% c("UNK", NA, "OTHER", "UNKNOWN"))
# fg_test_pfplot_all <- filter(pfplot_all, GFNARROWwalker %in% c("UNK", NA, "OTHER", "UNKNOWN"))
# fg_test_pfxy_all <- filter(pfxy_all, GFNARROWwalker %in% c("UNK", NA, "OTHER", "UNKNOWN"))
# 
# fg_sites_perccov_all <- unique(fg_test_perccov_all$SITE)
# fg_sites_pfplot_all <- unique(fg_test_pfplot_all$SITE)
# fg_sites_pfxy_all <- unique(fg_test_pfxy_all$SITE)
# 
# fg_sites <- c(fg_sites_perccov_all, fg_sites_pfplot_all, fg_sites_pfxy_all)
# 
# rm(fg_test_perccov_all, fg_test_pfplot_all, fg_test_pfxy_all, fg_sites_perccov_all,
#    fg_sites_pfplot_all, fg_sites_pfxy_all)

# Filter percentage cover dataframe
perccov_all_2 <- perccov_all %>% 
  filter(STATUS %in% status, TREATMENT %in% treatment, GFNARROWwalker %in% fg) %>% 
  unite(SiteSubsitePlotYear, c("SITE", "SUBSITE", "PLOT", "YEAR"), sep = ":", remove = FALSE) %>%
  unite(SiteSubsitePlot, c("SITE", "SUBSITE", "PLOT"), sep = ":", remove = FALSE) %>%
  unite(SiteSubsite, c("SITE", "SUBSITE"), sep = ":", remove = FALSE) %>%
  mutate(ABUNDANCE = ifelse(is.na(ABUNDANCE), 0, ABUNDANCE)) %>% # STEPSTONES:TUNDRA1 and :TUNDRA2 which are NA should be 0
  filter(SiteSubsite != "SADVENT:WET_PHOTO") # Bad subsite, removing to avoid cover calculation issues

# Filter point frame (no xy coords) dataframe
pfplot_all_2 <- pfplot_all %>% 
  filter(STATUS %in% status, TREATMENT %in% treatment, GFNARROWwalker %in% fg) %>% 
  unite(SiteSubsitePlotYear, c("SITE", "SUBSITE", "PLOT", "YEAR"), sep = ":", remove = FALSE) %>%
  unite(SiteSubsitePlot, c("SITE", "SUBSITE", "PLOT"), sep = ":", remove = FALSE) %>%
  unite(SiteSubsite, c("SITE", "SUBSITE"), sep = ":", remove = FALSE) %>%
  mutate(ABUNDANCE = ifelse(is.na(ABUNDANCE), 0, ABUNDANCE)) %>% # KANGERS should be 0 (investigated below)
  dplyr::select(-COVER_UNDERSTORY) # Remove superfluous column

# Investigate those with Abundance = NA
kanger.na <- pfplot_all_2 %>%
  filter(SiteSubsite %in% c("KANGER:BASHFUL", "KANGER:DOPEY", "KANGER:SNEEZY"))
    # When it's a 1 the info is in there, but there are no 0s. It's always the same species over the years so I think these are Abundance = 0

# Remove intermediate variables
rm(kanger.na)

# Filter point frame (xy coords) dataframe
pfxy_all_2 <- pfxy_all %>% 
  filter(STATUS %in% status, TREATMENT %in% treatment, GFNARROWwalker %in% fg) %>% 
  unite(SiteSubsitePlotYear, c("SITE", "SUBSITE", "PLOT", "YEAR"), sep = ":", remove = FALSE) %>%
  unite(SiteSubsitePlot, c("SITE", "SUBSITE", "PLOT"), sep = ":", remove = FALSE) %>%
  unite(SiteSubsite, c("SITE", "SUBSITE"), sep = ":", remove = FALSE) %>%
  mutate(ABUNDANCE = ifelse(is.na(ABUNDANCE), 1, ABUNDANCE)) # There was 1 NA abundance value for QHI:HE only, should be 1 instead as no 0s are recorded

# There are some summed plots that don't belong in XY
summed <- pfxy_all_2 %>% 
  filter(HIT == "sum") %>% # Identify summed plots
  dplyr::select(-X, -Y, -HIT) # Remove xy specific columns

# Remove from the xy database
pfxy_all_2 <- pfxy_all_2 %>% 
  filter(HIT != "sum" | is.na(HIT)) # Remove summed hits

# Add to summed database
pfplot_all_2 <- rbind(pfplot_all_2, summed)

# Remove unneeded dataframe
rm(summed)


# PERCENTAGE COVER CHECKS AND EDITS (perccover_all_3) ----

# Do the cover values add up to 100?
perccov_sum <- perccov_all_2 %>%
  filter(ValueType == "percent_cover") %>% 
  group_by(SiteSubsitePlotYear) %>%
  summarise(cover_sum = sum(ABUNDANCE)) %>% 
  ungroup()
    # Quite a lot of values over 100 so they need to be made proportional too so all values are comparable

# Confirm that 1 row = 1 species
perccov_test <- perccov_all_2 %>%
  group_by(SiteSubsitePlotYear) %>% 
  mutate(NumberRows = n(),
         NumberSpecies = length(unique(SPECIES_NAME))) %>%
  mutate(SameOrNot = ifelse(NumberRows == NumberSpecies, "Same", "Different")) %>%
  ungroup()

# Check if this is because of the missing species names or actually there are repeated species names (repeats!)
perccov_difs <- filter(perccov_test, SameOrNot == "Different")

# (2) Fix dataframe for merging: Add up values per species so we end up with only one row per species
perccov_dupl_species <- perccov_difs %>%
  filter(SiteSubsitePlotYear != "BARROW:CAREX_MOIST_MEADOW_MICROTOPO:BC02.5:1999") %>% # Remove site with duplicates
  group_by(SiteSubsitePlotYear, SPECIES_NAME) %>%
  mutate(AbundanceFixed = sum(ABUNDANCE)) %>%
  ungroup() %>%
  group_by(SiteSubsitePlotYear) %>%
  distinct(SPECIES_NAME, .keep_all = TRUE) %>% # Retain only one row per species
  ungroup() %>%
  mutate(ABUNDANCE = AbundanceFixed) %>% # Renaming abundance column with fixed abundance values
  dplyr::select(-AbundanceFixed)

# (3) Fix dataframe for merging: One site has duplicate values: all records have exactly the same values twice
perccov_dupl_barrow <- perccov_difs %>%
  filter(SiteSubsitePlotYear == "BARROW:CAREX_MOIST_MEADOW_MICROTOPO:BC02.5:1999") %>% 
  distinct(SPECIES_NAME, .keep_all = TRUE) # Creates dataframe with no duplicates

# Remove the 'different' values from the overall perccover_all_2 dataframe and bind with fixed records to generate one correct perccover_all dataframe
perccov_fixed <- perccov_test %>% 
  filter(SameOrNot != "Different") %>% 
  rbind(., perccov_dupl_species, perccov_dupl_barrow) %>% 
  dplyr::select(-NumberRows, -NumberSpecies, -SameOrNot) # Remove intermediate columns

# Convert all values to relative cover
perccov_all_3 <- perccov_fixed %>%
  group_by(SiteSubsitePlotYear) %>% 
  mutate(TotalAbundance = sum(ABUNDANCE)) %>%
  mutate(RelCover = (ABUNDANCE/TotalAbundance)*100) %>% # 5733 obs
  ungroup()

# Confirm that total cover values add up to 100 in every plotXyear (TRUE)
perccov_check <- perccov_all_3 %>%
  group_by(SiteSubsitePlotYear) %>% 
  mutate(TotalCover = sum(RelCover)) %>% 
  distinct(SiteSubsitePlotYear, .keep_all = TRUE) %>% 
  ungroup()


# POINT-FRAMING (SUMMED) CHECKS AND EDITS (pfplot_all_3) ----

# Confirm that 1 row = 1 species
pfplot_test <- pfplot_all_2 %>% 
  group_by(SiteSubsitePlotYear) %>% 
  mutate(NumberRows = n(),
         NumberSpecies = length(unique(SPECIES_NAME))) %>%
  mutate(SameOrNot = ifelse(NumberRows == NumberSpecies, "Same", "Different")) %>%
  ungroup()

# Check if this is because of the missing species names or actually there are repeated species names
pfplot_difs <- pfplot_test %>%
  filter(SameOrNot == "Different")

# Duplicates: exactly the same species and values on repeat
duplicate_sites <- c("ABISKO:PEATLAND:AA1:2000", "ABISKO:PEATLAND:AA1:2002", "ABISKO:PEATLAND:AA1:2004", "ABISKO:PEATLAND:AA1:2006", "ABISKO:PEATLAND:AA1:2008",
                     "ABISKO:PEATLAND:AA2:2000", "ABISKO:PEATLAND:AA2:2002", "ABISKO:PEATLAND:AA2:2004", "ABISKO:PEATLAND:AA2:2006", "ABISKO:PEATLAND:AA2:2008",
                     "ABISKO:PEATLAND:AA3:2000", "ABISKO:PEATLAND:AA3:2002", "ABISKO:PEATLAND:AA3:2004", "ABISKO:PEATLAND:AA3:2006", "ABISKO:PEATLAND:AA3:2008",
                     "ABISKO:PEATLAND:AA4:2000", "ABISKO:PEATLAND:AA4:2002", "ABISKO:PEATLAND:AA4:2004", "ABISKO:PEATLAND:AA4:2006", "ABISKO:PEATLAND:AA4:2008",
                     "ABISKO:PEATLAND:AA5:2000", "ABISKO:PEATLAND:AA5:2002", "ABISKO:PEATLAND:AA5:2004", "ABISKO:PEATLAND:AA5:2006", "ABISKO:PEATLAND:AA5:2008")

# (1) Fix dataframe for merging: Multiple plots have duplicate values, all records have exactly the same values twice
pfplot_dupl_sites <- pfplot_difs %>%
  filter(SiteSubsitePlotYear %in% duplicate_sites) %>%
  group_by(SiteSubsitePlotYear) %>%
  distinct(SPECIES_NAME, .keep_all = TRUE) %>% # Retain only one record from each duplicate
  ungroup()

# (2) Fix dataframe for merging: Add up values per species so we end up with only one row per species
pfplot_dupl_species <- pfplot_difs %>%
  filter(SiteSubsitePlotYear %notin% duplicate_sites) %>% # Remove sites with all duplicates
  group_by(SiteSubsitePlotYear, SPECIES_NAME) %>%
  mutate(AbundanceFixed = sum(ABUNDANCE)) %>%
  ungroup() %>%
  group_by(SiteSubsitePlotYear) %>%
  distinct(SPECIES_NAME, .keep_all = TRUE) %>% # Retain only one row per species
  ungroup() %>%
  mutate(ABUNDANCE = AbundanceFixed) %>% # Renaming abundance column with fixed abundance values
  dplyr::select(-AbundanceFixed)

# Remove the 'different' values from the overall perccover_all_2 dataframe and bind with fixed records to generate one correct perccover_all dataframe
pfplot_fixed <- pfplot_test %>% 
  filter(SameOrNot != "Different") %>% 
  rbind(., pfplot_dupl_sites, pfplot_dupl_species) %>% 
  dplyr::select(-NumberRows, -NumberSpecies, -SameOrNot) # Remove intermediate columns

# Convert all values to relative cover
pfplot_all_3 <- pfplot_fixed %>%
  group_by(SiteSubsitePlotYear) %>% 
  mutate(TotalAbundance = sum(ABUNDANCE)) %>%
  mutate(RelCover = (ABUNDANCE/TotalAbundance)*100) %>% # 15827 obs
  ungroup()

# Confirm that total cover values add up to 100 in every plotXyear
pfplot_check <- pfplot_all_3 %>%
  group_by(SiteSubsitePlotYear) %>% 
  mutate(TotalCover = sum(RelCover)) %>% 
  ungroup() %>% 
  distinct(SiteSubsitePlotYear, .keep_all = TRUE)


# POINT-FRAMING (XY) CHECKS AND EDITS (pfxy_all_3) ----

# There are XY data that only have one coordinate - filling in the NA cell so we avoid problems with cover calculation
pfxy_na_x <- pfxy_all_2 %>% 
  filter(is.na(X)) # No NAs

pfxy_na_y <- pfxy_all_2 %>% 
  filter(is.na(Y)) # 138814 NAs

# Replace Y coords that are NA by 0s, create a unique coordinate
pfxy_XY <- pfxy_all_2 %>%
  mutate(Y = ifelse(is.na(Y), 0, Y)) %>% 
  unite(XY, c("X", "Y"), sep = "_", remove = FALSE)

# Remove LATNJA:CAREX subsite as known to be duplicate
pfxy_latnja_carex <- pfxy_XY %>% 
  filter(SiteSubsite != "LATNJA:CAREX")

# STEP 1: Convert species abundance to presence/absence (2-D, ignoring multiple hits of same species at each xy coord, just 1)
pfxy_presence_absence <- pfxy_latnja_carex %>% 
  group_by(SiteSubsitePlotYear, XY) %>%
  distinct(SPECIES_NAME, .keep_all = TRUE) %>% 
  mutate(ABUNDANCE = ifelse(ABUNDANCE > 1, 1, ABUNDANCE)) %>%
  ungroup()

# STEP 2: Calculate unique species hits per plot and total unique species hits per plot
pfxy_presence_absence_2 <- pfxy_presence_absence %>% 
  group_by(SiteSubsitePlotYear, SPECIES_NAME) %>% 
  mutate(UniqueSpHitsPlot = n()) %>% 
  distinct(SiteSubsitePlotYear, SPECIES_NAME, .keep_all = TRUE) %>% 
  ungroup() %>%
  dplyr::select(-X, -Y, -XY, -HIT) %>% 
  group_by(SiteSubsitePlotYear) %>%
  mutate(TotalUniqueSpHitsPlot = sum(UniqueSpHitsPlot)) %>%
  ungroup()

# STEP 3: Calculate cover per species
pfxy_all_3 <- pfxy_presence_absence_2 %>% 
  mutate(RelCover = (UniqueSpHitsPlot / TotalUniqueSpHitsPlot) * 100) #56797

# Confirm that total cover values add up to 100 in every plotXyear
pfxy_check <- pfxy_all_3 %>%
  group_by(SiteSubsitePlotYear) %>% 
  mutate(TotalCover = sum(RelCover)) %>% 
  distinct(SiteSubsitePlotYear, .keep_all = TRUE) %>% 
  ungroup()


# COMBINE THE DATAFRAMES (itex.1) ----

# Ensure the columns match between the dataframes
perccov_all_F <- dplyr::select(perccov_all_3, -TotalAbundance, -ABUNDANCE)
pfplot_all_F <- dplyr::select(pfplot_all_3, -TotalAbundance, -ABUNDANCE)
pfxy_all_F <- dplyr::select(pfxy_all_3, -ABUNDANCE, -UniqueSpHitsPlot, -TotalUniqueSpHitsPlot)

# Combine into one ITEX dataframe & remove 'extra' plots with 0 cover (plots only have 0 or 100)
itex.1 <- rbind(perccov_all_F, pfplot_all_F, pfxy_all_F) %>% # 78357
  mutate(RelCover = if_else(is.nan(RelCover), 0, RelCover)) %>% # It's just 0/0
  group_by(SiteSubsitePlotYear) %>% 
  mutate(COVER = sum(RelCover)) %>% 
  ungroup() %>% 
  filter(COVER == 100) %>% 
  dplyr::select(-COVER)


# ITEX METADATA (itex.2) ----

# Metadata file edited manually so the subsites that didn't match up do (mainly just adding extra site name to subsite column so they match)
# itex.metadata.input <- read.csv("scripts/mgarciacriado/data/new_itex/TVC_SITE_SUBSITE_UPDATED2020_edited.csv")
itex.metadata.input <- read.csv("data/itex_june2022/input_itex_metadata.csv")

# Retain only relevant columns
itex.metadata <- itex.metadata.input %>%
  unite(SiteSubsite, c("SITE", "SUBSITE"), sep = ":", remove = FALSE) %>%
  dplyr::select(SiteSubsite, Moisture.category, LAT, LONG, ELEV, AZONE, SurveyedArea, PlotSize_m2) %>%
  rename(MOISTURE = Moisture.category)

# Standardise the names of LATNJA sites in the full ITEX database
itex.latnja.rename <- itex.1 %>% 
  mutate(SITE = str_replace_all(SITE, pattern = "LATNJAJAURE", replacement = "LATNJA"), # Standardising LATNJA names
         SiteSubsite = str_replace_all(SiteSubsite, pattern = "LATNJAJAURE", replacement = "LATNJA"),
         SiteSubsitePlot = str_replace_all(SiteSubsitePlot, pattern = "LATNJAJAURE", replacement = "LATNJA"),
         SiteSubsitePlotYear = str_replace_all(SiteSubsitePlotYear, pattern = "LATNJAJAURE", replacement = "LATNJA"))

# Remove Latnja mesic_meadow as likely inconsistent methods over time (CAREX removed in previous section as duplicate with this subsite also)
itex.latnja.remove <- itex.latnja.rename %>% 
  filter(SiteSubsite != "LATNJA:MESIC_MEADOW")

# Merge with composition data
itex.2 <- left_join(itex.latnja.remove, itex.metadata, by = "SiteSubsite")


# SITE CHECKS (itex.3) ----

# Identify all subsites not contiguously joined to the Arctic (below 60 degN)
nonarctic.sites <- itex.2 %>% 
  filter(LAT < 60) %>% 
  dplyr::select(SITE) %>% 
  unique()

# Create vector of these sites
nonarctic.sites <- nonarctic.sites$SITE

# Plots that had inconsistent surveyed areas or didn't identify species in the whole subsite
incorrect.subsites <- c("SADVENT:WET_PHOTO", "SADVENT:MES_PHOTO", "ALEXFIORD:LEVDOLOMITE",
                        "ALEXFIORD:LEVGRANITE", "SVERDRUP:SVERDRUP")

# Spot Check: there is one Alexfiord "total" plot?
itex.cass.alex <- filter(itex.2, SiteSubsite == "ALEXFIORD:CASSIOPE_COVER")
  # The only cassiope plot, doesn't appear to be summary of others, contains other species, probably fine to include

# Remove years from plot names (impacts change over time analyses when grouping)
itex.incorrect.plot.names <- itex.2 %>% 
  mutate(Year_in_plot = ifelse(str_detect(PLOT, pattern = as.character(YEAR)), "YES", "NO")) %>% # Detects if year is in plot name
  mutate(PLOT = ifelse(Year_in_plot == "YES", str_remove(PLOT, pattern = as.character(YEAR)), PLOT)) %>% # Removes year string from plot name
  mutate(PLOT = str_remove(PLOT, pattern = "^_+"), # Removes leading "_"
         PLOT = str_remove(PLOT, pattern = "_$")) %>%  # Removes trailing "_"
  unite(SiteSubsitePlotYear, c("SITE", "SUBSITE", "PLOT", "YEAR"), sep = ":", remove = FALSE) %>% # Replace as plot columns changed
  unite(SiteSubsitePlot, c("SITE", "SUBSITE", "PLOT"), sep = ":", remove = FALSE) %>% 
  dplyr::select(-Year_in_plot)

# Remove inconsistent sites
itex.3 <- itex.incorrect.plot.names %>% 
  filter(SITE %notin% nonarctic.sites,
         SiteSubsite %notin% incorrect.subsites) 


# SURVEYED AREA CHECKS (itex.4) ----

# Find non-numeric values for SurveyedArea
unique(itex.3$SurveyedArea)

# Create a vector of non-numeric SurveyedAreas
non.numeric.area <- c("", "Email", "mixed", "0.1089 or 0.2178")

# Get list of sites with non-numeric area
non.numeric.area.cut <- filter(itex.3, SurveyedArea %in% non.numeric.area) %>% arrange(SiteSubsite)
non.numeric.area.sites <- unique(non.numeric.area.cut$SiteSubsite)

# Replace missing Auyuittuq with 1 (as found in PlotSize column), and the doubled up Jameson Land site as per metadata calculations
itex.area <- itex.3 %>% 
  mutate(SurveyedArea = ifelse(SiteSubsite %in% c("AUYUITTUQ:OWL RIVER") & SurveyedArea %in% non.numeric.area, 1, SurveyedArea),
         SurveyedArea = ifelse(SiteSubsite %in% c("JAMESONLAND:TYSKIT"), 0.2178, SurveyedArea),
         SurveyedArea = ifelse(SiteSubsite %in% c("JAMESONLAND:TYSKIT") & str_detect(PLOT, paste(c("13", "CASEMP", "BBN"),collapse = '|')), 0.1089, SurveyedArea))
                                 
# Filter out areas with non-numeric or missing SurveyedAreas and retain plots 1m2 or less
itex.4 <- itex.area %>% 
  filter(SurveyedArea %notin% non.numeric.area) %>% # Remove non-numeric surveyed areas
  mutate(SurveyedArea = as.numeric(SurveyedArea)) %>% # Convert to numeric
  filter(SurveyedArea <= 1 | is.na(SurveyedArea)) # Retain plots only 1m2 or less (and the NAs)

# Check remaining values for SurveyedArea
unique(itex.4$SurveyedArea)


# MOISTURE CHECKS (itex.5) ----

# Did some subsites not get metadata?
unique(itex.4$MOISTURE) # "" but no NAs

# Replace missing MOISTURE information with NAs
itex.5 <- mutate(itex.4, MOISTURE = ifelse(MOISTURE == "", NA, MOISTURE))


# LAT & LONG CHECKS (itex.6) ----

# Unique values for LAT and LONG
unique(itex.5$LAT) # NAs
unique(itex.5$LONG) # NAs

# Identify sites with missing LAT information
missing.LAT <- filter(itex.5, is.na(LAT))
missing.LAT <- unique(missing.LAT$SiteSubsite)

# Identify sites with missing LONG information
missing.LONG <- filter(itex.5, is.na(LONG))
missing.LONG <- unique(missing.LONG$SiteSubsite)

# Missing LAT and LONG: "DISKO:DRYHEATH_FLUX" "DISKO:WETFEN_FLUX"
    # Missing info is not in the previous version of ITEX as these are new sites

# Identify sites known to have LAT and LONG that fall in the ocean (identified during climate extraction)
water.LAT.LONG <- itex.5 %>% 
  dplyr::select(SiteSubsite, YEAR, LAT, LONG) %>% 
  filter(LAT %in% c(74.28, 74.29)) %>% 
  distinct(SiteSubsite, .keep_all = TRUE)

water.LAT.LONG <- unique(water.LAT.LONG$SiteSubsite)

  # Incorrect LAT and LONG: 8 x "ZACKENGERG" SITES

# Manually input new LAT and LONG info for those missing or in water (all above 60 so not removed for being 'non-Arctic')
itex.6 <- itex.5 %>% 
  mutate(LAT = ifelse(SiteSubsite == "DISKO:DRYHEATH_FLUX", 69.27, LAT),
         LAT = ifelse(SiteSubsite == "DISKO:WETFEN_FLUX", 69.43, LAT),
         LAT = ifelse(SiteSubsite %in% water.LAT.LONG, 74.47427, LAT),
         LONG = ifelse(SiteSubsite == "DISKO:DRYHEATH_FLUX", -53.45, LONG),
         LONG = ifelse(SiteSubsite == "DISKO:WETFEN_FLUX", -53.78, LONG),
         LONG = ifelse(SiteSubsite %in% water.LAT.LONG,-20.52895, LONG),
         # Info. taken from Wikipedia - find a better source!
          # https://en.wikipedia.org/wiki/International_Tundra_Experiment#Disko,_Greenland
         LAT = ifelse(SiteSubsite == "NARSARSUAQ:HIGH_ELEVATION", 61.16, LAT),
         LONG = ifelse(SiteSubsite == "NARSARSUAQ:HIGH_ELEVATION", -45.40, LONG))
        # Info. taken from paper
          # https://cdnsciencepub.com/doi/10.1139/AS-2020-0041


# GEOGRAPHICAL DATA (itex.7) ----

  # Assign all sites to a Region dependent on their continent and glacial history during the LGM (~18,000 B.P.)

# Create a dataframe of coordinate pairs with dataset_sitename attached
coords.itex <- itex.6 %>% 
  dplyr::select(LONG, LAT) %>% 
  distinct(LONG, LAT, .keep_all = TRUE)

# Get country names for each pair of points
coords.country <- coords2country(coords.itex)

# Get continent names for each pair of points
coords.continent <- coords2continent(coords.itex)

# Join country and continent names to coordinate pairs dataframe
coords.c.c <- data.frame(coords.itex, coords.country, coords.continent) %>% 
  rename(Country = coords.country, 
         Continent = coords.continent) %>% 
  mutate(Country = as.character(Country),
         Continent = as.character(Continent))

# Fill in the missing country and continent data based on the LAT and LONG coords
coords.c.c.fix <- coords.c.c %>% 
  mutate(Country = case_when(LAT == 69.27000 & LONG == -53.450000 ~ "Greenland",
                             LAT == 69.57488 & LONG == -138.863470 ~ "Canada",
                             LAT == 69.57646 & LONG == -138.867740 ~ "Canada",
                             LAT == 69.38900 & LONG == -81.530000 ~ "Canada",
                             TRUE ~ Country),
         Continent = case_when(LAT == 69.27000 & LONG == -53.450000 ~ "Greenland",
                               LAT == 69.57488 & LONG == -138.863470 ~ "North America",
                               LAT == 69.57646 & LONG == -138.867740 ~ "North America",
                               LAT == 69.38900 & LONG == -81.530000 ~ "North America",
                               TRUE ~  Continent)) %>% 
  mutate(Continent = ifelse(Country == "Greenland", "Greenland", Continent))

# Create spatial points object of the unique coordinate pairs
lgm.coords <- coords.itex %>% 
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

# Write to file for map later
write.csv(coords.c.c.lgm, "data/output_01_coords_lgm.csv")

# Count how many locations were and weren't glaciated by continent
lgm.count <- coords.c.c.lgm %>% 
  group_by(Continent, glaciated) %>% 
  summarise(count = length(LAT)) %>% 
  ungroup() %>% 
  mutate(Continent = ifelse(Continent == "Europe", "Eurasia", Continent)) %>% 
  add_row(Continent = "Greenland", glaciated = "No", count = 0) %>% # No non-glaciated records from Greenland
  add_row(Continent = "Eurasia", glaciated = "No", count = 0) %>% # No non-glaciated records from Eurasia
  arrange(Continent)

# Create 'Region' column based on the continent and whether glaciated or not
coords.region <- coords.c.c.lgm %>% 
  mutate(Region = case_when(Continent %in% c("Europe", "Asia") & Country != "Iceland" ~ "Eurasia",
                            Continent == "Greenland" | Country == "Iceland" ~ "GreenIceland",
                            Continent == "North America" & glaciated == "Yes" ~ "North America-East",
                            Continent == "North America" & glaciated == "No" ~ "North America-West")) %>%
  dplyr::select(-c(Country, Continent, glaciated)) # No longer need country, continent or glaciated records


# Join Region information to the overall combined traits dataframe and add additional geographic information
itex.geographic <- left_join(itex.6, coords.region, by = c("LAT" = "LAT", "LONG" = "LONG")) %>% 
  mutate(AlpArc = ifelse(LAT >= 66.4 & ELEV >= 1000, "Arctic-Alpine",
                         ifelse(LAT >= 66.4, "Arctic",
                                ifelse(ELEV >= 1000, "Alpine", "Subarctic")))) %>% 
  mutate(lat_grid = plyr::round_any(LAT, 0.5, f = floor)) %>% # Gridcell information
  mutate(lon_grid = ifelse(LONG >0, plyr::round_any(LONG, 0.5, f = floor), 
                           plyr::round_any(LONG, 0.5, f = ceiling))) %>%
  mutate(gridcell = paste0(lat_grid, "_", lon_grid))
      # Three subarctic alpine sites retained as close to the Arctic

# Check all ITEX records now assigned to a Region and AlpArc
unique(itex.geographic$Region) # No NAs
unique(itex.geographic$AlpArc)

# Determine NA alparc
alparc.na <- filter(itex.geographic, is.na(AlpArc))
alparc.na.sites <- unique(alparc.na$SITE)
  # NAs at nine sites

# Modify the missing AlpArc records to finalise geographic information
itex.7 <- itex.geographic %>% 
  mutate(AlpArc = ifelse(SITE %in% c("DISKO", "NUUK", "NAKKALA", "BILLEFJORDEN",
                                     "PYRAMIDEN", "ADVENT", "JAMESONLAND", "IGLOOLIK") & is.na(AlpArc), "Arctic", AlpArc),
         AlpArc = ifelse(SITE %in% c("NARSARSUAQ") & is.na(AlpArc), "Subarctic", AlpArc))


# MORPHOSPECIES CHECKS (itex.7) ----

# Create a vector of unique itex species
itex.7.species <- unique(itex.7$SPECIES_NAME)
  # One NA species (removed below ITEX 8)
  # No missing species names

# Remove the plot(s) with NA (or algae) species in them (have to remove whole plot so cover = 100)
itex.7.na.species <- itex.7 %>% 
  mutate(NA_SPECIES = ifelse(is.na(SPECIES_NAME) | SPECIES_NAME == "XXXALGAE", 1, 0)) %>% # Assigns value of 1 to NA or algae species so total for plot isn't 0
  group_by(SiteSubsitePlotYear) %>% 
  mutate(PlotSumNA = sum(NA_SPECIES)) %>% 
  ungroup() %>% 
  filter(PlotSumNA == 0) %>% # Removed one plot
  dplyr::select(-c(NA_SPECIES, PlotSumNA))

# Standardise subspecies
morphospecies.1 <- itex.7.na.species %>% 
  mutate(SPECIES_NAME = trimws(SPECIES_NAME)) %>% 
  mutate(SPECIES_NAME = ifelse(SPECIES_NAME == "Ledum palustre subsp. groenlandicum", "Ledum palustre", SPECIES_NAME),
         SPECIES_NAME = ifelse(SPECIES_NAME == "Tephroseris integrifolia subsp. atropurpurea", "Tephroseris integrifolia", SPECIES_NAME),
         SPECIES_NAME = ifelse(SPECIES_NAME == "Silene uralensis subsp. apetala", "Silene uralensis", SPECIES_NAME),
         SPECIES_NAME = ifelse(SPECIES_NAME == "Luzula spicata/confusa", "XXXLuzula", SPECIES_NAME),
         SPECIES_NAME = ifelse(SPECIES_NAME == "Eriophorum scheuchzeri/chamissonis", "XXXEriophorum", SPECIES_NAME),
         SPECIES_NAME = ifelse(SPECIES_NAME == "Cardamine bellidifolia subsp. alpina", "Cardamine bellidifolia", SPECIES_NAME),
         SPECIES_NAME = ifelse(SPECIES_NAME == "Carex aquatilis var. minor", "Carex aquatilis", SPECIES_NAME),
         SPECIES_NAME = ifelse(SPECIES_NAME == "Eriophorum angustifolium subsp. triste", "Eriophorum angustifolium", SPECIES_NAME),
         SPECIES_NAME = ifelse(SPECIES_NAME == "Empetrum nigrum subsp. hermaphroditum", "Empetrum nigrum", SPECIES_NAME),
         SPECIES_NAME = ifelse(SPECIES_NAME == "Dryas integrifolia x octopetala", "XXXDryas", SPECIES_NAME),
         SPECIES_NAME = ifelse(SPECIES_NAME == "Betula pubescens var. pumila", "Betula pubescens", SPECIES_NAME),
         SPECIES_NAME = ifelse(SPECIES_NAME == "Rumex alpestris subsp. lapponicus", "Rumex alpestris", SPECIES_NAME),
         SPECIES_NAME = ifelse(SPECIES_NAME == "Salix arctica/arctophila", "XXXSalix", SPECIES_NAME),
         SPECIES_NAME = ifelse(SPECIES_NAME == "Salix lanata x reticulata", "XXXSalix", SPECIES_NAME),
         SPECIES_NAME = ifelse(SPECIES_NAME == "Anthoxanthum odoratum subsp. nipponicum", "Anthoxanthum odoratum", SPECIES_NAME),
         SPECIES_NAME = ifelse(SPECIES_NAME == "Carex microchaeta&rupestris", "XXXCarex", SPECIES_NAME),
         SPECIES_NAME = ifelse(SPECIES_NAME == "Deschampsia flexuosa/Juncus trifidus", "XXXGraminoid", SPECIES_NAME),
         SPECIES_NAME = ifelse(SPECIES_NAME == "Empetrum nigrum/Empetrum nigrum subsp. hermaphroditum", "Empetrum nigrum", SPECIES_NAME),
         SPECIES_NAME = ifelse(SPECIES_NAME == "Empetrum nigrum/Phyllodoce caerulea", "XXXEShrub", SPECIES_NAME))

# Identify SPECIES_NAMEs that are actually genus/family names (often only one word)
non.species <- morphospecies.1 %>% 
  dplyr::select(SPECIES_NAME) %>% 
  unique() %>% 
  arrange(SPECIES_NAME)

# Create a vector of SPECIES_NAMEs that are actually genus/family/gf names (NOTE the spaces)
g.f.morphospecies <- c("Alchemilla", "Anemone", "Antennaria", "Arnica", "Astragalus",
                       "Calamagrostis", "Cardamine", "Carex", "Deschampsia", "Cyperaceae",
                       "Deschampsia", "Draba", "Dryas", "Epilobium", "Eriophorum",
                       "Festuca", "Galium", "Gentiana", "Hepatica", "Juncus NA", "Luzula",
                       "Minuartia NA", "Oxytropis", "Pedicularis", "Petasites", "Poa", "Poaceae",
                       "Polemonium", "Polygonum", "Sagina NA", "Salix", "Saxifraga", "Senecio",
                       "Stellaria", "Tofieldia", "Viola", "Graminoid unknown")

# Transform into morphospecies
morphospecies.2 <- morphospecies.1 %>% 
  mutate(SPECIES_REF = SPECIES_NAME) %>% 
  mutate(SPECIES_NAME = ifelse(SPECIES_REF %in% g.f.morphospecies, str_remove_all(string = SPECIES_NAME, pattern = " "), SPECIES_NAME)) %>%
  mutate(SPECIES_NAME = ifelse(SPECIES_REF %in% g.f.morphospecies, str_remove_all(string = SPECIES_NAME, pattern = "NA"), SPECIES_NAME)) %>%
  mutate(SPECIES_NAME = ifelse(SPECIES_REF %in% g.f.morphospecies, str_remove_all(string = SPECIES_NAME, pattern = "unknown"), SPECIES_NAME)) %>%
  mutate(SPECIES_NAME = ifelse(SPECIES_REF %in% g.f.morphospecies, paste0("XXX", toupper(SPECIES_NAME), ":", SITE), SPECIES_NAME)) %>% 
  mutate(SPECIES_NAME = ifelse(str_detect(SPECIES_NAME, pattern = "XXX"), toupper(SPECIES_NAME), SPECIES_NAME)) %>% # Convert names to upper case for standardisation
  dplyr::select(-SPECIES_REF)

# Check the morphospecies in the dataset
morphospecies.only <- morphospecies.2 %>% 
  filter(str_detect(SPECIES_NAME, pattern = "XXX")) %>% 
  dplyr::select(SPECIES_NAME, SITE) %>% 
  unique()

# Apply manual corrections to incorrectly named morphospecies (at all levels)
itex.8 <- morphospecies.2 %>% 
  mutate(SPECIES_NAME = ifelse(SPECIES_NAME == "XXXASTERACEAE", paste0(SPECIES_NAME, ":", SITE), SPECIES_NAME),
         SPECIES_NAME = ifelse(SPECIES_NAME == "XXXCAREX", paste0(SPECIES_NAME, ":", SITE), SPECIES_NAME),
         SPECIES_NAME = ifelse(SPECIES_NAME == "XXXCAREX:SWEDEN", paste0("XXXCAREX", ":", SITE), SPECIES_NAME),
         SPECIES_NAME = ifelse(SPECIES_NAME == "XXXCOMP", paste0("XXXFORB", ":", SITE), SPECIES_NAME),
         SPECIES_NAME = ifelse(SPECIES_NAME == "XXXCOPTIDIUM:ADVENT2", "XXXCOPTIDIUM:ADVENT", SPECIES_NAME),
         SPECIES_NAME = ifelse(SPECIES_NAME == "XXXDRABA", paste0(SPECIES_NAME, ":", SITE), SPECIES_NAME),
         SPECIES_NAME = ifelse(SPECIES_NAME == "XXXDRABA:ADVENT2", "XXXDRABA:ADVENT", SPECIES_NAME),
         SPECIES_NAME = ifelse(SPECIES_NAME == "XXXDRABA:PYR_AWS" & SITE == "BILLEFJORDEN", "XXXDRABA:BILLEFJORDEN", SPECIES_NAME),
         SPECIES_NAME = ifelse(SPECIES_NAME == "XXXDRABA:PYR_AWS" & SITE == "PYRAMIDEN", "XXXDRABA:PYRAMIDEN", SPECIES_NAME),
         SPECIES_NAME = ifelse(SPECIES_NAME == "XXXDRYAS", paste0(SPECIES_NAME, ":", SITE), SPECIES_NAME),
         SPECIES_NAME = ifelse(SPECIES_NAME == "XXXEQUISETUM", paste0(SPECIES_NAME, ":", SITE), SPECIES_NAME),
         SPECIES_NAME = ifelse(SPECIES_NAME == "XXXERIOPHORUM", paste0(SPECIES_NAME, ":", SITE), SPECIES_NAME),
         SPECIES_NAME = ifelse(SPECIES_NAME == "XXXESHRUB", paste0(SPECIES_NAME, ":", SITE), SPECIES_NAME),
         SPECIES_NAME = ifelse(SPECIES_NAME == "XXXFORB", paste0(SPECIES_NAME, ":", SITE), SPECIES_NAME),
         SPECIES_NAME = ifelse(SPECIES_NAME == "XXXGRAM", paste0("XXXGRAMINOID", ":", SITE), SPECIES_NAME),
         SPECIES_NAME = ifelse(SPECIES_NAME == "XXXGRAMINOID", paste0(SPECIES_NAME, ":", SITE), SPECIES_NAME),
         SPECIES_NAME = ifelse(SPECIES_NAME == "XXXGRAMINOID UNKNOWN:THINGVELLIR", "XXXGRAMINOID:THINGVELLIR", SPECIES_NAME),
         SPECIES_NAME = ifelse(SPECIES_NAME == "XXXGRAMINU", paste0("XXXGRAMINOID", ":", SITE), SPECIES_NAME),
         SPECIES_NAME = ifelse(SPECIES_NAME == "XXXGRASS", paste0("XXXGRAMINOID", ":", SITE), SPECIES_NAME),
         SPECIES_NAME = ifelse(SPECIES_NAME == "XXXGRASSUNK", paste0("XXXGRAMINOID", ":", SITE), SPECIES_NAME),
         SPECIES_NAME = ifelse(SPECIES_NAME == "XXXHERB", paste0("XXXFORB", ":", SITE), SPECIES_NAME),
         SPECIES_NAME = ifelse(str_detect(SPECIES_NAME, pattern = "XXXHIERASPP:"), paste0("XXXHIERACIUM", ":", SITE), SPECIES_NAME),
         SPECIES_NAME = ifelse(SPECIES_NAME == "XXXJUNCUS NA:THUFUVER", "XXXJUNCUS:THUFUVER", SPECIES_NAME),
         SPECIES_NAME = ifelse(SPECIES_NAME == "XXXKOBRESIA", paste0(SPECIES_NAME, ":", SITE), SPECIES_NAME),
         SPECIES_NAME = ifelse(SPECIES_NAME == "XXXLEGSPP", paste0("XXXFORB", ":", SITE), SPECIES_NAME),
         SPECIES_NAME = ifelse(SPECIES_NAME == "XXXLUZULA", paste0(SPECIES_NAME, ":", SITE), SPECIES_NAME),
         SPECIES_NAME = ifelse(SPECIES_NAME == "XXXMINUARTIA NA:THUFUVER", "XXXMINUARTIA:THUFUVER", SPECIES_NAME),
         SPECIES_NAME = ifelse(SPECIES_NAME == "XXXOTHERASTER", paste0("XXXASTERACEAE", ":", SITE), SPECIES_NAME),
         SPECIES_NAME = ifelse(SPECIES_NAME == "XXXOTHERFORB", paste0("XXXFORB", ":", SITE), SPECIES_NAME),
         SPECIES_NAME = ifelse(SPECIES_NAME == "XXXOTHERGRAM", paste0("XXXGRAMINOID", ":", SITE), SPECIES_NAME),
         SPECIES_NAME = ifelse(SPECIES_NAME == "XXXOTHERHERB", paste0("XXXFORB", ":", SITE), SPECIES_NAME),
         SPECIES_NAME = ifelse(SPECIES_NAME == "XXXOXYTROPIS", paste0(SPECIES_NAME, ":", SITE), SPECIES_NAME),
         SPECIES_NAME = ifelse(SPECIES_NAME == "XXXPED", paste0("XXXPEDICULARIS", ":", SITE), SPECIES_NAME),
         SPECIES_NAME = ifelse(SPECIES_NAME == "XXXPEDICULARIS:ADVENT2", "XXXPEDICULARIS:ADVENT", SPECIES_NAME),
         SPECIES_NAME = ifelse(SPECIES_NAME == "XXXPOA", paste0(SPECIES_NAME, ":", SITE), SPECIES_NAME),
         SPECIES_NAME = ifelse(SPECIES_NAME == "XXXPOA:PYR_AWS" & SITE == "BILLEFJORDEN", "XXXPOA:BILLEFJORDEN", SPECIES_NAME),
         SPECIES_NAME = ifelse(SPECIES_NAME == "XXXPOA:PYR_AWS" & SITE == "PYRAMIDEN", "XXXPOA:PYRAMIDEN", SPECIES_NAME),
         SPECIES_NAME = ifelse(SPECIES_NAME == "XXXRANUNCULACEAE", paste0(SPECIES_NAME, ":", SITE), SPECIES_NAME),
         SPECIES_NAME = ifelse(SPECIES_NAME == "XXXRANUNCULUS:ADVENT2", "XXXRANUNCULUS:ADVENT", SPECIES_NAME),
         SPECIES_NAME = ifelse(SPECIES_NAME == "XXXRUSH", paste0("XXXGRAMINOID", ":", SITE), SPECIES_NAME),
         SPECIES_NAME = ifelse(SPECIES_NAME == "XXXSAGINA NA:THUFUVER", "XXXSAGINA:THUFUVER", SPECIES_NAME),
         SPECIES_NAME = ifelse(SPECIES_NAME == "XXXSAGINASPP:HOLTAVORDUHEIDI", "XXXSAGINA:HOLTAVORDUHEIDI", SPECIES_NAME),
         SPECIES_NAME = ifelse(SPECIES_NAME == "XXXSALIX", paste0(SPECIES_NAME, ":", SITE), SPECIES_NAME),
         SPECIES_NAME = ifelse(SPECIES_NAME == "XXXSAXIFRAGA", paste0(SPECIES_NAME, ":", SITE), SPECIES_NAME),
         SPECIES_NAME = ifelse(SPECIES_NAME == "XXXSAXSPP:HOLTAVORDUHEIDI", "XXXSAXIFRAGA:HOLTAVORDUHEIDI", SPECIES_NAME),
         SPECIES_NAME = ifelse(SPECIES_NAME == "XXXSAXSPP:OXNADALSHEIDI", "XXXSAXIFRAGA:OXNADALSHEIDI", SPECIES_NAME),
         SPECIES_NAME = ifelse(SPECIES_NAME == "XXXSEDGE", paste0(SPECIES_NAME, ":", SITE), SPECIES_NAME),
         SPECIES_NAME = ifelse(SPECIES_NAME == "XXXTARAXACUM:SWEDEN", paste0("XXXTARAXACUM", ":", SITE), SPECIES_NAME),
         SPECIES_NAME = ifelse(str_detect(SPECIES_NAME, pattern = "XXXTARAXSPP:"), paste0("XXXTARAXACUM", ":", SITE), SPECIES_NAME),
         SPECIES_NAME = ifelse(SPECIES_NAME == "XXXUNIFOR", paste0("XXXFORB", ":", SITE), SPECIES_NAME),
         SPECIES_NAME = ifelse(SPECIES_NAME == "XXXUNKDIC", paste0("XXXFORB", ":", SITE), SPECIES_NAME),
         SPECIES_NAME = ifelse(SPECIES_NAME == "XXXUNKFOR", paste0("XXXFORB", ":", SITE), SPECIES_NAME),
         SPECIES_NAME = ifelse(SPECIES_NAME == "XXXUNKGRA", paste0("XXXGRAMINOID", ":", SITE), SPECIES_NAME),
         SPECIES_NAME = ifelse(SPECIES_NAME == "XXXUNKGRAM", paste0("XXXGRAMINOID", ":", SITE), SPECIES_NAME),
         SPECIES_NAME = ifelse(SPECIES_NAME == "XXXUNKOXYTROPIS", paste0("XXXOXYTROPIS", ":", SITE), SPECIES_NAME),
         SPECIES_NAME = ifelse(SPECIES_NAME == "XXXWOODYDRYAS:SADVENT", "XXXDRYAS:SADVENT", SPECIES_NAME),
         SPECIES_NAME = ifelse(SPECIES_NAME == "XXXWOODYSALIX:SADVENT", "XXXSALIX:SADVENT", SPECIES_NAME),
         SPECIES_NAME = ifelse(SPECIES_NAME == "XXXORCHID:MODRUVELLIR", "XXXORCHIDACEAE:MODRUVELLIR", SPECIES_NAME))

# Check the corrections applied
morphospecies.check <-  itex.8 %>% 
  filter(str_detect(SPECIES_NAME, pattern = "XXX")) %>% 
  dplyr::select(SPECIES_NAME, SITE) %>% 
  unique()


# TAXONOMY CHECK & INCORRECT SPECIES ID (itex.9) ----

# Create vector of unique species in ITEX dataset
species.itex.8 <- unique(itex.8$SPECIES_NAME)

# # Run taxonomy checker on full dataset to check species names
# species.itex.checked <- TPL(species.itex.8)
# 
# # Write csv of taxonomy checker results
# write.csv(species.itex.checked, "data/output_01_tpl_species_itex.csv")

# Read csv of taxonomy checker results
species.itex.checked <- read.csv("data/output_01_tpl_species_itex.csv")

# Create a trimmed dataframe of the checked species
species.itex.checked.2 <- species.itex.checked %>%
  dplyr::select(Taxon, New.Genus, New.Species, New.Taxonomic.status, Family) %>%
  mutate(Name_TPL = paste(New.Genus, New.Species, sep = " ")) %>%
  relocate(Name_TPL, New.Genus, Family, .before = New.Taxonomic.status) %>%
  dplyr::select(-New.Species) %>% 
  rename(New.Species = Name_TPL, New.Family = Family) %>% 
  arrange(Taxon)

# Edit the XXX species in the dataframe
species.itex.checked.3 <- species.itex.checked.2 %>%
  mutate(New.Family = ifelse(str_detect(Taxon, pattern = "XXX"), NA, New.Family),
         New.Genus = ifelse(str_detect(Taxon, pattern = "XXX"), NA, New.Genus),
         New.Species = ifelse(str_detect(Taxon, pattern = "XXX"), Taxon, New.Species),
         New.Taxonomic.status = ifelse(str_detect(Taxon, pattern = "XXX"), "Morphospecies", New.Taxonomic.status))

# Records to amend (Blank and NA status): accept the unresolved records as should be consistent throughout TPL checks
species.itex.ammend <- species.itex.checked.3 %>% 
  filter(New.Taxonomic.status %in% c("", NA))

# Manually fix the broken records in the TPL dataframe (5)
species.itex.checked.4 <- species.itex.checked.3 %>% 
  mutate(New.Species = ifelse(Taxon %in% c("Aconitum delphiniifolium"), "Aconitum delphinifolium", New.Species), # Two spellings of same plant (TPL version)
         New.Family = ifelse(Taxon %in% c("Aconitum delphiniifolium"), "Ranunculaceae", New.Family),
         New.Taxonomic.status = ifelse(Taxon %in% c("Aconitum delphiniifolium"), "Fixed", New.Taxonomic.status),
         New.Species = ifelse(Taxon %in% c("Cardamine digitalis"), "Cardamine digitata", New.Species), # Mixes two genuses, only found at QHI (see below)
         New.Genus = ifelse(Taxon %in% c("Cardamine digitalis"), "Cardamine", New.Genus),
         New.Family = ifelse(Taxon %in% c("Cardamine digitalis"), "Brassicaceae", New.Family),
         New.Taxonomic.status = ifelse(Taxon %in% c("Cardamine digitalis"), "Morphospecies", New.Taxonomic.status),
         New.Species = ifelse(Taxon %in% c("Chamaeorchis alpina"), "Chamorchis alpina", New.Species), # Spelling error with genus
         New.Genus = ifelse(Taxon %in% c("Chamaeorchis alpina"), "Chamorchis", New.Genus),
         New.Family = ifelse(Taxon %in% c("Chamaeorchis alpina"), "Orchidaceae", New.Family),
         New.Taxonomic.status = ifelse(Taxon %in% c("Chamaeorchis alpina"), "Fixed", New.Taxonomic.status),
        New.Species = ifelse(Taxon %in% c("Valeriana stichensis"), "Valeriana sitchensis", New.Species), # Spelling error in species name
         New.Family = ifelse(Taxon %in% c("Valeriana stichensis"), "Caprifoliaceae", New.Family),
         New.Taxonomic.status = ifelse(Taxon %in% c("Valeriana stichensis"), "Fixed", New.Taxonomic.status))

# Rejoin new species information to full ITEX dataset
itex.9 <- left_join(itex.8, species.itex.checked.4, by = c("SPECIES_NAME" = "Taxon")) %>% 
  relocate(New.Species, New.Genus, New.Family, New.Taxonomic.status, .after = SPECIES_NAME)


# ASSIMILATE SPECIES NAMES (itex.10) ----

# For non-morphospecies utilise new provisional names/genuses/families, for morphospecies retain new name but retain genuses from previous
itex.10 <- itex.9 %>% 
  mutate(Prov.Species = New.Species,
         Prov.Genus = ifelse(!New.Taxonomic.status == "Morphospecies", New.Genus, GENUS),
         Prov.Family = ifelse(!New.Taxonomic.status == "Morphospecies", New.Family, NA)) %>% 
  relocate(Prov.Species, Prov.Genus, Prov.Family, .after = New.Taxonomic.status)


# FINAL GENUS, FAMILY CHECKED DATASET (itex.11) ----

# Check the non-morphospecies for genus and families
non.morphospecies.final <- itex.10 %>% 
  dplyr::select(New.Taxonomic.status, Prov.Species, Prov.Genus, Prov.Family) %>% 
  filter(!New.Taxonomic.status == "Morphospecies") %>% 
  dplyr::select(-New.Taxonomic.status) %>% 
  unique()

# Check the number of species matches number of unique (duplicates)
length(unique(non.morphospecies.final$Prov.Species)) # 374
length(non.morphospecies.final$Prov.Species) # 374 (NO duplicates)
  # All non-morphospecies have correct genus and family information

# Check the morphospecies for genus and family information
morphospecies.final.1 <- itex.10 %>% 
  dplyr::select(New.Taxonomic.status, Prov.Species, Prov.Genus, Prov.Family) %>% 
  filter(New.Taxonomic.status == "Morphospecies") %>% 
  dplyr::select(-New.Taxonomic.status) %>% 
  unique() %>% 
  arrange(Prov.Species)

# Extract list of families by genus from non-morphospecies
families <- non.morphospecies.final %>% 
  dplyr::select(Prov.Genus, Prov.Family) %>% 
  unique()

# Check the number of genuses matches number of unique (duplicates)
length(unique(families$Prov.Genus)) # 142
length(families$Prov.Genus) # 142 (NO duplicates)
  # All genuses correspond to just one family

# Join family information to morphospecies dataframe
morphospecies.final.2 <- morphospecies.final.1 %>% 
  dplyr::select(-Prov.Family) %>% 
  left_join(., families, by = c("Prov.Genus" = "Prov.Genus"))

# Amend missing genus information
morphospecies.final.3 <- morphospecies.final.2 %>%
  mutate(Prov.Genus = ifelse(str_detect(Prov.Species, pattern = paste(c("XXXFORB", "XXXGRAMINOID", "XXXGRASS",
                                                                        "XXXSEDGE", "XXXUNKNOWN", "XXXCYPERACEAE",
                                                                        "XXXPOACEAE", "XXXORCHIDACEAE", "XXXUNKNOWN"), collapse = "|")), NA, Prov.Genus),
         Prov.Genus = ifelse(str_detect(Prov.Species, pattern = c("XXXPOA:")), "Poa", Prov.Genus),
         Prov.Genus = ifelse(str_detect(Prov.Species, pattern = c("XXXCAREX:")), "Carex", Prov.Genus),
         Prov.Genus = ifelse(str_detect(Prov.Species, pattern = c("XXXDRABA:")), "Draba", Prov.Genus),
         Prov.Genus = ifelse(str_detect(Prov.Species, pattern = c("XXXLUZULA:")), "Luzula", Prov.Genus),
         Prov.Genus = ifelse(str_detect(Prov.Species, pattern = c("XXXEQUISETUM")), "Equisetum", Prov.Genus),
         Prov.Genus = ifelse(str_detect(Prov.Species, pattern = c("XXXHIERACIUM")), "Hieracium", Prov.Genus),
         Prov.Genus = ifelse(str_detect(Prov.Species, pattern = c("XXXLUZULA")), "Luzula", Prov.Genus),
         Prov.Genus = ifelse(str_detect(Prov.Species, pattern = c("XXXOXYTROPIS")), "Oxytropis", Prov.Genus),
         Prov.Genus = ifelse(str_detect(Prov.Species, pattern = c("XXXSAGINA")), "Sagina", Prov.Genus),
         Prov.Genus = ifelse(str_detect(Prov.Species, pattern = c("XXXSAXIFRAGA")), "Saxifraga", Prov.Genus),
         Prov.Genus = ifelse(str_detect(Prov.Species, pattern = c("XXXTARAXACUM")), "Taraxacum", Prov.Genus))

# Amend missing family information
morphospecies.final.4 <- morphospecies.final.3 %>%
  mutate(Prov.Family = ifelse(str_detect(Prov.Species, pattern = paste(c("XXXFORB", "XXXGRAMINOID", "XXXGRASS",
                                                                         "XXXSEDGE", "XXXUNKNOWN", "XXXESHRUB"), collapse = "|")), NA, Prov.Family),
         Prov.Family = ifelse(str_detect(Prov.Species, pattern = c("XXXPOLYGONUM:")), "Polygonaceae", Prov.Family),
         Prov.Family = ifelse(str_detect(Prov.Species, pattern = c("XXXSENECIO:")), "Compositae", Prov.Family),
         Prov.Family = ifelse(str_detect(Prov.Species, pattern = c("XXXVIOLA:")), "Violaceae", Prov.Family),
         Prov.Family = ifelse(str_detect(Prov.Species, pattern = c("XXXALCHEMILLA:")), "Rosaceae", Prov.Family),
         Prov.Family = ifelse(str_detect(Prov.Species, pattern = c("XXXCAREX:")), "Cyperaceae", Prov.Family),
         Prov.Family = ifelse(str_detect(Prov.Species, pattern = c("XXXPOA:")), "Poaceae", Prov.Family),
         Prov.Family = ifelse(str_detect(Prov.Species, pattern = c("XXXPOACEAE")), "Poaceae", Prov.Family),
         Prov.Family = ifelse(str_detect(Prov.Species, pattern = c("XXXDRABA:")), "Brassicaceae", Prov.Family),
         Prov.Family = ifelse(str_detect(Prov.Species, pattern = c("XXXLUZULA:")), "Juncaceae", Prov.Family),
         Prov.Family = ifelse(str_detect(Prov.Species, pattern = c("XXXRANUNCULACEAE:")), "Ranunculaceae", Prov.Family),
         Prov.Family = ifelse(str_detect(Prov.Species, pattern = c("XXXSAXIFRAGA:")), "Saxifragaceae", Prov.Family),
         Prov.Family = ifelse(str_detect(Prov.Species, pattern = c("XXXASTERACEAE:")), "Asteraceae", Prov.Family),
         Prov.Family = ifelse(str_detect(Prov.Species, pattern = c("XXXCYPERACEAE")), "Cyperaceae", Prov.Family),
         Prov.Family = ifelse(str_detect(Prov.Species, pattern = c("XXXEQUISETUM")), "Equisetaceae", Prov.Family),
         Prov.Family = ifelse(str_detect(Prov.Species, pattern = c("XXXHEPATICA")), "Ranunculaceae", Prov.Family),
         Prov.Family = ifelse(str_detect(Prov.Species, pattern = c("XXXHIERACIUM")), "Asteraceae", Prov.Family),
         Prov.Family = ifelse(str_detect(Prov.Species, pattern = c("XXXORCHIDACEAE")), "Orchidaceae", Prov.Family),
         Prov.Family = ifelse(str_detect(Prov.Species, pattern = c("XXXOXYTROPIS")), "Leguminosae", Prov.Family),
         Prov.Family = ifelse(str_detect(Prov.Species, pattern = c("XXXSAGINA")), "Caryophyllaceae", Prov.Family),
         Prov.Family = ifelse(str_detect(Prov.Species, pattern = c("XXXTARAXACUM")), "Asteraceae", Prov.Family))

# Final checks on morphospecies dataframe
morphospecies.final.5 <- morphospecies.final.4 %>% 
  distinct(Prov.Species, Prov.Genus, Prov.Family, .keep_all = TRUE) # Removes one now assimilated record

# Combine the final species dataframes
species.itex.final <- rbind(non.morphospecies.final, morphospecies.final.5) %>% 
  mutate(Final.Species = Prov.Species) %>% # Create final species column
  rename(Final.Genus = Prov.Genus, Final.Family = Prov.Family) %>% # Rename final genuses and families column
  relocate(Final.Species, Final.Genus, Final.Family, .after = Prov.Species) %>% 
  mutate(Final.Family = ifelse(Final.Family == "Compositae", "Asteraceae", Final.Family)) %>% 
  distinct(Prov.Species, Final.Species, Final.Genus, Final.Family, .keep_all = TRUE) # Remove duplicates

# Check each genus is only in one family
species.itex.final.check <- species.itex.final %>% 
  group_by(Final.Species) %>% 
  summarise(count = length(unique(Final.Genus))) %>% 
  ungroup()
    # All species have unique genus

# Check each genus is only in one family
unique(species.itex.final.check$count)

# Check each genus is only in one family
genus.itex.final.check <- species.itex.final %>% 
  group_by(Final.Genus) %>% 
  summarise(count = length(unique(Final.Family))) %>% 
  ungroup()
    # All genuses have unique family

# Join to main ITEX database and remove superfluous columns
itex.11 <- left_join(itex.10, species.itex.final, by = c("Prov.Species" = "Prov.Species")) %>% 
  dplyr::select(-c(GENUS, SPECIES_NAME, # Original columns
                   New.Species, New.Genus, New.Family, New.Taxonomic.status, # TPL columns
                   Prov.Species, Prov.Genus, Prov.Family)) %>% # First attempt at merging two
  relocate(Final.Species, Final.Genus, Final.Family, GFNARROWwalker, .after = ORIGINAL_NAME) %>% 
  rename(SPECIES = Final.Species, GENUS = Final.Genus, FAMILY = Final.Family, GROWTH_FORM = GFNARROWwalker)


# GROWTH FORM CHECKS (itex.12) ----

# Check values in GROWTH_FORM column
unique(itex.11$GROWTH_FORM)

# Get species for which growth forms is NA
gf.na <- filter(itex.11, is.na(GROWTH_FORM)) %>% 
  dplyr::select(SPECIES) %>% 
  unique()

# Create second column for more conglomerated growth forms and fixing the missing ones

gf.itex <- itex.11 %>% 
  mutate(FuncGroup = NA,
         FuncGroup = ifelse(GROWTH_FORM == "SEVER", "Evergreen_Shrub", FuncGroup),
         FuncGroup = ifelse(GROWTH_FORM == "SDECI", "Deciduous_Shrub", FuncGroup),
         FuncGroup = ifelse(GROWTH_FORM == "FORB", "Forb", FuncGroup),
         FuncGroup = ifelse(GROWTH_FORM %in% c("RUSH", "SEDGE", "GRASS", "GRAMINOIDU", "GRAMU"), "Graminoid", FuncGroup),
         FuncGroup = ifelse(SPECIES == "Carex maritima", "Graminoid", FuncGroup),
         FuncGroup = ifelse(SPECIES == "Bupleurum americanum", "Forb", FuncGroup), # https://plants.usda.gov/home/plantProfile?symbol=BUPLE
         FuncGroup = ifelse(SPECIES == "XXXEQUISETUM:ALEXFIORD", "Horsetail", FuncGroup)) %>% 
  relocate(FuncGroup, .after = GROWTH_FORM) %>% 
  filter(!FuncGroup %in% c("Horsetail")) # Remove horsetails as not one of the vascular plants we are investigating

# Now need to recalculate the Relative Cover again without the horsetails
gf.itex.cover.corr <- gf.itex %>% 
  group_by(SiteSubsitePlotYear) %>% 
  mutate(TotalRelCover = sum(RelCover)) %>%
  mutate(NewRelCover = (RelCover/TotalRelCover)*100) %>% # 5733 obs
  ungroup() %>% 
  mutate(RelCover = ifelse(TotalRelCover == 100, RelCover, NewRelCover)) %>% 
  dplyr::select(-c(NewRelCover, TotalRelCover))

# See if any genuses span FuncGroups (definitely can!), also picks up species with multiple GFs listed
gf.check <- gf.itex.cover.corr %>% 
  group_by(GENUS) %>% 
  mutate(count_fg = length(unique(FuncGroup))) %>% 
  ungroup() %>% 
  filter(count_fg > 1, !is.na(GENUS)) %>% 
  dplyr::select(SPECIES, GENUS, GROWTH_FORM, FuncGroup, count_fg) %>% 
  unique() %>% 
  arrange(GENUS)

# Manually correct incorrect records
itex.12 <- gf.itex.cover.corr %>% 
  mutate(FuncGroup = case_when(GENUS == "Luzula" ~ "Graminoid",
                               # SPECIES == "Comarum palustre" ~ "Forb",
                               GENUS == "Dryas" ~ "Evergreen_Shrub",
                               # SPECIES == "Linnaea borealis" ~ "Forb",
                               SPECIES == "Thymus praecox" ~ "Forb",
                               SPECIES == "Vaccinium vitis-idaea" ~ "Evergreen_Shrub", # Vaccinium can be both deciduous and evergreen
                               SPECIES == "Vaccinium uliginosum" ~ "Deciduous_Shrub",
                               SPECIES == "Vaccinium myrtillus" ~ "Deciduous_Shrub",
                               SPECIES == "Vaccinium microcarpum" ~ "Evergreen_Shrub",
                               SPECIES == "XXXFORB:QHI" ~ "Forb",
                               TRUE ~ FuncGroup))

# Check every species is only in one functional group
fg.check <- itex.12 %>% 
  group_by(SPECIES) %>% 
  summarise(FG_count = length(unique(FuncGroup))) %>% 
  ungroup() %>% 
  arrange(desc(FG_count))

# Visualise number
unique(fg.check$FG_count)


# REMOVE PLOTS WITH FAMILY OR HIGHER LEVEL MORPHOSPECIES (itex.13) ----

# How many unique plot:year combos currently retained
number.plots <- as.numeric(length(unique(itex.12$SiteSubsitePlotYear))) # 11,800

# All morphospecies at the family, fg or higher level are lacking genus information (checked), so mark those
remove.morphospecies <- itex.12 %>% 
  mutate(To_Remove = ifelse(is.na(GENUS), 1, 0)) %>% # If no NA values, all species recorded as 0
  group_by(SiteSubsitePlotYear) %>% 
  summarise(plot_sum = sum(To_Remove)) %>% # So when totalled = 0 so retain, if more than 0, clearly species needs removal
  ungroup() %>% 
  mutate(Remove = ifelse(plot_sum == 0, "No", "Yes")) %>% 
  dplyr::select(-plot_sum) # Creates key of plots to remove

# Number of plots removed (checks)
remove.number.plots <- remove.morphospecies %>% 
  group_by(Remove) %>% 
  summarise(total = length(SiteSubsitePlotYear)) %>% 
  ungroup()
    # Retain: 11,342; Remove: 458

# Join key to it dataframe and remove plots with family and higher level morphospecies
itex.13 <- left_join(itex.12, remove.morphospecies, by = c("SiteSubsitePlotYear" = "SiteSubsitePlotYear")) %>% 
  filter(!Remove == "Yes") %>% 
  dplyr::select(-Remove)

# Check correct number of plots retained
length(unique(itex.13$SiteSubsitePlotYear)) # 11,342 (correct!)

# Check no species with NA genus information retained (e.g. higher than family level morphospecies)
remove.morphospecies.check <- itex.13 %>% 
  filter(is.na(GENUS)) # 0 plots (correct!)


# PLOT OF MORPHOSPECIES PROPORTIONS ----
  
# Produce counts based on the dataframe counting plots to remove/keep (remove.number.plots)
remove.number.plots.2 <- remove.number.plots %>%
  mutate(Proportion = (total / number.plots) * 100) %>% 
  mutate(Proportion = round(Proportion, digits = 2)) %>% 
  arrange(Remove) %>% 
  mutate(csum = rev(cumsum(rev(Proportion))), # y-position of labels
         pos = Proportion/2 + lead(csum, 1),
         pos = if_else(is.na(pos), Proportion/2, pos))

# Create a pie chart of the proportion of plots containing morphospecies information ABOVE genus level
(pie.remove.morphospecies <- ggplot(remove.number.plots.2, aes(x = "", y = Proportion, fill = Remove)) +
    geom_bar(stat = "identity", width = 1, colour = "black") +
    scale_fill_manual(values = c("#F6D55C", "#ED553B")) +
    coord_polar("y", start = 0) +
    geom_label_repel(aes(y = pos, label = paste0(Proportion, "%")), size = 4.5, nudge_x = 1, show.legend = FALSE) +
    labs(title = "Proportion of Plot & Year Combinations containing Family or Higher Level Morphospecies",
         subtitle = "Combinations containing Family or Higher Level Morphospecies to be Removed",
         fill = "Contains Family or Higher Level Morphospecies") +
    theme_void() +
    theme(plot.title = element_text(face = "bold", size = 20, hjust = 0.5),
          plot.subtitle = element_text(size = 16, hjust = 0.5),
          legend.position = "bottom"))

# Save pie chart
ggsave(pie.remove.morphospecies, file = "figures/plot-year-xxx-proportion.png", width = 12, height = 10)


# REMOVE PLOTS WITH ONLY ONE SPECIES (breaks FD calculations) (itex.14) ----

# # Remove plot:year combos that contain just a single species
# itex.14 <- itex.13 %>%
#   group_by(SiteSubsitePlotYear) %>%
#   mutate(Species_Count = length(unique(SPECIES))) %>%
#   ungroup() %>%
#   filter(Species_Count > 1) %>% # Retains plot:year combos which contain more than one species
#   dplyr::select(-Species_Count)

# ALTERNATIVE: retain plots that have a single species (decided: 27/02/2024)
itex.14 <- itex.13


# CORRECT RELCOVER AND ENSURE SPECIES NOT DUPLICATED THROUGH RENAMING (itex.15) ----

# Remove records where relative cover == 0
itex.rel.cover <- itex.14 %>% 
  filter(RelCover != 0)

# Determine if any species occur twice in a plot due to renaming (occasionaly)
itex.species.dupl <- itex.rel.cover %>% 
  group_by(SiteSubsitePlotYear, SPECIES) %>% 
  summarise(Duplicate = length(RelCover)) %>% 
  ungroup()

# Combine records where species have same name but different records
itex.15 <- itex.rel.cover %>% 
  group_by(SiteSubsitePlotYear, SPECIES) %>% 
  mutate(RelCover_Corr = sum(RelCover)) %>% 
  ungroup() %>% 
  distinct(SiteSubsitePlotYear, SPECIES, RelCover_Corr, .keep_all = TRUE) %>% 
  dplyr::select(-RelCover) %>% 
  rename(RelCover = RelCover_Corr) %>% 
  relocate(RelCover, .after = FuncGroup)

# Check no species duplicates now
itex.species.dupl.2 <- itex.15 %>% 
  group_by(SiteSubsitePlotYear, SPECIES) %>% 
  summarise(Duplicate = length(RelCover)) %>% 
  ungroup()

# Check cover for all plots still equals 100
itex.species.dupl.cover <- itex.15 %>% 
  group_by(SiteSubsitePlotYear) %>% 
  summarise(TOTALCOVER = sum(RelCover)) %>% 
  ungroup()


# METHODS CHECKS (itex.16) ----

# Observe unique value types
unique(itex.15$ValueType)

# Standardise and create column of simplified methods
itex.16 <- itex.15 %>% 
  mutate(ValueType = ifelse(ValueType == "pf_all_xy", "pf_all_XY", ValueType), # Standardising capitalisation of XY
         ValueType = ifelse(ValueType == "pf_topbot_xy", "pf_topbot_XY", ValueType)) %>% 
  mutate(Method = case_when(ValueType == "percent_cover" ~ "Cover",
                            ValueType %in% c("pf_all_plot", "pf_topbot_plot") ~ "Point-framing (sum)",
                            ValueType %in% c("pf_all_XY", "pf_top_XY", "pf_topbot_XY") ~ "Point-framing (XY)"))


# FUNCTIONAL GROUP DOMINATION (itex.17) ----

# Percentage cover per functional group
fg.cover <- itex.16 %>% 
  group_by(SiteSubsitePlotYear, FuncGroup) %>% 
  mutate(FuncPlotCover = sum(RelCover)) %>%
  ungroup()

# Create a key for joining plot-dominating FG information to itex dataframe
fg.cover.key <- fg.cover %>% 
  dplyr::select(SiteSubsitePlotYear, FuncGroup, FuncPlotCover) %>% # Only columns we need for calculating plot-dominating FG
  distinct(SiteSubsitePlotYear, FuncGroup, .keep_all = TRUE) %>% # Keeps one record of each type of FG in each plot:year combo
  pivot_wider(names_from = FuncGroup, values_from = FuncPlotCover, values_fill = list(FuncPlotCover = 0)) %>% # 4643 plots (correct!)
  mutate(PlotDominatingFG = case_when(Evergreen_Shrub > 50 ~ "Evergreen_Shrub-Dominated",
                                      Deciduous_Shrub > 50 ~ "Deciduous_Shrub-Dominated",
                                      Graminoid > 50 ~ "Graminoid-Dominated",
                                      Forb > 50 ~ "Forb-Dominated",
                                      TRUE ~ "None")) %>%
  rename(EShrubCover = Evergreen_Shrub, DShrubCover = Deciduous_Shrub, GraminoidCover = Graminoid, ForbCover = Forb) %>%
  relocate(DShrubCover, .after = EShrubCover)

# Rejoin FG cover and dominating information to main ITEX dataframe
itex.17 <- left_join(itex.16, fg.cover.key, by = c("SiteSubsitePlotYear" = "SiteSubsitePlotYear"))


# YEAR VARIABLES (itex.18) ----

# Calculate columns for repeats and earliest/latest years
itex.18 <- itex.17 %>% 
  group_by(SiteSubsitePlot) %>%
  mutate(Repeats = length(unique(YEAR)),
         YEAR_latest = max(YEAR),
         YEAR_earliest = min(YEAR)) %>%
  ungroup() %>%
  relocate(YEAR_earliest, YEAR_latest, Repeats, .after = YEAR)


# FINAL STATUS, SURVEYAREA AND TREATMENT CHECKS ----

# Get a vector of SiteSubite combos with NA STATUS
na.status <- filter(itex.18, is.na(STATUS))
na.status <- unique(na.status$SITE)

# Get a vector of SiteSubite combos with NA TREATMENT (none)
na.treatment <- filter(itex.18, is.na(TREATMENT))
na.treatment <- unique(na.treatment$SITE)

# Get a vector of SiteSubite combos with NA SURVEYED AREA (none)
na.surveyarea <- filter(itex.18, is.na(SurveyedArea))
na.surveyarea <- unique(na.surveyarea$SITE)

# Get complete dataframe of sites with crucial missing info, and what they are missing
na.all <- c(na.status, na.treatment, na.surveyarea)
na.all <- data.frame(na.all) %>% 
  rename(SITE = na.all) %>% 
  mutate(STATUS_IS_NA = ifelse(SITE %in% na.status, "YES", "NO"),
         TREATMENT_IS_NA = ifelse(SITE %in% na.treatment, "YES", "NO"),
         SURVEYAREA_IS_NA = ifelse(SITE %in% na.surveyarea, "YES", "NO")) %>% 
  distinct(SITE, STATUS_IS_NA, TREATMENT_IS_NA, SURVEYAREA_IS_NA, .keep_all = TRUE) %>% 
  arrange(SITE)

# Export list of missing info
write.csv(na.all, file = "data/output_01_status_treatment_surveyarea_NA.csv")


# FINAL TIDY (itex.19) ----

itex.19 <- itex.18 %>%
  
  # Remove cols, only got LIVE and CTL plots remaining, TISSUE not used
  dplyr::select(-c(STATUS, TREATMENT, TISSUE, PlotSize_m2)) %>%
  
  # Reorder columns into logical order
  relocate(ValueType, Method, SurveyedArea, .after = RelCover) %>% 
  relocate(lat_grid, lon_grid, gridcell, Region, AlpArc, AZONE, ELEV, MOISTURE, .after = LONG) %>% 
  
  # Arrange by plot:year alphabetically to make checking trait calcuations easier
  arrange(SiteSubsitePlotYear) %>% 
  
  # Remove (86) rows where cover = 0 (don't need to remove full plot as not impacting other records)
  filter(RelCover > 0) %>%
  
  # Add rowID column for calculating trait values later
  mutate(ID = row_number()) %>% 
  relocate(ID, .before = SiteSubsitePlotYear)


# EXPORT CLEAN ITEX DATASET ----

# Save full dataset
write.csv(itex.19, "data/output_01_itex_full.csv")


# EXPORT LIST OF SUBSITES AND PLOTS ----

# Create table of subsites and number of plots in each
itex.subsites.plots <- itex.19 %>% 
  group_by(SiteSubsite) %>% 
  summarise(NumberPlots = length(unique(SiteSubsitePlot))) %>% 
  ungroup() %>% 
  arrange(SiteSubsite)

# Export table
write.csv(itex.subsites.plots, "data/output_01_itex_subsites_plots.csv")


# MAP OF ITEX SITES INCLUDED IN STUDY ----

# Load in dataset for making map
itex.19 <- read.csv("data/output_01_itex_full.csv")

# Load in coordinates dataframe
coords.c.c.lgm <- read.csv("data/output_01_coords_lgm.csv")

# Produce dataframe for making plot
itex.map <- itex.19 %>% 
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

# Join the glaciated data to these points
itex.map.gl <- left_join(itex.map, dplyr::select(coords.c.c.lgm, LAT, LONG, glaciated), by = c("LAT" = "LAT", "LONG" = "LONG"))
                                       
# Transform coordinates to UTM
itex.map.transform <- transform_coord(itex.map.gl, lon = "LONG", lat = "LAT", bind = TRUE)

# Map based on size of circles
(map.itex.sites.size <- basemap(limits = 55, grid.size = 0.05, grid.col = "#949494",
                           land.size = 0.05, land.col = "#dcdcdc", land.border.col = "#000000") +
    geom_point(data = itex.map.transform,
               aes(x = lon.proj, y = lat.proj, size = Number_Plots_Binned, fill = Region),
               alpha = 0.85, colour = "#000000", shape = 21) +
    scale_size_discrete(range = c(4,12,20,28,36)) +
    geom_label_repel(data = subset(itex.map.transform, SITE %in% c("Akureyri", "Alexfiord", "Utqiagvik", "Latnja", "Pyramiden", "Wolfcreek")),
                     aes(lon.proj, lat.proj, label = SITE), color = "black", box.padding = 2,
                     segment.color = "black", segment.size = 0.7, fill = "white", label.size = 0.4,  size = 4.5) +
    labs(title = "ITEX Site Locations",
         # subtitle = "Glaciated during the LGM (~18,000 B.P.) (Ehlers and Gibbard, 2004)",
         subtitle = "",
         size = "Number of Plots",
         fill = "Glaciated?") +
    theme_4_map()
  )

# Export map
ggsave(map.itex.sites.size, file = "figures/itex_sites_polarview_size.png", height = 8, width = 8)


# Map based on size of circles
(map.itex.sites.shape <- basemap(limits = 55, grid.size = 0.05, grid.col = "#949494",
                                land.size = 0.05, land.col = "#dcdcdc", land.border.col = "#000000") +
    geom_point(data = itex.map.transform,
               aes(x = lon.proj, y = lat.proj, shape = glaciated, fill = Number_Plots_Binned),
               alpha = 0.85, colour = "#000000", size = 4) +
    scale_shape_manual(values = c(24, 21)) +
    scale_fill_viridis(discrete = TRUE, option = "A", begin = 0.2, end = 1, direction = -1) +
    geom_label_repel(data = subset(itex.map.transform, SITE %in% c("Akureyri", "Alexfiord", "Utqiagvik", "Latnja", "Pyramiden", "Wolfcreek")),
                     aes(lon.proj, lat.proj, label = SITE), color = "black", box.padding = 2,
                     segment.color = "black", segment.size = 0.7, fill = "white", label.size = 0.4,  size = 4.5) +
    labs(title = "ITEX Site Locations",
         # subtitle = "Glaciated during the LGM (~18,000 B.P.) (Ehlers and Gibbard, 2004)",
         subtitle = "",
         fill = "Number of Plots",
         shape = "Glaciated?") +
    guides(fill = guide_legend(override.aes = list(shape = 21))) + # Fixes bug in the legend
    theme_4_map()
)

# Export map
ggsave(map.itex.sites.shape, file = "figures/itex_sites_polarview_shape.png", height = 8, width = 8)


# EXTRA: EXPORT UNIQUE MEASUREMENT LOCATIONS ----

# Dataset of coordinate pairs of unique ITEX measurement locations
itex.locations <- itex.19 %>% 
  dplyr::select(LONG, LAT) %>%
  drop_na() %>%
  distinct(LONG, LAT, .keep_all = TRUE)

# Export coordinates to csv with correct sites retained
write.csv(itex.locations, "data/output_01_itex_coords.csv")


# EXTRA: EXPORT SPECIES LISTS BY SUBSITE, GRIDCELL, REGION & FULL ITEX ----

# Removing XXX species as can't determine if same or not between subsites/gridcells/regions 

# Create am empty list of species at each SUBSITE
species.list.SUBSITE <- list()

# Create loop to extract species at each SUBSITE and append as vector to list
for (i in as.character(unique(itex.19$SUBSITE))){
  
  # Filter out XXX species, filter to the index SUBSITE and arrange species in alphabetical order
  itex.SUBSITE <- itex.19 %>%
    filter(!str_detect(SPECIES, pattern = "XXX")) %>% 
    filter(SUBSITE == i) %>% 
    arrange(SPECIES)
  
  # Assign the index subsite as the SUBSITE name
  name.SUBSITE <- i
  
  # Create a vector of the species names in that SUBSITE
  species.SUBSITE <- as.character(unique(itex.SUBSITE$SPECIES))
  
  # Append the SUBSITE name and species vector to the list
  species.list.SUBSITE[[name.SUBSITE]] <- species.SUBSITE
  
}

# Remove unwanted variables
rm(itex.SUBSITE, name.SUBSITE, species.SUBSITE, i)

# Export list as .RData object
save(species.list.SUBSITE, file = "data/output_01_species_list_SUBSITE.RData")


# Create am empty list of species at each gridcell
species.list.gridcell <- list()

# Create loop to extract species at each gridcell and append as vector to list
for (i in as.character(unique(itex.19$gridcell))){
  
  # Filter out XXX species, filter to the index gridcell and arrange species in alphabetical order
  itex.gridcell <- itex.19 %>%
    filter(!str_detect(SPECIES, pattern = "XXX")) %>% 
    filter(gridcell == i) %>% 
    arrange(SPECIES)
  
  # Assign the index subsite as the gridcell name
  name.gridcell <- i
  
  # Create a vector of the species names in that gridcell
  species.gridcell <- as.character(unique(itex.gridcell$SPECIES))
  
  # Append the gridcell name and species vector to the list
  species.list.gridcell[[name.gridcell]] <- species.gridcell
  
}

# Remove unwanted variables
rm(itex.gridcell, name.gridcell, species.gridcell, i)

# Export list as .RData object
save(species.list.gridcell, file = "data/output_01_species_list_gridcell.RData")


# Create am empty list of species at each Region
species.list.Region <- list()

# Create loop to extract species at each Region and append as vector to list
for (i in as.character(unique(itex.19$Region))){
  
  # Filter out XXX species, filter to the index Region and arrange species in alphabetical order
  itex.Region <- itex.19 %>%
    filter(!str_detect(SPECIES, pattern = "XXX")) %>% 
    filter(Region == i) %>% 
    arrange(SPECIES)
  
  # Assign the index subsite as the Region name
  name.Region <- i
  
  # Create a vector of the species names in that Region
  species.Region <- as.character(unique(itex.Region$SPECIES))
  
  # Append the Region name and species vector to the list
  species.list.Region[[name.Region]] <- species.Region
  
}

# Remove unwanted variables
rm(itex.Region, name.Region, species.Region, i)

# Export list as .RData object
save(species.list.Region, file = "data/output_01_species_list_Region.RData")


# Create an overall species list for the ITEX dataset
species.list.ITEX <- list()

# Filter out XXX species and arrange species in alphabetical order
itex.FULL <- itex.19 %>% 
  filter(!str_detect(SPECIES, pattern = "XXX")) %>%
  arrange(SPECIES)

# Append the unique non-XXX species from the ITEX dataset
species.list.ITEX[["ITEX"]] <- as.character(unique(itex.FULL$SPECIES))

# Remove unwanted variables
rm(itex.FULL)

# Export list as .RData object
save(species.list.ITEX, file = "data/output_01_species_list_ITEX.RData")

