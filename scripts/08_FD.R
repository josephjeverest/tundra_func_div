# 08 - Calculating FD indices
# Joseph Everest
# December 2021, adapted March 2022, July 2023


# LOAD PACKAGES & FUNCTIONS ----

# Load required packages
library(tidyverse)
library(devtools)
library(FD) # Version: FD_1.0-12 previously, now FD_1.0-12.3

# Load functions
source("scripts/08_FD_FUNCTION.R")


# **[CHANGE]** - DECIDE ON TRAIT INPUT, PCs & MODEL RUN SELECTION ----

# Determine whether to run FD calculations or just import
run.FD.calculations <- FALSE

# Determine whether to calculate ITEX input or just import
run.itex.input <- FALSE

# Use PCA or not (now deemed that PCA is incorrect; default = FALSE)
run.with.PCA <- FALSE

# # Determine number of PCs to use (if running with PCA)
# pc.count <- 3

# Determine 'm' input for FD function (if running without PCA)
m.input <- 3 # Options: 'min', 'max' or 3

    # Determine filepaths and PC counts depending on whether using PCA or not
    if (run.with.PCA == TRUE){
      
      # Generate filepath for saving outputs
      pc.filepath <- paste0("_PCA_pc", pc.count)
      
    } else { # End of run.with.PCA == TRUE if

        # Generate filepath for saving outputs
        pc.filepath <- paste0("_nonPCA_", m.input)
      
    } # End of run.with.PCA == TRUE else


# TRAIT INPUT & DATA PREP. (trait.input) ----

# Load in trait values by species
if (run.with.PCA == TRUE){
  
  # Load in PCA trait data
  traits <- read.csv(paste0("data/output_07_traits", pc.filepath, ".csv")) %>%
    dplyr::select(-X) %>% 
    rename(ID = rowname)
  
} else { # End of if loop
 
  # Load in non-PCA trait data
  traits <- read.csv("data/output_06_traits_final.csv") %>%
    dplyr::select(-X)
  
} # End of else loop


# ITEX INPUT & DATA PREP. (itex.input) ----

# Load in ITEX data
itex <- read.csv("data/output_04_itex_all.csv") %>% 
  dplyr::select(-X)

# Run loop for if need to remove plots based on PCAs
if (run.with.PCA == TRUE){
  
  # Determine number of species per plot and remove any with equal or less species than traits
  itex.sr <- itex %>% 
    group_by(SiteSubsitePlotYear) %>% 
    mutate(SR = length(unique(SPECIES))) %>% 
    ungroup()
  
  # Export list of plots with their corresponding species richness
  write.csv(dplyr::select(itex.sr, SiteSubsitePlotYear, SiteSubsitePlot, SR),
            file = "data/output_08_itex_sr_list.csv",
            row.names = FALSE)
  
  # Remove any plots with equal or less species than traits
  itex.trimmed <- itex.sr %>% 
    filter(SR >= pc.count + 1) %>% 
    dplyr::select(-SR)
  
} else { # End of if loop
  
  # Just carry over itex input to itex
  itex.trimmed <- itex
  
} # End of else loop


# Get vector of retained IDs for trimming trait data
itex.IDs <- unique(itex.trimmed$ID)


# Generate wide form of the ITEX data with just plot_ID and species' abundances
if (run.itex.input == TRUE){
  
  # Generate itex input (long run time)
  itex.input <- itex.trimmed %>%
    dplyr::select(SiteSubsitePlotYear, ID, RelCover) %>%
    pivot_wider(names_from = ID, values_from = RelCover) %>%
    replace(is.na(.), 0) %>%
    dplyr::select(-SiteSubsitePlotYear) # Rows in trait.input must match columns in itex.input (32,218)
  
  # Save itex.input (long loading time)
  save(itex.input, file = paste0("data/output_08_itex_input", pc.filepath, ".RData"))
  
} else { # End of if loop
  
  # Load in itex.input file
  get(load(paste0("data/output_08_itex_input", pc.filepath, ".RData")))
  
} # End of else loop


# TRIM TRAIT DATA AND PREPARE ----

# Trim based on IDs retained in ITEX
traits.trimmed <- traits %>% 
  filter(ID %in% itex.IDs)

# Make the record ID into row name
trait.input <- data.frame(traits.trimmed, row.names = 1)


# FD: CALCULATION (FRic, FEve & SR) (fd._) ----

# Create vector of rownames in itex.input
itex.rownames <- rownames(itex.input)

# Create broken-down vectors to run loops on as fails when trying to run 5000+ iterations at once
itex_rownames_01 <- itex.rownames[1:1000]
itex_rownames_02 <- itex.rownames[1001:2000]
itex_rownames_03 <- itex.rownames[2001:3000]
itex_rownames_04 <- itex.rownames[3001:4000]
itex_rownames_05 <- itex.rownames[4001:5000]
itex_rownames_06 <- itex.rownames[5001:6000]
itex_rownames_07 <- itex.rownames[6001:length(itex.rownames)] # Total number of rows in itex.input
# itex_rownames_test <- itex.rownames[6000:6010] # TEST VARIABLE

# Run FD calculation loops on each subset of ITEX plots
if (run.FD.calculations == TRUE){
  
  # Run FD calculation
  FD.calc(itex_rownames_01)
  FD.calc(itex_rownames_02)
  FD.calc(itex_rownames_03)
  FD.calc(itex_rownames_04)
  FD.calc(itex_rownames_05)
  FD.calc(itex_rownames_06)
  FD.calc(itex_rownames_07)
  # FD.calc(itex_rownames_test) # TEST VARIABLE
  
}

# Remove unwanted itex_rownames_... vectors
rm(itex_rownames_01, itex_rownames_02, itex_rownames_03, itex_rownames_04,
   itex_rownames_05, itex_rownames_06, itex_rownames_07)

# Load in the seven separate FD loop outputs
FD.output.01 <- read.csv(paste0("data/", list.files("data/", pattern = c(paste0("output_08_fd_output_", "01", pc.filepath)))))
FD.output.02 <- read.csv(paste0("data/", list.files("data/", pattern = c(paste0("output_08_fd_output_", "02", pc.filepath)))))
FD.output.03 <- read.csv(paste0("data/", list.files("data/", pattern = c(paste0("output_08_fd_output_", "03", pc.filepath)))))
FD.output.04 <- read.csv(paste0("data/", list.files("data/", pattern = c(paste0("output_08_fd_output_", "04", pc.filepath)))))
FD.output.05 <- read.csv(paste0("data/", list.files("data/", pattern = c(paste0("output_08_fd_output_", "05", pc.filepath)))))
FD.output.06 <- read.csv(paste0("data/", list.files("data/", pattern = c(paste0("output_08_fd_output_", "06", pc.filepath)))))
FD.output.07 <- read.csv(paste0("data/", list.files("data/", pattern = c(paste0("output_08_fd_output_", "07", pc.filepath)))))
# FD.output.test <- read.csv(paste0("data/", list.files("data/", pattern = c(paste0("output_08_fd_output_", "test", pc.filepath)))))

# Join the seven separate outputs into one final dataframe
FD.combined <- rbind(FD.output.01, FD.output.02, FD.output.03, FD.output.04,
                     FD.output.05, FD.output.06, FD.output.07)

# If using m = 3 criterion, set the FRic for plots with three species to NA
if (m.input == 3){
  
  # Change FRic for when SR == 3
  FD.combined <- FD.combined %>% 
    mutate(FRic = ifelse(nbsp == 3, NA, FRic))
  
}

# Save combined output as one dataframe
write.csv(FD.combined, paste0("data/output_08_fd_output_raw", pc.filepath, ".csv"),
          row.names = FALSE)
          
# Remove unwanted variables
rm(FD.output.01, FD.output.02, FD.output.03, FD.output.04,
   FD.output.05, FD.output.06, FD.output.07)


# FD: JOIN FD OUTPUT TO ITEX (combo.1) ----

# Add row_ID to FD outputs for joining back to ITEX plot data AND remove unwanted variables
FD.final <- FD.combined %>% 
  dplyr::select(FRic, FEve, nbsp, FDis, qual.FRic, Warning) %>% 
  mutate(row_ID = rownames(.)) %>% 
  mutate(row_ID = as.numeric(row_ID)) %>% 
  relocate(row_ID, .before = FRic)

# Generate dataframe of unique ITEX sites with plot-based information (lose species-level)
itex.plots <- itex %>% 
  group_by(SiteSubsitePlotYear) %>% 
  mutate(Manual_SR = length(unique(SPECIES))) %>% # Determine manual SR for checking matching to FD output
  ungroup() %>% 
  dplyr::select(-c(ID, ORIGINAL_NAME, SPECIES, GENUS, FAMILY, GROWTH_FORM, RelCover, FuncGroup)) %>% 
  distinct(SiteSubsitePlotYear, .keep_all = TRUE) %>% 
  mutate(row_ID = row_number()) %>% 
  mutate(row_ID = as.numeric(row_ID)) %>% 
  relocate(row_ID, SiteSubsitePlotYear, SiteSubsitePlot, SiteSubsite, SITE, SUBSITE, PLOT, .before = YEAR)
  
# Join FD output data to ITEX plot data
combo.initial <- left_join(itex.plots, FD.final, by = c("row_ID" = "row_ID"))

# Check SR in original dataset and calculated dataset is correct
combo.test <- combo.initial %>% 
  mutate(SR_check = ifelse(nbsp == Manual_SR, TRUE, FALSE)) # Should all be true

# All values in the SR_check column should be TRUE
combo.test.output <- unique(combo.test$SR_check)

# If correct, can output without Manual_SR
if (combo.test.output == TRUE){
  
  # Create combo.1 output
  combo.1 <- combo.initial %>% 
    dplyr::select(-Manual_SR)
  
}


# REMOVE RESURVEYED PLOTS WITH DIFFERING NAMES ACROSS DIFFERENT YEARS ----

# There are plots where a = start year, whilst b & c = subsequent years in the same plots
combo.letters <- combo.1 %>% 
  mutate(LetterInPlot = ifelse(str_detect(PLOT, "[abc]$"), "YES", "NO"))

# Create a database that DOES NOT have letters in plots
combo.letters.no <- combo.letters %>% 
  filter(LetterInPlot == "NO")

# Create a database that DOES have letters in plots
combo.letters.yes <- combo.letters %>% 
  filter(LetterInPlot == "YES") %>% 
  mutate(PLOT = ifelse(LetterInPlot == "YES", str_remove(PLOT, "[abc]$"), PLOT)) %>% # Remove letter
  mutate(SiteSubsitePlotYear = paste0(SITE, ":", SUBSITE, ":", PLOT, ":", YEAR), # Reproduce without letter in plot/year name
         SiteSubsitePlot = paste0(SITE, ":", SUBSITE, ":", PLOT)) %>% 
  group_by(PLOT) %>% 
  filter(YEAR == max(YEAR)) %>%
  ungroup()

# Merge the two dataframes back together
combo.2 <- rbind(combo.letters.no, combo.letters.yes) %>% 
  arrange(row_ID) %>% 
  dplyr::select(-LetterInPlot)

# Remove unwanted objects
rm(combo.letters, combo.letters.no, combo.letters.yes)


# REMOVE PLOTS WITH THE 'ZERO DISTANCE(S)' WARNING ----

# Check what warnings are present in the warning column
unique.warnings <- unique(combo.2$Warning)

# Create the dataframe without warning plots
combo.3 <- combo.2 %>% 
  mutate(RETAIN = FALSE,
         RETAIN = ifelse(Warning == "no non-missing arguments to min; returning Inf", TRUE, RETAIN),
         RETAIN = ifelse(is.na(Warning), TRUE, RETAIN)) %>% 
  filter(RETAIN == TRUE) %>% 
  dplyr::select(-c(Warning, RETAIN))
  

# RUN STATISTICS ON PLOTS WITH THE 'ZERO DISTANCE(S)' WARNING ----

# Create dataframe of the warning plots
combo.warnings <- combo.2 %>% 
  filter(Warning == "Zero distance(s)")

# Number of plots
number.warnings <- nrow(combo.warnings)

# Percentage of plots
percentage.warnings <- (nrow(combo.warnings) / nrow(combo.2)) * 100

# Average number of species in removed plots
sr.warnings <- (sum(combo.warnings$nbsp) / nrow(combo.warnings))

# Average FRic value for removed plots
FRic.warnings <- (sum(combo.warnings$FRic) / nrow(combo.warnings))


# GENERATE PLOT LEVEL DATA (combo.all/combo.latest) ----

# Create column for repeat surveys and latest year at each plot
combo.all <- combo.3 %>% 
  relocate(FRic, FEve, nbsp, FDis, qual.FRic, .after = Region) %>% 
  rename(SR = nbsp)

# Create dataframe for latest years only
combo.latest <- combo.all %>% 
  filter(YEAR == YEAR_latest)


# EXPORT DATAFRAMES ----

# Export dataframe for all years
write.csv(combo.all, paste0("data/output_08_fd_output_combined", pc.filepath, "_all_years.csv"),
          row.names = FALSE)

# Export dataframe for latest years
write.csv(combo.latest, paste0("data/output_08_fd_output_combined", pc.filepath, "_latest_years.csv"),
          row.names = FALSE)
