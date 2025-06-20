# 08 - Spatial Models
# Joseph Everest
# December 2021, adapted February 2022, November 2022, July 2023, February 2024, January 2025


# LOAD PACKAGES ----

# Load (and install) required packages
library(tidyverse)
library(mgcv)
library(brms)
library(lme4)
library(ggeffects)
library(gridExtra)
library (tidybayes)
# devtools::install_github('m-clark/lazerhawk')
library(lazerhawk) # brms() summary tables to df
library(viridis)
library(broom)


# LOAD THEMES AND CUSTOM FUNCTIONS ----

# Load ggplot themes from separate source script
source("scripts/00_ggplot_themes.R")

# Load custom spatial functions
source("scripts/09_models_FUNCTION.R")


# **[CHANGE]** - DETERMINE INPUTS FOR TRAIT SELECTION, VARIABLES, COLOUR PALETTES & DISTRIBUTIONS ----

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

# Determine whether to check if convergence was true
check.convergence <- TRUE

# Determine which variables we want to generate slopes over time for
temporal.variables <- c("FRic", "FEve", "SR", "FDis", "ShrubCover", "GraminoidCover", "ForbCover")

# Determine colour palette inputs for continuous variables
FRic.colour <- "plasma"
FEve.colour <- "viridis"
SR.colour <- "mako"
FDis.colour <- "plasma"
ShrubCover.colour <- "plasma"
GraminoidCover.colour <- "viridis"
ForbCover.colour <- "mako"

# Determine colour palette inputs for single variables
FRic_slopes.colour <- "#B53977"
FEve_slopes.colour <- "#53C367"
SR_slopes.colour <- "#44B1A2"
FDis_slopes.colour <- "#B53977"
ShrubCover_slopes.colour <- "#B53977"
GraminoidCover_slopes.colour <- "#53C367"
ForbCover_slopes.colour <- "#44B1A2"

# Determine model distribution inputs for spatial variables
FRic.distribution <- "lognormal"
FEve.distribution <- "gaussian"
SR.distribution <- "negbinomial"
FDis.distribution <- "gaussian"

# Determine model distribution inputs for temporal variables
FRic_slopes.distribution <- "gaussian"
FEve_slopes.distribution <- "gaussian"
SR_slopes.distribution <- "gaussian"
FDis_slopes.distribution <- "gaussian"
ShrubCover_slopes.distribution <- "gaussian"
GraminoidCover_slopes.distribution <- "gaussian"
ForbCover_slopes.distribution <- "gaussian"


# IMPORT & MANIPULATE DATA - SPATIAL ----

# Import the ITEX and FD data for the latest year at each plot ONLY
combo.latest <- read.csv(paste0("data/output_08_fd_output_combined", pc.filepath, "_latest_years.csv")) %>% 
  # filter(FRic > 0) %>% # Remove rows where FRic == 0
  mutate(ShrubCover = EShrubCover + DShrubCover) %>% # Create variable for combined shrub cover
  dplyr::select(-EShrubCover, -DShrubCover) %>% 
  relocate(ShrubCover, .before = GraminoidCover) %>% 
  mutate(PlotDominatingFG = case_when(ShrubCover > 50 ~ "Shrub-Dominated", # Need to modify the PlotDominatingFG column as have now combined shrub cover
                                      GraminoidCover > 50 ~ "Graminoid-Dominated",
                                      ForbCover > 50 ~ "Forb-Dominated",
                                      TRUE ~ "None")) %>% 
  mutate(MOISTURE = ifelse(SiteSubsite %in% c("LATNJA:DRY_MEADOW", "LATNJA:TUSSOCK_TUNDRA"), "MOIST", MOISTURE))

# Run extra if statement to ensure plots have more species than traits
if (run.with.PCA == TRUE){
  
  # Run filter
  combo.latest <- filter(combo.latest, SR >= pc.count + 1)
  
}

# IMPORT & MANIPULATE DATA - TEMPORAL ----

# Import the ITEX and FD data for all years at each plot
combo.all.full <- read.csv(paste0("data/output_08_fd_output_combined", pc.filepath, "_all_years.csv")) %>% 
  # filter(FRic > 0, !is.na(FEve)) %>% # Remove rows where FRic == 0
  mutate(ShrubCover = EShrubCover + DShrubCover) %>% # Create variable for combined shrub cover
  dplyr::select(-EShrubCover, -DShrubCover) %>% 
  relocate(ShrubCover, .before = GraminoidCover) %>% 
  mutate(PlotDominatingFG = case_when(ShrubCover > 50 ~ "Shrub-Dominated", # Need to modify the PlotDominatingFG column as have now combined shrub cover
                                      GraminoidCover > 50 ~ "Graminoid-Dominated",
                                      ForbCover > 50 ~ "Forb-Dominated",
                                      TRUE ~ "None")) %>% 
  mutate(MOISTURE = ifelse(SiteSubsite %in% c("LATNJA:DRY_MEADOW", "LATNJA:TUSSOCK_TUNDRA"), "MOIST", MOISTURE))

# Run if statement to remove all plots that have at least one year with under the required SR
if (run.with.PCA == TRUE){
  
  # Run filter
  combo.all.cut <- filter(combo.all.full, SR >= pc.count + 1) # Run extra if statement to ensure plots have more species than traits

  # Import itex species richness list per plot
  itex.richness <- read.csv("data/output_08_itex_sr_list.csv")
  
  # Generate metric for whether less than the PC count
  itex.richness.remove <-  itex.richness %>% 
    mutate(REMOVE.count = ifelse(SR == pc.count, 1, 0)) %>% 
    group_by(SiteSubsitePlot) %>% 
    summarise(REMOVE = sum(REMOVE.count)) %>% 
    ungroup() %>% 
    mutate(REMOVE = ifelse(REMOVE > 0, TRUE, FALSE))
  
  # Join to temporal data and remove all plots that have at least one year with under the required SR
  combo.all <- combo.all.cut %>% 
    left_join(., itex.richness.remove, by = c("SiteSubsitePlot" = "SiteSubsitePlot")) %>% 
    filter(REMOVE != TRUE) %>% # Lose approximately 600 plots
    dplyr::select(-REMOVE)
  
} else { # End of if statement
  
  # Just retain combo.all.five.years as combo.all
  combo.all <- combo.all.full
  
}


# Create dataframe to append slopes to
slopes.combined <- data.frame()

# Run loop to generate slopes for each of the selected variables
for (i in temporal.variables){
  
  # Produce dataframe for running models on
  slopes.input.data.full <- combo.all %>% 
    dplyr::select(i, YEAR, SiteSubsitePlot) %>% 
    rename(x_variable = i) %>% # Rename first column to x_variable
    filter(!is.na(x_variable))
  
  # Remove instances where less than five years of data recorded at that plot
  slopes.input.data <- slopes.input.data.full %>% 
    group_by(SiteSubsitePlot) %>% 
    mutate(YEAR_earliest = min(YEAR),
           YEAR_latest = max(YEAR)) %>% 
    ungroup() %>% 
    mutate(timeframe = YEAR_latest - YEAR_earliest) %>% 
    filter(timeframe >= 5) %>% 
    dplyr::select(-c(YEAR_earliest, YEAR_latest, timeframe))
  
  # Run function to generate linear slopes
  slopes <- linear.slopes(i, slopes.input.data)
  
  # Append slopes to overall output
  slopes.combined <- rbind(slopes.combined, slopes)
  
}

# Create ITEX dataset with just one row per SiteSubsitePlot
combo.plots <- combo.all %>% 
  filter(Repeats > 1) %>% # To keep only the plots with time series
  dplyr::select(SiteSubsitePlot, SITE, SUBSITE, PLOT, MOISTURE, LAT, LONG, ELEV, PlotDominatingFG,
                SurveyedArea, AlpArc, lat_grid, lon_grid, gridcell, Region, TempAvSum) %>% # Unique plot-based info. (not time sensitive)
  unique() # Retains one record of each unique plot

# Convert dataframe to wide format for running bayesian models on and join to plot info.
slopes.input <- slopes.combined %>% 
  pivot_wider(names_from = "Slope_Type", values_from = "Slope") %>% 
  left_join(combo.plots, ., by = c("SiteSubsitePlot" = "SiteSubsitePlot"))


# ADD CLIMATE CHANGE DATA ----

# Import climate change data from Mariana's paper
climate.change <- read.csv("data/22clim_slopes.csv") %>% 
  dplyr::select(-X)

# Modify climate change data to leave only values of change per coordinate pair
climate.change.cut <- climate.change %>% 
  dplyr::select(LAT, LONG, WarmQSlope, PrecSlope) %>% 
  unique()

# Join to overall dataframe
slopes.input <- left_join(slopes.input, climate.change.cut, by = c("LAT" = "LAT", "LONG" = "LONG"))


# (Q1) - GEOGRAPHICAL MODELS ----

# Run FRic, FEve and SR models vs LAT
bayesian.spatial.continuous(run.FRic = FALSE,
                            censored = "Yes",
                            run.FEve = FALSE, 
                            quadratic = "No",
                            run.SR = FALSE,
                            run.FDis = FALSE,
                            x.var = "LAT",
                            data = combo.latest,
                            colour.by = "SR",
                            colour.discrete = FALSE,
                            chains = 2,
                            cores = 4,
                            warmup = 500,
                            iterations = 2000,
                            delta = 0.95,
                            treedepth = 13)

# Generate plots for the manuscript
plot.latitude(censored = "Yes")
plot.precip(censored = "Yes")
plot.moisture()
plot.region()
plot.latitude.time()
plot.temperature.time()
plot.precip.change()

# Run FRic, FEve and SR models vs Region
bayesian.spatial.categoric(run.FRic = FALSE, 
                           censored = "Yes",
                           run.FEve = FALSE, 
                           run.SR = FALSE, 
                           run.FDis = FALSE,
                           x.var = "Region",
                           data = combo.latest,
                           chains = 2,
                           cores = 4,
                           warmup = 500,
                           iterations = 2000,
                           delta = 0.95,
                           treedepth = 13)


# (Q2) - CLIMATE MODELS ----

# Run FRic, FEve and SR models vs TempAvSum (mean summer temperature)
bayesian.spatial.continuous(run.FRic = FALSE, 
                            censored = "Yes",
                            run.FEve = FALSE, 
                            quadratic = "No",
                            run.SR = FALSE, 
                            run.FDis = FALSE,
                            x.var = "TempAvSum",
                            data = combo.latest,
                            colour.by = "SR",
                            colour.discrete = FALSE,
                            chains = 2,
                            cores = 4,
                            warmup = 500,
                            iterations = 2000,
                            delta = 0.95,
                            treedepth = 13)

# Generate plots for the manuscript
plot.warmest.quarter(censored = "Yes")

# Run FRic, FEve and SR models vs PrecipAnn (annual precipitation)
bayesian.spatial.continuous(run.FRic = FALSE, 
                            censored = "Yes",
                            run.FEve = FALSE, 
                            quadratic = "No",
                            run.SR = FALSE, 
                            run.FDis = FALSE,
                            x.var = "PrecipAnn",
                            data = combo.latest,
                            colour.by = "PlotDominatingFG",
                            colour.discrete = TRUE,
                            chains = 2,
                            cores = 4,
                            warmup = 500,
                            iterations = 2000,
                            delta = 0.95,
                            treedepth = 13)

# Run FRic, FEve and SR models vs MOISTURE (soil moisture)
bayesian.spatial.categoric(run.FRic = FALSE, 
                           censored = "Yes",
                           run.FEve = FALSE, 
                           run.SR = FALSE, 
                           run.FDis = FALSE,
                           x.var = "MOISTURE",
                           data = combo.latest,
                           chains = 2,
                           cores = 4,
                           warmup = 500,
                           iterations = 2000,
                           delta = 0.95,
                           treedepth = 13)


# (Q3) - FUNCTIONAL GROUP MODELS ----

# Run FRic, FEve and SR models vs ShrubCover (percentage of plot shrub-covered)
bayesian.spatial.continuous(run.FRic = FALSE, 
                            censored = "Yes",
                            run.FEve = FALSE, 
                            quadratic = "Yes",
                            run.SR = FALSE, 
                            run.FDis = FALSE,
                            x.var = "ShrubCover",
                            data = combo.latest,
                            colour.by = "SR",
                            colour.discrete = FALSE,
                            chains = 2,
                            cores = 4,
                            warmup = 500,
                            iterations = 2000,
                            delta = 0.95,
                            treedepth = 17)

# Run quadratic FRic, FEve and SR models vs ShrubCover
bayesian.spatial.quadratic(run.FRic = FALSE,
                           censored = "No",
                           run.FEve = FALSE,
                           run.SR = FALSE,
                           run.FDis = FALSE,
                           x.var = "ShrubCover",
                           data = combo.latest,
                           colour.by = "SR",
                           colour.discrete = FALSE,
                           chains = 2,
                           cores = 4,
                           warmup = 500,
                           iterations = 2000,
                           delta = 0.95,
                           treedepth = 13)

# Run FRic, FEve and SR models vs ForbCover (percentage of plot forb-covered)
bayesian.spatial.continuous(run.FRic = FALSE, 
                            censored = "Yes",
                            run.FEve = FALSE, 
                            quadratic = "Yes",
                            run.SR = FALSE, 
                            run.FDis = FALSE,
                            x.var = "ForbCover",
                            data = combo.latest,
                            colour.by = "SR",
                            colour.discrete = FALSE,
                            chains = 2,
                            cores = 4,
                            warmup = 500,
                            iterations = 2000,
                            delta = 0.95,
                            treedepth = 17)

# Run quadratic FRic, FEve and SR models vs ForbCover
bayesian.spatial.quadratic(run.FRic = FALSE,
                           censored = "No",
                           run.FEve = FALSE,
                           run.SR = FALSE,
                           run.FDis = FALSE,
                           x.var = "ForbCover",
                           data = combo.latest,
                           colour.by = "SR",
                           colour.discrete = FALSE,
                           chains = 2,
                           cores = 4,
                           warmup = 500,
                           iterations = 2000,
                           delta = 0.95,
                           treedepth = 13)

# Run FRic, FEve and SR models vs GraminoidCover (percentage of plot graminoid-covered)
bayesian.spatial.continuous(run.FRic = FALSE, 
                            censored = "Yes",
                            run.FEve = FALSE, 
                            quadratic = "Yes",
                            run.SR = FALSE, 
                            run.FDis = FALSE,
                            x.var = "GraminoidCover",
                            data = combo.latest,
                            colour.by = "SR",
                            colour.discrete = FALSE,
                            chains = 2,
                            cores = 4,
                            warmup = 500,
                            iterations = 2000,
                            delta = 0.95,
                            treedepth = 17)

# Run quadratic FRic, FEve and SR models vs GraminoidCover
bayesian.spatial.quadratic(run.FRic = FALSE,
                           censored = "No",
                           run.FEve = FALSE,
                           run.SR = FALSE,
                           run.FDis = FALSE,
                           x.var = "GraminoidCover",
                           data = combo.latest,
                           colour.by = "SR",
                           colour.discrete = FALSE,
                           chains = 2,
                           cores = 4,
                           warmup = 500,
                           iterations = 2000,
                           delta = 0.95,
                           treedepth = 13)

# Run function to generate combined output plots for all three cover types
plot.combined.cover(censored.FRic = "Yes",quadratic.FEve = "No")
plot.combined.cover.quadratic(censored = "No")


# Run FRic, FEve and SR models vs PlotDominatingFG (functional group dominating plot)
bayesian.spatial.categoric(run.FRic = FALSE, 
                           censored = "Yes",
                           run.FEve = FALSE, 
                           run.SR = FALSE, 
                           run.FDis = FALSE,
                           x.var = "PlotDominatingFG",
                           data = combo.latest,
                           chains = 2,
                           cores = 4,
                           warmup = 500,
                           iterations = 2000,
                           delta = 0.95,
                           treedepth = 13)


# (Q4) - TEMPORAL MODELS ----

# Run FRic, FEve, SR, FDis, ShrubCover, GraminoidCover & ForbCover models vs time
bayesian.temporal.change(run.FRic = FALSE, 
                         run.FEve = FALSE, 
                         run.SR = FALSE, 
                         run.FDis = FALSE,
                         run.ShrubCover = FALSE, 
                         run.GraminoidCover = FALSE, 
                         run.ForbCover = FALSE, 
                         data = slopes.input,
                         chains = 2,
                         cores = 4,
                         warmup = 500,
                         iterations = 2000,
                         delta = 0.95,
                         treedepth = 13)

# Generate manuscript outputs
plot.change.histograms()

# Run FRic, FEve and SR slopes vs LAT (latitude)
bayesian.temporal.continuous(run.FRic = FALSE, 
                             run.FEve = FALSE, 
                             run.SR = FALSE, 
                             run.FDis = FALSE,
                             x.var = "LAT",
                             data = slopes.input,
                             colour.by = "PlotDominatingFG", # Must be categoric, e.g. "Region" or "PlotDominatingFG"
                             chains = 2,
                             cores = 4,
                             warmup = 500,
                             iterations = 2000,
                             delta = 0.95,
                             treedepth = 13)

# Run FRic, FEve and SR slopes vs TempAvSum (mean summer temperature)
bayesian.temporal.continuous(run.FRic = FALSE, 
                             run.FEve = FALSE, 
                             run.SR = FALSE, 
                             run.FDis = FALSE,
                             x.var = "TempAvSum",
                             data = slopes.input,
                             colour.by = "PlotDominatingFG", # Must be categoric, e.g. "Region" or "PlotDominatingFG"
                             chains = 2,
                             cores = 4,
                             warmup = 500,
                             iterations = 2000,
                             delta = 0.95,
                             treedepth = 13)

# Run FRic, FEve and SR slopes vs WarmQSlope (change in summer temperature)
bayesian.temporal.continuous(run.FRic = FALSE, 
                             run.FEve = FALSE, 
                             run.SR = FALSE, 
                             run.FDis = FALSE,
                             x.var = "WarmQSlope",
                             data = slopes.input,
                             colour.by = "PlotDominatingFG", # Must be categoric, e.g. "Region" or "PlotDominatingFG"
                             chains = 2,
                             cores = 4,
                             warmup = 500,
                             iterations = 2000,
                             delta = 0.95,
                             treedepth = 13)

# Plot manuscript plots for temperature change vs functional diversity change
plot.temp.change()

# Run FRic, FEve and SR slopes vs PrecSlope (change in annual precipitation)
bayesian.temporal.continuous(run.FRic = FALSE, 
                             run.FEve = TRUE, 
                             run.SR = TRUE, 
                             run.FDis = TRUE,
                             x.var = "PrecSlope",
                             data = slopes.input,
                             colour.by = "PlotDominatingFG", # Must be categoric, e.g. "Region" or "PlotDominatingFG"
                             chains = 2,
                             cores = 4,
                             warmup = 500,
                             iterations = 2000,
                             delta = 0.95,
                             treedepth = 13)

# Run FRic, FEve and SR slopes vs ShrubCover_slopes (change in shrub cover)
bayesian.temporal.continuous(run.FRic = FALSE, 
                             run.FEve = FALSE, 
                             run.SR = FALSE, 
                             run.FDis = FALSE,
                             x.var = "ShrubCover_slopes",
                             data = slopes.input,
                             colour.by = "PlotDominatingFG", # Must be categoric, e.g. "Region" or "PlotDominatingFG"
                             chains = 2,
                             cores = 4,
                             warmup = 500,
                             iterations = 2000,
                             delta = 0.95,
                             treedepth = 13)

# Run FRic, FEve and SR slopes vs GraminoidCover_slopes (change in graminoid cover)
bayesian.temporal.continuous(run.FRic = FALSE,
                             run.FEve = FALSE,
                             run.SR = FALSE,
                             run.FDis = FALSE,
                             x.var = "GraminoidCover_slopes",
                             data = slopes.input,
                             colour.by = "PlotDominatingFG", # Must be categoric, e.g. "Region" or "PlotDominatingFG"
                             chains = 2,
                             cores = 4,
                             warmup = 500,
                             iterations = 2000,
                             delta = 0.95,
                             treedepth = 13)

# Run FRic, FEve and SR slopes vs ForbCover_slopes (change in forb cover)
bayesian.temporal.continuous(run.FRic = FALSE,
                             run.FEve = FALSE,
                             run.SR = FALSE,
                             run.FDis = FALSE,
                             x.var = "ForbCover_slopes",
                             data = slopes.input,
                             colour.by = "PlotDominatingFG", # Must be categoric, e.g. "Region" or "PlotDominatingFG"
                             chains = 2,
                             cores = 4,
                             warmup = 500,
                             iterations = 2000,
                             delta = 0.95,
                             treedepth = 13)

# Plot the combined cover result
plot.combined.cover.change()


# (Q5) - METRIC MODEL ----

# Run SR vs FD metrics models
bayesian.metric.comparison(run.FRic = FALSE,
                           quadratic = "Yes", # Or "No"
                           run.FEve = FALSE,
                           run.FDis = FALSE,
                           x.var = "SR",
                           data = combo.latest,
                           chains = 2,
                           cores = 4,
                           warmup = 500,
                           iterations = 2000,
                           delta = 0.95,
                           treedepth = 13)

# Run SR vs FD GAMs
GAM.metric.comparison(run = TRUE)


# CONVERGENCE ISSUES AND MODEL RESULTS ----

# Run function to check for model convergence
bayesian.results.convergence.warnings()

# Load in convergence output
results.convergence.warnings <- read.csv(paste0("data/model_outputs_new/",
                                        "convergence_summaries", pc.filepath, ".csv"))

# Tidy model outputs to only retain ones we want
results.output <- results.convergence.warnings %>% 
  dplyr::select(-c(warning, error, rhat, problematic_rhat_count)) %>%
  mutate(significant = ifelse(significant == "Yes", "Y", "N")) %>% 
  mutate(model = ifelse(model == "m_FRic_lognormal_quadratic_SR_nonPCA_min", "m_FRic_quadratic_SR", model)) %>% 
  mutate(model = str_remove(model, pattern = "_nonPCA_min"),
         model = str_remove(model, pattern = "m_"),
         model = str_remove(model, pattern = "_gaussian")) %>% 
  filter(!str_detect(model, pattern = "FEve_ForbCover"),
         !str_detect(model, pattern = "FEve_GraminoidCover"),
         !str_detect(model, pattern = "FEve_ShrubCover"),
         !str_detect(model, pattern = "lognormal")) %>% 
  mutate(Structure = ifelse(str_detect(model, pattern = "quadratic"), "Quadratic", "Linear"),
         Structure = ifelse(str_detect(model, pattern = "censored"), "Censored", Structure),
         model = str_remove(model, pattern = "_quadratic"),
         model = str_remove(model, pattern = "_censored")) %>% 
  
  mutate(model.split = model) %>% 
  separate(model.split, into = c("a", "b", "c", "d"), sep = "_") %>%
  mutate(y.var = NA, # Determine y variables
         y.var = ifelse(is.na(c) & b != "slopes", a, y.var),
         y.var = ifelse(b == "slopes", paste0(a, " ", b), y.var),
         y.var = ifelse(b == "lognormal", a, y.var)) %>% 
  mutate(x.var = NA, # Detemine x variables
         x.var = ifelse(is.na(c) & b != "slopes", b, x.var),
         x.var = ifelse(is.na(c) & b == "slopes", "1", x.var),
         x.var = ifelse(!is.na(c) & is.na(d), c, x.var),
         x.var = ifelse(!is.na(d), paste0(c, " ", d), x.var)) %>% 
  dplyr::select(-c(a, b, c, d, model, significant)) %>% 
  
  mutate(Covariate = ifelse(Covariate == "REPLACE", NA, Covariate),
         Covariate = ifelse(x.var == "1", NA, Covariate),
         x.var.type = ifelse(x.var == "1", NA, x.var.type),
         x.var.type = ifelse(x.var.type == "continuous", "Continuous", x.var.type)) %>% 
  
  mutate(Retain = ifelse(Structure == "Quadratic" & Covariate == "Reference", FALSE, TRUE)) %>%
  filter(Retain == TRUE) %>%
  dplyr::select(-Retain) %>% 
  
  mutate(Covariate = ifelse(Covariate == "centred_", "1st Order", Covariate),
         Covariate = ifelse(Covariate == "Icentred_E2", "2nd Order", Covariate),
         Covariate = ifelse(Covariate == "centred_SR", "1st Order", Covariate),
         Covariate = ifelse(Covariate == "Icentred_SRE2", "2nd Order", Covariate),
         
         
         Covariate = ifelse(Covariate == "ShrubMDominated", "Shrub Dominated", Covariate),
         Covariate = ifelse(Covariate == "GraminoidMDominated", "Graminoid Dominated", Covariate),
         Covariate = ifelse(Covariate == "ForbMDominated", "Forb Dominated", Covariate),
         Covariate = ifelse(Covariate == "NorthAmericaMEast", "N.America-E", Covariate),
         Covariate = ifelse(Covariate == "NorthAmericaMWest", "N.America-W", Covariate),
         Covariate = ifelse(Covariate == "WET", "Wet", Covariate),
         Covariate = ifelse(Covariate == "MOIST", "Moist", Covariate)) %>%
  
  mutate(x.var = str_replace_all(x.var, pattern = "Cover slopes", replacement = " change"),
         x.var = str_replace_all(x.var, pattern = "Cover", replacement = "s"),
         y.var = str_replace_all(y.var, pattern = "Cover slopes", replacement = " change"),
         y.var = str_replace_all(y.var, pattern = "slopes", replacement = "change"),
         x.var = ifelse(x.var == "PlotDominatingFG", "Dominant FG", x.var),
         x.var = ifelse(x.var == "Graminoid change", "Gram. change", x.var),
         y.var = ifelse(y.var == "Graminoid change", "Gram. change", y.var),
         Covariate = str_remove(Covariate, pattern = " Dominated")) %>% 
  
  # str_replace(x.var, pattern = "Cover slopes", replacement = " change") %>% 
  # str_replace(x.var, pattern = "Cover", replacement = "s") %>% 
  
  rename(X = x.var,
         "X Type" = x.var.type,
         Y = y.var,
         Estimate = estimate,
         "Lower CI (95%)" = lower.CI,
         "Upper CI (95%)" = upper.CI) %>% 
  
  relocate(Y, X, "X Type", Covariate, Structure, .after = ) %>% 
  arrange(X, Y)

# # Tidy model results output
# results.output <- results.convergence.warnings %>% 
#   mutate(PCA = ifelse(str_detect(model, "PCA_pc"), TRUE, FALSE), # Create PCA columns
#          NumberOfPCs = ifelse(str_detect(model, "pc3"), 3, NA),
#          NumberOfPCs = ifelse(str_detect(model, "pc4"), 4, NumberOfPCs)) %>%
#   mutate(CriteriaOfPCoA = ifelse(str_detect(model, "nonPCA"), "max (incorrect)", NA),
#          CriteriaOfPCoA = ifelse(str_detect(model, "nonPCA_min"), "min (correct)", CriteriaOfPCoA)) %>% 
#   mutate(model = str_remove(model, "m_"), # Trim the model name to variables
#          # model = str_remove(model, paste0("_PCA_pc", pc.count, ".RData")),
#          model = str_remove(model, paste0(pc.filepath, ".RData"))) %>%
#   mutate(model.split = model) %>% # Determine individual items
#   separate(model.split, into = c("a", "b", "c", "d"), sep = "_") %>%
#   mutate(y.var = NA, # Determine y variables
#          y.var = ifelse(is.na(c) & b != "slopes", a, y.var),
#          y.var = ifelse(b == "slopes", paste0(a, " ", b), y.var),
#          y.var = ifelse(b == "lognormal", a, y.var)) %>% 
#   mutate(x.var = NA, # Detemine x variables
#          x.var = ifelse(is.na(c) & b != "slopes", b, x.var),
#          x.var = ifelse(is.na(c) & b == "slopes", "1", x.var),
#          x.var = ifelse(!is.na(c) & is.na(d), c, x.var),
#          x.var = ifelse(!is.na(d), paste0(c, " ", d), x.var)) %>% #,
#   mutate(y.distribution = ifelse(b == "lognormal", b, NA), # Determine FRic distributions
#          y.distribution = ifelse(b == "gamma", b, y.distribution)) %>% 
#   relocate(x.var, x.var.type, y.var, y.distribution, PCA, NumberOfPCs, CriteriaOfPCoA, estimate, lower.CI,
#            upper.CI, significant, warning, error, rhat, problematic_rhat_count, .before = ) %>% # Tidy up dataframe
#   arrange(y.var, x.var) %>% 
#   dplyr::select(-c(model, a, b, c, d))

# Ouptut the results to .csv
write.csv(results.output, file = paste0("data/model_outputs_new/",
                                        "results_output", pc.filepath, ".csv"), row.names = FALSE)


# COMBINE THE CONVERGENCE OUTPUTS INTO ONE ----

# Remove results output file so don't add to it already
unlink("data/model_outputs_new/results_output_full*")

# Generate vector of all the convergence outputs
filepaths.convergence <- paste0("data/model_outputs_new/",
                                list.files("data/model_outputs_new/",
                                           pattern = "results_output*"))

# Create overall output dataframe
results.output.full <- data.frame()

# Run loop to append all outputs to a full output
for (i in filepaths.convergence){
  
  # Load in output
  results.output.indiv <- read.csv(paste0(i))
  
  # Append to overall dataframe
  results.output.full <- rbind(results.output.full, results.output.indiv)
  
  # Arrange the output by columns
  results.output.full <- results.output.full %>% 
    arrange(PCA, CriteriaOfPCoA, NumberOfPCs, y.var, x.var)
  
} # End of for loop

# Save the output as a single .csv
write.csv(results.output.full,
          file = "data/model_outputs_new/results_output_full.csv",
          row.names = FALSE)


# CHECK FOR CORRELATIONS BETWEEN PCs ----

# Only run if PCAs being used
if (run.with.PCA == TRUE){
  
  # Run bayesian model to compare FRic/FEve values when calculated with different numbers of PCs
  bayesian.pc.count.comparison(run.FRic = FALSE,
                               run.FEve = FALSE,
                               pc.count.1 = 3,
                               pc.count.2 = 4,
                               colour.by = "PlotDominatingFG", # Must be discrete
                               chains = 2,
                               cores = 4,
                               warmup = 500,
                               iterations = 2000,
                               delta = 0.95,
                               treedepth = 13)
  
  # Run linear model to compare FRic/FEve values when calculated with different numbers of PCs
  linear.pc.count.comparison(run.model = FALSE,
                             pc.count.1 = 3,
                             pc.count.2 = 4,
                             species.count = 2, # 1 or 2 more species than traits
                             colour.by = "SR",
                             colour.discrete = FALSE)
  
}
