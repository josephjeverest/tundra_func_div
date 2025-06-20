# 07 - Trait data reduction and plot removal
# Joseph Everest
# January 2024


# LOAD PACKAGES ----

# Load the required packages
library(tidyverse)


# DETERMINE NUMBER OF PCs TO RETAIN ----

# Determine number of PCs
pc.count <- 3


# IMPORT DATA ----

# Load in trait values by species
traits <- read.csv("data/output_06_traits_final.csv") %>% 
  dplyr::select(-X)


# RUN PCA ON DATA ----

# Restructure the trait data to run dimensionality reduction
traits.restructure <- traits %>% 
  column_to_rownames(var = "ID")

# Run PCA on the restructured trait data
traits.PCA <- prcomp(traits.restructure, scale = TRUE)

# Reverse sign of results as eigenvectors negative in R by default
traits.PCA$x <- -1*traits.PCA$x

# Get results of the spectra PCA
traits.PCA.results <- data.frame(traits.PCA$x)

# Calculate variance explained by each principle component
traits.variance <- traits.PCA$sdev^2 / sum(traits.PCA$sdev^2)

# Retain only the required PCs
traits.PCA.selected <- traits.PCA.results[, 1:pc.count]

# Restructure data back to the original format for mantel test calculations
traits.PCA.restructure <- traits.PCA.selected %>%
  rownames_to_column(.) 


# EXPORT PCA TRAIT DATA ----

# Export PCAed trait data to .csv
write.csv(traits.PCA.restructure, file = paste0("data/output_07_traits_PCA_pc", pc.count, ".csv"))
