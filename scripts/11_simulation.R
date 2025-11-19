# 11 - Simulation Study
# December 2022, adapted November 2025


# LOAD PACKAGES, THEMES, FUNCTIONS & TRAIT SELECTION ----

# Load packages
library(tidyverse)
library(gridExtra)
library(devtools)
library(FD)

# Load themes and functions
source("scripts/00_ggplot_themes.R")
source("scripts/08_FD_FUNCTION.R")


# Use PCA or not (now deemed that PCA is incorrect; default = FALSE)
run.with.PCA <- FALSE

# Determine number of PCs to use (if running with PCA)
pc.count <- 3

# Determine 'm' input for FD function (if running without PCA) - for "max" use ""
m.input <- "min"

# Determine filepaths and PC counts depending on whether using PCA or not
if (run.with.PCA == TRUE){
  
  # Generate filepath for saving outputs
  pc.filepath <- paste0("_PCA_pc", pc.count)
  
} else { # End of run.with.PCA == TRUE if
  
  # Now run loop for if max or min
  if (m.input == "min"){
    
    # Generate filepath for saving outputs
    pc.filepath <- paste0("_nonPCA_", m.input)
    
  } else { # End of m.input == "min" if
    
    # Generate filepath for saving outputs
    pc.filepath <- "_nonPCA"
    
  } # End of m.input == "min" else
  
} # End of run.with.PCA == TRUE else


# IMPORT DATA ----

# Import manually created simulation dataset(s)
simul.ALL <- read.csv("data/input_simulation_plots_FEve_3.csv")

# MUTATE DATA ----

# Add an ID column
ID.ALL <- simul.ALL %>% 
  mutate(ID = row.names(.)) %>% 
  relocate(ID, .before =)

# Create simulation cover input
cover.ALL <- ID.ALL %>% 
  dplyr::select(PLOT, ID, COVER) %>% 
  pivot_wider(names_from = "ID", values_from = "COVER") %>% 
  replace(is.na(.), 0) %>% 
  dplyr::select(-PLOT)

# Create simulation trait input
trait.ALL <- ID.ALL %>%
  dplyr::select(TRAIT1, TRAIT2, TRAIT3)


# RUN FD CALCULATIONS ----

# Create a vector of simulation rownames
rownames.ALL <- rownames(cover.ALL)

# Rename ALL input objects to match that of ITEX/traits
itex.input <- cover.ALL
trait.input <- trait.ALL
itex_rownames_simulation_ALL <- rownames.ALL

# Run FD calculations on simulated FRic data
FD.calc(itex_rownames_simulation_ALL)

output.ALL <- read.csv("data/output_08_fd_output_simulation_ALL_nonPCA_min.csv") %>%
  mutate(PLOT = rownames(.)) %>% 
  relocate(PLOT, .before = ) %>% 
  dplyr::select(PLOT, FRic, FEve, FDis) %>% 
  mutate(PLOT = as.numeric(PLOT))

# Plot type
plot.type <- simul.FEve %>% 
  dplyr::select(PLOT, DESCRIPTION) %>% 
  unique()

# Join and tidy outputs
output.tidy <- output.ALL %>% 
  left_join(., plot.type, by = c("PLOT" = "PLOT")) %>% 
  pivot_longer(names_to = "Metric", values_to = "Value", cols = c("FRic", "FEve", "FDis")) %>% 
  mutate(PLOT = ifelse(PLOT %in% c(1, 3, 5), "Before", "After"),
         Value = round(Value, digits = 2)) %>% 
  pivot_wider(names_from = PLOT, values_from = Value) %>% 
  mutate(Difference = ((After - Before)/Before)*100,
         Difference = round(Difference, digits = 2)) %>% 
  rename(`Plot 1` = Before, `Plot 2` = After) %>% 
  dplyr::select(-DESCRIPTION)

# Export output
write.csv(output.tidy, file = "data/output_08_fd_output_simulation_ALL_figures.csv", row.names = FALSE)
 

# PRODUCE PLOTS ----

# Create vector of six plot numbers per simulation
simulated.plots <- c(1, 2, 3, 4, 5, 6)

# Create temporary data point for bodging the scale
scale.bodge <- data.frame("SPECIES" = "A", "TRAIT_VALUE" = 200)

# Create empty list to output simulated FEve plots into
output.plots.ALL <- list()

# Run loop over each plot to create plots
for (i in simulated.plots){
  
  # Create FEve plots
  (plot.ALL <- ggplot(filter(simul.ALL, PLOT == i), aes(group = 1)) +
     
     geom_point(data = scale.bodge, aes(x = SPECIES, y = 0.5*TRAIT_VALUE), colour = "#FFFFFF") +

     geom_bar(aes(x = SPECIES, y = COVER), stat = "identity", colour = "#000000", alpha = 0.65,
              fill = ifelse(i == 1, "#E685C3",
                            ifelse(i == 2, "#F5D7EA",
                                   ifelse(i == 3, "#67B867",
                                          ifelse(i == 4, "#D1F0D1",
                                                 ifelse(i == 5, "#1C86EE", "#B1CCE8")))))) +
                            
     geom_line(aes(x = SPECIES, y = 0.5*TRAIT1), stat = "identity", colour = "#C7C7C7", size = 1.5) +
     geom_line(aes(x = SPECIES, y = 0.5*TRAIT2), stat = "identity", colour = "#969696", size = 1.5) +
     geom_line(aes(x = SPECIES, y = 0.5*TRAIT3), stat = "identity", colour = "#6B6B6B", size = 1.5) +
     
     geom_point(aes(x = SPECIES, y = 0.5*TRAIT1), colour = "#000000", fill = "#000000", size = 3, shape = 21) +
     geom_point(aes(x = SPECIES, y = 0.5*TRAIT2), colour = "#000000", fill = "#000000", size = 3, shape = 22) +
     geom_point(aes(x = SPECIES, y = 0.5*TRAIT3), colour = "#000000", fill = "#000000", size = 3, shape = 24) +
     
     scale_y_continuous(sec.axis = sec_axis(~.*2, name = "Trait Units\n")) +
     
     labs(title = paste0("Plot ", i),
          subtitle = paste0(ifelse(i %in% c(1,2), "Changing Species Richness",
                                   ifelse(i %in% c(3,4), "Changing Trait Values", "Changing Cover"))),
          x = "\nSpecies",
          y = "Cover (%)\n") +
     theme_1())

  # Append plot to list
  output.plots.ALL[[paste0("Plot: ", i)]] <- plot.ALL
  
}

# Create a panel of the outputted plots
panel.ALL <- arrangeGrob(output.plots.ALL[[1]], output.plots.ALL[[2]], output.plots.ALL[[3]], 
                         output.plots.ALL[[4]], output.plots.ALL[[5]], output.plots.ALL[[6]], 
                         ncol = 2, widths = c(1,1), heights = c(1,1,1), top = "")

# Export the panel to image
ggsave(panel.ALL, file = "figures/outputs_new/simulation_ALL.png",
       width = 20, height = 25)

