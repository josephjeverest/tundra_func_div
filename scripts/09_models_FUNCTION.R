# 09 - FUNCTIONS: Spatial Models
# Joseph Everest
# December 2021, adapted February 2022, November 2022, July 2023, February 2024, January 2025


# FUNCTION: SPATIAL BAYESIAN MODELS FOR FD WITH A CONTINUOUS X-VARIABLE ----

bayesian.spatial.continuous <- function(run.FRic, censored, run.FEve, quadratic, run.SR, run.FDis, x.var, data, colour.by,
                                        colour.discrete, chains, cores, warmup, iterations, delta, treedepth){
  # if (1 = 1){ # testing why breaks
  # FRic ----
  
  # Produce dataframe for running models on
  input.data <- data %>% 
    dplyr::select(x.var, FRic, FEve, SR, FDis, SurveyedArea, gridcell, SUBSITE) %>% 
    rename(x_variable = x.var) # Rename first column to x_variable
  
  # Tidy import data for running censored models
  input.data.2 <- input.data %>% 
    mutate(RETAIN = ifelse(SR >= 4 & is.na(FRic), FALSE, TRUE)) %>% # Remove erroneous NA value for SR == 4 
    filter(RETAIN == TRUE) %>% 
    dplyr::select(-RETAIN) %>% 
    mutate(log_FRic = log(FRic)) # For gaussian distribution
  
  # Determine the median value of FRic for when SR < 4 (i.e. FRic is NA)
  median.FRic.4 <- median(filter(input.data.2, SR == 4)$log_FRic)
  
  # Create columns for modelling
  input.data.3 <- input.data.2 %>% 
    mutate(FRic2 = ifelse(SR >= 1 & SR < 4, 0.000000000000000001, log_FRic), # Create interval censored columns
           FRic3 = ifelse(SR >= 1 & SR < 4, median.FRic.4, log_FRic),
           cen2 = ifelse(SR >= 1 & SR < 4, "interval", "none")) %>% 
    mutate(FRic_na = ifelse(cen2 == "none", log_FRic, NA))
  
  # Modify combo.latest for colouring points
  combo.colour <- combo.latest %>% 
    rename(x_variable = x.var) %>% # Rename first column to x_variable
    rename(colour.variable = colour.by)
  
  # Determine which version of FRic to run/import (censored or not)
  if (run.FRic == TRUE & censored == "No"){
    
    # Test check
    cat("running fric \n")
    
    # Function for running the Bayesian model
    FRic.mod <- brm(FRic ~ x_variable + log(SurveyedArea) + (1 | gridcell / SUBSITE),
                    data = input.data, family = FRic.distribution, chains = chains,
                    warmup = warmup, iter = iterations, cores = cores,
                    control = list(adapt_delta = delta,
                                   max_treedepth = treedepth))
    
    # Export model output
    save(FRic.mod, file = paste0("data/model_outputs_new/m_FRic_",
                                 FRic.distribution, "_", x.var, pc.filepath, ".RData"))
    
  } # End of run and not censored
  
  if (run.FRic == FALSE & censored == "No"){
    
    # Load in model output
    FRic.mod <- get(load(paste0("data/model_outputs_new/m_FRic_",
                                FRic.distribution, "_", x.var, pc.filepath, ".RData")))
    
  } # End of import and not censored
  
  if (run.FRic == TRUE & censored == "Yes"){
    
    # Define the stanvars
    mean_FRic <- mean(input.data.3$FRic_na, na.rm = TRUE)
    sd_FRic <- sd(input.data.3$FRic_na, na.rm = TRUE)
    stanvars <- stanvar(mean_FRic, name = "mean_FRic") + stanvar(sd_FRic,   name = "sd_FRic")
    
    # Run the censored brms() model
    FRic.mod <- brm(FRic2 | cens(cen2, FRic3) ~ x_variable + log(SurveyedArea) + (1 | gridcell / SUBSITE),
                    data = input.data.3, family = "gaussian",
                    chains = chains, warmup = warmup, iter = iterations,
                    cores = cores, control = list(adapt_delta = delta,
                                                  max_treedepth = treedepth))
    
    # Export model output
    save(FRic.mod, file = paste0("data/model_outputs_new/m_FRic_",
                                 "censored_", x.var, pc.filepath, ".RData"))
    
  } # End of run and censored
  
  if (run.FRic == FALSE & censored == "Yes"){
    
    # Load in model output
    FRic.mod <- get(load(paste0("data/model_outputs_new/m_FRic_",
                                "censored_", x.var, pc.filepath, ".RData")))
    
  } # End of import and censored
  

  # Process predictions differently dependent on whether censored or not
  if (censored == "No"){
    
    # Convert model output into a dataframe (with 4.d.p.)
    FRic.df <- brms_SummaryTable(FRic.mod, formatOptions = list(digits = 4, nsmall = 4), round = 4)
    
    # Extract the confidence intervals as a list for use in the plotting
    FRic.df.ci <- FRic.df %>% 
      filter(Covariate %in% c("x_variable")) %>% 
      dplyr::select("l-95% CI", "u-95% CI") %>% 
      pivot_longer(cols = 1:2, names_to = "Interval", values_to = "Value")
    
    # Save the confidence intervals as a list and remove the intermediate dataframe
    FRic.ci <- as.list(FRic.df.ci$Value) # Adding values to the list
    names(FRic.ci) <- FRic.df.ci$Interval # Adding names to the values
    rm(FRic.df, FRic.df.ci) # Remove unnecessary dataframes of summary and CIs
    
    # Range to predict variable
    range.to.predict <- paste0("x_variable [", min(input.data$x_variable), ":", max(input.data$x_variable), "]")
    
    # Extract the prediction data frame
    FRic.pred <- ggpredict(FRic.mod, terms = range.to.predict, back.transform = TRUE)
    
  } else { # End of censored == "No"
    
    # Convert model output into a dataframe (with 4.d.p.)
    FRic.df <- brms_SummaryTable(FRic.mod, formatOptions = list(digits = 4, nsmall = 4), round = 4)
    
    # Extract the confidence intervals as a list for use in the plotting
    FRic.df.ci <- FRic.df %>% 
      filter(Covariate %in% c("x_variable")) %>%
      # rename(lowerCI = "l-95% CI",
      #        upperCI = "u-95% CI") %>% 
      # mutate(Estimate = round(exp(as.numeric(Estimate)), digits = 4),
      #        lowerCI = round(exp(as.numeric(lowerCI)), digits = 4),
      #        upperCI = round(exp(as.numeric(upperCI)), digits = 4)) %>% 
      # rename("l-95% CI" = lowerCI,
      #        "u-95% CI" = upperCI) %>% 
      dplyr::select("l-95% CI", "u-95% CI") %>% 
      pivot_longer(cols = 1:2, names_to = "Interval", values_to = "Value")
    
    # Save the confidence intervals as a list and remove the intermediate dataframe
    FRic.ci <- as.list(FRic.df.ci$Value) # Adding values to the list
    names(FRic.ci) <- FRic.df.ci$Interval # Adding names to the values
    rm(FRic.df, FRic.df.ci) # Remove unnecessary dataframes of summary and CIs
    
    # Range to predict variable
    range.to.predict <- paste0("x_variable [", min(input.data.3$x_variable), ":", max(input.data.3$x_variable), "]")
    
    # Extract the prediction data frame and exponentiate the outputs
    FRic.pred <- ggpredict(FRic.mod, terms = range.to.predict, back.transform = FALSE) %>% 
      mutate(conf.low = exp(conf.low),
             conf.high = exp(conf.high),
             predicted = exp(predicted))
    
  } # End of censored == "Yes
  
  # Extract dataframe with only the x-variable for plotting
  data.x.var.only <- dplyr::select(data, x.var)
  
  # Add in if loop for if the colour variable is also SR
  if (colour.by == "FRic"){
    
    # Add column back into combo.latest
    combo.colour <- combo.colour %>% 
      mutate(FRic = colour.variable)
    
  }
  
  # Plot model outputs
  (FRic.plot <- ggplot(FRic.pred) +
      geom_point(data = combo.colour, aes(x = data.x.var.only[,1], y = FRic, fill = colour.variable), # Adds original FRic vs LAT data points and colours by region
                 colour = c("#000000"), alpha = 0.5, shape = 21, size = 3) +
      geom_line(aes(x = x, y = predicted), colour = "grey10") + # Adds line for predicted values of FRic vs LAT
      geom_ribbon(aes(x = x, ymin = conf.low, ymax = conf.high),
                  fill = "grey10", alpha = 0.27, colour = "#000000", linetype = 2) + # Adds c.intervals for predictions as ribbon
      scale_fill_viridis(option = FRic.colour, begin = 0.2, end = 1, direction = -1, discrete = colour.discrete,
                         guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")) +
      labs(y = "FRic \n",
           x = paste0("\n ", x.var),
           fill = paste0(colour.by, "\n"),
           title = paste0("CIs: ", FRic.ci$"l-95% CI", " to ", FRic.ci$"u-95% CI"),
           subtitle = paste(ifelse((FRic.ci$"l-95% CI" < 0 & FRic.ci$"u-95% CI" < 0) | # Automatically asigns significance in subtitle based on CIs
                                     (FRic.ci$"l-95% CI" > 0 & FRic.ci$"u-95% CI" > 0),
                                   "Significant: CIs do NOT span 0", "Not Significant: CIs span 0"))) +
      theme_1() +
      theme(legend.position = "right",
            plot.subtitle = element_text(colour = ifelse((FRic.ci$"l-95% CI" < 0 & FRic.ci$"u-95% CI" < 0) | # Automatically colours subtitle based on significance from CIs
                                                           (FRic.ci$"l-95% CI" > 0 & FRic.ci$"u-95% CI" > 0),
                                                         paste("#006400"), paste("#8b0000")))))
  # FEve ----
  
  # Produce dataframe for running models on
  input.data <- data %>% 
    dplyr::select(x.var, FRic, FEve, SR, FDis, SurveyedArea, gridcell, SUBSITE) %>% 
    rename(x_variable = x.var) 
  
  # Determine mean value for centering
  mean_x_variable <- mean(input.data$x_variable)

  # Center the x_variable
  input.data.centred <- input.data %>% 
    mutate(centred_x_variable = x_variable - mean_x_variable) %>% 
    data.frame(.)
  
  # Determine which version of FRic to run/import (censored or not)
  if (run.FEve == TRUE & quadratic == "No"){
    
    # Function for running the Bayesian model
    FEve.mod <- brm(FEve ~ x_variable + log(SurveyedArea) + (1 | gridcell / SUBSITE),
                    data = input.data, family = FEve.distribution, chains = chains,
                    warmup = warmup, iter = iterations, cores = cores,
                    control = list(adapt_delta = delta,
                                   max_treedepth = treedepth))
    
    # Export model output
    save(FEve.mod, file = paste0("data/model_outputs_new/m_FEve_",
                                 x.var, pc.filepath, ".RData"))
    
  } # End of run and not censored
  
  if (run.FEve == FALSE & quadratic == "No"){
    
    # Load in model output
    FEve.mod <- get(load(paste0("data/model_outputs_new/m_FEve_",
                                x.var, pc.filepath, ".RData")))
    
  } # End of import and not censored
  
  if (run.FEve == TRUE & quadratic == "Yes"){
    
    # Function for running the Bayesian model
    FEve.mod <- brm(FEve ~ centred_x_variable + I(centred_x_variable^2) + log(SurveyedArea) + (1 | gridcell / SUBSITE),
                    data = input.data.centred, family = FEve.distribution, chains = chains,
                    warmup = warmup, iter = iterations, cores = cores,
                    control = list(adapt_delta = delta,
                                   max_treedepth = treedepth))
    
    # Export model output
    save(FEve.mod, file = paste0("data/model_outputs_new/m_FEve_",
                                 "quadratic_", x.var, pc.filepath, ".RData"))
    
  } # End of run and censored
  
  if (run.FEve == FALSE & quadratic == "Yes"){
    
    # Export model output
    FEve.mod <- get(load(paste0("data/model_outputs_new/m_FEve_",
                                "quadratic_", x.var, pc.filepath, ".RData")))
    
  } # End of import and censored
  
  # Extract dataframe with only the x-variable for plotting
  data.x.var.only <- dplyr::select(data, x.var)
  
  # Add in if loop for if the colour variable is also SR
  if (colour.by == "FEve"){
    
    # Add column back into combo.latest
    combo.colour <- combo.colour %>% 
      mutate(FEve = colour.variable)
    
  }
  
  # Run a second if loop for generating CIs and predictions dependent on whether quadratic or not
  if (quadratic == "No"){
    
    # Convert model output into a dataframe (with 4.d.p.)
    FEve.df <- brms_SummaryTable(FEve.mod, formatOptions = list(digits = 4, nsmall = 4), round = 4)
    
    # Extract the confidence intervals as a list for use in the plotting
    FEve.df.ci <- FEve.df %>% 
      filter(Covariate %in% c("x_variable")) %>% 
      dplyr::select("l-95% CI", "u-95% CI") %>% 
      pivot_longer(cols = 1:2, names_to = "Interval", values_to = "Value")
    
    # Save the confidence intervals as a list and remove the intermediate dataframe
    FEve.ci <- as.list(FEve.df.ci$Value) # Adding values to the list
    names(FEve.ci) <- FEve.df.ci$Interval # Adding names to the values
    rm(FEve.df, FEve.df.ci) # Remove unnecessary dataframes of summary and CIs  
    
    # Range to predict variable
    range.to.predict <- paste0("x_variable [", min(input.data$x_variable), ":", max(input.data$x_variable), "]")
    
    # Extract the prediction data frame
    FEve.pred <- ggpredict(FEve.mod, terms = range.to.predict)
    
    # Plot model outputs
    (FEve.plot <- ggplot(FEve.pred) +
        geom_point(data = combo.colour, aes(x = data.x.var.only[,1], y = FEve, fill = colour.variable), # Adds original FEve vs LAT data points and colours by region
                   colour = c("#000000"), alpha = 0.5, shape = 21, size = 3) +
        geom_line(aes(x = x, y = predicted), colour = "grey10") + # Adds line for predicted values of FEve vs LAT
        geom_ribbon(aes(x = x, ymin = conf.low, ymax = conf.high),
                    fill = "grey10", alpha = 0.27, colour = "#000000", linetype = 2) + # Adds c.intervals for predictions as ribbon
        scale_fill_viridis(option = FEve.colour, begin = 0.2, end = 1, direction = -1, discrete = colour.discrete,
                           guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")) +
        labs(y = "FEve \n",
             x = paste0("\n ", x.var),
             fill = paste0(colour.by, "\n"),
             title = paste0("CIs: ", FEve.ci$"l-95% CI", " to ", FEve.ci$"u-95% CI"),
             subtitle = paste(ifelse((FEve.ci$"l-95% CI" < 0 & FEve.ci$"u-95% CI" < 0) | # Automatically asigns significance in subtitle based on CIs
                                       (FEve.ci$"l-95% CI" > 0 & FEve.ci$"u-95% CI" > 0),
                                     "Significant: CIs do NOT span 0", "Not Significant: CIs span 0"))) +
        theme_1() +
        theme(legend.position = "right",
              plot.subtitle = element_text(colour = ifelse((FEve.ci$"l-95% CI" < 0 & FEve.ci$"u-95% CI" < 0) | # Automatically colours subtitle based on significance from CIs
                                                             (FEve.ci$"l-95% CI" > 0 & FEve.ci$"u-95% CI" > 0),
                                                           paste("#006400"), paste("#8b0000")))))
    
  } else { # End of non-quadratic if statement
    
    # Convert model output into a dataframe (with 4.d.p.)
    FEve.df <- brms_SummaryTable(FEve.mod, formatOptions = list(digits = 4, nsmall = 4), round = 4)
    
    # Extract the confidence intervals as a list for use in the plotting
    FEve.df.ci <- FEve.df %>% 
      filter(Covariate %in% c("centred_x_variable")) %>% 
      dplyr::select("l-95% CI", "u-95% CI") %>% 
      pivot_longer(cols = 1:2, names_to = "Interval", values_to = "Value")
    
    # Save the confidence intervals as a list and remove the intermediate dataframe
    FEve.ci <- as.list(FEve.df.ci$Value) # Adding values to the list
    names(FEve.ci) <- FEve.df.ci$Interval # Adding names to the values
    rm(FEve.df, FEve.df.ci) # Remove unnecessary dataframes of summary and CIs  
    
    # Create template predictions dataframe
    FEve.pred.df = data.frame(centred_x_variable = seq(0 - mean_x_variable, 100 - mean_x_variable), SurveyedArea = 1)
    
    # Add predictions to template datframe
    FEve.pred = add_epred_draws(FEve.pred.df, FEve.mod, re_formula = NA) %>% 
      mutate(x_variable = centred_x_variable + mean_x_variable) # Uncentre the x variable
    
    # Plot model outputs
    (FEve.plot <- ggplot(data = FEve.pred, aes(x = x_variable, y = .epred)) +
        geom_point(data = input.data, aes(x = x_variable, y = FEve, fill = colour.by)) +
        stat_lineribbon(aes(y = .epred), .width = 0.95, fill = "grey10",
                        alpha = 0.27, colour = "#000000", linetype = 2) +
        # scale_fill_viridis(option = FEve.colour, begin = 0.2, end = 1, direction = -1, discrete = colour.discrete,
        #                    guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")) +
        labs(y = "FEve \n",
             x = paste0("\n ", x.var),
             fill = paste0(colour.by, "\n"),
             title = paste0("CIs: ", FEve.ci$"l-95% CI", " to ", FEve.ci$"u-95% CI"),
             subtitle = paste(ifelse((FEve.ci$"l-95% CI" < 0 & FEve.ci$"u-95% CI" < 0) | # Automatically asigns significance in subtitle based on CIs
                                       (FEve.ci$"l-95% CI" > 0 & FEve.ci$"u-95% CI" > 0),
                                     "Significant: CIs do NOT span 0", "Not Significant: CIs span 0"))) +
        theme_1() +
        theme(legend.position = "right",
              plot.subtitle = element_text(colour = ifelse((FEve.ci$"l-95% CI" < 0 & FEve.ci$"u-95% CI" < 0) | # Automatically colours subtitle based on significance from CIs
                                                             (FEve.ci$"l-95% CI" > 0 & FEve.ci$"u-95% CI" > 0),
                                                           paste("#006400"), paste("#8b0000")))))
    
  } # End of quadratic else statement
  
  # SR ----
  
  # REMAKE the dataframe for running models on SR
  input.data <- data %>% 
    dplyr::select(x.var, FRic, FEve, SR, SurveyedArea, gridcell, SUBSITE) %>% 
    rename(x_variable = x.var) # Rename first column to x_variable
  
  # If you want to run model for first time
  if (run.SR == TRUE){
    
    # Function for running the Bayesian model
    SR.mod <- brm(SR ~ x_variable + log(SurveyedArea) + (1 | gridcell / SUBSITE),
                  data = input.data, family = SR.distribution, chains = chains,
                  warmup = warmup, iter = iterations, cores = cores,
                  control = list(adapt_delta = delta,
                                 max_treedepth = treedepth))
    
    # Export model output
    save(SR.mod, file = paste0("data/model_outputs_new/m_SR_",
                               x.var, pc.filepath, ".RData"))
    
  } else {
    
    # Load in model output
    SR.mod <- get(load(paste0("data/model_outputs_new/m_SR_",
                              x.var, pc.filepath, ".RData")))
    
  }
  
  # Convert model output into a dataframe (with 4.d.p.)
  SR.df <- brms_SummaryTable(SR.mod, formatOptions = list(digits = 4, nsmall = 4), round = 4)
  
  # Extract the confidence intervals as a list for use in the plotting
  SR.df.ci <- SR.df %>% 
    filter(Covariate %in% c("x_variable")) %>% 
    dplyr::select("l-95% CI", "u-95% CI") %>% 
    pivot_longer(cols = 1:2, names_to = "Interval", values_to = "Value")
  
  # Save the confidence intervals as a list and remove the intermediate dataframe
  SR.ci <- as.list(SR.df.ci$Value) # Adding values to the list
  names(SR.ci) <- SR.df.ci$Interval # Adding names to the values
  rm(SR.df, SR.df.ci) # Remove unnecessary dataframes of summary and CIs  
  
  # Range to predict variable
  range.to.predict <- paste0("x_variable [", min(input.data$x_variable), ":", max(input.data$x_variable), "]")
  
  # Extract the prediction data frame
  SR.pred <- ggpredict(SR.mod, terms = range.to.predict,
                       # type = "random", condition = c(gridcell = 1266), allow_new_levels = TRUE,
                       back.transform = FALSE)
  
  # Add in if loop for if the colour variable is also SR
  if (colour.by == "SR"){
    
    # Add column back into combo.latest
    combo.colour <- combo.colour %>% 
      mutate(SR = colour.variable)
    
  }
  
  # Plot model outputs
  (SR.plot <- ggplot(SR.pred) +
      geom_point(data = combo.colour, aes(x = x_variable, y = SR, fill = colour.variable), # Adds original SR vs LAT data points and colours by region
                 colour = c("#000000"), alpha = 0.5, shape = 21, size = 3) +
      geom_line(aes(x = x, y = predicted), colour = "grey10") + # Adds line for predicted values of SR vs LAT
      geom_ribbon(aes(x = x, ymin = conf.low, ymax = conf.high),
                  fill = "grey10", alpha = 0.27, colour = "#000000", linetype = 2) + # Adds c.intervals for predictions as ribbon
      scale_fill_viridis(option = SR.colour, begin = 0.2, end = 1, direction = -1, discrete = colour.discrete,
                         guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")) +
      labs(y = "SR \n",
           x = paste0("\n ", x.var),
           fill = paste0(colour.by, "\n"),
           title = paste0("CIs: ", SR.ci$"l-95% CI", " to ", SR.ci$"u-95% CI"),
           subtitle = paste(ifelse((SR.ci$"l-95% CI" < 0 & SR.ci$"u-95% CI" < 0) | # Automatically asigns significance in subtitle based on CIs
                                     (SR.ci$"l-95% CI" > 0 & SR.ci$"u-95% CI" > 0),
                                   "Significant: CIs do NOT span 0", "Not Significant: CIs span 0"))) +
      theme_1() +
      theme(legend.position = "right",
            plot.subtitle = element_text(colour = ifelse((SR.ci$"l-95% CI" < 0 & SR.ci$"u-95% CI" < 0) | # Automatically colours subtitle based on significance from CIs
                                                           (SR.ci$"l-95% CI" > 0 & SR.ci$"u-95% CI" > 0),
                                                         paste("#006400"), paste("#8b0000")))))
  
  # FDis ----
  
  # REMAKE the dataframe for running models on FDis
  input.data <- data %>% 
    dplyr::select(x.var, FRic, FEve, SR, FDis, SurveyedArea, gridcell, SUBSITE) %>% 
    rename(x_variable = x.var) # Rename first column to x_variable
  
  # If you want to run model for first time
  if (run.FDis == TRUE){
    
    # Function for running the Bayesian model
    FDis.mod <- brm(FDis ~ x_variable + log(SurveyedArea) + (1 | gridcell / SUBSITE),
                    data = input.data, family = FDis.distribution, chains = chains,
                    warmup = warmup, iter = iterations, cores = cores,
                    control = list(adapt_delta = delta,
                                   max_treedepth = treedepth))
    
    # Export model output
    save(FDis.mod, file = paste0("data/model_outputs_new/m_FDis_",
                               x.var, pc.filepath, ".RData"))
    
  } else {
    
    # Load in model output
    FDis.mod <- get(load(paste0("data/model_outputs_new/m_FDis_",
                                x.var, pc.filepath, ".RData")))
    
  }
  
  # Convert model output into a dataframe (with 4.d.p.)
  FDis.df <- brms_SummaryTable(FDis.mod, formatOptions = list(digits = 4, nsmall = 4), round = 4)
  
  # Extract the confidence intervals as a list for use in the plotting
  FDis.df.ci <- FDis.df %>% 
    filter(Covariate %in% c("x_variable")) %>%
    dplyr::select("l-95% CI", "u-95% CI") %>%
    pivot_longer(cols = 1:2, names_to = "Interval", values_to = "Value")
  
  # Save the confidence intervals as a list and remove the intermediate dataframe
  FDis.ci <- as.list(FDis.df.ci$Value) # Adding values to the list
  names(FDis.ci) <- FDis.df.ci$Interval # Adding names to the values
  rm(FDis.df, FDis.df.ci) # Remove unnecessary dataframes of summary and CIs  
  
  # Range to predict variable
  range.to.predict <- paste0("x_variable [", min(input.data$x_variable), ":", max(input.data$x_variable), "]")
  
  # Extract the prediction data frame
  FDis.pred <- ggpredict(FDis.mod, terms = range.to.predict, back.transform = FALSE)
  
  # Add in if loop for if the colour variable is also SR
  if (colour.by == "FDis"){
    
    # Add column back into combo.latest
    combo.colour <- combo.colour %>% 
      mutate(FDis = colour.variable)
    
  }
  
  # Plot model outputs
  (FDis.plot <- ggplot(FDis.pred) +
      geom_point(data = combo.colour, aes(x = x_variable, y = FDis, fill = colour.variable), # Adds original SR vs LAT data points and colours by region
                 colour = c("#000000"), alpha = 0.5, shape = 21, size = 3) +
      geom_line(aes(x = x, y = predicted), colour = "grey10") + # Adds line for predicted values of SR vs LAT
      geom_ribbon(aes(x = x, ymin = conf.low, ymax = conf.high),
                  fill = "grey10", alpha = 0.27, colour = "#000000", linetype = 2) + # Adds c.intervals for predictions as ribbon
      scale_fill_viridis(option = FDis.colour, begin = 0.2, end = 1, direction = -1, discrete = colour.discrete,
                         guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")) +
      labs(y = "FDis \n",
           x = paste0("\n ", x.var),
           fill = paste0(colour.by, "\n"),
           title = paste0("CIs: ", FDis.ci$"l-95% CI", " to ", FDis.ci$"u-95% CI"),
           subtitle = paste(ifelse((FDis.ci$"l-95% CI" < 0 & FDis.ci$"u-95% CI" < 0) | # Automatically asigns significance in subtitle based on CIs
                                     (FDis.ci$"l-95% CI" > 0 & FDis.ci$"u-95% CI" > 0),
                                   "Significant: CIs do NOT span 0", "Not Significant: CIs span 0"))) +
      theme_1() +
      theme(legend.position = "right",
            plot.subtitle = element_text(colour = ifelse((FDis.ci$"l-95% CI" < 0 & FDis.ci$"u-95% CI" < 0) | # Automatically colours subtitle based on significance from CIs
                                                           (FDis.ci$"l-95% CI" > 0 & FDis.ci$"u-95% CI" > 0),
                                                         paste("#006400"), paste("#8b0000")))))
  
  # Create a panel of all three plots
  combined.panel <- grid.arrange(FRic.plot, FEve.plot, SR.plot, FDis.plot, ncol = 2)
  
  # Add in if statements for saving the output
  if (censored == "No" & quadratic == "No"){
    
    # Export panel
    ggsave(combined.panel, filename = paste0("figures/outputs_new/combined_",
                                             x.var, pc.filepath, ".png"), width = 15, height = 14)
    
  }
  
  # Add in if statements for saving the output
  if (censored == "No" & quadratic == "Yes"){
    
    # Export panel
    ggsave(combined.panel, filename = paste0("figures/outputs_new/combined_",
                                             x.var, "_quadratic", pc.filepath, ".png"), width = 15, height = 14)
    
  }
  
  # Add in if statements for saving the output
  if (censored == "Yes" & quadratic == "No"){
    
    # Export panel
    ggsave(combined.panel, filename = paste0("figures/outputs_new/combined_",
                                             x.var, "_censored", pc.filepath, ".png"), width = 15, height = 14)
    
  }
  
  # Add in if statements for saving the output
  if (censored == "Yes" & quadratic == "Yes"){
    
    # Export panel
    ggsave(combined.panel, filename = paste0("figures/outputs_new/combined_",
                                             x.var, "_censored_quadratic", pc.filepath, ".png"), width = 15, height = 14)
    
  }
  
}


# FUNCTION: ALL QUADRATIC SPATIAL BAYESIAN MODELS FOR FD WITH A CONTINUOUS X-VARIABLE ----

bayesian.spatial.quadratic <- function(run.FRic, censored, run.FEve, run.SR, run.FDis, x.var, data, colour.by,
                                       colour.discrete, chains, cores, warmup, iterations, delta, treedepth){
  
  # FRic ----
  
  if (censored == "Yes"){
    
    # Produce dataframe for running models on
    input.data <- data %>% 
      dplyr::select(x.var, FRic, FEve, SR, FDis, SurveyedArea, gridcell, SUBSITE) %>% 
      rename(x_variable = x.var) # Rename first column to x_variable
    
    # Tidy import data for running censored models
    input.data.2 <- input.data %>% 
      mutate(RETAIN = ifelse(SR >= 4 & is.na(FRic), FALSE, TRUE)) %>% # Remove erroneous NA value for SR == 4 
      filter(RETAIN == TRUE) %>% 
      dplyr::select(-RETAIN) %>% 
      mutate(log_FRic = log(FRic)) # For gaussian distribution
    
    # Determine the median value of FRic for when SR < 4 (i.e. FRic is NA)
    median.FRic.4 <- median(filter(input.data.2, SR == 4)$log_FRic)
    
    # Create columns for modelling
    input.data.3 <- input.data.2 %>% 
      mutate(FRic2 = ifelse(SR >= 1 & SR < 4, 0.000000000000000001, log_FRic), # Create interval censored columns
             FRic3 = ifelse(SR >= 1 & SR < 4, median.FRic.4, log_FRic),
             cen2 = ifelse(SR >= 1 & SR < 4, "interval", "none")) %>% 
      mutate(FRic_na = ifelse(cen2 == "none", log_FRic, NA))
    
    # Determine mean value for centering
    mean_x_variable <- mean(input.data.3$x_variable)
    
    # Center the x_variable
    input.data.centred <- input.data.3 %>% 
      mutate(centred_x_variable = x_variable - mean_x_variable) %>% 
      data.frame(.)
    
    if (run.FRic == TRUE){
      
      # Define the stanvars
      mean_FRic <- mean(input.data.centred$FRic_na, na.rm = TRUE)
      sd_FRic <- sd(input.data.centred$FRic_na, na.rm = TRUE)
      stanvars <- stanvar(input.data.centred, name = "mean_FRic") + stanvar(sd_FRic,   name = "sd_FRic")
      
      # Run the censored brms() model
      FRic.mod <- brm(FRic2 | cens(cen2, FRic3) ~ centred_x_variable + I(centred_x_variable^2) + log(SurveyedArea) + (1 | gridcell / SUBSITE),
                      data = input.data.centred, family = "gaussian",
                      chains = chains, warmup = warmup, iter = iterations,
                      cores = cores, control = list(adapt_delta = delta,
                                                    max_treedepth = treedepth))
      
      # Export model output
      save(FRic.mod, file = paste0("data/model_outputs_new/m_FRic_",
                                   "censored_quadratic_", x.var, pc.filepath, ".RData"))
      
    } else {
      
      # Load in model output
      FRic.mod <- get(load(paste0("data/model_outputs_new/m_FRic_",
                                  "censored_quadratic_", x.var, pc.filepath, ".RData")))
      
    } # End of import and censored
    
    # Convert model output into a dataframe (with 4.d.p.)
    FRic.df <- brms_SummaryTable(FRic.mod, formatOptions = list(digits = 4, nsmall = 4), round = 4)
    
    # Extract the confidence intervals as a list for use in the plotting
    FRic.df.ci <- FRic.df %>% 
      filter(Covariate %in% c("Icentred_x_variableE2")) %>% 
      dplyr::select("l-95% CI", "u-95% CI") %>% 
      pivot_longer(cols = 1:2, names_to = "Interval", values_to = "Value")
    
    # Save the confidence intervals as a list and remove the intermediate dataframe
    FRic.ci <- as.list(FRic.df.ci$Value) # Adding values to the list
    names(FRic.ci) <- FRic.df.ci$Interval # Adding names to the values
    rm(FRic.df, FRic.df.ci) # Remove unnecessary dataframes of summary and CIs  
    
    # Create template predictions dataframe
    FRic.pred.df = data.frame(centred_x_variable = seq(0 - mean_x_variable, 100 - mean_x_variable), SurveyedArea = 1)
    
    # Add predictions to template datframe
    FRic.pred = add_epred_draws(FRic.pred.df, FRic.mod, re_formula = NA) %>% 
      mutate(x_variable = centred_x_variable + mean_x_variable, # Uncentre the x variable
             ".epred" = exp(".epred"))
    
  } else { # End of if censored == TRUE, start of if censored == FALSE
    
    # Produce dataframe for running models on
    input.data <- data %>% 
      dplyr::select(x.var, FRic, FEve, SR, FDis, SurveyedArea, gridcell, SUBSITE) %>% 
      rename(x_variable = x.var) 
    
    # Determine mean value for centering
    mean_x_variable <- mean(input.data$x_variable)
    
    # Center the x_variable
    input.data.centred <- input.data %>% 
      mutate(centred_x_variable = x_variable - mean_x_variable) %>% 
      data.frame(.)
    
    
    if (run.FRic == TRUE){
      
      # Function for running the Bayesian model
      FRic.mod <- brm(FRic ~ centred_x_variable + I(centred_x_variable^2) + log(SurveyedArea) + (1 | gridcell / SUBSITE),
                      data = input.data.centred, family = FRic.distribution, chains = chains,
                      warmup = warmup, iter = iterations, cores = cores,
                      control = list(adapt_delta = delta,
                                     max_treedepth = treedepth))
      
      # Export model output
      save(FRic.mod, file = paste0("data/model_outputs_new/m_FRic_",
                                   "quadratic_", x.var, pc.filepath, ".RData"))
      
    } else {
      
      # Export model output
      FRic.mod <- get(load(paste0("data/model_outputs_new/m_FRic_",
                                  "quadratic_", x.var, pc.filepath, ".RData")))
      
    } # End of import and censored
    
    summary(FRic.mod)
    
    # Convert model output into a dataframe (with 4.d.p.)
    FRic.df <- brms_SummaryTable(FRic.mod, formatOptions = list(digits = 5, nsmall = 5), round = 5)
    
    # Extract the confidence intervals as a list for use in the plotting
    FRic.df.ci <- FRic.df %>% 
      filter(Covariate %in% c("Icentred_x_variableE2")) %>% 
      dplyr::select("l-95% CI", "u-95% CI") %>% 
      pivot_longer(cols = 1:2, names_to = "Interval", values_to = "Value")
    
    # Save the confidence intervals as a list and remove the intermediate dataframe
    FRic.ci <- as.list(FRic.df.ci$Value) # Adding values to the list
    names(FRic.ci) <- FRic.df.ci$Interval # Adding names to the values
    rm(FRic.df, FRic.df.ci) # Remove unnecessary dataframes of summary and CIs  
    
    # Create template predictions dataframe
    FRic.pred.df = data.frame(centred_x_variable = seq(0 - mean_x_variable, 100 - mean_x_variable), SurveyedArea = 1)
    
    # Add predictions to template datframe
    FRic.pred = add_epred_draws(FRic.pred.df, FRic.mod, re_formula = NA) %>% 
      mutate(x_variable = centred_x_variable + mean_x_variable) # Uncentre the x variable
    
  } # End of else censored == FALSE
  
  # Plot model outputs
  (FRic.plot <- ggplot(data = FRic.pred, aes(x = x_variable, y = .epred)) +
      geom_point(data = input.data, aes(x = x_variable, y = FRic, fill = "#53C367"),
                 shape = 21) +
      stat_lineribbon(aes(y = .epred), .width = 0.95, fill = "#53C367",
                      alpha = 0.27, colour = "#000000", linetype = 2) +
      labs(y = "FRic \n",
           x = paste0("\n ", x.var),
           fill = paste0(colour.by, "\n"),
           title = paste0("CIs: ", FRic.ci$"l-95% CI", " to ", FRic.ci$"u-95% CI"),
           subtitle = paste(ifelse((FRic.ci$"l-95% CI" < 0 & FRic.ci$"u-95% CI" < 0) | # Automatically asigns significance in subtitle based on CIs
                                     (FRic.ci$"l-95% CI" > 0 & FRic.ci$"u-95% CI" > 0),
                                   "Significant: CIs do NOT span 0", "Not Significant: CIs span 0"))) +
      theme_1() +
      theme(legend.position = "right",
            plot.subtitle = element_text(colour = ifelse((FRic.ci$"l-95% CI" < 0 & FRic.ci$"u-95% CI" < 0) | # Automatically colours subtitle based on significance from CIs
                                                           (FRic.ci$"l-95% CI" > 0 & FRic.ci$"u-95% CI" > 0),
                                                         paste("#006400"), paste("#8b0000")))))
  
  # FEve ----
  
  # Produce dataframe for running models on
  input.data <- data %>% 
    dplyr::select(x.var, FRic, FEve, SR, FDis, SurveyedArea, gridcell, SUBSITE) %>% 
    rename(x_variable = x.var) 
  
  # Determine mean value for centering
  mean_x_variable <- mean(input.data$x_variable)
  
  # Center the x_variable
  input.data.centred <- input.data %>% 
    mutate(centred_x_variable = x_variable - mean_x_variable) %>% 
    data.frame(.)
  
  
  if (run.FEve == TRUE){
    
    # Function for running the Bayesian model
    FEve.mod <- brm(FEve ~ centred_x_variable + I(centred_x_variable^2) + log(SurveyedArea) + (1 | gridcell / SUBSITE),
                    data = input.data.centred, family = FEve.distribution, chains = chains,
                    warmup = warmup, iter = iterations, cores = cores,
                    control = list(adapt_delta = delta,
                                   max_treedepth = treedepth))
    
    # Export model output
    save(FEve.mod, file = paste0("data/model_outputs_new/m_FEve_",
                                 "quadratic_", x.var, pc.filepath, ".RData"))
    
  } else {
    
    # Export model output
    FEve.mod <- get(load(paste0("data/model_outputs_new/m_FEve_",
                                "quadratic_", x.var, pc.filepath, ".RData")))
    
  } # End of import and censored
  
  # Convert model output into a dataframe (with 4.d.p.)
  FEve.df <- brms_SummaryTable(FEve.mod, formatOptions = list(digits = 5, nsmall = 5), round = 5)
  
  # Extract the confidence intervals as a list for use in the plotting
  FEve.df.ci <- FEve.df %>% 
    filter(Covariate %in% c("Icentred_x_variableE2")) %>% 
    dplyr::select("l-95% CI", "u-95% CI") %>% 
    pivot_longer(cols = 1:2, names_to = "Interval", values_to = "Value")
  
  # Save the confidence intervals as a list and remove the intermediate dataframe
  FEve.ci <- as.list(FEve.df.ci$Value) # Adding values to the list
  names(FEve.ci) <- FEve.df.ci$Interval # Adding names to the values
  rm(FEve.df, FEve.df.ci) # Remove unnecessary dataframes of summary and CIs  
  
  # Create template predictions dataframe
  FEve.pred.df = data.frame(centred_x_variable = seq(0 - mean_x_variable, 100 - mean_x_variable), SurveyedArea = 1)
  
  # Add predictions to template datframe
  FEve.pred = add_epred_draws(FEve.pred.df, FEve.mod, re_formula = NA) %>% 
    mutate(x_variable = centred_x_variable + mean_x_variable) # Uncentre the x variable
  
  # Plot model outputs
  (FEve.plot <- ggplot(data = FEve.pred, aes(x = x_variable, y = .epred)) +
      geom_point(data = input.data, aes(x = x_variable, y = FEve, fill = "#53C367"),
                 shape = 21) +
      stat_lineribbon(aes(y = .epred), .width = 0.95, fill = "#53C367",
                      alpha = 0.27, colour = "#000000", linetype = 2) +
      labs(y = "FEve \n",
           x = paste0("\n ", x.var),
           fill = paste0(colour.by, "\n"),
           title = paste0("CIs: ", FEve.ci$"l-95% CI", " to ", FEve.ci$"u-95% CI"),
           subtitle = paste(ifelse((FEve.ci$"l-95% CI" < 0 & FEve.ci$"u-95% CI" < 0) | # Automatically asigns significance in subtitle based on CIs
                                     (FEve.ci$"l-95% CI" > 0 & FEve.ci$"u-95% CI" > 0),
                                   "Significant: CIs do NOT span 0", "Not Significant: CIs span 0"))) +
      theme_1() +
      theme(legend.position = "right",
            plot.subtitle = element_text(colour = ifelse((FEve.ci$"l-95% CI" < 0 & FEve.ci$"u-95% CI" < 0) | # Automatically colours subtitle based on significance from CIs
                                                           (FEve.ci$"l-95% CI" > 0 & FEve.ci$"u-95% CI" > 0),
                                                         paste("#006400"), paste("#8b0000")))))
  
  # SR ----
  
  # Produce dataframe for running models on
  input.data <- data %>% 
    dplyr::select(x.var, FRic, FEve, SR, FDis, SurveyedArea, gridcell, SUBSITE) %>% 
    rename(x_variable = x.var) 
  
  # Determine mean value for centering
  mean_x_variable <- mean(input.data$x_variable)
  
  # Center the x_variable
  input.data.centred <- input.data %>% 
    mutate(centred_x_variable = x_variable - mean_x_variable) %>% 
    data.frame(.)
  
  
  if (run.SR == TRUE){
    
    # Function for running the Bayesian model
    SR.mod <- brm(SR ~ centred_x_variable + I(centred_x_variable^2) + log(SurveyedArea) + (1 | gridcell / SUBSITE),
                  data = input.data.centred, family = SR.distribution, chains = chains,
                  warmup = warmup, iter = iterations, cores = cores,
                  control = list(adapt_delta = delta,
                                 max_treedepth = treedepth))
    
    # Export model output
    save(SR.mod, file = paste0("data/model_outputs_new/m_SR_",
                               "quadratic_", x.var, pc.filepath, ".RData"))
    
  } else {
    
    # Export model output
    SR.mod <- get(load(paste0("data/model_outputs_new/m_SR_",
                              "quadratic_", x.var, pc.filepath, ".RData")))
    
  } # End of import and censored
  
  # Convert model output into a dataframe (with 4.d.p.)
  SR.df <- brms_SummaryTable(SR.mod, formatOptions = list(digits = 5, nsmall = 5), round = 5)
  
  # Extract the confidence intervals as a list for use in the plotting
  SR.df.ci <- SR.df %>% 
    filter(Covariate %in% c("Icentred_x_variableE2")) %>% 
    dplyr::select("l-95% CI", "u-95% CI") %>% 
    pivot_longer(cols = 1:2, names_to = "Interval", values_to = "Value")
  
  # Save the confidence intervals as a list and remove the intermediate dataframe
  SR.ci <- as.list(SR.df.ci$Value) # Adding values to the list
  names(SR.ci) <- SR.df.ci$Interval # Adding names to the values
  rm(SR.df, SR.df.ci) # Remove unnecessary dataframes of summary and CIs  
  
  # Create template predictions dataframe
  SR.pred.df = data.frame(centred_x_variable = seq(0 - mean_x_variable, 100 - mean_x_variable), SurveyedArea = 1)
  
  # Add predictions to template datframe
  SR.pred = add_epred_draws(SR.pred.df, SR.mod, re_formula = NA) %>% 
    mutate(x_variable = centred_x_variable + mean_x_variable) # Uncentre the x variable
  
  # Plot model outputs
  (SR.plot <- ggplot(data = SR.pred, aes(x = x_variable, y = .epred)) +
      geom_point(data = input.data, aes(x = x_variable, y = SR, fill = "#44B1A2")) +
      stat_lineribbon(aes(y = .epred), .width = 0.95, fill = "#44B1A2",
                      alpha = 0.27, colour = "#000000", linetype = 2) +
      labs(y = "SR \n",
           x = paste0("\n ", x.var),
           fill = paste0(colour.by, "\n"),
           title = paste0("CIs: ", SR.ci$"l-95% CI", " to ", SR.ci$"u-95% CI"),
           subtitle = paste(ifelse((SR.ci$"l-95% CI" < 0 & SR.ci$"u-95% CI" < 0) | # Automatically asigns significance in subtitle based on CIs
                                     (SR.ci$"l-95% CI" > 0 & SR.ci$"u-95% CI" > 0),
                                   "Significant: CIs do NOT span 0", "Not Significant: CIs span 0"))) +
      theme_1() +
      theme(legend.position = "right",
            plot.subtitle = element_text(colour = ifelse((SR.ci$"l-95% CI" < 0 & SR.ci$"u-95% CI" < 0) | # Automatically colours subtitle based on significance from CIs
                                                           (SR.ci$"l-95% CI" > 0 & SR.ci$"u-95% CI" > 0),
                                                         paste("#006400"), paste("#8b0000")))))
  
  # FDis ----
  
  # Produce dataframe for running models on
  input.data <- data %>% 
    dplyr::select(x.var, FRic, FEve, SR, FDis, SurveyedArea, gridcell, SUBSITE) %>% 
    rename(x_variable = x.var) 
  
  # Determine mean value for centering
  mean_x_variable <- mean(input.data$x_variable)
  
  # Center the x_variable
  input.data.centred <- input.data %>% 
    mutate(centred_x_variable = x_variable - mean_x_variable) %>% 
    data.frame(.)
  
  
  if (run.FDis == TRUE){
    
    # Function for running the Bayesian model
    FDis.mod <- brm(FDis ~ centred_x_variable + I(centred_x_variable^2) + log(SurveyedArea) + (1 | gridcell / SUBSITE),
                    data = input.data.centred, family = FDis.distribution, chains = chains,
                    warmup = warmup, iter = iterations, cores = cores,
                    control = list(adapt_delta = delta,
                                   max_treedepth = treedepth))
    
    # Export model output
    save(FDis.mod, file = paste0("data/model_outputs_new/m_FDis_",
                                 "quadratic_", x.var, pc.filepath, ".RData"))
    
  } else {
    
    # Export model output
    FDis.mod <- get(load(paste0("data/model_outputs_new/m_FDis_",
                                "quadratic_", x.var, pc.filepath, ".RData")))
    
  } # End of import and censored
  
  # Convert model output into a dataframe (with 4.d.p.)
  FDis.df <- brms_SummaryTable(FDis.mod, formatOptions = list(digits = 5, nsmall = 5), round = 5)
  
  # Extract the confidence intervals as a list for use in the plotting
  FDis.df.ci <- FDis.df %>% 
    filter(Covariate %in% c("Icentred_x_variableE2")) %>% 
    dplyr::select("l-95% CI", "u-95% CI") %>% 
    pivot_longer(cols = 1:2, names_to = "Interval", values_to = "Value")
  
  # Save the confidence intervals as a list and remove the intermediate dataframe
  FDis.ci <- as.list(FDis.df.ci$Value) # Adding values to the list
  names(FDis.ci) <- FDis.df.ci$Interval # Adding names to the values
  rm(FDis.df, FDis.df.ci) # Remove unnecessary dataframes of summary and CIs  
  
  # Create template predictions dataframe
  FDis.pred.df = data.frame(centred_x_variable = seq(0 - mean_x_variable, 100 - mean_x_variable), SurveyedArea = 1)
  
  # Add predictions to template datframe
  FDis.pred = add_epred_draws(FDis.pred.df, FDis.mod, re_formula = NA) %>% 
    mutate(x_variable = centred_x_variable + mean_x_variable) # Uncentre the x variable
  
  # Plot model outputs
  (FDis.plot <- ggplot(data = FDis.pred, aes(x = x_variable, y = .epred)) +
      geom_point(data = input.data, aes(x = x_variable, y = FDis, fill = "#FF7F00")) +
      stat_lineribbon(aes(y = .epred), .width = 0.95, fill = "#FF7F00",
                      alpha = 0.27, colour = "#000000", linetype = 2) +
      # scale_fill_viridis(option = FDis.colour, begin = 0.2, end = 1, direction = -1, discrete = colour.discrete,
      #                    guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")) +
      labs(y = "FDis \n",
           x = paste0("\n ", x.var),
           fill = paste0(colour.by, "\n"),
           title = paste0("CIs: ", FDis.ci$"l-95% CI", " to ", FDis.ci$"u-95% CI"),
           subtitle = paste(ifelse((FDis.ci$"l-95% CI" < 0 & FDis.ci$"u-95% CI" < 0) | # Automatically asigns significance in subtitle based on CIs
                                     (FDis.ci$"l-95% CI" > 0 & FDis.ci$"u-95% CI" > 0),
                                   "Significant: CIs do NOT span 0", "Not Significant: CIs span 0"))) +
      theme_1() +
      theme(legend.position = "right",
            plot.subtitle = element_text(colour = ifelse((FDis.ci$"l-95% CI" < 0 & FDis.ci$"u-95% CI" < 0) | # Automatically colours subtitle based on significance from CIs
                                                           (FDis.ci$"l-95% CI" > 0 & FDis.ci$"u-95% CI" > 0),
                                                         paste("#006400"), paste("#8b0000")))))
  
  # Panel ----
  
  # Create a panel of all four plots
  combined.panel <- grid.arrange(FRic.plot, FEve.plot, SR.plot, FDis.plot, ncol = 2)
  
  # Add in if statements for saving the output
  if (censored == "No") {
    
    # Export panel
    ggsave(combined.panel, filename = paste0("figures/outputs_new/combined_",
                                             x.var, "_ALL_quadratic", pc.filepath,".png"), width = 15, height = 14)
    
  } else {
    
    # Export panel
    ggsave(combined.panel, filename = paste0("figures/outputs_new/combined_",
                                             x.var, "_censored_ALL_quadratic", pc.filepath,".png"), width = 15, height = 14)
    
  }
  
} # End of function


# FUNCTION: SPATIAL BAYESIAN MODELS FOR FD WITH A CATEGORIC X-VARIABLE ----

bayesian.spatial.categoric <- function(run.FRic, censored, run.FEve, run.SR, run.FDis, x.var, data,
                                       chains, cores, warmup, iterations, delta, treedepth){

  # FRic ----
  
  # Produce dataframe for running models on
  input.data <- data %>% 
    dplyr::select(x.var, FRic, FEve, SR, FDis, SurveyedArea, gridcell, SUBSITE) %>% 
    rename(x_variable = x.var) # Rename first column to x_variable
  
  # Tidy import data for running censored models
  input.data.2 <- input.data %>% 
    mutate(RETAIN = ifelse(SR >= 4 & is.na(FRic), FALSE, TRUE)) %>% # Remove erroneous NA value for SR == 4 
    filter(RETAIN == TRUE) %>% 
    dplyr::select(-RETAIN) %>% 
    mutate(log_FRic = log(FRic)) # For gaussian distribution
  
  # Determine the median value of FRic for when SR < 4 (i.e. FRic is NA)
  median.FRic.4 <- median(filter(input.data.2, SR == 4)$log_FRic)
  
  # Create columns for modelling
  input.data.3 <- input.data.2 %>% 
    mutate(FRic2 = ifelse(SR >= 1 & SR < 4, 0.000000000000000001, log_FRic), # Create interval censored columns
           FRic3 = ifelse(SR >= 1 & SR < 4, median.FRic.4, log_FRic),
           cen2 = ifelse(SR >= 1 & SR < 4, "interval", "none")) %>% 
    mutate(FRic_na = ifelse(cen2 == "none", log_FRic, NA))
  
  # Determine which version of FRic to run/import (censored or not)
  if (run.FRic == TRUE & censored == "No"){
    
    # Function for running the Bayesian model
    FRic.mod <- brm(FRic ~ x_variable + log(SurveyedArea) + (1 | gridcell / SUBSITE),
                    data = input.data, family = FRic.distribution, chains = chains,
                    warmup = warmup, iter = iterations, cores = cores,
                    control = list(adapt_delta = delta,
                                   max_treedepth = treedepth))
    
    # Export model output
    save(FRic.mod, file = paste0("data/model_outputs_new/m_FRic_",
                                 FRic.distribution, "_", x.var, pc.filepath, ".RData"))
    
  } # End of run and not censored
  
  if (run.FRic == FALSE & censored == "No"){
    
    # Load in model output
    FRic.mod <- get(load(paste0("data/model_outputs_new/m_FRic_",
                                FRic.distribution, "_", x.var, pc.filepath, ".RData")))
    
  } # End of import and not censored
  
  if (run.FRic == TRUE & censored == "Yes"){
    
    # Define the stanvars
    mean_FRic <- mean(input.data.3$FRic_na, na.rm = TRUE)
    sd_FRic <- sd(input.data.3$FRic_na, na.rm = TRUE)
    stanvars <- stanvar(mean_FRic, name = "mean_FRic") + stanvar(sd_FRic,   name = "sd_FRic")
    
    # Run the censored brms() model
    FRic.mod <- brm(FRic2 | cens(cen2, FRic3) ~ x_variable + log(SurveyedArea) + (1 | gridcell / SUBSITE),
                    data = input.data.3, family = "gaussian",
                    chains = chains, warmup = warmup, iter = iterations,
                    cores = cores, control = list(adapt_delta = delta,
                                                  max_treedepth = treedepth))
    
    # Export model output
    save(FRic.mod, file = paste0("data/model_outputs_new/m_FRic_",
                                 "censored_", x.var, pc.filepath, ".RData"))
    
  } # End of run and censored
  
  if (run.FRic == FALSE & censored == "Yes"){
    
    # Load in model output
    FRic.mod <- get(load(paste0("data/model_outputs_new/m_FRic_",
                                "censored_", x.var, pc.filepath, ".RData")))
    
  } # End of import and censored
  
  # Process predictions differently dependent on whether censored or not
  if (censored == "No"){
    
    # Extract the prediction data frame
    FRic.pred <- ggpredict(FRic.mod, terms = "x_variable",
                           # type = "random", condition = c(gridcell = 1266), allow_new_levels = TRUE,
                           back.transform = TRUE)
    
  } else { # End of censored == "No"
    
    # Extract the prediction data frame and exponentiate the outputs
    FRic.pred <- ggpredict(FRic.mod, terms = "x_variable",
                           # type = "random", condition = c(gridcell = 1266), allow_new_levels = TRUE,
                           back.transform = FALSE) %>% 
      mutate(predicted = exp(predicted),
             conf.low = exp(conf.low),
             conf.high = exp(conf.high))
    
  } # End of censored == "Yes
  
  # Plot the model outputs
  (FRic.plot <- ggplot(FRic.pred) +
      geom_errorbar(aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high), width = 0.2) +
      geom_point(aes(x = x, y = predicted, fill = x), stat = "identity", colour = c("#000000"), alpha = 0.7, shape = 21, size = 5) +
      scale_fill_viridis(option = FRic.colour, begin = 0.2, end = 1, direction = -1, discrete = TRUE) +
      labs(y = "FRic \n",
           x = "\n Region") +
      coord_flip() +
      theme_1() +
      theme(legend.position = "none",
            plot.title = element_text(size = 12),
            plot.subtitle = element_text(size = 6)))
  
  # FEve ----
  
  # If you want to run model for first time
  if (run.FEve == TRUE){
    
    # Function for running the Bayesian model
    FEve.mod <- brm(FEve ~ x_variable + log(SurveyedArea) + (1 | gridcell / SUBSITE),
                    data = input.data, family = FEve.distribution, chains = chains,
                    warmup = warmup, iter = iterations, cores = cores,
                    control = list(adapt_delta = delta,
                                   max_treedepth = treedepth))
    
    # Export model output
    save(FEve.mod, file = paste0("data/model_outputs_new/m_FEve_",
                                 x.var, pc.filepath, ".RData"))
    
  } else {
    
    # Load in model output
    FEve.mod <- get(load(paste0("data/model_outputs_new/m_FEve_",
                                x.var, pc.filepath, ".RData")))
    
  }
  
  
  # Extract the prediction data frame
  FEve.pred <- ggpredict(FEve.mod, terms = "x_variable", back.transform = FALSE)

  # Plot the model outputs
  (FEve.plot <- ggplot(FEve.pred) +
      geom_errorbar(aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high), width = 0.2) +
      geom_point(aes(x = x, y = predicted, fill = x), stat = "identity", colour = c("#000000"), alpha = 0.7, shape = 21, size = 5) +
      scale_fill_viridis(option = FEve.colour, begin = 0.2, end = 1, direction = -1, discrete = TRUE) +
      labs(y = "FEve \n",
           x = "\n Region") +
      coord_flip() +
      theme_1() +
      theme(legend.position = "none",
            plot.title = element_text(size = 12),
            plot.subtitle = element_text(size = 6)))
  
  # SR ----
  
  # If you want to run model for first time
  if (run.SR == TRUE){
    
    # Function for running the Bayesian model
    SR.mod <- brm(SR ~ x_variable + log(SurveyedArea) + (1 | gridcell / SUBSITE),
                  data = input.data, family = SR.distribution, chains = chains,
                  warmup = warmup, iter = iterations, cores = cores,
                  control = list(adapt_delta = delta,
                                 max_treedepth = treedepth))
    
    # Export model output
    save(SR.mod, file = paste0("data/model_outputs_new/m_SR_",
                               x.var, pc.filepath, ".RData"))
    
  } else {
    
    # Load in model output
    SR.mod <- get(load(paste0("data/model_outputs_new/m_SR_",
                              x.var, pc.filepath, ".RData")))
    
  }
  
  
  # Extract the prediction data frame
  SR.pred <- ggpredict(SR.mod, terms = "x_variable", back.transform = FALSE)

  # Plot the model outputs
  (SR.plot <- ggplot(SR.pred) +
      geom_errorbar(aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high), width = 0.2) +
      geom_point(aes(x = x, y = predicted, fill = x), stat = "identity", colour = c("#000000"), alpha = 0.7, shape = 21, size = 5) +
      scale_fill_viridis(option = SR.colour, begin = 0.2, end = 1, direction = -1, discrete = TRUE) +
      labs(y = "SR \n",
           x = "\n Region") +
      coord_flip() +
      theme_1() +
      theme(legend.position = "none",
            plot.title = element_text(size = 12),
            plot.subtitle = element_text(size = 6)))
  
  # FDis ----
  
  # If you want to run model for first time
  if (run.FDis == TRUE){
    
    # Function for running the Bayesian model
    FDis.mod <- brm(FDis ~ x_variable + log(SurveyedArea) + (1 | gridcell / SUBSITE),
                    data = input.data, family = FDis.distribution, chains = chains,
                    warmup = warmup, iter = iterations, cores = cores,
                    control = list(adapt_delta = delta,
                                   max_treedepth = treedepth))
    
    # Export model output
    save(FDis.mod, file = paste0("data/model_outputs_new/m_FDis_",
                                 x.var, pc.filepath, ".RData"))
    
  } else {
    
    # Load in model output
    FDis.mod <- get(load(paste0("data/model_outputs_new/m_FDis_",
                                x.var, pc.filepath, ".RData")))
    
  }
  
  
  # Extract the prediction data frame
  FDis.pred <- ggpredict(FDis.mod, terms = "x_variable", back.transform = FALSE)
  
  # Plot the model outputs
  (FDis.plot <- ggplot(FDis.pred) +
      geom_errorbar(aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high), width = 0.2) +
      geom_point(aes(x = x, y = predicted, fill = x), stat = "identity", colour = c("#000000"), alpha = 0.7, shape = 21, size = 5) +
      scale_fill_viridis(option = FDis.colour, begin = 0.2, end = 1, direction = -1, discrete = TRUE) +
      labs(y = "FDis \n",
           x = "\n Region") +
      coord_flip() +
      theme_1() +
      theme(legend.position = "none",
            plot.title = element_text(size = 12),
            plot.subtitle = element_text(size = 6)))
  
  
  # Create a panel of all three plots
  combined.panel <- grid.arrange(FRic.plot, FEve.plot, SR.plot, FDis.plot, ncol = 2)
  
  # Export panel based on whether censored or not
  if (censored == "No"){
    
    # Export panel
    ggsave(combined.panel, filename = paste0("figures/outputs_new/combined_",
                                             x.var, pc.filepath, ".png"), width = 15, height = 14)
    
  } # End of censored == "No"
  
  if (censored == "Yes"){
    
    # Export panel
    ggsave(combined.panel, filename = paste0("figures/outputs_new/combined_",
                                             x.var, "_censored", pc.filepath, ".png"), width = 15, height = 14)
    
  } # End of censored == "Yes"
    
}


# FUNCTION: GENERATE LINEAR SLOPES PER PLOT ----

linear.slopes <- function(x.var, data){

  # Generate dataframe of slopes
  slopes <- data %>% 
    nest_by(SiteSubsitePlot) %>%
    mutate(mod = list(lm(x_variable ~ YEAR, data = data))) %>% 
    summarise(tidy(mod)) %>% 
    ungroup() %>% 
    filter(term == "YEAR") %>% 
    dplyr::select(SiteSubsitePlot, estimate) %>%
    mutate(Slope_Type = paste0(x.var, "_slopes")) %>% 
    rename(Slope = estimate)
  
}


# FUNCTION: TEMPORAL BAYESIAN MODELS FOR FD & COVER CHANGE OVER TIME ----

bayesian.temporal.change <- function(run.FRic, run.FEve, run.SR, run.FDis, run.ShrubCover,
                                     run.GraminoidCover, run.ForbCover, data, chains,
                                     cores, warmup, iterations, delta, treedepth){

  # FRic ----
  
  # If you want to run model for first time
  if (run.FRic == TRUE){
    
    # Function for running the Bayesian model
    FRic.mod <- brm(FRic_slopes ~ 1 + log(SurveyedArea) + (1 | gridcell / SUBSITE), data = data,
                    family = FRic_slopes.distribution, chains = chains,
                    warmup = warmup, iter = iterations, cores = cores,
                    control = list(adapt_delta = delta,
                                   max_treedepth = treedepth))
    
    # Export model output
    save(FRic.mod, file = paste0("data/model_outputs_new/m_FRic_slopes",
                                 pc.filepath, ".RData"))
    
  } else {
    
    # Load in model output
    FRic.mod <- get(load(paste0("data/model_outputs_new/m_FRic_slopes",
                                pc.filepath, ".RData")))
    
  }
  
  # Save model output as a dataframe
  FRic.df <- brms_SummaryTable(FRic.mod, formatOptions = list(digits = 4, nsmall = 4), round = 4)
  
  # Create intermediate dataframe to calculate mean FRic slope
  FRic_slopes.df <- data %>% 
    filter(!is.na(FRic_slopes))
  
  # Calculate mean FRic slope value
  FRic_slopes.mean <- mean(FRic_slopes.df$FRic_slopes)
  
  # Determine bin width for plotting
  FRic_slopes.bin.width <- (max(filter(data, !is.na(FRic_slopes))$FRic_slopes) - min(filter(data, !is.na(FRic_slopes))$FRic_slopes))/40
  
  # Plot model outputs
  (FRic.plot <- ggplot(data) +
      geom_histogram(aes(x = FRic_slopes), stat = "bin", binwidth = FRic_slopes.bin.width,
                     fill = "gray", alpha = 1) +
      geom_density(aes(x = FRic_slopes), stat = "bin", binwidth = FRic_slopes.bin.width,
                   fill = FRic_slopes.colour, colour = FRic_slopes.colour, alpha = 0.4, size = 0.5) +
      geom_vline(aes(xintercept = FRic_slopes.mean), colour = "#000000", linetype = "dashed", size = 1) +
      labs(y = "Number of Plots \n",
           x = "\n Change in FRic (per year)",
           title = paste0("Mean FRic Change: ", FRic.df$Estimate, " per year", "\n",
                          "95% CIs: [", FRic.df$`l-95% CI`, ", ", FRic.df$`u-95% CI`, "]"),
           subtitle = paste0(ifelse((FRic.df$"l-95% CI" < 0 & FRic.df$"u-95% CI" < 0) | # Automatically asigns significance in subtitle based on CIs
                                      (FRic.df$"l-95% CI" > 0 & FRic.df$"u-95% CI" > 0),
                                    "Significant: CIs do NOT span 0", "Not Significant: CIs span 0"))) +
      theme_1() +
      theme(plot.subtitle = element_text(colour = ifelse((FRic.df$"l-95% CI" < 0 & FRic.df$"u-95% CI" < 0) | # Automatically colours subtitle based on significance from CIs
                                                           (FRic.df$"l-95% CI" > 0 & FRic.df$"u-95% CI" > 0),
                                                         paste("#006400"), paste("#8b0000"))))
      
      )
  
  # FEve ----
  
  # If you want to run model for first time
  if (run.FEve == TRUE){
    
    # Function for running the Bayesian model
    FEve.mod <- brm(FEve_slopes ~ 1 + log(SurveyedArea) + (1 | gridcell / SUBSITE), data = data,
                    family = FEve_slopes.distribution, chains = chains,
                    warmup = warmup, iter = iterations, cores = cores,
                    control = list(adapt_delta = delta,
                                   max_treedepth = treedepth))
    
    # Export model output
    save(FEve.mod, file = paste0("data/model_outputs_new/m_FEve_slopes",
                                 pc.filepath, ".RData"))
    
  } else {
    
    # Load in model output
    FEve.mod <- get(load(paste0("data/model_outputs_new/m_FEve_slopes",
                                pc.filepath, ".RData")))
    
  }
  
  # Save model output as a dataframe
  FEve.df <- brms_SummaryTable(FEve.mod, formatOptions = list(digits = 4, nsmall = 4), round = 4)
  
  # Create intermediate dataframe to calculate mean FEve slope
  FEve_slopes.df <- data %>% 
    filter(!is.na(FEve_slopes))
  
  # Calculate mean FEve slope value
  FEve_slopes.mean <- mean(FEve_slopes.df$FEve_slopes)
  
  # Determine bin width for plotting
  FEve_slopes.bin.width <- (max(filter(data, !is.na(FEve_slopes))$FEve_slopes) - min(filter(data, !is.na(FEve_slopes))$FEve_slopes))/40
  
  # Plot model outputs
  (FEve.plot <- ggplot(data) +
      geom_histogram(aes(x = FEve_slopes), stat = "bin", binwidth = FEve_slopes.bin.width,
                     fill = "gray", alpha = 1) +
      geom_density(aes(x = FEve_slopes), stat = "bin", binwidth = FEve_slopes.bin.width,
                   fill = FEve_slopes.colour, colour = FEve_slopes.colour, alpha = 0.4, size = 0.5) +
      geom_vline(aes(xintercept = FEve_slopes.mean), colour = "#000000", linetype = "dashed", size = 1) +
      labs(y = "Number of Plots \n",
           x = "\n Change in FEve (per year)",
           title = paste0("Mean FEve Change: ", FEve.df$Estimate, " per year", "\n",
                          "95% CIs: [", FEve.df$`l-95% CI`, ", ", FEve.df$`u-95% CI`, "]"),
           subtitle = paste0(ifelse((FEve.df$"l-95% CI" < 0 & FEve.df$"u-95% CI" < 0) | # Automatically asigns significance in subtitle based on CIs
                                      (FEve.df$"l-95% CI" > 0 & FEve.df$"u-95% CI" > 0),
                                    "Significant: CIs do NOT span 0", "Not Significant: CIs span 0"))) +
      theme_1() +
      theme(plot.subtitle = element_text(colour = ifelse((FEve.df$"l-95% CI" < 0 & FEve.df$"u-95% CI" < 0) | # Automatically colours subtitle based on significance from CIs
                                                           (FEve.df$"l-95% CI" > 0 & FEve.df$"u-95% CI" > 0),
                                                         paste("#006400"), paste("#8b0000"))))
    
  )
  
  # SR ----
  
  # If you want to run model for first time
  if (run.SR == TRUE){
    
    # Function for running the Bayesian model
    SR.mod <- brm(SR_slopes ~ 1 + log(SurveyedArea) + (1 | gridcell / SUBSITE), data = data,
                  family = SR_slopes.distribution, chains = chains,
                  warmup = warmup, iter = iterations, cores = cores,
                  control = list(adapt_delta = delta,
                                 max_treedepth = treedepth))
    
    # Export model output
    save(SR.mod, file = paste0("data/model_outputs_new/m_SR_slopes",
                               pc.filepath, ".RData"))
    
  } else {
    
    # Load in model output
    SR.mod <- get(load(paste0("data/model_outputs_new/m_SR_slopes",
                              pc.filepath, ".RData")))
    
  }
  
  # Save model output as a dataframe
  SR.df <- brms_SummaryTable(SR.mod, formatOptions = list(digits = 4, nsmall = 4), round = 4)
  
  # Create intermediate dataframe to calculate mean SR slope
  SR_slopes.df <- data %>% 
    filter(!is.na(SR_slopes))
  
  # Calculate mean SR slope value
  SR_slopes.mean <- mean(SR_slopes.df$SR_slopes)
  
  # Determine bin width for plotting
  SR_slopes.bin.width <- (max(filter(data, !is.na(SR_slopes))$SR_slopes) - min(filter(data, !is.na(SR_slopes))$SR_slopes))/40
  
  # Plot model outputs
  (SR.plot <- ggplot(data) +
      geom_histogram(aes(x = SR_slopes), stat = "bin", binwidth = SR_slopes.bin.width,
                     fill = "gray", alpha = 1) +
      geom_density(aes(x = SR_slopes), stat = "bin", binwidth = SR_slopes.bin.width,
                   fill = SR_slopes.colour, colour = SR_slopes.colour, alpha = 0.4, size = 0.5) +
      geom_vline(aes(xintercept = SR_slopes.mean), colour = "#000000", linetype = "dashed", size = 1) +
      labs(y = "Number of Plots \n",
           x = "\n Change in SR (per year)",
           title = paste0("Mean SR Change: ", SR.df$Estimate, " per year", "\n",
                          "95% CIs: [", SR.df$`l-95% CI`, ", ", SR.df$`u-95% CI`, "]"),
           subtitle = paste0(ifelse((SR.df$"l-95% CI" < 0 & SR.df$"u-95% CI" < 0) | # Automatically asigns significance in subtitle based on CIs
                                      (SR.df$"l-95% CI" > 0 & SR.df$"u-95% CI" > 0),
                                    "Significant: CIs do NOT span 0", "Not Significant: CIs span 0"))) +
      theme_1() +
      theme(plot.subtitle = element_text(colour = ifelse((SR.df$"l-95% CI" < 0 & SR.df$"u-95% CI" < 0) | # Automatically colours subtitle based on significance from CIs
                                                           (SR.df$"l-95% CI" > 0 & SR.df$"u-95% CI" > 0),
                                                         paste("#006400"), paste("#8b0000"))))
    
  )
  
  # FDis ----
  
  # If you want to run model for first time
  if (run.FDis == TRUE){
    
    # Function for running the Bayesian model
    FDis.mod <- brm(FDis_slopes ~ 1 + log(SurveyedArea) + (1 | gridcell / SUBSITE), data = data,
                    family = FDis_slopes.distribution, chains = chains,
                    warmup = warmup, iter = iterations, cores = cores,
                    control = list(adapt_delta = delta,
                                   max_treedepth = treedepth))
    
    # Export model output
    save(FDis.mod, file = paste0("data/model_outputs_new/m_FDis_slopes",
                                 pc.filepath, ".RData"))
    
  } else {
    
    # Load in model output
    FDis.mod <- get(load(paste0("data/model_outputs_new/m_FDis_slopes",
                                pc.filepath, ".RData")))
    
  }
  
  # Save model output as a dataframe
  FDis.df <- brms_SummaryTable(FDis.mod, formatOptions = list(digits = 4, nsmall = 4), round = 4)
  
  # Create intermediate dataframe to calculate mean SR slope
  FDis_slopes.df <- data %>% 
    filter(!is.na(FDis_slopes))
  
  # Calculate mean SR slope value
  FDis_slopes.mean <- mean(FDis_slopes.df$FDis_slopes)
  
  # Determine bin width for plotting
  FDis_slopes.bin.width <- (max(filter(data, !is.na(FDis_slopes))$FDis_slopes) - min(filter(data, !is.na(FDis_slopes))$FDis_slopes))/40
  
  # Plot model outputs
  (FDis.plot <- ggplot(data) +
      geom_histogram(aes(x = FDis_slopes), stat = "bin", binwidth = FDis_slopes.bin.width,
                     fill = "gray", alpha = 1) +
      geom_density(aes(x = FDis_slopes), stat = "bin", binwidth = FDis_slopes.bin.width,
                   fill = FDis_slopes.colour, colour = FDis_slopes.colour, alpha = 0.4, size = 0.5) +
      geom_vline(aes(xintercept = FDis_slopes.mean), colour = "#000000", linetype = "dashed", size = 1) +
      labs(y = "Number of Plots \n",
           x = "\n Change in FDis (per year)",
           title = paste0("Mean FDis Change: ", FDis.df$Estimate, " per year", "\n",
                          "95% CIs: [", FDis.df$`l-95% CI`, ", ", FDis.df$`u-95% CI`, "]"),
           subtitle = paste0(ifelse((FDis.df$"l-95% CI" < 0 & FDis.df$"u-95% CI" < 0) | # Automatically asigns significance in subtitle based on CIs
                                      (FDis.df$"l-95% CI" > 0 & FDis.df$"u-95% CI" > 0),
                                    "Significant: CIs do NOT span 0", "Not Significant: CIs span 0"))) +
      theme_1() +
      theme(plot.subtitle = element_text(colour = ifelse((FDis.df$"l-95% CI" < 0 & FDis.df$"u-95% CI" < 0) | # Automatically colours subtitle based on significance from CIs
                                                           (FDis.df$"l-95% CI" > 0 & FDis.df$"u-95% CI" > 0),
                                                         paste("#006400"), paste("#8b0000"))))
    
  )
  
  # Create a panel of all three plots
  combined.panel.FD <- grid.arrange(FRic.plot, FEve.plot, SR.plot, FDis.plot, ncol = 2)
  
  # Export panel
  ggsave(combined.panel.FD, filename = paste0("figures/outputs_new/combined_FD_change",
                                              pc.filepath, ".png"), width = 15, height = 14)
  
  # Shrubs ----
  
  # If you want to run model for first time
  if (run.ShrubCover == TRUE){
    
    # Function for running the Bayesian model
    ShrubCover.mod <- brm(ShrubCover_slopes ~ 1 + log(SurveyedArea) + (1 | gridcell / SUBSITE), data = data,
                          family = ShrubCover_slopes.distribution, chains = chains,
                          warmup = warmup, iter = iterations, cores = cores,
                          control = list(adapt_delta = delta,
                                         max_treedepth = treedepth))
    
    # Export model output
    save(ShrubCover.mod, file = paste0("data/model_outputs_new/m_ShrubCover_slopes",
                                       pc.filepath, ".RData"))
    
  } else {
    
    # Load in model output
    ShrubCover.mod <- get(load(paste0("data/model_outputs_new/m_ShrubCover_slopes",
                                      pc.filepath, ".RData")))
    
  }
  
  # Save model output as a dataframe
  ShrubCover.df <- brms_SummaryTable(ShrubCover.mod, formatOptions = list(digits = 4, nsmall = 4), round = 4)
  
  # Create intermediate dataframe to calculate mean ShrubCover slope
  ShrubCover_slopes.df <- data %>% 
    filter(!is.na(ShrubCover_slopes))
  
  # Calculate mean ShrubCover slope value
  ShrubCover_slopes.mean <- mean(ShrubCover_slopes.df$ShrubCover_slopes)
  
  # Determine bin width for plotting
  ShrubCover_slopes.bin.width <- (max(filter(data, !is.na(ShrubCover_slopes))$ShrubCover_slopes) - min(filter(data, !is.na(ShrubCover_slopes))$ShrubCover_slopes))/40
  
  # Plot model outputs
  (ShrubCover.plot <- ggplot(data) +
      geom_histogram(aes(x = ShrubCover_slopes), stat = "bin", binwidth = ShrubCover_slopes.bin.width,
                     fill = "gray", alpha = 1) +
      geom_density(aes(x = ShrubCover_slopes), stat = "bin", binwidth = ShrubCover_slopes.bin.width,
                   fill = ShrubCover_slopes.colour, colour = ShrubCover_slopes.colour, alpha = 0.4, size = 0.5) +
      geom_vline(aes(xintercept = ShrubCover_slopes.mean), colour = "#000000", linetype = "dashed", size = 1) +
      labs(y = "Number of Plots \n",
           x = "\n Change in ShrubCover (per year)",
           title = paste0("Mean ShrubCover Change: ", ShrubCover.df$Estimate, " per year", "\n",
                          "95% CIs: [", ShrubCover.df$`l-95% CI`, ", ", ShrubCover.df$`u-95% CI`, "]"),
           subtitle = paste0(ifelse((ShrubCover.df$"l-95% CI" < 0 & ShrubCover.df$"u-95% CI" < 0) | # Automatically asigns significance in subtitle based on CIs
                                      (ShrubCover.df$"l-95% CI" > 0 & ShrubCover.df$"u-95% CI" > 0),
                                    "Significant: CIs do NOT span 0", "Not Significant: CIs span 0"))) +
      theme_1() +
      theme(plot.subtitle = element_text(colour = ifelse((ShrubCover.df$"l-95% CI" < 0 & ShrubCover.df$"u-95% CI" < 0) | # Automatically colours subtitle based on significance from CIs
                                                           (ShrubCover.df$"l-95% CI" > 0 & ShrubCover.df$"u-95% CI" > 0),
                                                         paste("#006400"), paste("#8b0000"))))
    
  )
  
  # Graminoids ----
  
  # If you want to run model for first time
  if (run.GraminoidCover == TRUE){
    
    # Function for running the Bayesian model
    GraminoidCover.mod <- brm(GraminoidCover_slopes ~ 1 + log(SurveyedArea) + (1 | gridcell / SUBSITE), data = data,
                              family = GraminoidCover_slopes.distribution, chains = chains,
                              warmup = warmup, iter = iterations, cores = cores,
                              control = list(adapt_delta = delta,
                                             max_treedepth = treedepth))
    
    # Export model output
    save(GraminoidCover.mod, file = paste0("data/model_outputs_new/m_GraminoidCover_slopes",
                                           pc.filepath, ".RData"))
    
  } else {
    
    # Load in model output
    GraminoidCover.mod <- get(load(paste0("data/model_outputs_new/m_GraminoidCover_slopes",
                                          pc.filepath, ".RData")))
    
  }
  
  # Save model output as a dataframe
  GraminoidCover.df <- brms_SummaryTable(GraminoidCover.mod, formatOptions = list(digits = 4, nsmall = 4), round = 4)
  
  # Create intermediate dataframe to calculate mean GraminoidCover slope
  GraminoidCover_slopes.df <- data %>% 
    filter(!is.na(GraminoidCover_slopes))
  
  # Calculate mean GraminoidCover slope value
  GraminoidCover_slopes.mean <- mean(GraminoidCover_slopes.df$GraminoidCover_slopes)
  
  # Determine bin width for plotting
  GraminoidCover_slopes.bin.width <- (max(filter(data, !is.na(GraminoidCover_slopes))$GraminoidCover_slopes) - min(filter(data, !is.na(GraminoidCover_slopes))$GraminoidCover_slopes))/40
  
  # Plot model outputs
  (GraminoidCover.plot <- ggplot(data) +
      geom_histogram(aes(x = GraminoidCover_slopes), stat = "bin", binwidth = GraminoidCover_slopes.bin.width,
                     fill = "gray", alpha = 1) +
      geom_density(aes(x = GraminoidCover_slopes), stat = "bin", binwidth = GraminoidCover_slopes.bin.width,
                   fill = GraminoidCover_slopes.colour, colour = GraminoidCover_slopes.colour, alpha = 0.4, size = 0.5) +
      geom_vline(aes(xintercept = GraminoidCover_slopes.mean), colour = "#000000", linetype = "dashed", size = 1) +
      labs(y = "Number of Plots \n",
           x = "\n Change in GraminoidCover (per year)",
           title = paste0("Mean GraminoidCover Change: ", GraminoidCover.df$Estimate, " per year", "\n",
                          "95% CIs: [", GraminoidCover.df$`l-95% CI`, ", ", GraminoidCover.df$`u-95% CI`, "]"),
           subtitle = paste0(ifelse((GraminoidCover.df$"l-95% CI" < 0 & GraminoidCover.df$"u-95% CI" < 0) | # Automatically asigns significance in subtitle based on CIs
                                      (GraminoidCover.df$"l-95% CI" > 0 & GraminoidCover.df$"u-95% CI" > 0),
                                    "Significant: CIs do NOT span 0", "Not Significant: CIs span 0"))) +
      theme_1() +
      theme(plot.subtitle = element_text(colour = ifelse((GraminoidCover.df$"l-95% CI" < 0 & GraminoidCover.df$"u-95% CI" < 0) | # Automatically colours subtitle based on significance from CIs
                                                           (GraminoidCover.df$"l-95% CI" > 0 & GraminoidCover.df$"u-95% CI" > 0),
                                                         paste("#006400"), paste("#8b0000"))))
    
  )
  
  # Forbs ----
  
  # If you want to run model for first time
  if (run.ForbCover == TRUE){
    
    # Function for running the Bayesian model
    ForbCover.mod <- brm(ForbCover_slopes ~ 1 + log(SurveyedArea) + (1 | gridcell / SUBSITE), data = data,
                         family = ForbCover_slopes.distribution, chains = chains,
                         warmup = warmup, iter = iterations, cores = cores,
                         control = list(adapt_delta = delta,
                                        max_treedepth = treedepth))
    
    # Export model output
    save(ForbCover.mod, file = paste0("data/model_outputs_new/m_ForbCover_slopes",
                                      pc.filepath, ".RData"))
    
  } else {
    
    # Load in model output
    ForbCover.mod <- get(load(paste0("data/model_outputs_new/m_ForbCover_slopes",
                                     pc.filepath, ".RData")))
    
  }
  
  # Save model output as a dataframe
  ForbCover.df <- brms_SummaryTable(ForbCover.mod, formatOptions = list(digits = 4, nsmall = 4), round = 4)
  
  # Create intermediate dataframe to calculate mean ForbCover slope
  ForbCover_slopes.df <- data %>% 
    filter(!is.na(ForbCover_slopes))
  
  # Calculate mean ForbCover slope value
  ForbCover_slopes.mean <- mean(ForbCover_slopes.df$ForbCover_slopes)
  
  # Determine bin width for plotting
  ForbCover_slopes.bin.width <- (max(filter(data, !is.na(ForbCover_slopes))$ForbCover_slopes) - min(filter(data, !is.na(ForbCover_slopes))$ForbCover_slopes))/40
  
  # Plot model outputs
  (ForbCover.plot <- ggplot(data) +
      geom_histogram(aes(x = ForbCover_slopes), stat = "bin", binwidth = ForbCover_slopes.bin.width,
                     fill = "gray", alpha = 1) +
      geom_density(aes(x = ForbCover_slopes), stat = "bin", binwidth = ForbCover_slopes.bin.width,
                   fill = ForbCover_slopes.colour, colour = ForbCover_slopes.colour, alpha = 0.4, size = 0.5) +
      geom_vline(aes(xintercept = ForbCover_slopes.mean), colour = "#000000", linetype = "dashed", size = 1) +
      labs(y = "Number of Plots \n",
           x = "\n Change in ForbCover (per year)",
           title = paste0("Mean ForbCover Change: ", ForbCover.df$Estimate, " per year", "\n",
                          "95% CIs: [", ForbCover.df$`l-95% CI`, ", ", ForbCover.df$`u-95% CI`, "]"),
           subtitle = paste0(ifelse((ForbCover.df$"l-95% CI" < 0 & ForbCover.df$"u-95% CI" < 0) | # Automatically asigns significance in subtitle based on CIs
                                      (ForbCover.df$"l-95% CI" > 0 & ForbCover.df$"u-95% CI" > 0),
                                    "Significant: CIs do NOT span 0", "Not Significant: CIs span 0"))) +
      theme_1() +
      theme(plot.subtitle = element_text(colour = ifelse((ForbCover.df$"l-95% CI" < 0 & ForbCover.df$"u-95% CI" < 0) | # Automatically colours subtitle based on significance from CIs
                                                           (ForbCover.df$"l-95% CI" > 0 & ForbCover.df$"u-95% CI" > 0),
                                                         paste("#006400"), paste("#8b0000"))))
    
  )
  
  # Create a panel of all three plots
  combined.panel.Cover <- grid.arrange(ShrubCover.plot, GraminoidCover.plot, ForbCover.plot, ncol = 3)
  
  # Export panel
  ggsave(combined.panel.Cover, filename = paste0("figures/outputs_new/combined_Cover_change",
                                                 pc.filepath, ".png"), width = 23, height = 7)
  
}


# FUNCTION: TEMPORAL BAYESIAN MODELS FOR FD WITH A CONTINUOUS X-VARIABLE ----

bayesian.temporal.continuous <- function(run.FRic, run.FEve, run.SR, run.FDis, x.var, data, colour.by,
                                         chains, cores, warmup, iterations, delta, treedepth){
  
  # Produce dataframe for running models on
  input.data <- data %>% 
    dplyr::select(x.var, FRic_slopes, FEve_slopes, SR_slopes, FDis_slopes,
                  SurveyedArea, gridcell, SUBSITE, PlotDominatingFG) %>% 
    rename(x_variable = x.var) # Rename first column to x_variable
  
  # Modify data for colouring points
  data <- data %>% 
    rename(colour.variable = colour.by)
  
  # FRic ----
  
  # If you want to run model for first time
  if (run.FRic == TRUE){
    
    # Function for running the Bayesian model
    FRic.mod <- brm(FRic_slopes ~ x_variable + log(SurveyedArea) + (1 | gridcell / SUBSITE),
                    data = input.data, family = FRic_slopes.distribution, chains = chains,
                    warmup = warmup, iter = iterations, cores = cores,
                    control = list(adapt_delta = delta,
                                   max_treedepth = treedepth))
    
    # Export model output
    save(FRic.mod, file = paste0("data/model_outputs_new/m_FRic_slopes_",
                                 x.var, pc.filepath, ".RData"))
    
  } else {
    
    # Load in model output
    FRic.mod <- get(load(paste0("data/model_outputs_new/m_FRic_slopes_",
                                x.var, pc.filepath, ".RData")))
    
  }
  
  # Convert model output into a dataframe (with 4.d.p.)
  FRic.df <- brms_SummaryTable(FRic.mod, formatOptions = list(digits = 5, nsmall = 5), round = 5)
  
  # Extract the confidence intervals as a list for use in the plotting
  FRic.df.ci <- FRic.df %>% 
    filter(Covariate %in% c("x_variable")) %>% 
    dplyr::select("l-95% CI", "u-95% CI") %>% 
    pivot_longer(cols = 1:2, names_to = "Interval", values_to = "Value")
  
  # Save the confidence intervals as a list and remove the intermediate dataframe
  FRic.ci <- as.list(FRic.df.ci$Value) # Adding values to the list
  names(FRic.ci) <- FRic.df.ci$Interval # Adding names to the values
  rm(FRic.df, FRic.df.ci) # Remove unnecessary dataframes of summary and CIs  
  
  # Range to predict variable
  range.to.predict <- paste0("x_variable [", min(input.data$x_variable, na.rm = TRUE), ":", max(input.data$x_variable, na.rm = TRUE), 
                             " by =", (max(input.data$x_variable, na.rm = TRUE) - min(input.data$x_variable, na.rm = TRUE))/50, "]")

  # Extract the prediction data frame
  FRic.pred <- ggpredict(FRic.mod, terms = range.to.predict)
  
  # Extract dataframe with only the x-variable for plotting
  data.x.var.only <- dplyr::select(data, x.var)
  
  # Plot model outputs
  (FRic.plot <- ggplot(FRic.pred) +
      geom_line(aes(x = x, y = predicted), colour = "grey10") + # Adds line for predicted values of FRic vs LAT
      geom_ribbon(aes(x = x, ymin = conf.low, ymax = conf.high), fill = "grey10", alpha = 0.2) + # Adds c.intervals for predictions as ribbon
      geom_point(data = data, aes(x = data.x.var.only[,1], y = FRic_slopes, fill = colour.variable), # Adds original FRic vs LAT data points and colours by region
                 colour = c("#000000"), alpha = 0.5, shape = 21, size = 3) +
      scale_fill_viridis(option = FRic.colour, begin = 0.2, end = 1, direction = -1, discrete = TRUE) +
      labs(y = "FRic_slopes \n",
           x = paste0("\n ", x.var),
           fill = paste0(colour.by),
           title = paste0("CIs: ", FRic.ci$"l-95% CI", " to ", FRic.ci$"u-95% CI"),
           subtitle = paste(ifelse((FRic.ci$"l-95% CI" < 0 & FRic.ci$"u-95% CI" < 0) | # Automatically asigns significance in subtitle based on CIs
                                     (FRic.ci$"l-95% CI" > 0 & FRic.ci$"u-95% CI" > 0),
                                   "Significant: CIs do NOT span 0", "Not Significant: CIs span 0"))) +
      theme_1() +
      theme(legend.position = "bottom",
            plot.subtitle = element_text(colour = ifelse((FRic.ci$"l-95% CI" < 0 & FRic.ci$"u-95% CI" < 0) | # Automatically colours subtitle based on significance from CIs
                                                           (FRic.ci$"l-95% CI" > 0 & FRic.ci$"u-95% CI" > 0),
                                                         paste("#006400"), paste("#8b0000")))))

  # FEve ----
  
  # If you want to run model for first time
  if (run.FEve == TRUE){
    
    # Function for running the Bayesian model
    FEve.mod <- brm(FEve_slopes ~ x_variable + log(SurveyedArea) + (1 | gridcell / SUBSITE),
                    data = input.data, family = FEve_slopes.distribution, chains = chains,
                    warmup = warmup, iter = iterations, cores = cores,
                    control = list(adapt_delta = delta,
                                   max_treedepth = treedepth))
    
    # Export model output
    save(FEve.mod, file = paste0("data/model_outputs_new/m_FEve_slopes_",
                                 x.var, pc.filepath, ".RData"))
    
  } else {
    
    # Load in model output
    FEve.mod <- get(load(paste0("data/model_outputs_new/m_FEve_slopes_",
                                x.var, pc.filepath, ".RData")))
    
  }
  
  # Convert model output into a dataframe (with 4.d.p.)
  FEve.df <- brms_SummaryTable(FEve.mod, formatOptions = list(digits = 5, nsmall = 5), round = 5)
  
  # Extract the confidence intervals as a list for use in the plotting
  FEve.df.ci <- FEve.df %>% 
    filter(Covariate %in% c("x_variable")) %>% 
    dplyr::select("l-95% CI", "u-95% CI") %>% 
    pivot_longer(cols = 1:2, names_to = "Interval", values_to = "Value")
  
  # Save the confidence intervals as a list and remove the intermediate dataframe
  FEve.ci <- as.list(FEve.df.ci$Value) # Adding values to the list
  names(FEve.ci) <- FEve.df.ci$Interval # Adding names to the values
  rm(FEve.df, FEve.df.ci) # Remove unnecessary dataframes of summary and CIs  
  
  # Range to predict variable
  range.to.predict <- paste0("x_variable [", min(input.data$x_variable, na.rm = TRUE), ":", max(input.data$x_variable, na.rm = TRUE), 
                             " by =", (max(input.data$x_variable, na.rm = TRUE) - min(input.data$x_variable, na.rm = TRUE))/50, "]")
  
  # Extract the prediction data frame
  FEve.pred <- ggpredict(FEve.mod, terms = range.to.predict, back.transform = FALSE)
  
  # Extract dataframe with only the x-variable for plotting
  data.x.var.only <- dplyr::select(data, x.var)
  
  # Plot model outputs
  (FEve.plot <- ggplot(FEve.pred) +
      geom_line(aes(x = x, y = predicted), colour = "grey10") + # Adds line for predicted values of FEve vs LAT
      geom_ribbon(aes(x = x, ymin = conf.low, ymax = conf.high), fill = "grey10", alpha = 0.2) + # Adds c.intervals for predictions as ribbon
      geom_point(data = data, aes(x = data.x.var.only[,1], y = FEve_slopes, fill = colour.variable), # Adds original FEve vs LAT data points and colours by region
                 colour = c("#000000"), alpha = 0.5, shape = 21, size = 3) +
      scale_fill_viridis(option = FEve.colour, begin = 0.2, end = 1, direction = -1, discrete = TRUE) +
      labs(y = "FEve_slopes \n",
           x = paste0("\n ", x.var),
           fill = paste0(colour.by),
           title = paste0("CIs: ", FEve.ci$"l-95% CI", " to ", FEve.ci$"u-95% CI"),
           subtitle = paste(ifelse((FEve.ci$"l-95% CI" < 0 & FEve.ci$"u-95% CI" < 0) | # Automatically asigns significance in subtitle based on CIs
                                     (FEve.ci$"l-95% CI" > 0 & FEve.ci$"u-95% CI" > 0),
                                   "Significant: CIs do NOT span 0", "Not Significant: CIs span 0"))) +
      theme_1() +
      theme(legend.position = "bottom",
            plot.subtitle = element_text(colour = ifelse((FEve.ci$"l-95% CI" < 0 & FEve.ci$"u-95% CI" < 0) | # Automatically colours subtitle based on significance from CIs
                                                           (FEve.ci$"l-95% CI" > 0 & FEve.ci$"u-95% CI" > 0),
                                                         paste("#006400"), paste("#8b0000")))))

  # SR ----
  
  # If you want to run model for first time
  if (run.SR == TRUE){
    
    # Function for running the Bayesian model
    SR.mod <- brm(SR_slopes ~ x_variable + log(SurveyedArea) + (1 | gridcell / SUBSITE),
                  data = input.data, family = SR_slopes.distribution, chains = chains,
                  warmup = warmup, iter = iterations, cores = cores,
                  control = list(adapt_delta = delta,
                                 max_treedepth = treedepth))
    
    # Export model output
    save(SR.mod, file = paste0("data/model_outputs_new/m_SR_slopes_",
                               x.var, pc.filepath, ".RData"))
    
  } else {
    
    # Load in model output
    SR.mod <- get(load(paste0("data/model_outputs_new/m_SR_slopes_",
                              x.var, pc.filepath, ".RData")))
    
  }
  
  # Convert model output into a dataframe (with 4.d.p.)
  SR.df <- brms_SummaryTable(SR.mod, formatOptions = list(digits = 5, nsmall = 5), round = 5)
  
  # Extract the confidence intervals as a list for use in the plotting
  SR.df.ci <- SR.df %>% 
    filter(Covariate %in% c("x_variable")) %>% 
    dplyr::select("l-95% CI", "u-95% CI") %>% 
    pivot_longer(cols = 1:2, names_to = "Interval", values_to = "Value")
  
  # Save the confidence intervals as a list and remove the intermediate dataframe
  SR.ci <- as.list(SR.df.ci$Value) # Adding values to the list
  names(SR.ci) <- SR.df.ci$Interval # Adding names to the values
  rm(SR.df, SR.df.ci) # Remove unnecessary dataframes of summary and CIs  
  
  # Range to predict variable
  range.to.predict <- paste0("x_variable [", min(input.data$x_variable, na.rm = TRUE), ":", max(input.data$x_variable, na.rm = TRUE), 
                             " by =", (max(input.data$x_variable, na.rm = TRUE) - min(input.data$x_variable, na.rm = TRUE))/50, "]")
  
  # Extract the prediction data frame
  SR.pred <- ggpredict(SR.mod, terms = range.to.predict, back.transform = FALSE)
  
  # Extract dataframe with only the x-variable for plotting
  data.x.var.only <- dplyr::select(data, x.var)
  
  # Plot model outputs
  (SR.plot <- ggplot(SR.pred) +
      geom_line(aes(x = x, y = predicted), colour = "grey10") + # Adds line for predicted values of SR vs LAT
      geom_ribbon(aes(x = x, ymin = conf.low, ymax = conf.high), fill = "grey10", alpha = 0.2) + # Adds c.intervals for predictions as ribbon
      geom_point(data = data, aes(x = data.x.var.only[,1], y = SR_slopes, fill = colour.variable), # Adds original SR vs LAT data points and colours by region
                 colour = c("#000000"), alpha = 0.5, shape = 21, size = 3) +
      scale_fill_viridis(option = SR.colour, begin = 0.2, end = 1, direction = -1, discrete = TRUE) +
      labs(y = "SR_slopes \n",
           x = paste0("\n ", x.var),
           fill = paste0(colour.by),
           title = paste0("CIs: ", SR.ci$"l-95% CI", " to ", SR.ci$"u-95% CI"),
           subtitle = paste(ifelse((SR.ci$"l-95% CI" < 0 & SR.ci$"u-95% CI" < 0) | # Automatically asigns significance in subtitle based on CIs
                                     (SR.ci$"l-95% CI" > 0 & SR.ci$"u-95% CI" > 0),
                                   "Significant: CIs do NOT span 0", "Not Significant: CIs span 0"))) +
      theme_1() +
      theme(legend.position = "bottom",
            plot.subtitle = element_text(colour = ifelse((SR.ci$"l-95% CI" < 0 & SR.ci$"u-95% CI" < 0) | # Automatically colours subtitle based on significance from CIs
                                                           (SR.ci$"l-95% CI" > 0 & SR.ci$"u-95% CI" > 0),
                                                         paste("#006400"), paste("#8b0000")))))
  
  # FDis ----
  
  # If you want to run model for first time
  if (run.FDis == TRUE){
    
    # Function for running the Bayesian model
    FDis.mod <- brm(SR_slopes ~ x_variable + log(SurveyedArea) + (1 | gridcell / SUBSITE),
                    data = input.data, family = FDis_slopes.distribution, chains = chains,
                    warmup = warmup, iter = iterations, cores = cores,
                    control = list(adapt_delta = delta,
                                   max_treedepth = treedepth))
    
    # Export model output
    save(FDis.mod, file = paste0("data/model_outputs_new/m_FDis_slopes_",
                                 x.var, pc.filepath, ".RData"))
    
  } else {
    
    # Load in model output
    FDis.mod <- get(load(paste0("data/model_outputs_new/m_FDis_slopes_",
                                x.var, pc.filepath, ".RData")))
    
  }
  
  # Convert model output into a dataframe (with 4.d.p.)
  FDis.df <- brms_SummaryTable(FDis.mod, formatOptions = list(digits = 5, nsmall = 5), round = 5)
  
  # Extract the confidence intervals as a list for use in the plotting
  FDis.df.ci <- FDis.df %>% 
    filter(Covariate %in% c("x_variable")) %>% 
    dplyr::select("l-95% CI", "u-95% CI") %>% 
    pivot_longer(cols = 1:2, names_to = "Interval", values_to = "Value")
  
  # Save the confidence intervals as a list and remove the intermediate dataframe
  FDis.ci <- as.list(FDis.df.ci$Value) # Adding values to the list
  names(FDis.ci) <- FDis.df.ci$Interval # Adding names to the values
  rm(FDis.df, FDis.df.ci) # Remove unnecessary dataframes of summary and CIs  
  
  # Range to predict variable
  range.to.predict <- paste0("x_variable [", min(input.data$x_variable, na.rm = TRUE), ":", max(input.data$x_variable, na.rm = TRUE), 
                             " by =", (max(input.data$x_variable, na.rm = TRUE) - min(input.data$x_variable, na.rm = TRUE))/50, "]")
  
  # Extract the prediction data frame
  FDis.pred <- ggpredict(FDis.mod, terms = range.to.predict, back.transform = FALSE)
  
  # Extract dataframe with only the x-variable for plotting
  data.x.var.only <- dplyr::select(data, x.var)
  
  # Plot model outputs
  (FDis.plot <- ggplot(FDis.pred) +
      geom_line(aes(x = x, y = predicted), colour = "grey10") + # Adds line for predicted values of SR vs LAT
      geom_ribbon(aes(x = x, ymin = conf.low, ymax = conf.high), fill = "grey10", alpha = 0.2) + # Adds c.intervals for predictions as ribbon
      geom_point(data = data, aes(x = data.x.var.only[,1], y = FDis_slopes, fill = colour.variable), # Adds original SR vs LAT data points and colours by region
                 colour = c("#000000"), alpha = 0.5, shape = 21, size = 3) +
      scale_fill_viridis(option = FDis.colour, begin = 0.2, end = 1, direction = -1, discrete = TRUE) +
      labs(y = "FDis_slopes \n",
           x = paste0("\n ", x.var),
           fill = paste0(colour.by),
           title = paste0("CIs: ", FDis.ci$"l-95% CI", " to ", FDis.ci$"u-95% CI"),
           subtitle = paste(ifelse((FDis.ci$"l-95% CI" < 0 & FDis.ci$"u-95% CI" < 0) | # Automatically asigns significance in subtitle based on CIs
                                     (FDis.ci$"l-95% CI" > 0 & FDis.ci$"u-95% CI" > 0),
                                   "Significant: CIs do NOT span 0", "Not Significant: CIs span 0"))) +
      theme_1() +
      theme(legend.position = "bottom",
            plot.subtitle = element_text(colour = ifelse((FDis.ci$"l-95% CI" < 0 & FDis.ci$"u-95% CI" < 0) | # Automatically colours subtitle based on significance from CIs
                                                           (FDis.ci$"l-95% CI" > 0 & FDis.ci$"u-95% CI" > 0),
                                                         paste("#006400"), paste("#8b0000")))))
  
  
  # Create a panel of all three plots
  combined.panel <- grid.arrange(FRic.plot, FEve.plot, SR.plot, FDis.plot, ncol = 2)
  
  # Export panel
  ggsave(combined.panel, filename = paste0("figures/outputs_new/combined_slopes_",
                                           x.var, pc.filepath, ".png"), width = 15, height = 14)
  
}


# FUNCTION: SR vs FD METRICS ----

bayesian.metric.comparison <- function(run.FRic, quadratic, run.FEve, run.FDis, x.var, data,
                                       chains, cores, warmup, iterations, delta, treedepth){
  
  # Modify combo.latest for colouring points
  combo.latest <- combo.latest %>% 
    mutate(PlotDominatingFG = str_remove_all(PlotDominatingFG, pattern = "-Dominated"))
  
  
  # FRic ----
  
  # Run if loop for if linear
  if (quadratic == "No"){
    
    # Generate input dataframe
    input.data <- data %>% 
      dplyr::select(x.var, FRic, FEve, FDis, SurveyedArea, gridcell, SUBSITE) %>% 
      rename(x_variable = x.var) # Rename first column to x_variable
    
    # Run model in if loop
    if (run.model == TRUE){
      
      # Function for running the Bayesian model
      SR.FRic.mod <- brm(FRic ~ x_variable + log(SurveyedArea) + (1 | gridcell / SUBSITE),
                         data = input.data, family = FRic.distribution, chains = chains,
                         warmup = warmup, iter = iterations, cores = cores,
                         control = list(adapt_delta = delta,
                                        max_treedepth = treedepth))
      
      # Export model output
      save(SR.FRic.mod, file = paste0("data/model_outputs_new/m_FRic_",
                                      FRic.distribution, "_SR", pc.filepath, ".RData"))
      
    } else {
      
      # Load in model output
      SR.FRic.mod <- get(load(paste0("data/model_outputs_new/m_FRic_",
                                     FRic.distribution, "_SR", pc.filepath, ".RData")))
      
    } # End of if else loop
    
    # Convert model output into a dataframe (with 4.d.p.)
    SR.FRic.df <- brms_SummaryTable(SR.FRic.mod, formatOptions = list(digits = 4, nsmall = 4), round = 4)
    
    # Extract the confidence intervals as a list for use in the plotting
    SR.FRic.df.ci <- SR.FRic.df %>% 
      filter(Covariate %in% c("x_variable")) %>% 
      dplyr::select("l-95% CI", "u-95% CI") %>% 
      pivot_longer(cols = 1:2, names_to = "Interval", values_to = "Value")
    
    # Save the confidence intervals as a list and remove the intermediate dataframe
    SR.FRic.ci <- as.list(SR.FRic.df.ci$Value) # Adding values to the list
    names(SR.FRic.ci) <- SR.FRic.df.ci$Interval # Adding names to the values
    rm(SR.FRic.df, SR.FRic.df.ci) # Remove unnecessary dataframes of summary and CIs  
    
    # Extract the prediction data frame
    SR.FRic.pred <- ggpredict(SR.FRic.mod, terms = "x_variable [4:20]",  back.transform = TRUE)
    
    # Extract dataframe with only the x-variable for plotting
    data.x.var.only <- dplyr::select(data, SR)
    
    # Plot model outputs
    (SR.FRic.plot <- ggplot(SR.FRic.pred) +
        geom_line(aes(x = x, y = predicted), colour = "grey10") + # Adds line for predicted values of FRic vs LAT
        geom_ribbon(aes(x = x, ymin = conf.low, ymax = conf.high), fill = "grey10", alpha = 0.2) + # Adds c.intervals for predictions as ribbon
        geom_point(data = combo.latest, aes(x = data.x.var.only[,1], y = FRic, fill = colour.variable), # Adds original FRic vs LAT data points and colours by region
                   colour = c("#000000"), alpha = 0.5, shape = 21, size = 3) +
        # scale_fill_viridis(option = FRic.colour, begin = 0.2, end = 1, direction = -1, discrete = colour.discrete,
        #                    guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")) +
        scale_fill_gradient2(low = "#D02090",
                             # mid = "#228B22",
                             mid = "#FFFFFF",
                             high = "#1E90FF", midpoint = 0.5,
                             guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")) +
        scale_y_continuous(limits = c(0, max(combo.latest$FRic)*1.05)) +
        labs(y = "FRic \n",
             x = "\n SR",
             fill = paste0(colour.by),
             title = paste0("CIs: ", SR.FRic.ci$"l-95% CI", " to ", SR.FRic.ci$"u-95% CI"),
             subtitle = paste(ifelse((SR.FRic.ci$"l-95% CI" < 0 & SR.FRic.ci$"u-95% CI" < 0) | # Automatically asigns significance in subtitle based on CIs
                                       (SR.FRic.ci$"l-95% CI" > 0 & SR.FRic.ci$"u-95% CI" > 0),
                                     "Significant: CIs do NOT span 0", "Not Significant: CIs span 0"))) +
        theme_1() +
        theme(legend.position = "right",
              plot.subtitle = element_text(colour = ifelse((SR.FRic.ci$"l-95% CI" < 0 & SR.FRic.ci$"u-95% CI" < 0) | # Automatically colours subtitle based on significance from CIs
                                                             (SR.FRic.ci$"l-95% CI" > 0 & SR.FRic.ci$"u-95% CI" > 0),
                                                           paste("#006400"), paste("#8b0000")))))
    
    # Export plot
    ggsave(SR.FRic.plot, filename = paste0("figures/outputs_new/SR_vs_FRic",
                                           pc.filepath, ".png"), width = 8, height = 7)
    
  } # End of not quadratic if loop
  
  # Run if loop for is quadratic
  if (quadratic == "Yes"){
    
    # Generate input dataframe
    input.data <- data %>% 
      dplyr::select(x.var, FRic, FEve, FDis, SurveyedArea, gridcell, SUBSITE) %>% 
      rename(SR = x.var) # Rename first column to x_variable
    
    # Determine mean value for centering
    mean.SR <- mean(input.data$SR, na.rm = TRUE)
    
    # Center the x_variable
    input.data <- input.data %>% 
      mutate(centred_SR = SR - mean.SR)
    
    # Run model in if loop
    if (run.FRic == TRUE){
      
      # Function for running the Bayesian model
      SR.FRic.mod <- brm(FRic ~ centred_SR + I(centred_SR^2) + (1 | gridcell / SUBSITE),
                         data = input.data, family = FRic.distribution, chains = chains,
                         warmup = warmup, iter = iterations, cores = cores,
                         control = list(adapt_delta = delta,
                                        max_treedepth = treedepth))
      
      # Export model output
      save(SR.FRic.mod, file = paste0("data/model_outputs_new/m_FRic_",
                                      FRic.distribution, "_quadratic_SR", pc.filepath, ".RData"))
      
    } else {
      
      # Load in model output
      SR.FRic.mod <- get(load(paste0("data/model_outputs_new/m_FRic_",
                                     FRic.distribution, "_quadratic_SR", pc.filepath, ".RData")))
      
    } # End of if else loop
    
    # Check summary of model
    summary(SR.FRic.mod)
    
    # Convert model output into a dataframe (with 4.d.p.)
    SR.FRic.df <- brms_SummaryTable(SR.FRic.mod, formatOptions = list(digits = 4, nsmall = 4), round = 4)
    
    # Extract the confidence intervals as a list for use in the plotting
    SR.FRic.df.ci <- SR.FRic.df %>% 
      filter(Covariate %in% c("centred_SR")) %>% 
      dplyr::select("l-95% CI", "u-95% CI") %>% 
      pivot_longer(cols = 1:2, names_to = "Interval", values_to = "Value")
    
    # Save the confidence intervals as a list and remove the intermediate dataframe
    SR.FRic.ci <- as.list(SR.FRic.df.ci$Value) # Adding values to the list
    names(SR.FRic.ci) <- SR.FRic.df.ci$Interval # Adding names to the values
    rm(SR.FRic.df, SR.FRic.df.ci) # Remove unnecessary dataframes of summary and CIs  
    
    # Create template predictions dataframe
    SR.FRic.pred.df = data.frame(centred_SR = seq(0 - mean.SR, 100 - mean.SR))
    
    # Add predictions to template datframe
    SR.FRic.pred = add_epred_draws(SR.FRic.pred.df, SR.FRic.mod, re_formula = NA) %>% 
      mutate(SR = centred_SR + mean.SR) # Uncentre the x variable
    
    # Plot model outputs
    (SR.FRic.plot <- ggplot(data = SR.FRic.pred, aes(x = SR, y = .epred)) +
        
        geom_jitter(data = input.data, aes(x = SR, y = FRic), fill = "#D02090", colour = "#000000", shape = 21, alpha = 0.2, size = 2) +
        
        stat_lineribbon(aes(y = .epred), .width = 0.95, fill = "#D02090", alpha = 0.35, colour = "#D02090", linewidth = 1.1,
                        linetype = ifelse((SR.FRic.ci$"l-95% CI" < 0 & SR.FRic.ci$"u-95% CI" < 0) | # Automatically assigns significance in subtitle based on CIs
                                            (SR.FRic.ci$"l-95% CI" > 0 & SR.FRic.ci$"u-95% CI" > 0), 1, 2)) +
        
        stat_lineribbon(aes(y = .epred), .width = 0.95, fill = NA, colour = "#D02090", linewidth = 1.1,
                        linetype = ifelse((SR.FRic.ci$"l-95% CI" < 0 & SR.FRic.ci$"u-95% CI" < 0) | # Automatically assigns significance in subtitle based on CIs
                                            (SR.FRic.ci$"l-95% CI" > 0 & SR.FRic.ci$"u-95% CI" > 0), 1, 2)) +
        
        # geom_hex(data = combo.latest, aes(x = SR, y = FRic, fill = ), bins = 18) +
        # scale_fill_gradient(low = "#FFE1FF", high = "#D02090") +
        
        scale_x_continuous(limits = c(0,max(combo.latest$SR, na.rm = TRUE))) +
        scale_y_continuous(limits = c(0, 61)) +
        labs(y = "FRic \n",
             x = "SR \n",
             fill = "FEve \n",
             title = paste0("CIs: ", SR.FRic.ci$"l-95% CI", " to ", SR.FRic.ci$"u-95% CI"),
             subtitle = paste(ifelse((SR.FRic.ci$"l-95% CI" < 0 & SR.FRic.ci$"u-95% CI" < 0) | # Automatically asigns significance in subtitle based on CIs
                                       (SR.FRic.ci$"l-95% CI" > 0 & SR.FRic.ci$"u-95% CI" > 0),
                                     "Significant: CIs do NOT span 0", "Not Significant: CIs span 0"))) +
        theme_1() +
        theme(legend.position = "right",
              plot.subtitle = element_text(colour = ifelse((SR.FRic.ci$"l-95% CI" < 0 & SR.FRic.ci$"u-95% CI" < 0) | # Automatically colours subtitle based on significance from CIs
                                                             (SR.FRic.ci$"l-95% CI" > 0 & SR.FRic.ci$"u-95% CI" > 0),
                                                           paste("#006400"), paste("#8b0000")))))
    
    
    # FEve ----
    
      # Generate input dataframe
      input.data <- data %>%
        dplyr::select(x.var, FRic, FEve, FDis, SurveyedArea, gridcell, SUBSITE) %>%
        rename(x_variable = x.var) # Rename first column to x_variable
    
    # Run model in if loop
    if (run.FEve == TRUE){
      
      # Function for running the Bayesian model
      SR.FEve.mod <- brm(FEve ~ x_variable + log(SurveyedArea) + (1 | gridcell / SUBSITE),
                         data = input.data, family = FEve.distribution, chains = chains,
                         warmup = warmup, iter = iterations, cores = cores,
                         control = list(adapt_delta = delta,
                                        max_treedepth = treedepth))
      
      # Export model output
      save(SR.FEve.mod, file = paste0("data/model_outputs_new/m_FEve_",
                                      FEve.distribution, "_SR", pc.filepath, ".RData"))
      
    } else {
      
      # Load in model output
      SR.FEve.mod <- get(load(paste0("data/model_outputs_new/m_FEve_",
                                     FEve.distribution, "_SR", pc.filepath, ".RData")))
      
    } # End of if else loop
    
    # Convert model output into a dataframe (with 4.d.p.)
    SR.FEve.df <- brms_SummaryTable(SR.FEve.mod, formatOptions = list(digits = 4, nsmall = 4), round = 4)
    
    # Extract the confidence intervals as a list for use in the plotting
    SR.FEve.df.ci <- SR.FEve.df %>% 
      filter(Covariate %in% c("x_variable")) %>% 
      dplyr::select("l-95% CI", "u-95% CI") %>% 
      pivot_longer(cols = 1:2, names_to = "Interval", values_to = "Value")
    
    # Save the confidence intervals as a list and remove the intermediate dataframe
    SR.FEve.ci <- as.list(SR.FEve.df.ci$Value) # Adding values to the list
    names(SR.FEve.ci) <- SR.FEve.df.ci$Interval # Adding names to the values
    rm(SR.FEve.df, SR.FEve.df.ci) # Remove unnecessary dataframes of summary and CIs  
    
    # Range to predict variable
    range.to.predict <- paste0("x_variable [", min(input.data$x_variable), ":", max(input.data$x_variable), "]")
    
    # Extract the prediction data frame
    SR.FEve.pred <- ggpredict(SR.FEve.mod, terms = range.to.predict)

    # Plot model outputs
    (SR.FEve.plot <- ggplot(SR.FEve.pred) +
        
        geom_jitter(data = combo.latest, aes(x = SR, y = FEve), fill = "#228B22", colour = "#000000", alpha = 0.2, shape = 21, size = 3) +
        
        geom_line(aes(x = x, y = predicted), colour = "#228B22", linewidth = 1.1, # Adds line for predicted values
                  linetype = ifelse((SR.FEve.ci$"l-95% CI" < 0 & SR.FEve.ci$"u-95% CI" < 0) | # Automatically assigns significance in subtitle based on CIs
                                      (SR.FEve.ci$"l-95% CI" > 0 & SR.FEve.ci$"u-95% CI" > 0), 1, 2)) + 
        geom_ribbon(aes(x = x, ymin = conf.low, ymax = conf.high), fill = "#228B22", alpha = 0.35) + # Adds c.intervals for predictions as ribbon
        
        scale_y_continuous(limits = c(0, max(combo.latest$FEve)*1.05)) +
        labs(y = "FEve \n",
             x = "\n SR",
             fill = paste0(colour.by),
             title = paste0("CIs: ", SR.FEve.ci$"l-95% CI", " to ", SR.FEve.ci$"u-95% CI"),
             subtitle = paste(ifelse((SR.FEve.ci$"l-95% CI" < 0 & SR.FEve.ci$"u-95% CI" < 0) | # Automatically asigns significance in subtitle based on CIs
                                       (SR.FEve.ci$"l-95% CI" > 0 & SR.FEve.ci$"u-95% CI" > 0),
                                     "Significant: CIs do NOT span 0", "Not Significant: CIs span 0"))) +
        theme_1() +
        theme(legend.position = "right",
              plot.subtitle = element_text(colour = ifelse((SR.FEve.ci$"l-95% CI" < 0 & SR.FEve.ci$"u-95% CI" < 0) | # Automatically colours subtitle based on significance from CIs
                                                             (SR.FEve.ci$"l-95% CI" > 0 & SR.FEve.ci$"u-95% CI" > 0),
                                                           paste("#006400"), paste("#8b0000")))))
    
    
    # FDis ----
    
    # Run model in if loop
    if (run.FDis == TRUE){
      
      # Function for running the Bayesian model
      SR.FDis.mod <- brm(FDis ~ x_variable + log(SurveyedArea) + (1 | gridcell / SUBSITE),
                         data = input.data, family = FDis.distribution, chains = chains,
                         warmup = warmup, iter = iterations, cores = cores,
                         control = list(adapt_delta = delta,
                                        max_treedepth = treedepth))
      
      # Export model output
      save(SR.FDis.mod, file = paste0("data/model_outputs_new/m_FDis_",
                                      FDis.distribution, "_SR", pc.filepath, ".RData"))
      
    } else {
      
      # Load in model output
      SR.FDis.mod <- get(load(paste0("data/model_outputs_new/m_FDis_",
                                     FDis.distribution, "_SR", pc.filepath, ".RData")))
      
    } # End of if else loop
    
    # Convert model output into a dataframe (with 4.d.p.)
    SR.FDis.df <- brms_SummaryTable(SR.FDis.mod, formatOptions = list(digits = 4, nsmall = 4), round = 4)
    
    # Extract the confidence intervals as a list for use in the plotting
    SR.FDis.df.ci <- SR.FDis.df %>% 
      filter(Covariate %in% c("x_variable")) %>% 
      dplyr::select("l-95% CI", "u-95% CI") %>% 
      pivot_longer(cols = 1:2, names_to = "Interval", values_to = "Value")
    
    # Save the confidence intervals as a list and remove the intermediate dataframe
    SR.FDis.ci <- as.list(SR.FDis.df.ci$Value) # Adding values to the list
    names(SR.FDis.ci) <- SR.FDis.df.ci$Interval # Adding names to the values
    rm(SR.FDis.df, SR.FDis.df.ci) # Remove unnecessary dataframes of summary and CIs  
    
    # Extract the prediction data frame
    SR.FDis.pred <- ggpredict(SR.FDis.mod, terms = "x_variable [4:20]",  back.transform = TRUE)
    
    # Extract dataframe with only the x-variable for plotting
    data.x.var.only <- dplyr::select(data, SR)
    
    # Plot model outputs
    (SR.FDis.plot <- ggplot(SR.FDis.pred) +
        
        geom_jitter(data = combo.latest, aes(x = SR, y = FDis), fill = "#EE7600", colour = "#000000", alpha = 0.2, shape = 21, size = 3) +
        
        geom_line(aes(x = x, y = predicted), colour = "#EE7600", linewidth = 1.1, # Adds line for predicted values
                  linetype = ifelse((SR.FDis.ci$"l-95% CI" < 0 & SR.FDis.ci$"u-95% CI" < 0) | # Automatically assigns significance in subtitle based on CIs
                                      (SR.FDis.ci$"l-95% CI" > 0 & SR.FDis.ci$"u-95% CI" > 0), 1, 2)) + 
        geom_ribbon(aes(x = x, ymin = conf.low, ymax = conf.high), fill = "#EE7600", alpha = 0.35) + # Adds c.intervals for predictions as ribbon
        
        scale_y_continuous(limits = c(0, max(combo.latest$FDis)*1.05)) +
        labs(y = "FDis \n",
             x = "\n SR",
             fill = paste0(colour.by),
             title = paste0("CIs: ", SR.FDis.ci$"l-95% CI", " to ", SR.FDis.ci$"u-95% CI"),
             subtitle = paste(ifelse((SR.FDis.ci$"l-95% CI" < 0 & SR.FDis.ci$"u-95% CI" < 0) | # Automatically asigns significance in subtitle based on CIs
                                       (SR.FDis.ci$"l-95% CI" > 0 & SR.FDis.ci$"u-95% CI" > 0),
                                     "Significant: CIs do NOT span 0", "Not Significant: CIs span 0"))) +
        theme_1() +
        theme(legend.position = "right",
              plot.subtitle = element_text(colour = ifelse((SR.FDis.ci$"l-95% CI" < 0 & SR.FDis.ci$"u-95% CI" < 0) | # Automatically colours subtitle based on significance from CIs
                                                             (SR.FDis.ci$"l-95% CI" > 0 & SR.FDis.ci$"u-95% CI" > 0),
                                                           paste("#006400"), paste("#8b0000")))))
    
    # Panel ----
    
    # Create panel of all three plots
    SR.panel <- grid.arrange(SR.FRic.plot, SR.FEve.plot, SR.FDis.plot, ncol = 3)
    
    # Export plot
    ggsave(SR.panel, filename = paste0("figures/outputs_new/SR_vs_FD_metrics",
                                       pc.filepath, ".png"), width = 17.25, height = 5.5)
    
  } # End of quadratic if loop
  
} # End of function


# FUNCTION: SR vs FD GAMs ----

GAM.metric.comparison <- function(run){
  
  if(run == TRUE){
    
    # Modify combo.latest for colouring points
    combo.latest <- combo.latest %>% 
      mutate(PlotDominatingFG = str_remove_all(PlotDominatingFG, pattern = "-Dominated"))
    
    
    # FRic ----
    
    # Run the GAM for SR vs FRic
    SR.FRic.GAM <- mgcv::gam(FRic ~ s(SR, k = 2), data = combo.latest)
    
    # Save outputs as a list
    SR.FRic.GAM.list <- as.list(summary(SR.FRic.GAM))
    
    # Generate output as dataframe
    SR.FRic.GAM.output <- data.frame(SR.FRic.GAM.list[1:4])
    
    # Extract the significance value of the GAM
    SR.FRic.GAM.signif <- as.numeric(SR.FRic.GAM.output[1,4])
    
    # Range to predict variable
    range.to.predict <- paste0("SR [", min(combo.latest$SR), ":", max(combo.latest$SR), "]")
    
    # Extract preditions from GAM
    SR.FRic.GAM.predictions <- data.frame(ggpredict(SR.FRic.GAM, terms = range.to.predict))
    
    # Plot the outputs
    (SR.FRic.GAM.plot <- ggplot() +
        geom_point(data = combo.latest, aes(x = SR, y = FRic), colour = "#000000", fill = "#D02090", shape = 21, size = 2, alpha = 0.2) +
        geom_ribbon(data = SR.FRic.GAM.predictions, aes(x = x, ymin = conf.low, ymax = conf.high),
                    fill = "#D02090", alpha = 0.2, colour = "#D02090", linewidth = 0.05) + # Adds c.intervals for predictions as ribbon
        geom_line(data = SR.FRic.GAM.predictions, aes(x = x, y = predicted), colour = "#000000", linewidth = 1.1, # Adds line for predicted values of FRic vs LAT
                  linetype = ifelse(SR.FRic.GAM.signif >= 0.05, 2, 1)) +
        labs(y = "Functional Richness\n",
             x = "Species Richness",
             title = paste0("GAM: p = ", SR.FRic.GAM.signif),
             subtitle = paste(ifelse(SR.FRic.GAM.signif >= 0.05, "Not Significant: p >= 0.05", "Significant: p < 0.05"))) +
        theme_1() +
        theme(legend.position = "none",
              plot.subtitle = element_text(colour = ifelse(SR.FRic.GAM.signif < 0.05,
                                                           paste("#006400"), paste("#8b0000")))))
    
    
    # FEve ----
    
    # Run the GAM for SR vs FEve
    SR.FEve.GAM <- mgcv::gam(FEve ~ s(SR, k = 2), data = combo.latest)
    
    # Save outputs as a list
    SR.FEve.GAM.list <- as.list(summary(SR.FEve.GAM))
    
    # Generate output as dataframe
    SR.FEve.GAM.output <- data.frame(SR.FEve.GAM.list[1:4])
    
    # Extract the significance value of the GAM
    SR.FEve.GAM.signif <- as.numeric(SR.FEve.GAM.output[1,4])
    
    # Range to predict variable
    range.to.predict <- paste0("SR [", min(combo.latest$SR), ":", max(combo.latest$SR), "]")
    
    # Extract preditions from GAM
    SR.FEve.GAM.predictions <- data.frame(ggpredict(SR.FEve.GAM, terms = range.to.predict))
    
    # Plot the outputs
    (SR.FEve.GAM.plot <- ggplot() +
        geom_point(data = combo.latest, aes(x = SR, y = FEve), colour = "#000000", fill = "#228B22", shape = 21, size = 2, alpha = 0.2) +
        geom_ribbon(data = SR.FEve.GAM.predictions, aes(x = x, ymin = conf.low, ymax = conf.high),
                    fill = "#228B22", alpha = 0.2, colour = "#228B22", linewidth = 0.05) + # Adds c.intervals for predictions as ribbon
        geom_line(data = SR.FEve.GAM.predictions, aes(x = x, y = predicted), colour = "#000000", linewidth = 1.1, # Adds line for predicted values of FEve vs LAT
                  linetype = ifelse(SR.FEve.GAM.signif >= 0.05, 2, 1)) +
        labs(y = "Functional Evenness\n",
             x = "Species Richness",
             title = paste0("GAM: p = ", SR.FEve.GAM.signif),
             subtitle = paste(ifelse(SR.FEve.GAM.signif >= 0.05, "Not Significant: p >= 0.05", "Significant: p < 0.05"))) +
        theme_1() +
        theme(legend.position = "none",
              plot.subtitle = element_text(colour = ifelse(SR.FEve.GAM.signif < 0.05,
                                                           paste("#006400"), paste("#8b0000")))))
    
    
    # FDis ----
    
    # Run the GAM for SR vs FDis
    SR.FDis.GAM <- mgcv::gam(FDis ~ s(SR, k = 2), data = combo.latest)
    
    # Save outputs as a list
    SR.FDis.GAM.list <- as.list(summary(SR.FDis.GAM))
    
    # Generate output as dataframe
    SR.FDis.GAM.output <- data.frame(SR.FDis.GAM.list[1:4])
    
    # Extract the significance value of the GAM
    SR.FDis.GAM.signif <- as.numeric(SR.FDis.GAM.output[1,4])
    
    # Range to predict variable
    range.to.predict <- paste0("SR [", min(combo.latest$SR), ":", max(combo.latest$SR), "]")
    
    # Extract preditions from GAM
    SR.FDis.GAM.predictions <- data.frame(ggpredict(SR.FDis.GAM, terms = range.to.predict))
    
    # Plot the outputs
    (SR.FDis.GAM.plot <- ggplot() +
        geom_point(data = combo.latest, aes(x = SR, y = FDis), colour = "#000000", fill = "#EE7600", shape = 21, size = 2, alpha = 0.2) +
        geom_ribbon(data = SR.FDis.GAM.predictions, aes(x = x, ymin = conf.low, ymax = conf.high),
                    fill = "#EE7600", alpha = 0.2, colour = "#EE7600", linewidth = 0.05) + # Adds c.intervals for predictions as ribbon
        geom_line(data = SR.FDis.GAM.predictions, aes(x = x, y = predicted), colour = "#000000", linewidth = 1.1, # Adds line for predicted values of FDis vs LAT
                  linetype = ifelse(SR.FDis.GAM.signif >= 0.05, 2, 1)) +
        labs(y = "Functional Dispersion\n",
             x = "Species Richness",
             title = paste0("GAM: p = ", SR.FDis.GAM.signif),
             subtitle = paste(ifelse(SR.FDis.GAM.signif >= 0.05, "Not Significant: p >= 0.05", "Significant: p < 0.05"))) +
        theme_1() +
        theme(legend.position = "none",
              plot.subtitle = element_text(colour = ifelse(SR.FDis.GAM.signif < 0.05,
                                                           paste("#006400"), paste("#8b0000")))))
    
    
    # Panel ----
    
    # Create panel
    SR.GAM.panel <- grid.arrange(SR.FRic.GAM.plot, SR.FEve.GAM.plot, SR.FDis.GAM.plot, ncol = 3)
    
    # Export the panel
    ggsave(SR.GAM.panel, filename = paste0("figures/outputs_new/SR_vs_FD_GAMS",
                                           pc.filepath, ".png"), width = 17.25, height = 5.5)    
    
  } # End of if run == TRUE
  
} # End of function


# FUNCTION: WRITE FUNCTION TO CHECK CONVERGENCE OF EACH MODEL ----

bayesian.results.convergence.warnings <- function(){
  
  # Run this function if check.convergence == TRUE
  if (check.convergence == TRUE){
    
    # Load bespoke tryCatch function into the local environment
    # Found at: https://stackoverflow.com/questions/4948361/how-do-i-save-warnings-and-errors-as-output-from-a-function
    myTryCatch <- function(expr) {
      warn <- err <- NULL
      value <- withCallingHandlers(
        tryCatch(expr, error=function(e) {
          err <<- e
          NULL
        }), warning=function(w) {
          warn <<- w
          invokeRestart("muffleWarning")
        })
      list(value=value, warning=warn, error=err)
    }
    
    # Determine all the model filenames in the 'model_outputs_new' folder
    model.filepaths <- paste0("data/model_outputs_new/",
                              list.files("data/model_outputs_new/",
                                         pattern = paste0("*", pc.filepath, ".RData")))
    
    # Generate blank output dataframe
    convergence.warnings <- data.frame()
    
    # Add in model counter for checking
    model.counter <- 0
    
    # Run loops to check results and convergence for each metric and each variable
    for (i in model.filepaths){
      
      # Update model counter
      model.counter <- model.counter + 1
      
      # Load in the model output
      model <- get(load(paste0(i)))
      
      # Convert model output into a dataframe (with 4.d.p.)
      model.df <- brms_SummaryTable(model, formatOptions = list(digits = 4, nsmall = 4), round = 4)
      
      # Trimmed dataframe
      model.df.cut <- model.df %>% 
        mutate(Retain = ifelse(str_detect(Covariate, "x_variable"), TRUE, FALSE),
               Retain = ifelse(str_detect(Covariate, "centred_SR"), TRUE, Retain),
               Retain = ifelse(str_detect(Covariate, "Intercept"), TRUE, Retain)) %>% 
        filter(Retain == TRUE) %>% 
        rename(lower.CI = "l-95% CI",
               upper.CI = "u-95% CI") %>% 
        mutate(Estimate = as.numeric(Estimate),
               Est.Error = as.numeric(Est.Error),
               lower.CI = as.numeric(lower.CI),
               upper.CI = as.numeric(upper.CI))
      
      # Run loops depending on whether continuous x-variable or y-variable
      if (nrow(model.df.cut) == 2){
        
        # Trim to just the x-variable
        model.df.cut <- model.df.cut %>% 
          mutate(Retain = ifelse(str_detect(Covariate, "x_variable"), TRUE, FALSE)) %>% 
          filter(Retain == TRUE) %>%
          dplyr::select(-Retain)
        
        # Determine the x variable type, estimate, lower CI, upper CI, and whether CIs span zero
        x.var.type <- "continuous"
        estimate <- model.df.cut$Estimate
        lower.interval <- model.df.cut$lower.CI
        upper.interval <- model.df.cut$upper.CI
        CI.span <- ifelse((lower.interval < 0 & upper.interval < 0) | 
                            (lower.interval > 0 & upper.interval > 0), "Yes", "No")
        
        # Get model summmary within tryCatch
        model.summary <- myTryCatch(summary(model))
        
        # Extract the rhat values from the model output
        rhat.full <- data.frame(rhat(model)) %>% 
          rename(rhat = rhat.model.)
        
        # Determine whether rhat values are problematic (e.g. +- 0.05 from 1)
        rhat.summary <- rhat.full %>% 
          mutate(problematic = ifelse(rhat >= 1.05 | rhat <= 0.95, TRUE, FALSE))
        
        # Extract useful statistics regarding rhats
        rhat.mean <- as.numeric(round(mean(rhat.full$rhat), digits = 2))
        rhat.problematic.count <- as.numeric(nrow(filter(rhat.summary, problematic == TRUE)))
        
        # Extract just model name
        model.name <- str_remove(i, pattern = "data/model_outputs_new/")
        model.name <- str_remove(model.name, pattern = ".RData")
        
        # Add warnings and model name to dataframe
        convergence.summary <- data.frame("model" = paste0(model.name),
                                          "x.var.type" = x.var.type,
                                          "Covariate" = "REPLACE",
                                          "estimate" = estimate,
                                          "lower.CI" = lower.interval,
                                          "upper.CI" = upper.interval,
                                          "significant" = CI.span,
                                          "warning" = ifelse(is.null(model.summary$warning$message[1]),
                                                             NA, as.character(model.summary$warning$message[1])),
                                          "error" = ifelse(is.null(model.summary$error$message[1]),
                                                           NA, as.character(model.summary$error$message[1])),
                                          "rhat" = rhat.mean,
                                          "problematic_rhat_count" = rhat.problematic.count)
        
        # Append convergence summary to the overall convergence output
        convergence.warnings <- rbind(convergence.warnings, convergence.summary)
        
      } else { # End of continuous (e.g. 1 row)
        
        # Extract just model name
        model.name <- str_remove(i, pattern = "data/model_outputs_new/")
        model.name <- str_remove(model.name, pattern = ".RData")
        
        # Modify table for inclusion in outputs
        convergence.summary <- model.df.cut %>% 
          dplyr::select(-Retain, -Est.Error) %>% 
          mutate(model = model.name,
                 x.var.type = "Categoric",
                 Covariate = ifelse(str_detect(Covariate, pattern = "Intercept"), "Reference", Covariate),
                 Covariate = str_remove(Covariate, pattern = "x_variable"),
                 significant = ifelse((lower.interval < 0 & upper.interval < 0) | 
                                        (lower.interval > 0 & upper.interval > 0), "Yes", "No"),
                 warning = NA, error = NA, rhat = NA, problematic_rhat_count = NA) %>% 
          rename(estimate = Estimate) %>% 
          relocate(model, x.var.type, .after = ) 
        
        # Append convergence summary to the overall convergence output
        convergence.warnings <- rbind(convergence.warnings, convergence.summary)
        
      } # End of if else loop for categoric (e.g. more than 1 row)
      
    } # End of variables loop
    
    # Export the convergence output
    write.csv(convergence.warnings, file = paste0("data/model_outputs_new/",
                                                  "convergence_summaries", pc.filepath, ".csv"), row.names = FALSE)
    
  } # End of if statement for check.convergence == TRUE
  
} # End of function

# bayesian.results.convergence.warnings <- function(){
#   
#   # Run this function if check.convergence == TRUE
#   if (check.convergence == TRUE){
#     
#     # Load bespoke tryCatch function into the local environment
#     # Found at: https://stackoverflow.com/questions/4948361/how-do-i-save-warnings-and-errors-as-output-from-a-function
#     myTryCatch <- function(expr) {
#       warn <- err <- NULL
#       value <- withCallingHandlers(
#         tryCatch(expr, error=function(e) {
#           err <<- e
#           NULL
#         }), warning=function(w) {
#           warn <<- w
#           invokeRestart("muffleWarning")
#         })
#       list(value=value, warning=warn, error=err)
#     }
#     
#     # Determine all the model filenames in the 'model_outputs_new' folder
#     model.filepaths <- paste0("data/model_outputs_new/",
#                         list.files("data/model_outputs_new/",
#                                    pattern = paste0("*", pc.filepath, ".RData")))
#     
#     # Generate blank output dataframe
#     convergence.warnings <- data.frame()
#     
#     # Add in model counter for checking
#     model.counter <- 0
#     
#     # Run loops to check results and convergence for each metric and each variable
#     for (i in model.filepaths){
#         
#       # Update model counter
#       model.counter <- model.counter + 1
#       
#       # Load in the model output
#       model <- get(load(paste0(i)))
#       
#       # Convert model output into a dataframe (with 4.d.p.)
#       model.df <- brms_SummaryTable(model, formatOptions = list(digits = 4, nsmall = 4), round = 4)
#       
#       # Trimmed dataframe
#       model.df.cut <- model.df %>% 
#         mutate(Retain = ifelse(str_detect(Covariate, "x_variable"), TRUE, FALSE)) %>% 
#         filter(Retain == TRUE) %>% 
#         rename(lower.CI = "l-95% CI",
#                upper.CI = "u-95% CI") %>% 
#         mutate(Estimate = as.numeric(Estimate),
#                Est.Error = as.numeric(Est.Error),
#                lower.CI = as.numeric(lower.CI),
#                upper.CI = as.numeric(upper.CI))
#       
#       # Run loops depending on whether continuous x-variable or y-variable
#       if (nrow(model.df.cut) == 1){
#         
#         # Determine the x variable type, estimate, lower CI, upper CI, and whether CIs span zero
#         x.var.type <- "continuous"
#         estimate <- model.df.cut$Estimate
#         lower.interval <- model.df.cut$lower.CI
#         upper.interval <- model.df.cut$upper.CI
#         CI.span <- ifelse((lower.interval < 0 & upper.interval < 0) | 
#                             (lower.interval > 0 & upper.interval > 0), "Yes", "No")
#         
#       } else { # End of continuous (e.g. 1 row)
#         
#         # Determine the x variable type, estimate, lower CI, upper CI, and whether CIs span zero
#         x.var.type <- "categoric"
#         estimate <- NA
#         lower.interval <- NA
#         upper.interval <- NA
# 
#         # Generate a determining if any variables in x significant (0 = NO, 1 = YES)
#         model.df.CIs <- model.df.cut %>% 
#           mutate(Significant = 0,
#                  Significant = ifelse(lower.CI < 0 & upper.CI < 0, 1, Significant),
#                  Significant = ifelse(lower.CI > 0 & upper.CI > 0, 1, Significant))
# 
#         # Determine if any CIs significant
#         CI.span <- ifelse(sum(model.df.CIs$Significant) > 0, "Yes", "No")
# 
#       } # End of if else loop for categoric (e.g. more than 1 row)
#         
#       
#       # Get model summmary within tryCatch
#       model.summary <- myTryCatch(summary(model))
#       
#       # Extract the rhat values from the model output
#       rhat.full <- data.frame(rhat(model)) %>% 
#         rename(rhat = rhat.model.)
#       
#       # Determine whether rhat values are problematic (e.g. +- 0.05 from 1)
#       rhat.summary <- rhat.full %>% 
#         mutate(problematic = ifelse(rhat >= 1.05 | rhat <= 0.95, TRUE, FALSE))
#       
#       # Extract useful statistics regarding rhats
#       rhat.mean <- as.numeric(round(mean(rhat.full$rhat), digits = 2))
#       rhat.problematic.count <- as.numeric(nrow(filter(rhat.summary, problematic == TRUE)))
#       
#       # Extract just model name
#       model.name <- str_remove(i, pattern = "data/model_outputs_new/")
#       
#       # Add warnings and model name to dataframe
#       convergence.summary <- data.frame("model" = paste0(model.name),
#                                         "x.var.type" = x.var.type,
#                                         "estimate" = estimate,
#                                         "lower.CI" = lower.interval,
#                                         "upper.CI" = upper.interval,
#                                         "significant" = CI.span,
#                                         "warning" = ifelse(is.null(model.summary$warning$message[1]),
#                                                            NA, as.character(model.summary$warning$message[1])),
#                                         "error" = ifelse(is.null(model.summary$error$message[1]),
#                                                          NA, as.character(model.summary$error$message[1])),
#                                         "rhat" = rhat.mean,
#                                         "problematic_rhat_count" = rhat.problematic.count)
#       
#       # Append convergence summary to the overall convergence output
#       convergence.warnings <- rbind(convergence.warnings, convergence.summary)
#       
#     } # End of variables loop
#     
#     # Export the convergence output
#     write.csv(convergence.warnings, file = paste0("data/model_outputs_new/",
#                                                   "convergence_summaries", pc.filepath, ".csv"), row.names = FALSE)
#     
#   } # End of if statement for check.convergence == TRUE
#   
# } # End of function


# FUNCTION: LINEAR METRIC COMPARISON ----

bayesian.pc.count.comparison <- function(run.FRic, run.FEve, pc.count.1, pc.count.2, colour.by,
                                         chains, cores, warmup, iterations, delta, treedepth){
  
  # Import FD outputs for each PC count
  FD.output.1 <- read.csv(paste0("data/output_08_fd_output_combined",
                                 pc.filepath, "_all_years.csv"))
  FD.output.2 <- read.csv(paste0("data/output_08_fd_output_combined",
                                 pc.filepath, "_all_years.csv"))
  
  # Modify dataframes
  FD.tidy.1 <- FD.output.1 %>%
    dplyr::select(SiteSubsitePlotYear, FRic, FEve, SR) %>% # Leave SR to check same in both
    rename(FRic.1 = FRic, FEve.1 = FEve)
  FD.tidy.2 <- FD.output.2 %>%
    dplyr::select(SiteSubsitePlotYear, FRic, FEve, SR, SurveyedArea,
                  gridcell, SUBSITE, PlotDominatingFG) %>% 
    rename(FRic.2 = FRic, FEve.2 = FEve)
  
  # Join together the two dataframes
  FD.tidy.complete <- left_join(FD.tidy.2, FD.tidy.1,
                                by = c("SR" = "SR", "SiteSubsitePlotYear" ="SiteSubsitePlotYear")) %>% 
    dplyr::select(-SR) %>% 
    relocate(SiteSubsitePlotYear, FRic.1, FRic.2, FEve.1, FEve.2, .after = ) %>% 
    mutate(PlotDominatingFG = ifelse(PlotDominatingFG == "Evergreen_Shrub-Dominated",
                                     "E.Shrub", PlotDominatingFG),
           PlotDominatingFG = ifelse(PlotDominatingFG == "Deciduous_Shrub-Dominated",
                                     "D.Shrub", PlotDominatingFG),
           PlotDominatingFG = ifelse(PlotDominatingFG == "Graminoid-Dominated",
                                     "Graminoid", PlotDominatingFG),
           PlotDominatingFG = ifelse(PlotDominatingFG == "Forb-Dominated",
                                     "Forb", PlotDominatingFG)) %>% 
    rename(colour.variable = colour.by)
  
  if (run.FRic == TRUE){
    
    # Function for running the Bayesian model
    FRic.mod <- brm(FRic.2 ~ FRic.1 + log(SurveyedArea) + (1 | gridcell / SUBSITE),
                    data = FD.tidy.complete, family = FRic.distribution, chains = chains,
                    warmup = warmup, iter = iterations, cores = cores,
                    control = list(adapt_delta = delta,
                                   max_treedepth = treedepth))
    
    # Export model output
    save(FRic.mod, file = paste0("data/model_outputs_new/m_FRic_",
                                 "PC_comparison_pc", pc.count.1, "_vs_pc", pc.count.2, ".RData"))
                                 

  } else {
    
    # Load in model output
    FRic.mod <- get(load(paste0("data/model_outputs_new/m_FRic_",
                                "PC_comparison_pc", pc.count.1, "_vs_pc", pc.count.2, ".RData")))
    
  } # End of FRic if loop
  
  # Convert model output into a dataframe (with 4.d.p.)
  FRic.df <- brms_SummaryTable(FRic.mod, formatOptions = list(digits = 4, nsmall = 4), round = 4)
  
  # Extract the confidence intervals as a list for use in the plotting
  FRic.df.ci <- FRic.df %>% 
    filter(Covariate %in% c("FRic.1")) %>% 
    dplyr::select("l-95% CI", "u-95% CI") %>% 
    pivot_longer(cols = 1:2, names_to = "Interval", values_to = "Value")
  
  # Save the confidence intervals as a list and remove the intermediate dataframe
  FRic.ci <- as.list(FRic.df.ci$Value) # Adding values to the list
  names(FRic.ci) <- FRic.df.ci$Interval # Adding names to the values
  rm(FRic.df, FRic.df.ci) # Remove unnecessary dataframes of summary and CIs
  
  # Extract the prediction data frame
  FRic.pred <- ggpredict(FRic.mod, terms = "FRic.1 [sample = 30]",
                         # type = "random", condition = c(gridcell = 1266), allow_new_levels = TRUE,
                         back.transform = TRUE)
  
  # Extract dataframe with only the x-variable for plotting
  data.x.var.only <- dplyr::select(FD.tidy.complete, FRic.1)
  
  # Plot the modelled comparisons against one another
  (pc.comparison.FRic <- ggplot() +
      geom_point(data = FD.tidy.complete, aes(x = FRic.1, y = FRic.2, fill = colour.variable),
        colour = "#000000", alpha = 0.5, shape = 21, size = 2) +
      
      geom_line(data = FRic.pred, aes(x = x, y = predicted), colour = "grey10") + # Adds line for predicted values of FRic vs LAT
      geom_ribbon(data = FRic.pred, aes(x = x, ymin = conf.low, ymax = conf.high),
                  fill = "grey10", alpha = 0.2) + # Adds c.intervals for predictions as ribbon
      scale_fill_viridis(option = FRic.colour, begin = 0.2, end = 1, direction = -1, discrete = TRUE) +
      labs(y = paste0("FRic PC", pc.count.2, "\n"),
           x = paste0("\n FRic PC", pc.count.1),
           title = paste0("CIs: ", FRic.ci$"l-95% CI", " to ", FRic.ci$"u-95% CI"),
           subtitle = paste(ifelse((FRic.ci$"l-95% CI" < 0 & FRic.ci$"u-95% CI" < 0) | # Automatically asigns significance in subtitle based on CIs
                                     (FRic.ci$"l-95% CI" > 0 & FRic.ci$"u-95% CI" > 0),
                                   "Significant: CIs do NOT span 0", "Not Significant: CIs span 0")),
           fill = "Dominant \n Functional \n Group") +
      theme_1() +
      theme(legend.position = "right",
            plot.subtitle = element_text(colour = ifelse((FRic.ci$"l-95% CI" < 0 & FRic.ci$"u-95% CI" < 0) | # Automatically colours subtitle based on significance from CIs
                                                           (FRic.ci$"l-95% CI" > 0 & FRic.ci$"u-95% CI" > 0),
                                                         paste("#006400"), paste("#8b0000")))))
  
  if (run.FEve == TRUE){
    
    # Function for running the Bayesian model
    FEve.mod <- brm(FEve.2 ~ FEve.1 + log(SurveyedArea) + (1 | gridcell / SUBSITE),
                    data = FD.tidy.complete, family = FEve.distribution, chains = chains,
                    warmup = warmup, iter = iterations, cores = cores,
                    control = list(adapt_delta = delta,
                                   max_treedepth = treedepth))
    
    # Export model output
    save(FEve.mod, file = paste0("data/model_outputs_new/m_FEve_",
                                 "PC_comparison_pc", pc.count.1, "_vs_pc", pc.count.2, ".RData"))
    
    
  } else {
    
    # Load in model output
    FEve.mod <- get(load(paste0("data/model_outputs_new/m_FEve_",
                                "PC_comparison_pc", pc.count.1, "_vs_pc", pc.count.2, ".RData")))
    
  } # End of FEve if loop
  
  # Convert model output into a dataframe (with 4.d.p.)
  FEve.df <- brms_SummaryTable(FEve.mod, formatOptions = list(digits = 4, nsmall = 4), round = 4)
  
  # Extract the confidence intervals as a list for use in the plotting
  FEve.df.ci <- FEve.df %>% 
    filter(Covariate %in% c("FEve.1")) %>% 
    dplyr::select("l-95% CI", "u-95% CI") %>% 
    pivot_longer(cols = 1:2, names_to = "Interval", values_to = "Value")
  
  # Save the confidence intervals as a list and remove the intermediate dataframe
  FEve.ci <- as.list(FEve.df.ci$Value) # Adding values to the list
  names(FEve.ci) <- FEve.df.ci$Interval # Adding names to the values
  rm(FEve.df, FEve.df.ci) # Remove unnecessary dataframes of summary and CIs
  
  # Extract the prediction data frame
  FEve.pred <- ggpredict(FEve.mod, terms = "FEve.1 [sample = 30]",
                         # type = "random", condition = c(gridcell = 1266), allow_new_levels = TRUE,
                         back.transform = FALSE)
  
  # Extract dataframe with only the x-variable for plotting
  data.x.var.only <- dplyr::select(FD.tidy.complete, FEve.1)
  
  # Plot the modelled comparisons against one another
  (pc.comparison.FEve <- ggplot() +
      geom_point(data = FD.tidy.complete, aes(x = FEve.1, y = FEve.2, fill = colour.variable),
                 colour = "#000000", alpha = 0.5, shape = 21, size = 2) +
      
      geom_line(data = FEve.pred, aes(x = x, y = predicted), colour = "grey10") + # Adds line for predicted values of FEve vs LAT
      geom_ribbon(data = FEve.pred, aes(x = x, ymin = conf.low, ymax = conf.high),
                  fill = "grey10", alpha = 0.2) + # Adds c.intervals for predictions as ribbon
      scale_fill_viridis(option = FEve.colour, begin = 0.2, end = 1, direction = -1, discrete = TRUE) +
      labs(y = paste0("FEve PC", pc.count.2, "\n"),
           x = paste0("\n FEve PC", pc.count.1),
           title = paste0("CIs: ", FEve.ci$"l-95% CI", " to ", FEve.ci$"u-95% CI"),
           subtitle = paste(ifelse((FEve.ci$"l-95% CI" < 0 & FEve.ci$"u-95% CI" < 0) | # Automatically asigns significance in subtitle based on CIs
                                     (FEve.ci$"l-95% CI" > 0 & FEve.ci$"u-95% CI" > 0),
                                   "Significant: CIs do NOT span 0", "Not Significant: CIs span 0")),
           fill = "Dominant \n Functional \n Group") +
      theme_1() +
      theme(legend.position = "right",
            plot.subtitle = element_text(colour = ifelse((FEve.ci$"l-95% CI" < 0 & FEve.ci$"u-95% CI" < 0) | # Automatically colours subtitle based on significance from CIs
                                                           (FEve.ci$"l-95% CI" > 0 & FEve.ci$"u-95% CI" > 0),
                                                         paste("#006400"), paste("#8b0000")))))
  
  # Create a panel of all three plots
  combined.panel <- grid.arrange(pc.comparison.FRic, pc.comparison.FEve, ncol = 2)
  
  # Export panel
  ggsave(combined.panel, filename = paste0("figures/outputs_new/combined_",
                                           "FD_bayesian_comparison_pc", pc.count.1, "_vs_", pc.count.2, ".png"),
         width = 16, height = 7)

  } # End of function


# FUNCTION: LINEAR METRIC COMPARISON ----
linear.pc.count.comparison <- function(run.model, pc.count.1, pc.count.2,
                                       species.count, colour.by, colour.discrete){
  
  # Determine whether to run the model
  if (run.model == TRUE){
    
    # Import FD outputs for each PC count
    FD.output.1 <- read.csv(paste0("data/output_08_fd_output_combined_",
                                   "PCA_pc", pc.count.1, "_all_years.csv"))
    FD.output.2 <- read.csv(paste0("data/output_08_fd_output_combined_",
                                   "PCA_pc", pc.count.2, "_all_years.csv"))
    
    # Modify dataframes
    FD.tidy.1 <- FD.output.1 %>%
      dplyr::select(SiteSubsitePlotYear, FRic, FEve, SR) %>% # Leave SR to check same in both
      rename(FRic.1 = FRic, FEve.1 = FEve) %>% 
      filter(SR >= pc.count.1 + species.count)
    FD.tidy.2 <- FD.output.2 %>%
      dplyr::select(SiteSubsitePlotYear, FRic, FEve, SR, SurveyedArea,
                    gridcell, SUBSITE, PlotDominatingFG) %>% 
      rename(FRic.2 = FRic, FEve.2 = FEve) %>% 
      filter(SR >= pc.count.2 + species.count)
    
    # Join together the two dataframes
    FD.tidy.complete <- left_join(FD.tidy.2, FD.tidy.1,
                                  by = c("SR" = "SR", "SiteSubsitePlotYear" ="SiteSubsitePlotYear")) %>% 
      relocate(SiteSubsitePlotYear, FRic.1, FRic.2, FEve.1, FEve.2, .after = ) %>% 
      mutate(PlotDominatingFG = ifelse(PlotDominatingFG == "Evergreen_Shrub-Dominated",
                                       "E.Shrub", PlotDominatingFG),
             PlotDominatingFG = ifelse(PlotDominatingFG == "Deciduous_Shrub-Dominated",
                                       "D.Shrub", PlotDominatingFG),
             PlotDominatingFG = ifelse(PlotDominatingFG == "Graminoid-Dominated",
                                       "Graminoid", PlotDominatingFG),
             PlotDominatingFG = ifelse(PlotDominatingFG == "Forb-Dominated",
                                       "Forb", PlotDominatingFG))
    
    # Function for running the linear model
    FRic.mod <- lmer(FRic.2 ~ FRic.1 + (1 | gridcell / SUBSITE), data = FD.tidy.complete)
    
    # Get model summary as dataframe
    FRic.summary <- data.frame(summary(FRic.mod)$coefficients)
    
    # Determine estimate and standard error values
    FRic.estimate <- FRic.summary[2,1]
    FRic.stderror <- FRic.summary[2,2]
    
    # Extract the prediction data frame
    FRic.pred <- ggpredict(FRic.mod, terms = "FRic.1 [0:19, sample = 30]",
                           back.transform = FALSE)
    
    # Extract dataframe with only the x-variable for plotting
    data.x.var.only <- dplyr::select(FD.tidy.complete, FRic.1)
    
    # Modify data for colouring points
    FD.tidy.plot <- FD.tidy.complete %>% 
      rename(colour.variable = colour.by)
    
    # Plot the modelled comparisons against one another
    (pc.comparison.FRic <- ggplot() +
        geom_point(data = FD.tidy.plot, aes(x = FRic.1, y = FRic.2, fill = colour.variable),
                   colour = "#000000", alpha = 0.5, shape = 21, size = 2) +
        geom_line(data = FRic.pred, aes(x = x, y = predicted), colour = "grey10") + # Adds line for predicted values of FRic vs LAT
        geom_ribbon(data = FRic.pred, aes(x = x, ymin = predicted - std.error,
                                          ymax = predicted + std.error),
                    fill = "grey10", alpha = 0.2) + # Adds c.intervals for predictions as ribbon
        scale_fill_viridis(option = FRic.colour, begin = 0.2, end = 1, direction = -1, discrete = colour.discrete) +
        labs(y = paste0("FRic PC", pc.count.2, "\n"),
             x = paste0("\n FRic PC", pc.count.1),
             title = paste0("Estimate: ", round(FRic.estimate, digits = 3)),
             subtitle = paste0("Standard Error: ", round(FRic.stderror, digits = 3)),
             caption = paste0("Species Count = PCs + ", as.character(species.count)),
             fill = paste0(colour.by)) +
        theme_1() +
        theme(legend.position = "right"))
    
    # Function for running the linear model
    FEve.mod <- lmer(FEve.2 ~ FEve.1 + (1 | gridcell / SUBSITE), data = FD.tidy.complete)
    
    # Get model summary as dataframe
    FEve.summary <- data.frame(summary(FEve.mod)$coefficients)
    
    # Determine estimate and standard error values
    FEve.estimate <- FEve.summary[2,1]
    FEve.stderror <- FEve.summary[2,2]
    
    # Extract the prediction data frame
    FEve.pred <- ggpredict(FEve.mod, terms = "FEve.1 [0:1, sample = 30]",
                           back.transform = FALSE)
    
    # Extract dataframe with only the x-variable for plotting
    data.x.var.only <- dplyr::select(FD.tidy.complete, FEve.1)
    
    # Modify data for colouring points
    FD.tidy.plot <- FD.tidy.complete %>% 
      rename(colour.variable = colour.by)
    
    # Plot the modelled comparisons against one another
    (pc.comparison.FEve <- ggplot() +
        geom_point(data = FD.tidy.plot, aes(x = FEve.1, y = FEve.2, fill = colour.variable),
                   colour = "#000000", alpha = 0.5, shape = 21, size = 2) +
        geom_line(data = FEve.pred, aes(x = x, y = predicted), colour = "grey10") + # Adds line for predicted values of FEve vs LAT
        geom_ribbon(data = FEve.pred, aes(x = x, ymin = predicted - std.error,
                                          ymax = predicted + std.error),
                    fill = "grey10", alpha = 0.2) + # Adds c.intervals for predictions as ribbon
        scale_fill_viridis(option = FEve.colour, begin = 0.2, end = 1, direction = -1, discrete = colour.discrete) +
        labs(y = paste0("FEve PC", pc.count.2, "\n"),
             x = paste0("\n FEve PC", pc.count.1),
             title = paste0("Estimate: ", round(FEve.estimate, digits = 3)),
             subtitle = paste0("Standard Error: ", round(FEve.stderror, digits = 3)),
             caption = paste0("Species Count = PCs + ", as.character(species.count)),
             fill = paste0(colour.by)) +
        theme_1() +
        theme(legend.position = "right"))
    
    # Create a panel of all three plots
    combined.panel <- grid.arrange(pc.comparison.FRic, pc.comparison.FEve, ncol = 2)
    
    # Export panel
    ggsave(combined.panel, filename = paste0("figures/outputs_new/combined_",
                                             "FD_linear_comparison_pc", pc.count.1, "_vs_", pc.count.2, ".png"),
           width = 16, height = 7)
    
  } # End of if loop
  
} # End of function


# FIGURE 3b: GENERATE LATITUDE PLOTS WITH COLOUR SCHEME ----

plot.latitude <- function(censored){
  
  # Determine which version of FRic to import (censored or not)
  
  if (censored == "No"){
    
    # Load in model output
    FRic.mod <- get(load(paste0("data/model_outputs_new/m_FRic_",
                                FRic.distribution, "_LAT", pc.filepath, ".RData")))
    
    # Extract the prediction data frame
    FRic.pred <- ggpredict(FRic.mod, terms = "x_variable [60:85, sample = 30]", back.transform = TRUE)
    
    # Convert model output into a dataframe (with 4.d.p.)
    FRic.df <- brms_SummaryTable(FRic.mod, formatOptions = list(digits = 4, nsmall = 4), round = 4)
    
    # Extract the confidence intervals as a list for use in the plotting
    FRic.df.ci <- FRic.df %>% 
      filter(Covariate %in% c("x_variable")) %>% 
      dplyr::select("l-95% CI", "u-95% CI") %>% 
      pivot_longer(cols = 1:2, names_to = "Interval", values_to = "Value")
    
  } # End of import and not censored
  
  if (censored == "Yes"){
    
    # Load in model output
    FRic.mod <- get(load(paste0("data/model_outputs_new/m_FRic_",
                                "censored_LAT", pc.filepath, ".RData")))
    
    # Extract the prediction data frame and exponentiate the outputs
    FRic.pred <- ggpredict(FRic.mod, terms = "x_variable [60:85, sample = 30]", back.transform = FALSE) %>% 
      mutate(conf.low = exp(conf.low),
             conf.high = exp(conf.high),
             predicted = exp(predicted))
    
    # Convert model output into a dataframe (with 4.d.p.)
    FRic.df <- brms_SummaryTable(FRic.mod, formatOptions = list(digits = 4, nsmall = 4), round = 4)
    
    # Extract the confidence intervals as a list for use in the plotting
    FRic.df.ci <- FRic.df %>% 
      filter(Covariate %in% c("x_variable")) %>%
      dplyr::select("l-95% CI", "u-95% CI") %>% 
      pivot_longer(cols = 1:2, names_to = "Interval", values_to = "Value")
    
  } # End of import and censored
  
  # Save the confidence intervals as a list and remove the intermediate dataframe
  FRic.ci <- as.list(FRic.df.ci$Value) # Adding values to the list
  names(FRic.ci) <- FRic.df.ci$Interval # Adding names to the values
  rm(FRic.df, FRic.df.ci) # Remove unnecessary dataframes of summary and CIs
  
  # Plot model outputs
  (FRic.plot <- ggplot(FRic.pred) +
      
      # Scatterplot
      # geom_point(data = combo.latest, aes(x = LAT, y = FRic), # Adds original FRic vs LAT data points and colours by region
      #            fill = "#D02090", colour = c("#000000"), alpha = 0.2, shape = 21, size = 3) +
      # geom_ribbon(aes(x = x, ymin = conf.low, ymax = conf.high),
      #             fill = "#D02090", alpha = 0.27, colour = "#D02090", linewdith = 0.1) + # Adds c.intervals for predictions as ribbon
      # geom_line(aes(x = x, y = predicted), colour = "#D02090", linewidth = 1.1,
      #           linetype = ifelse((FRic.ci$"l-95% CI" < 0 & FRic.ci$"u-95% CI" < 0) | # Automatically asigns significance in subtitle based on CIs
                                    # (FRic.ci$"l-95% CI" > 0 & FRic.ci$"u-95% CI" > 0), 1, 2)) + # Adds line for predicted values of FRic vs LAT
      
      # Hexplot
      geom_hex(data = combo.latest, aes(x = LAT, y = FRic), bins = 18) +
            scale_fill_gradient(low = "#FFE1FF", high = "#D02090") +
      geom_ribbon(aes(x = x, ymin = conf.low, ymax = conf.high),
                  fill = "#D02090", alpha = 0.2, colour = "#D02090", linewdith = 0.05) + # Adds c.intervals for predictions as ribbon
      geom_line(aes(x = x, y = predicted), colour = "#000000", linewidth = 1.1,
                linetype = ifelse((FRic.ci$"l-95% CI" < 0 & FRic.ci$"u-95% CI" < 0) | # Automatically asigns significance in subtitle based on CIs
                                    (FRic.ci$"l-95% CI" > 0 & FRic.ci$"u-95% CI" > 0), 1, 2)) + # Adds line for predicted values of FRic vs LAT
      
      labs(y = "Functional Richness\n",
           x = paste0("\nLatitude (", "\u00B0", "N)"),
           title = paste0("CIs: ", FRic.ci$"l-95% CI", " to ", FRic.ci$"u-95% CI"),
           subtitle = paste(ifelse((FRic.ci$"l-95% CI" < 0 & FRic.ci$"u-95% CI" < 0) | # Automatically asigns significance in subtitle based on CIs
                                     (FRic.ci$"l-95% CI" > 0 & FRic.ci$"u-95% CI" > 0),
                                   "Significant: CIs do NOT span 0", "Not Significant: CIs span 0"))) +
      theme_1() +
      theme(legend.position = "none",
            plot.subtitle = element_text(colour = ifelse((FRic.ci$"l-95% CI" < 0 & FRic.ci$"u-95% CI" < 0) | # Automatically colours subtitle based on significance from CIs
                                                           (FRic.ci$"l-95% CI" > 0 & FRic.ci$"u-95% CI" > 0),
                                                         paste("#006400"), paste("#8b0000")))))
  
  
  # Load in model output
  FEve.mod <- get(load(paste0("data/model_outputs_new/m_FEve_",
                              "LAT", pc.filepath, ".RData")))
  
  # Convert model output into a dataframe (with 4.d.p.)
  FEve.df <- brms_SummaryTable(FEve.mod, formatOptions = list(digits = 4, nsmall = 4), round = 4)
  
  # Extract the confidence intervals as a list for use in the plotting
  FEve.df.ci <- FEve.df %>% 
    filter(Covariate %in% c("x_variable")) %>% 
    dplyr::select("l-95% CI", "u-95% CI") %>% 
    pivot_longer(cols = 1:2, names_to = "Interval", values_to = "Value")
  
  # Save the confidence intervals as a list and remove the intermediate dataframe
  FEve.ci <- as.list(FEve.df.ci$Value) # Adding values to the list
  names(FEve.ci) <- FEve.df.ci$Interval # Adding names to the values
  rm(FEve.df, FEve.df.ci) # Remove unnecessary dataframes of summary and CIs  
  
  # Extract the prediction data frame
  FEve.pred <- ggpredict(FEve.mod, terms = "x_variable [60:85, sample = 30]", back.transform = FALSE)
  
  # Generate plot
  (FEve.plot <- ggplot(FEve.pred) +
      
      
      # # SCatterplot
      # geom_point(data = combo.latest, aes(x = LAT, y = FEve), # Adds original FEve vs LAT data points and colours by region
      #            fill = "#228B22", colour = c("#000000"), alpha = 0.2, shape = 21, size = 3) +
      # geom_ribbon(aes(x = x, ymin = conf.low, ymax = conf.high),
      #             fill = "#228B22", alpha = 0.27, colour = "#228B22", linewdith = 0.1) + # Adds c.intervals for predictions as ribbon
      # geom_line(aes(x = x, y = predicted), colour = "#228B22", linewidth = 1.1,
      #           linetype = ifelse((FEve.ci$"l-95% CI" < 0 & FEve.ci$"u-95% CI" < 0) | # Automatically asigns significance in subtitle based on CIs
      #                               (FEve.ci$"l-95% CI" > 0 & FEve.ci$"u-95% CI" > 0), 1, 2)) + # Adds line for predicted values of FEve vs LAT

      # Hexplot
      geom_hex(data = combo.latest, aes(x = LAT, y = FEve), bins = 18) +
      scale_fill_gradient(low = "#CFF099", high = "#228B22") +
      geom_ribbon(aes(x = x, ymin = conf.low, ymax = conf.high),
                  fill = "#228B22", alpha = 0.25, colour = "#228B22", linewdith = 0.05) + # Adds c.intervals for predictions as ribbon
      geom_line(aes(x = x, y = predicted), colour = "#000000", linewidth = 1.1,
                linetype = ifelse((FEve.ci$"l-95% CI" < 0 & FEve.ci$"u-95% CI" < 0) | # Automatically asigns significance in subtitle based on CIs
                                    (FEve.ci$"l-95% CI" > 0 & FEve.ci$"u-95% CI" > 0), 1, 2)) + # Adds line for predicted values of FRic vs LAT
      
      labs(y = "Functional Evenness\n",
           x = paste0("\nLatitude (", "\u00B0", "N)"),
           title = paste0("CIs: ", FEve.ci$"l-95% CI", " to ", FEve.ci$"u-95% CI"),
           subtitle = paste(ifelse((FEve.ci$"l-95% CI" < 0 & FEve.ci$"u-95% CI" < 0) | # Automatically asigns significance in subtitle based on CIs
                                     (FEve.ci$"l-95% CI" > 0 & FEve.ci$"u-95% CI" > 0),
                                   "Significant: CIs do NOT span 0", "Not Significant: CIs span 0"))) +
      theme_1() +
      theme(legend.position = "none",
            plot.subtitle = element_text(colour = ifelse((FEve.ci$"l-95% CI" < 0 & FEve.ci$"u-95% CI" < 0) | # Automatically colours subtitle based on significance from CIs
                                                           (FEve.ci$"l-95% CI" > 0 & FEve.ci$"u-95% CI" > 0),
                                                         paste("#006400"), paste("#8b0000")))))
  
  
  # Load in model output
  SR.mod <- get(load(paste0("data/model_outputs_new/m_SR_",
                            "LAT", pc.filepath, ".RData")))
  
  # Convert model output into a dataframe (with 4.d.p.)
  SR.df <- brms_SummaryTable(SR.mod, formatOptions = list(digits = 4, nsmall = 4), round = 4)
  
  # Extract the confidence intervals as a list for use in the plotting
  SR.df.ci <- SR.df %>% 
    filter(Covariate %in% c("x_variable")) %>% 
    dplyr::select("l-95% CI", "u-95% CI") %>% 
    pivot_longer(cols = 1:2, names_to = "Interval", values_to = "Value")
  
  # Save the confidence intervals as a list and remove the intermediate dataframe
  SR.ci <- as.list(SR.df.ci$Value) # Adding values to the list
  names(SR.ci) <- SR.df.ci$Interval # Adding names to the values
  rm(SR.df, SR.df.ci) # Remove unnecessary dataframes of summary and CIs  
  
  # Extract the prediction data frame
  SR.pred <- ggpredict(SR.mod, terms = "x_variable [60:85, sample = 30]", back.transform = FALSE)
  
  # Generate plot
  (SR.plot <- ggplot(SR.pred) +
      
      # Scatterplot
      # geom_point(data = combo.latest, aes(x = LAT, y = SR), # Adds original SR vs LAT data points and colours by region
      #            fill = "#1C86EE", colour = c("#000000"), alpha = 0.2, shape = 21, size = 3) +
      # geom_ribbon(aes(x = x, ymin = conf.low, ymax = conf.high),
      #             fill = "#1C86EE", alpha = 0.27, colour = "#1C86EE", linewdith = 0.1) + # Adds c.intervals for predictions as ribbon
      # geom_line(aes(x = x, y = predicted), colour = "#1C86EE", linewidth = 1.1,
      #           linetype = ifelse((SR.ci$"l-95% CI" < 0 & SR.ci$"u-95% CI" < 0) | # Automatically asigns significance in subtitle based on CIs
      #                               (SR.ci$"l-95% CI" > 0 & SR.ci$"u-95% CI" > 0), 1, 2)) + # Adds line for predicted values of SR vs LAT
      
      # Hexplot
      geom_hex(data = combo.latest, aes(x = LAT, y = SR), bins = 18) +
      scale_fill_gradient(low = "#ABE0EB", high = "#1C86EE") +
      geom_ribbon(aes(x = x, ymin = conf.low, ymax = conf.high),
                  fill = "#1C86EE", alpha = 0.25, colour = "#1C86EE", linewdith = 0.05) + # Adds c.intervals for predictions as ribbon
      geom_line(aes(x = x, y = predicted), colour = "#000000", linewidth = 1.1,
                linetype = ifelse((SR.ci$"l-95% CI" < 0 & SR.ci$"u-95% CI" < 0) | # Automatically asigns significance in subtitle based on CIs
                                    (SR.ci$"l-95% CI" > 0 & SR.ci$"u-95% CI" > 0), 1, 2)) + # Adds line for predicted values of SR vs LAT
      
      labs(y = "Species Richness\n",
           x = paste0("\nLatitude (", "\u00B0", "N)"),
           title = paste0("CIs: ", SR.ci$"l-95% CI", " to ", SR.ci$"u-95% CI"),
           subtitle = paste(ifelse((SR.ci$"l-95% CI" < 0 & SR.ci$"u-95% CI" < 0) | # Automatically asigns significance in subtitle based on CIs
                                     (SR.ci$"l-95% CI" > 0 & SR.ci$"u-95% CI" > 0),
                                   "Significant: CIs do NOT span 0", "Not Significant: CIs span 0"))) +
      theme_1() +
      theme(legend.position = "none",
            plot.subtitle = element_text(colour = ifelse((SR.ci$"l-95% CI" < 0 & SR.ci$"u-95% CI" < 0) | # Automatically colours subtitle based on significance from CIs
                                                           (SR.ci$"l-95% CI" > 0 & SR.ci$"u-95% CI" > 0),
                                                         paste("#006400"), paste("#8b0000")))))
  
  
  # Load in model output
  FDis.mod <- get(load(paste0("data/model_outputs_new/m_FDis_",
                              "LAT", pc.filepath, ".RData")))
  
  # Convert model output into a dataframe (with 4.d.p.)
  FDis.df <- brms_SummaryTable(FDis.mod, formatOptions = list(digits = 4, nsmall = 4), round = 4)
  
  # Extract the confidence intervals as a list for use in the plotting
  FDis.df.ci <- FDis.df %>% 
    filter(Covariate %in% c("x_variable")) %>% 
    dplyr::select("l-95% CI", "u-95% CI") %>% 
    pivot_longer(cols = 1:2, names_to = "Interval", values_to = "Value")
  
  # Save the confidence intervals as a list and remove the intermediate dataframe
  FDis.ci <- as.list(FDis.df.ci$Value) # Adding values to the list
  names(FDis.ci) <- FDis.df.ci$Interval # Adding names to the values
  rm(FDis.df, FDis.df.ci) # Remove unnecessary dataframes of summary and CIs  
  
  # Extract the prediction data frame
  FDis.pred <- ggpredict(FDis.mod, terms = "x_variable [60:85, sample = 30]", back.transform = FALSE)
  
  # Generate plot
  (FDis.plot <- ggplot(FDis.pred) +
      
      
      # # SCatterplot
      # geom_point(data = combo.latest, aes(x = LAT, y = FDis), # Adds original FDis vs LAT data points and colours by region
      #            fill = "#228B22", colour = c("#000000"), alpha = 0.2, shape = 21, size = 3) +
      # geom_ribbon(aes(x = x, ymin = conf.low, ymax = conf.high),
      #             fill = "#228B22", alpha = 0.27, colour = "#228B22", linewdith = 0.1) + # Adds c.intervals for predictions as ribbon
      # geom_line(aes(x = x, y = predicted), colour = "#228B22", linewidth = 1.1,
      #           linetype = ifelse((FDis.ci$"l-95% CI" < 0 & FDis.ci$"u-95% CI" < 0) | # Automatically asigns significance in subtitle based on CIs
      #                               (FDis.ci$"l-95% CI" > 0 & FDis.ci$"u-95% CI" > 0), 1, 2)) + # Adds line for predicted values of FDis vs LAT
      
    # Hexplot
    geom_hex(data = combo.latest, aes(x = LAT, y = FDis), bins = 18) +
      scale_fill_gradient(low = "#FFCF91", high = "#EE7600") +
      geom_ribbon(aes(x = x, ymin = conf.low, ymax = conf.high),
                  fill = "#EE7600", alpha = 0.25, colour = "#EE7600", linewdith = 0.05) + # Adds c.intervals for predictions as ribbon
      geom_line(aes(x = x, y = predicted), colour = "#000000", linewidth = 1.1,
                linetype = ifelse((FDis.ci$"l-95% CI" < 0 & FDis.ci$"u-95% CI" < 0) | # Automatically asigns significance in subtitle based on CIs
                                    (FDis.ci$"l-95% CI" > 0 & FDis.ci$"u-95% CI" > 0), 1, 2)) + # Adds line for predicted values of FRic vs LAT
      
      labs(y = "Functional Dispersion\n",
           x = paste0("\nLatitude (", "\u00B0", "N)"),
           title = paste0("CIs: ", FDis.ci$"l-95% CI", " to ", FDis.ci$"u-95% CI"),
           subtitle = paste(ifelse((FDis.ci$"l-95% CI" < 0 & FDis.ci$"u-95% CI" < 0) | # Automatically asigns significance in subtitle based on CIs
                                     (FDis.ci$"l-95% CI" > 0 & FDis.ci$"u-95% CI" > 0),
                                   "Significant: CIs do NOT span 0", "Not Significant: CIs span 0"))) +
      theme_1() +
      theme(legend.position = "none",
            plot.subtitle = element_text(colour = ifelse((FDis.ci$"l-95% CI" < 0 & FDis.ci$"u-95% CI" < 0) | # Automatically colours subtitle based on significance from CIs
                                                           (FDis.ci$"l-95% CI" > 0 & FDis.ci$"u-95% CI" > 0),
                                                         paste("#006400"), paste("#8b0000")))))
  
  
  # Create a panel of all three plots
  combined.panel <- grid.arrange(FRic.plot, FEve.plot, FDis.plot, SR.plot, ncol = 4)
  
  # Add in if statements for saving the output
  if (censored == "No"){
    
    # Export panel
    ggsave(combined.panel, filename = paste0("figures/outputs_new/manuscript_",
                                             "LAT", pc.filepath, ".png"), width = 18, height = 5)
    
  }
  
  # Add in if statements for saving the output
  if (censored == "Yes"){
    
    # Export panel
    ggsave(combined.panel, filename = paste0("figures/outputs_new/manuscript_",
                                             "LAT_censored", pc.filepath, ".png"), width = 18, height = 5)
    
  }
  
}


# FIGURE 4a: GENERATE WARMEST QUARTER PLOTS WITH COLOUR SCHEME ----

plot.warmest.quarter <- function(censored){
  
  # Determine which version of FRic to import (censored or not)

  if (censored == "No"){
    
    # Load in model output
    FRic.mod <- get(load(paste0("data/model_outputs_new/m_FRic_",
                                FRic.distribution, "_TempAvSum", pc.filepath, ".RData")))
    
    # Extract the prediction data frame
    FRic.pred <- ggpredict(FRic.mod, terms = "x_variable [0.9:14.5, sample = 30]", back.transform = TRUE)
    
    # Convert model output into a dataframe (with 4.d.p.)
    FRic.df <- brms_SummaryTable(FRic.mod, formatOptions = list(digits = 4, nsmall = 4), round = 4)
    
    # Extract the confidence intervals as a list for use in the plotting
    FRic.df.ci <- FRic.df %>% 
      filter(Covariate %in% c("x_variable")) %>% 
      dplyr::select("l-95% CI", "u-95% CI") %>% 
      pivot_longer(cols = 1:2, names_to = "Interval", values_to = "Value")
    
  } # End of import and not censored
  
  if (censored == "Yes"){
    
    # Load in model output
    FRic.mod <- get(load(paste0("data/model_outputs_new/m_FRic_",
                                "censored_TempAvSum", pc.filepath, ".RData")))
    
    # Extract the prediction data frame and exponentiate the outputs
    FRic.pred <- ggpredict(FRic.mod, terms = "x_variable [0.9:14.5, sample = 30]", back.transform = FALSE) %>% 
      mutate(conf.low = exp(conf.low),
             conf.high = exp(conf.high),
             predicted = exp(predicted))
    
    # Convert model output into a dataframe (with 4.d.p.)
    FRic.df <- brms_SummaryTable(FRic.mod, formatOptions = list(digits = 4, nsmall = 4), round = 4)
    
    # Extract the confidence intervals as a list for use in the plotting
    FRic.df.ci <- FRic.df %>% 
      filter(Covariate %in% c("x_variable")) %>%
      dplyr::select("l-95% CI", "u-95% CI") %>% 
      pivot_longer(cols = 1:2, names_to = "Interval", values_to = "Value")
    
  } # End of import and censored
  
  # Save the confidence intervals as a list and remove the intermediate dataframe
  FRic.ci <- as.list(FRic.df.ci$Value) # Adding values to the list
  names(FRic.ci) <- FRic.df.ci$Interval # Adding names to the values
  rm(FRic.df, FRic.df.ci) # Remove unnecessary dataframes of summary and CIs
  
  # Plot model outputs
  (FRic.plot <- ggplot(FRic.pred) +
      
      # Scatterplot
      # geom_point(data = combo.latest, aes(x = TempAvSum, y = FRic), # Adds original FRic vs LAT data points and colours by region
      #            fill = "#D02090", colour = c("#000000"), alpha = 0.2, shape = 21, size = 3) +
      # geom_ribbon(aes(x = x, ymin = conf.low, ymax = conf.high),
      #             fill = "#D02090", alpha = 0.27, colour = "#D02090", linewdith = 0.1) + # Adds c.intervals for predictions as ribbon
      # geom_line(aes(x = x, y = predicted), colour = "#D02090", linewidth = 1.1,
      #           linetype = ifelse((FRic.ci$"l-95% CI" < 0 & FRic.ci$"u-95% CI" < 0) | # Automatically asigns significance in subtitle based on CIs
      #                               (FRic.ci$"l-95% CI" > 0 & FRic.ci$"u-95% CI" > 0), 1, 2)) + # Adds line for predicted values of FRic vs LAT
      
      # Hexplot
      geom_hex(data = combo.latest, aes(x = TempAvSum, y = FRic), bins = 18) +
      scale_fill_gradient(low = "#FFE1FF", high = "#D02090") +
      geom_ribbon(aes(x = x, ymin = conf.low, ymax = conf.high),
                  fill = "#D02090", alpha = 0.2, colour = "#D02090", linewdith = 0.05) + # Adds c.intervals for predictions as ribbon
      geom_line(aes(x = x, y = predicted), colour = "#000000", linewidth = 1.1,
                linetype = ifelse((FRic.ci$"l-95% CI" < 0 & FRic.ci$"u-95% CI" < 0) | # Automatically asigns significance in subtitle based on CIs
                                    (FRic.ci$"l-95% CI" > 0 & FRic.ci$"u-95% CI" > 0), 1, 2)) + # Adds line for predicted values of FRic vs LAT
      
      
      labs(y = "Functional Richness\n",
           x = paste0("\nTemperature (", "\u00B0", "C)"),
           title = paste0("CIs: ", FRic.ci$"l-95% CI", " to ", FRic.ci$"u-95% CI"),
           subtitle = paste(ifelse((FRic.ci$"l-95% CI" < 0 & FRic.ci$"u-95% CI" < 0) | # Automatically asigns significance in subtitle based on CIs
                                     (FRic.ci$"l-95% CI" > 0 & FRic.ci$"u-95% CI" > 0),
                                   "Significant: CIs do NOT span 0", "Not Significant: CIs span 0"))) +
      theme_1() +
      theme(legend.position = "none",
            plot.subtitle = element_text(colour = ifelse((FRic.ci$"l-95% CI" < 0 & FRic.ci$"u-95% CI" < 0) | # Automatically colours subtitle based on significance from CIs
                                                           (FRic.ci$"l-95% CI" > 0 & FRic.ci$"u-95% CI" > 0),
                                                         paste("#006400"), paste("#8b0000")))))
  
  
    # Load in model output
    FEve.mod <- get(load(paste0("data/model_outputs_new/m_FEve_",
                                "TempAvSum", pc.filepath, ".RData")))
    
    # Convert model output into a dataframe (with 4.d.p.)
    FEve.df <- brms_SummaryTable(FEve.mod, formatOptions = list(digits = 4, nsmall = 4), round = 4)
    
    # Extract the confidence intervals as a list for use in the plotting
    FEve.df.ci <- FEve.df %>% 
      filter(Covariate %in% c("x_variable")) %>% 
      dplyr::select("l-95% CI", "u-95% CI") %>% 
      pivot_longer(cols = 1:2, names_to = "Interval", values_to = "Value")
    
    # Save the confidence intervals as a list and remove the intermediate dataframe
    FEve.ci <- as.list(FEve.df.ci$Value) # Adding values to the list
    names(FEve.ci) <- FEve.df.ci$Interval # Adding names to the values
    rm(FEve.df, FEve.df.ci) # Remove unnecessary dataframes of summary and CIs  
    
    # Extract the prediction data frame
    FEve.pred <- ggpredict(FEve.mod, terms = "x_variable [0.9:14.5, sample = 30]", back.transform = FALSE)
    
    # Generate plot
    (FEve.plot <- ggplot(FEve.pred) +
        
        # Scatterplot
        # geom_point(data = combo.latest, aes(x = TempAvSum, y = FEve), # Adds original FEve vs LAT data points and colours by region
        #            fill = "#228B22", colour = c("#000000"), alpha = 0.2, shape = 21, size = 3) +
        # geom_ribbon(aes(x = x, ymin = conf.low, ymax = conf.high),
        #             fill = "#228B22", alpha = 0.27, colour = "#228B22", linewdith = 0.1) + # Adds c.intervals for predictions as ribbon
        # geom_line(aes(x = x, y = predicted), colour = "#228B22", linewidth = 1.1,
        #           linetype = ifelse((FEve.ci$"l-95% CI" < 0 & FEve.ci$"u-95% CI" < 0) | # Automatically asigns significance in subtitle based on CIs
        #                               (FEve.ci$"l-95% CI" > 0 & FEve.ci$"u-95% CI" > 0), 1, 2)) + # Adds line for predicted values of FEve vs LAT
        
        # Hexplot
        geom_hex(data = combo.latest, aes(x = TempAvSum, y = FEve), bins = 18) +
        scale_fill_gradient(low = "#CFF099", high = "#228B22") +
        geom_ribbon(aes(x = x, ymin = conf.low, ymax = conf.high),
                    fill = "#228B22", alpha = 0.25, colour = "#228B22", linewdith = 0.05) + # Adds c.intervals for predictions as ribbon
        geom_line(aes(x = x, y = predicted), colour = "#000000", linewidth = 1.1,
                  linetype = ifelse((FEve.ci$"l-95% CI" < 0 & FEve.ci$"u-95% CI" < 0) | # Automatically asigns significance in subtitle based on CIs
                                      (FEve.ci$"l-95% CI" > 0 & FEve.ci$"u-95% CI" > 0), 1, 2)) + # Adds line for predicted values of FRic vs LAT
        
        labs(y = "Functional Evenness\n",
             x = paste0("\nTemperature (", "\u00B0", "C)"),
             title = paste0("CIs: ", FEve.ci$"l-95% CI", " to ", FEve.ci$"u-95% CI"),
             subtitle = paste(ifelse((FEve.ci$"l-95% CI" < 0 & FEve.ci$"u-95% CI" < 0) | # Automatically asigns significance in subtitle based on CIs
                                       (FEve.ci$"l-95% CI" > 0 & FEve.ci$"u-95% CI" > 0),
                                     "Significant: CIs do NOT span 0", "Not Significant: CIs span 0"))) +
        theme_1() +
        theme(legend.position = "none",
              plot.subtitle = element_text(colour = ifelse((FEve.ci$"l-95% CI" < 0 & FEve.ci$"u-95% CI" < 0) | # Automatically colours subtitle based on significance from CIs
                                                             (FEve.ci$"l-95% CI" > 0 & FEve.ci$"u-95% CI" > 0),
                                                           paste("#006400"), paste("#8b0000")))))

    
    # Load in model output
    SR.mod <- get(load(paste0("data/model_outputs_new/m_SR_",
                              "TempAvSum", pc.filepath, ".RData")))
    
    # Convert model output into a dataframe (with 4.d.p.)
    SR.df <- brms_SummaryTable(SR.mod, formatOptions = list(digits = 4, nsmall = 4), round = 4)
    
    # Extract the confidence intervals as a list for use in the plotting
    SR.df.ci <- SR.df %>% 
      filter(Covariate %in% c("x_variable")) %>% 
      dplyr::select("l-95% CI", "u-95% CI") %>% 
      pivot_longer(cols = 1:2, names_to = "Interval", values_to = "Value")
    
    # Save the confidence intervals as a list and remove the intermediate dataframe
    SR.ci <- as.list(SR.df.ci$Value) # Adding values to the list
    names(SR.ci) <- SR.df.ci$Interval # Adding names to the values
    rm(SR.df, SR.df.ci) # Remove unnecessary dataframes of summary and CIs  
    
    # Extract the prediction data frame
    SR.pred <- ggpredict(SR.mod, terms = "x_variable [0.9:14.5, sample = 30]", back.transform = FALSE)
    
    # Generate plot
    (SR.plot <- ggplot(SR.pred) +
        
        # SCatterplot
        # geom_point(data = combo.latest, aes(x = TempAvSum, y = SR), # Adds original SR vs LAT data points and colours by region
        #            fill = "#1C86EE", colour = c("#000000"), alpha = 0.2, shape = 21, size = 3) +
        # geom_ribbon(aes(x = x, ymin = conf.low, ymax = conf.high),
        #             fill = "#1C86EE", alpha = 0.27, colour = "#1C86EE", linewdith = 0.1) + # Adds c.intervals for predictions as ribbon
        # geom_line(aes(x = x, y = predicted), colour = "#1C86EE", linewidth = 1.1,
        #           linetype = ifelse((SR.ci$"l-95% CI" < 0 & SR.ci$"u-95% CI" < 0) | # Automatically asigns significance in subtitle based on CIs
        #                               (SR.ci$"l-95% CI" > 0 & SR.ci$"u-95% CI" > 0), 1, 2)) + # Adds line for predicted values of SR vs LAT
        
        # Hexplot
        geom_hex(data = combo.latest, aes(x = TempAvSum, y = SR), bins = 18) +
        scale_fill_gradient(low = "#ABE0EB", high = "#1C86EE") +
        geom_ribbon(aes(x = x, ymin = conf.low, ymax = conf.high),
                    fill = "#1C86EE", alpha = 0.25, colour = "#1C86EE", linewdith = 0.05) + # Adds c.intervals for predictions as ribbon
        geom_line(aes(x = x, y = predicted), colour = "#000000", linewidth = 1.1,
                  linetype = ifelse((SR.ci$"l-95% CI" < 0 & SR.ci$"u-95% CI" < 0) | # Automatically asigns significance in subtitle based on CIs
                                      (SR.ci$"l-95% CI" > 0 & SR.ci$"u-95% CI" > 0), 1, 2)) + # Adds line for predicted values of SR vs LAT
        
        labs(y = "Species Richness\n",
             x = paste0("\nTemperature (", "\u00B0", "C)"),
             title = paste0("CIs: ", SR.ci$"l-95% CI", " to ", SR.ci$"u-95% CI"),
             subtitle = paste(ifelse((SR.ci$"l-95% CI" < 0 & SR.ci$"u-95% CI" < 0) | # Automatically asigns significance in subtitle based on CIs
                                       (SR.ci$"l-95% CI" > 0 & SR.ci$"u-95% CI" > 0),
                                     "Significant: CIs do NOT span 0", "Not Significant: CIs span 0"))) +
        theme_1() +
        theme(legend.position = "none",
              plot.subtitle = element_text(colour = ifelse((SR.ci$"l-95% CI" < 0 & SR.ci$"u-95% CI" < 0) | # Automatically colours subtitle based on significance from CIs
                                                             (SR.ci$"l-95% CI" > 0 & SR.ci$"u-95% CI" > 0),
                                                           paste("#006400"), paste("#8b0000")))))
    
    
    # Load in model output
    FDis.mod <- get(load(paste0("data/model_outputs_new/m_FDis_",
                                "TempAvSum", pc.filepath, ".RData")))
    
    # Convert model output into a dataframe (with 4.d.p.)
    FDis.df <- brms_SummaryTable(FDis.mod, formatOptions = list(digits = 4, nsmall = 4), round = 4)
    
    # Extract the confidence intervals as a list for use in the plotting
    FDis.df.ci <- FDis.df %>% 
      filter(Covariate %in% c("x_variable")) %>% 
      dplyr::select("l-95% CI", "u-95% CI") %>% 
      pivot_longer(cols = 1:2, names_to = "Interval", values_to = "Value")
    
    # Save the confidence intervals as a list and remove the intermediate dataframe
    FDis.ci <- as.list(FDis.df.ci$Value) # Adding values to the list
    names(FDis.ci) <- FDis.df.ci$Interval # Adding names to the values
    rm(FDis.df, FDis.df.ci) # Remove unnecessary dataframes of summary and CIs  
    
    # Extract the prediction data frame
    FDis.pred <- ggpredict(FDis.mod, terms = "x_variable [0.9:14.5, sample = 30]", back.transform = FALSE)
    
    # Generate plot
    (FDis.plot <- ggplot(FDis.pred) +
        
        # SCatterplot
        # geom_point(data = combo.latest, aes(x = TempAvSum, y = FDis), # Adds original FDis vs LAT data points and colours by region
        #            fill = "#1C86EE", colour = c("#000000"), alpha = 0.2, shape = 21, size = 3) +
        # geom_ribbon(aes(x = x, ymin = conf.low, ymax = conf.high),
        #             fill = "#1C86EE", alpha = 0.27, colour = "#1C86EE", linewdith = 0.1) + # Adds c.intervals for predictions as ribbon
        # geom_line(aes(x = x, y = predicted), colour = "#1C86EE", linewidth = 1.1,
        #           linetype = ifelse((FDis.ci$"l-95% CI" < 0 & FDis.ci$"u-95% CI" < 0) | # Automatically asigns significance in subtitle based on CIs
        #                               (FDis.ci$"l-95% CI" > 0 & FDis.ci$"u-95% CI" > 0), 1, 2)) + # Adds line for predicted values of FDis vs LAT
        
      # Hexplot
      geom_hex(data = combo.latest, aes(x = TempAvSum, y = FDis), bins = 18) +
        scale_fill_gradient(low = "#FFCF91", high = "#EE7600") +
        geom_ribbon(aes(x = x, ymin = conf.low, ymax = conf.high),
                    fill = "#EE7600", alpha = 0.25, colour = "#EE7600", linewdith = 0.05) + # Adds c.intervals for predictions as ribbon
        geom_line(aes(x = x, y = predicted), colour = "#000000", linewidth = 1.1,
                  linetype = ifelse((FDis.ci$"l-95% CI" < 0 & FDis.ci$"u-95% CI" < 0) | # Automatically asigns significance in subtitle based on CIs
                                      (FDis.ci$"l-95% CI" > 0 & FDis.ci$"u-95% CI" > 0), 1, 2)) + # Adds line for predicted values of FRic vs LAT
        
        labs(y = "Functional Dispersion\n",
             x = paste0("\nTemperature (", "\u00B0", "C)"),
             title = paste0("CIs: ", FDis.ci$"l-95% CI", " to ", FDis.ci$"u-95% CI"),
             subtitle = paste(ifelse((FDis.ci$"l-95% CI" < 0 & FDis.ci$"u-95% CI" < 0) | # Automatically asigns significance in subtitle based on CIs
                                       (FDis.ci$"l-95% CI" > 0 & FDis.ci$"u-95% CI" > 0),
                                     "Significant: CIs do NOT span 0", "Not Significant: CIs span 0"))) +
        theme_1() +
        theme(legend.position = "none",
              plot.subtitle = element_text(colour = ifelse((FDis.ci$"l-95% CI" < 0 & FDis.ci$"u-95% CI" < 0) | # Automatically colours subtitle based on significance from CIs
                                                             (FDis.ci$"l-95% CI" > 0 & FDis.ci$"u-95% CI" > 0),
                                                           paste("#006400"), paste("#8b0000")))))
  
  
  # Create a panel of all three plots
  combined.panel <- grid.arrange(FRic.plot, FEve.plot, FDis.plot, SR.plot, ncol = 4)
  
  # Add in if statements for saving the output
  if (censored == "No"){
    
    # Export panel
    ggsave(combined.panel, filename = paste0("figures/outputs_new/manuscript_",
                                             "TempAvSum", pc.filepath, ".png"), width = 18, height = 5)
    
  }
  
  # Add in if statements for saving the output
  if (censored == "Yes"){
    
    # Export panel
    ggsave(combined.panel, filename = paste0("figures/outputs_new/manuscript_",
                                             "TempAvSum_censored", pc.filepath, ".png"), width = 18, height = 5)
    
  }

}


# FIGURE 4b: GENERATE COMBINED PLOTS FOR COVER VS METRICS ----

plot.combined.cover <- function(censored.FRic, quadratic.FEve){
  
  # Load in model outputs for FRic and process the predictions
  if (censored.FRic == "No"){
    
    # Load in model outputs
    FRic.shrub.mod <- get(load(paste0("data/model_outputs_new/m_FRic_",
                                      FRic.distribution, "_ShrubCover", pc.filepath, ".RData")))
    FRic.forb.mod <- get(load(paste0("data/model_outputs_new/m_FRic_",
                                     FRic.distribution, "_ForbCover", pc.filepath, ".RData")))
    FRic.gram.mod <- get(load(paste0("data/model_outputs_new/m_FRic_",
                                     FRic.distribution, "_GraminoidCover", pc.filepath, ".RData")))
    
    # Extract the prediction data frame
    FRic.shrub.pred <- ggpredict(FRic.shrub.mod, terms = "x_variable [0:100, sample = 30]", back.transform = TRUE) %>% 
      mutate(Functional_Group = "Shrubs")
    FRic.forb.pred <- ggpredict(FRic.forb.mod, terms = "x_variable [0:100, sample = 30]", back.transform = TRUE) %>% 
      mutate(Functional_Group = "Forbs")
    FRic.gram.pred <- ggpredict(FRic.gram.mod, terms = "x_variable [0:100, sample = 30]", back.transform = TRUE) %>% 
      mutate(Functional_Group = "Graminoids")
    
  } else { # End of if not censored
    
    # Load in model output
    FRic.shrub.mod <- get(load(paste0("data/model_outputs_new/m_FRic_",
                                      "censored_ShrubCover", pc.filepath, ".RData")))
    FRic.forb.mod <- get(load(paste0("data/model_outputs_new/m_FRic_",
                                     "censored_ForbCover", pc.filepath, ".RData")))
    FRic.gram.mod <- get(load(paste0("data/model_outputs_new/m_FRic_",
                                     "censored_GraminoidCover", pc.filepath, ".RData")))
    
    # Extract the prediction data frame and exponentiate the outputs
    FRic.shrub.pred <- ggpredict(FRic.shrub.mod, terms = "x_variable [0:100, sample = 30]", back.transform = TRUE) %>% 
      mutate(predicted = exp(predicted), conf.low = exp(conf.low), conf.high = exp(conf.high)) %>% 
      mutate(Functional_Group = "Shrubs")
    FRic.forb.pred <- ggpredict(FRic.forb.mod, terms = "x_variable [0:100, sample = 30]", back.transform = TRUE) %>% 
      mutate(predicted = exp(predicted), conf.low = exp(conf.low), conf.high = exp(conf.high)) %>% 
      mutate(Functional_Group = "Forbs")
    FRic.gram.pred <- ggpredict(FRic.gram.mod, terms = "x_variable [0:100, sample = 30]", back.transform = TRUE) %>% 
      mutate(predicted = exp(predicted), conf.low = exp(conf.low), conf.high = exp(conf.high)) %>% 
      mutate(Functional_Group = "Graminoids")
    
  } # End of else is censored
  
  # Combine the predictions outputs into one dataframe and rename variables
  FRic.predictions <- rbind(FRic.shrub.pred, FRic.forb.pred, FRic.gram.pred) %>% 
    rename(Cover = x) %>% 
    dplyr::select(-group)
  
  # Convert model output into a dataframe (with 4.d.p.)
  FRic.shrub.df <- brms_SummaryTable(FRic.shrub.mod, formatOptions = list(digits = 4, nsmall = 4), round = 4)
  FRic.forb.df <- brms_SummaryTable(FRic.forb.mod, formatOptions = list(digits = 4, nsmall = 4), round = 4)
  FRic.gram.df <- brms_SummaryTable(FRic.gram.mod, formatOptions = list(digits = 4, nsmall = 4), round = 4)
  
  # Extract the confidence intervals as a list for use in the plotting
  FRic.shrub.df.ci <- FRic.shrub.df %>% 
    filter(Covariate %in% c("x_variable")) %>% 
    dplyr::select("l-95% CI", "u-95% CI") %>% 
    pivot_longer(cols = 1:2, names_to = "Interval", values_to = "Value")
  FRic.forb.df.ci <- FRic.forb.df %>% 
    filter(Covariate %in% c("x_variable")) %>% 
    dplyr::select("l-95% CI", "u-95% CI") %>% 
    pivot_longer(cols = 1:2, names_to = "Interval", values_to = "Value")
  FRic.gram.df.ci <- FRic.gram.df %>% 
    filter(Covariate %in% c("x_variable")) %>% 
    dplyr::select("l-95% CI", "u-95% CI") %>% 
    pivot_longer(cols = 1:2, names_to = "Interval", values_to = "Value")
  
  # Save the confidence intervals as a list and remove the intermediate dataframe
  FRic.shrub.ci <- as.list(FRic.shrub.df.ci$Value) # Adding values to the list
  names(FRic.shrub.ci) <- FRic.shrub.df.ci$Interval # Adding names to the values
  
  FRic.forb.ci <- as.list(FRic.forb.df.ci$Value) # Adding values to the list
  names(FRic.forb.ci) <- FRic.forb.df.ci$Interval # Adding names to the values
  
  FRic.gram.ci <- as.list(FRic.gram.df.ci$Value) # Adding values to the list
  names(FRic.gram.ci) <- FRic.gram.df.ci$Interval # Adding names to the values
  
  # Create a ggplot() item combining the three outputs
  (FRic.plot <- ggplot(data = FRic.predictions) +
      geom_ribbon(aes(x = Cover, ymin = conf.low, ymax = conf.high, fill = Functional_Group),
                  alpha = 0.2, linetype = 2) +
      scale_fill_manual(values = c("#D02090", "#228B22", "#1C86EE")) +
      
      geom_line(data = filter(FRic.predictions, Functional_Group == "Forbs"),
                aes(x = Cover, y = predicted), colour = "#D02090",
                linetype = ifelse((FRic.forb.ci$"l-95% CI" < 0 & FRic.forb.ci$"u-95% CI" < 0) | # Automatically asigns significance in subtitle based on CIs
                                    (FRic.forb.ci$"l-95% CI" > 0 & FRic.forb.ci$"u-95% CI" > 0), 1, 2)) +
      
      geom_line(data = filter(FRic.predictions, Functional_Group == "Shrubs"),
                aes(x = Cover, y = predicted), colour = "#1C86EE",
                linetype = ifelse((FRic.shrub.ci$"l-95% CI" < 0 & FRic.shrub.ci$"u-95% CI" < 0) | # Automatically asigns significance in subtitle based on CIs
                                    (FRic.shrub.ci$"l-95% CI" > 0 & FRic.shrub.ci$"u-95% CI" > 0), 1, 2)) +
      
      geom_line(data = filter(FRic.predictions, Functional_Group == "Graminoids"),
                aes(x = Cover, y = predicted), colour = "#228B22",
                linetype = ifelse((FRic.gram.ci$"l-95% CI" < 0 & FRic.gram.ci$"u-95% CI" < 0) | # Automatically asigns significance in subtitle based on CIs
                                    (FRic.gram.ci$"l-95% CI" > 0 & FRic.gram.ci$"u-95% CI" > 0), 1, 2)) +
      
      labs(y = "FRic \n",
           x = "\nCover (%)",
           fill = "Functional\nGroup",
           colour = "Functional\nGroup") +
      theme_1() +
      theme(legend.position = "none",
            legend.key = element_rect(colour = "#000000")))
  
  
  # Run if else criteria for importing quadratic vs non-quadratic outputs for cover vs FEve  
  if (quadratic.FEve == "No"){
    
    # Load in model outputs for FEve and process the predictions
    FEve.shrub.mod <- get(load(paste0("data/model_outputs_new/m_FEve_",
                                      "ShrubCover", pc.filepath, ".RData")))
    FEve.forb.mod <- get(load(paste0("data/model_outputs_new/m_FEve_",
                                     "ForbCover", pc.filepath, ".RData")))
    FEve.gram.mod <- get(load(paste0("data/model_outputs_new/m_FEve_",
                                     "GraminoidCover", pc.filepath, ".RData")))
    
    # Extract the prediction data frame
    FEve.shrub.pred <- ggpredict(FEve.shrub.mod, terms = "x_variable [0:100, sample = 30]", back.transform = FALSE) %>% 
      mutate(Functional_Group = "Shrubs")
    FEve.forb.pred <- ggpredict(FEve.forb.mod, terms = "x_variable [0:100, sample = 30]", back.transform = FALSE) %>% 
      mutate(Functional_Group = "Forbs")
    FEve.gram.pred <- ggpredict(FEve.gram.mod, terms = "x_variable [0:100, sample = 30]", back.transform = FALSE) %>% 
      mutate(Functional_Group = "Graminoids")
    
    # Combine the predictions outputs into one dataframe and rename variables
    FEve.predictions <- rbind(FEve.shrub.pred, FEve.forb.pred, FEve.gram.pred) %>% 
      rename(Cover = x) %>% 
      dplyr::select(-group)
    
    # Convert model output into a dataframe (with 4.d.p.)
    FEve.shrub.df <- brms_SummaryTable(FEve.shrub.mod, formatOptions = list(digits = 4, nsmall = 4), round = 4)
    FEve.forb.df <- brms_SummaryTable(FEve.forb.mod, formatOptions = list(digits = 4, nsmall = 4), round = 4)
    FEve.gram.df <- brms_SummaryTable(FEve.gram.mod, formatOptions = list(digits = 4, nsmall = 4), round = 4)
    
    # Extract the confidence intervals as a list for use in the plotting
    FEve.shrub.df.ci <- FEve.shrub.df %>% 
      filter(Covariate %in% c("x_variable")) %>% 
      dplyr::select("l-95% CI", "u-95% CI") %>% 
      pivot_longer(cols = 1:2, names_to = "Interval", values_to = "Value")
    FEve.forb.df.ci <- FEve.forb.df %>% 
      filter(Covariate %in% c("x_variable")) %>% 
      dplyr::select("l-95% CI", "u-95% CI") %>% 
      pivot_longer(cols = 1:2, names_to = "Interval", values_to = "Value")
    FEve.gram.df.ci <- FEve.gram.df %>% 
      filter(Covariate %in% c("x_variable")) %>% 
      dplyr::select("l-95% CI", "u-95% CI") %>% 
      pivot_longer(cols = 1:2, names_to = "Interval", values_to = "Value")
    
    # Save the confidence intervals as a list and remove the intermediate dataframe
    FEve.shrub.ci <- as.list(FEve.shrub.df.ci$Value) # Adding values to the list
    names(FEve.shrub.ci) <- FEve.shrub.df.ci$Interval # Adding names to the values
    
    FEve.forb.ci <- as.list(FEve.forb.df.ci$Value) # Adding values to the list
    names(FEve.forb.ci) <- FEve.forb.df.ci$Interval # Adding names to the values
    
    FEve.gram.ci <- as.list(FEve.gram.df.ci$Value) # Adding values to the list
    names(FEve.gram.ci) <- FEve.gram.df.ci$Interval # Adding names to the values
    
    # Create a ggplot() item combining the three outputs
    (FEve.plot <- ggplot(data = FEve.predictions) +
        geom_ribbon(aes(x = Cover, ymin = conf.low, ymax = conf.high, fill = Functional_Group),
                    alpha = 0.2, linetype = 2) +
        scale_fill_manual(values = c("#D02090", "#228B22", "#1C86EE")) +
        
        geom_line(data = filter(FEve.predictions, Functional_Group == "Forbs"),
                  aes(x = Cover, y = predicted), colour = "#D02090",
                  linetype = ifelse((FEve.forb.ci$"l-95% CI" < 0 & FEve.forb.ci$"u-95% CI" < 0) | # Automatically asigns significance in subtitle based on CIs
                                      (FEve.forb.ci$"l-95% CI" > 0 & FEve.forb.ci$"u-95% CI" > 0), 1, 2)) +
        
        geom_line(data = filter(FEve.predictions, Functional_Group == "Shrubs"),
                  aes(x = Cover, y = predicted), colour = "#1C86EE",
                  linetype = ifelse((FEve.shrub.ci$"l-95% CI" < 0 & FEve.shrub.ci$"u-95% CI" < 0) | # Automatically asigns significance in subtitle based on CIs
                                      (FEve.shrub.ci$"l-95% CI" > 0 & FEve.shrub.ci$"u-95% CI" > 0), 1, 2)) +
        
        geom_line(data = filter(FEve.predictions, Functional_Group == "Graminoids"),
                  aes(x = Cover, y = predicted), colour = "#228B22",
                  linetype = ifelse((FEve.gram.ci$"l-95% CI" < 0 & FEve.gram.ci$"u-95% CI" < 0) | # Automatically asigns significance in subtitle based on CIs
                                      (FEve.gram.ci$"l-95% CI" > 0 & FEve.gram.ci$"u-95% CI" > 0), 1, 2)) +
        
        labs(y = "FEve \n",
             x = "\nCover (%)",
             fill = "Functional\nGroup",
             colour = "Functional\nGroup") +
        theme_1() +
        theme(legend.position = "none",
              legend.key = element_rect(colour = "#000000")))
    
  } else { # End of non-quadratic if statement
    
    # Load in model outputs for FEve and process the predictions
    FEve.shrub.mod <- get(load(paste0("data/model_outputs_new/m_FEve_",
                                      "quadratic_ShrubCover", pc.filepath, ".RData")))
    FEve.forb.mod <- get(load(paste0("data/model_outputs_new/m_FEve_",
                                     "quadratic_ForbCover", pc.filepath, ".RData")))
    FEve.gram.mod <- get(load(paste0("data/model_outputs_new/m_FEve_",
                                     "quadratic_GraminoidCover", pc.filepath, ".RData")))
    
    # Determine mean value for centering
    mean.x_variable.shrub <- mean(combo.latest$ShrubCover)
    mean.x_variable.forb <- mean(combo.latest$ForbCover)
    mean.x_variable.gram <- mean(combo.latest$GraminoidCover)
    
    # Create template predictions dataframe
    FEve.shrub.pred.df = data.frame(centred_x_variable = seq(0 - mean.x_variable.shrub, 100 - mean.x_variable.shrub), SurveyedArea = 1)
    FEve.forb.pred.df = data.frame(centred_x_variable = seq(0 - mean.x_variable.forb, 100 - mean.x_variable.forb), SurveyedArea = 1)
    FEve.gram.pred.df = data.frame(centred_x_variable = seq(0 - mean.x_variable.gram, 100 - mean.x_variable.gram), SurveyedArea = 1)
    
    # Add predictions to template datframe
    FEve.shrub.pred = add_epred_draws(FEve.shrub.pred.df, FEve.shrub.mod, re_formula = NA) %>% 
      mutate(x_variable = centred_x_variable + mean.x_variable.shrub) %>% # Uncentre the x variable
      mutate(Functional_Group = "Shrubs")
    FEve.forb.pred = add_epred_draws(FEve.forb.pred.df, FEve.forb.mod, re_formula = NA) %>% 
      mutate(x_variable = centred_x_variable + mean.x_variable.forb) %>% # Uncentre the x variable
      mutate(Functional_Group = "Forbs")
    FEve.gram.pred = add_epred_draws(FEve.gram.pred.df, FEve.gram.mod, re_formula = NA) %>% 
      mutate(x_variable = centred_x_variable + mean.x_variable.gram) %>% # Uncentre the x variable
      mutate(Functional_Group = "Graminoids")
    
    # Combine the predictions outputs into one dataframe and rename variables
    FEve.predictions <- rbind(FEve.shrub.pred, FEve.forb.pred, FEve.gram.pred)
    
    # Convert model output into a dataframe (with 4.d.p.)
    FEve.shrub.df <- brms_SummaryTable(FEve.shrub.mod, formatOptions = list(digits = 4, nsmall = 4), round = 4)
    FEve.forb.df <- brms_SummaryTable(FEve.forb.mod, formatOptions = list(digits = 4, nsmall = 4), round = 4)
    FEve.gram.df <- brms_SummaryTable(FEve.gram.mod, formatOptions = list(digits = 4, nsmall = 4), round = 4)
    
    # Extract the confidence intervals as a list for use in the plotting
    FEve.shrub.df.ci <- FEve.shrub.df %>% 
      filter(Covariate %in% c("Icentred_x_variableE2")) %>% 
      dplyr::select("l-95% CI", "u-95% CI") %>% 
      pivot_longer(cols = 1:2, names_to = "Interval", values_to = "Value")
    FEve.forb.df.ci <- FEve.forb.df %>% 
      filter(Covariate %in% c("Icentred_x_variableE2")) %>% 
      dplyr::select("l-95% CI", "u-95% CI") %>% 
      pivot_longer(cols = 1:2, names_to = "Interval", values_to = "Value")
    FEve.gram.df.ci <- FEve.gram.df %>% 
      filter(Covariate %in% c("Icentred_x_variableE2")) %>% 
      dplyr::select("l-95% CI", "u-95% CI") %>% 
      pivot_longer(cols = 1:2, names_to = "Interval", values_to = "Value")
    
    # Save the confidence intervals as a list and remove the intermediate dataframe
    FEve.shrub.ci <- as.list(FEve.shrub.df.ci$Value) # Adding values to the list
    names(FEve.shrub.ci) <- FEve.shrub.df.ci$Interval # Adding names to the values
    
    FEve.forb.ci <- as.list(FEve.forb.df.ci$Value) # Adding values to the list
    names(FEve.forb.ci) <- FEve.forb.df.ci$Interval # Adding names to the values
    
    FEve.gram.ci <- as.list(FEve.gram.df.ci$Value) # Adding values to the list
    names(FEve.gram.ci) <- FEve.gram.df.ci$Interval # Adding names to the values
    
    # Create a ggplot() item combining the three outputs
    (FEve.plot <- ggplot(data = FEve.predictions, aes(x = x_variable, y = .epred)) +
        
        stat_lineribbon(data = filter(FEve.predictions, Functional_Group == "Forbs"),
                        aes(y = .epred), .width = 0.95, alpha = 0.2, colour = "#D02090", fill = "#D02090", linewidth = 0.1,
                        linetype = ifelse((FEve.forb.ci$"l-95% CI" < 0 & FEve.forb.ci$"u-95% CI" < 0) | # Automatically asigns significance in subtitle based on CIs
                                            (FEve.forb.ci$"l-95% CI" > 0 & FEve.forb.ci$"u-95% CI" > 0), 1, 2)) +
        
        stat_lineribbon(data = filter(FEve.predictions, Functional_Group == "Shrubs"),
                        aes(y = .epred), .width = 0.95, alpha = 0.2, colour = "#1C86EE", fill = "#1C86EE", linewidth = 0.1,
                        linetype = ifelse((FEve.shrub.ci$"l-95% CI" < 0 & FEve.shrub.ci$"u-95% CI" < 0) | # Automatically asigns significance in subtitle based on CIs
                                            (FEve.shrub.ci$"l-95% CI" > 0 & FEve.shrub.ci$"u-95% CI" > 0), 1, 2)) +
        
        stat_lineribbon(data = filter(FEve.predictions, Functional_Group == "Graminoids"),
                        aes(y = .epred), .width = 0.95, alpha = 0.2, colour = "#228B22", fill = "#228B22", linewidth = 0.1,
                        linetype = ifelse((FEve.gram.ci$"l-95% CI" < 0 & FEve.gram.ci$"u-95% CI" < 0) | # Automatically asigns significance in subtitle based on CIs
                                            (FEve.gram.ci$"l-95% CI" > 0 & FEve.gram.ci$"u-95% CI" > 0), 1, 2)) +
        
        labs(y = "FEve \n",
             x = "\nCover (%)",
             fill = "Functional\nGroup",
             colour = "Functional\nGroup") +
        theme_1() +
        theme(legend.position = "none",
              legend.key = element_rect(colour = "#000000")))
    
  } # End of quadratic else statement
  
  
  # Load in model outputs for SR and process the predictions
  SR.shrub.mod <- get(load(paste0("data/model_outputs_new/m_SR_",
                                  "ShrubCover", pc.filepath, ".RData")))
  SR.forb.mod <- get(load(paste0("data/model_outputs_new/m_SR_",
                                 "ForbCover", pc.filepath, ".RData")))
  SR.gram.mod <- get(load(paste0("data/model_outputs_new/m_SR_",
                                 "GraminoidCover", pc.filepath, ".RData")))
  
  # Extract the prediction data frame
  SR.shrub.pred <- ggpredict(SR.shrub.mod, terms = "x_variable [0:100, sample = 30]", back.transform = FALSE) %>% 
    mutate(Functional_Group = "Shrubs")
  SR.forb.pred <- ggpredict(SR.forb.mod, terms = "x_variable [0:100, sample = 30]", back.transform = FALSE) %>% 
    mutate(Functional_Group = "Forbs")
  SR.gram.pred <- ggpredict(SR.gram.mod, terms = "x_variable [0:100, sample = 30]", back.transform = FALSE) %>% 
    mutate(Functional_Group = "Graminoids")
  
  # Combine the predictions outputs into one dataframe and rename variables
  SR.predictions <- rbind(SR.shrub.pred, SR.forb.pred, SR.gram.pred) %>% 
    rename(Cover = x) %>% 
    dplyr::select(-group)
  
  # Convert model output into a dataframe (with 4.d.p.)
  SR.shrub.df <- brms_SummaryTable(SR.shrub.mod, formatOptions = list(digits = 4, nsmall = 4), round = 4)
  SR.forb.df <- brms_SummaryTable(SR.forb.mod, formatOptions = list(digits = 4, nsmall = 4), round = 4)
  SR.gram.df <- brms_SummaryTable(SR.gram.mod, formatOptions = list(digits = 4, nsmall = 4), round = 4)
  
  # Extract the confidence intervals as a list for use in the plotting
  SR.shrub.df.ci <- SR.shrub.df %>% 
    filter(Covariate %in% c("x_variable")) %>% 
    dplyr::select("l-95% CI", "u-95% CI") %>% 
    pivot_longer(cols = 1:2, names_to = "Interval", values_to = "Value")
  SR.forb.df.ci <- SR.forb.df %>% 
    filter(Covariate %in% c("x_variable")) %>% 
    dplyr::select("l-95% CI", "u-95% CI") %>% 
    pivot_longer(cols = 1:2, names_to = "Interval", values_to = "Value")
  SR.gram.df.ci <- SR.gram.df %>% 
    filter(Covariate %in% c("x_variable")) %>% 
    dplyr::select("l-95% CI", "u-95% CI") %>% 
    pivot_longer(cols = 1:2, names_to = "Interval", values_to = "Value")
  
  # Save the confidence intervals as a list and remove the intermediate dataframe
  SR.shrub.ci <- as.list(SR.shrub.df.ci$Value) # Adding values to the list
  names(SR.shrub.ci) <- SR.shrub.df.ci$Interval # Adding names to the values
  
  SR.forb.ci <- as.list(SR.forb.df.ci$Value) # Adding values to the list
  names(SR.forb.ci) <- SR.forb.df.ci$Interval # Adding names to the values
  
  SR.gram.ci <- as.list(SR.gram.df.ci$Value) # Adding values to the list
  names(SR.gram.ci) <- SR.gram.df.ci$Interval # Adding names to the values
  
  # Create a ggplot() item combining the three outputs
  (SR.plot <- ggplot(data = SR.predictions) +
      geom_ribbon(aes(x = Cover, ymin = conf.low, ymax = conf.high, fill = Functional_Group),
                  alpha = 0.2, linetype = 2) +
      scale_fill_manual(values = c("#D02090", "#228B22", "#1C86EE")) +
      
      geom_line(data = filter(SR.predictions, Functional_Group == "Forbs"),
                aes(x = Cover, y = predicted), colour = "#D02090",
                linetype = ifelse((SR.forb.ci$"l-95% CI" < 0 & SR.forb.ci$"u-95% CI" < 0) | # Automatically asigns significance in subtitle based on CIs
                                    (SR.forb.ci$"l-95% CI" > 0 & SR.forb.ci$"u-95% CI" > 0), 1, 2)) +
      
      geom_line(data = filter(SR.predictions, Functional_Group == "Shrubs"),
                aes(x = Cover, y = predicted), colour = "#1C86EE",
                linetype = ifelse((SR.shrub.ci$"l-95% CI" < 0 & SR.shrub.ci$"u-95% CI" < 0) | # Automatically asigns significance in subtitle based on CIs
                                    (SR.shrub.ci$"l-95% CI" > 0 & SR.shrub.ci$"u-95% CI" > 0), 1, 2)) +
      
      geom_line(data = filter(SR.predictions, Functional_Group == "Graminoids"),
                aes(x = Cover, y = predicted), colour = "#228B22",
                linetype = ifelse((SR.gram.ci$"l-95% CI" < 0 & SR.gram.ci$"u-95% CI" < 0) | # Automatically asigns significance in subtitle based on CIs
                                    (SR.gram.ci$"l-95% CI" > 0 & SR.gram.ci$"u-95% CI" > 0), 1, 2)) +
      
      labs(y = "SR \n",
           x = "\nCover (%)",
           fill = "Functional\nGroup",
           colour = "Functional\nGroup") +
      theme_1() +
      theme(legend.position = "none",
            legend.key = element_rect(colour = "#000000")))
  
  # Create a panel of all three plots
  combined.panel <- grid.arrange(FRic.plot, FEve.plot, SR.plot, ncol = 3, widths = c(1,1,1))
  
  
  # Add in if statements for saving the output
  if (censored.FRic == "No" & quadratic.FEve == "No"){
    
    # Export panel
    ggsave(combined.panel, filename = paste0("figures/outputs_new/manuscript_",
                                             "all_cover", pc.filepath, ".png"), width = 23, height = 7)
    
  }
  
  # Add in if statements for saving the output
  if (censored.FRic == "No" & quadratic.FEve == "Yes"){
    
    # Export panel
    ggsave(combined.panel, filename = paste0("figures/outputs_new/manuscript_",
                                             "all_cover_quadratic", pc.filepath, ".png"), width = 23, height = 7)
    
  }
  
  # Add in if statements for saving the output
  if (censored.FRic == "Yes" & quadratic.FEve == "No"){
    
    # Export panel
    ggsave(combined.panel, filename = paste0("figures/outputs_new/manuscript_",
                                             "all_cover_censored", pc.filepath, ".png"), width = 23, height = 7)
    
  }
  
  # Add in if statements for saving the output
  if (censored.FRic == "Yes" & quadratic.FEve == "Yes"){
    
    # Export panel
    ggsave(combined.panel, filename = paste0("figures/outputs_new/manuscript_",
                                             "all_cover_censored_quadratic", pc.filepath, ".png"), width = 23, height = 7)
    
  }
  
} # End of function


# FIGURE 4b (QUAD): GENERATE COMBINED QUADRATIC PLOTS FOR COVER VS METRICS ----

plot.combined.cover.quadratic <- function(censored){
  
  
  # FRic ----
  
  # Load in model outputs depending on whether censored or not
  if (censored == "Yes"){
    
    # Load in model outputs for FEve and process the predictions
    FRic.shrub.mod <- get(load(paste0("data/model_outputs_new/m_FRic_",
                                      "censored_quadratic_ShrubCover", pc.filepath, ".RData")))
    FRic.forb.mod <- get(load(paste0("data/model_outputs_new/m_FRic_",
                                     "censored_quadratic_ForbCover", pc.filepath, ".RData")))
    FRic.gram.mod <- get(load(paste0("data/model_outputs_new/m_FRic_",
                                     "censored_quadratic_GraminoidCover", pc.filepath, ".RData")))
    
  } else {
    
    # Load in model outputs for FEve and process the predictions
    FRic.shrub.mod <- get(load(paste0("data/model_outputs_new/m_FRic_",
                                      "quadratic_ShrubCover", pc.filepath, ".RData")))
    FRic.forb.mod <- get(load(paste0("data/model_outputs_new/m_FRic_",
                                     "quadratic_ForbCover", pc.filepath, ".RData")))
    FRic.gram.mod <- get(load(paste0("data/model_outputs_new/m_FRic_",
                                     "quadratic_GraminoidCover", pc.filepath, ".RData")))
    
  }
  
  # Determine mean value for centering
  mean.x_variable.shrub <- mean(combo.latest$ShrubCover)
  mean.x_variable.forb <- mean(combo.latest$ForbCover)
  mean.x_variable.gram <- mean(combo.latest$GraminoidCover)
  
  # Create template predictions dataframe
  FRic.shrub.pred.df = data.frame(centred_x_variable = seq(0 - mean.x_variable.shrub, 100 - mean.x_variable.shrub), SurveyedArea = 1)
  FRic.forb.pred.df = data.frame(centred_x_variable = seq(0 - mean.x_variable.forb, 100 - mean.x_variable.forb), SurveyedArea = 1)
  FRic.gram.pred.df = data.frame(centred_x_variable = seq(0 - mean.x_variable.gram, 100 - mean.x_variable.gram), SurveyedArea = 1)
  
  # Add predictions to template datframe
  FRic.shrub.pred = add_epred_draws(FRic.shrub.pred.df, FRic.shrub.mod, re_formula = NA) %>% 
    mutate(x_variable = centred_x_variable + mean.x_variable.shrub) %>% # Uncentre the x variable
    mutate(Functional_Group = "Shrubs")
  FRic.forb.pred = add_epred_draws(FRic.forb.pred.df, FRic.forb.mod, re_formula = NA) %>% 
    mutate(x_variable = centred_x_variable + mean.x_variable.forb) %>% # Uncentre the x variable
    mutate(Functional_Group = "Forbs")
  FRic.gram.pred = add_epred_draws(FRic.gram.pred.df, FRic.gram.mod, re_formula = NA) %>% 
    mutate(x_variable = centred_x_variable + mean.x_variable.gram) %>% # Uncentre the x variable
    mutate(Functional_Group = "Graminoids")
  
  # Combine the predictions outputs into one dataframe and rename variables
  FRic.predictions <- rbind(FRic.shrub.pred, FRic.forb.pred, FRic.gram.pred)
  
  # Convert model output into a dataframe (with 4.d.p.)
  FRic.shrub.df <- brms_SummaryTable(FRic.shrub.mod, formatOptions = list(digits = 5, nsmall = 4), round = 5)
  FRic.forb.df <- brms_SummaryTable(FRic.forb.mod, formatOptions = list(digits = 5, nsmall = 4), round = 5)
  FRic.gram.df <- brms_SummaryTable(FRic.gram.mod, formatOptions = list(digits = 5, nsmall = 4), round = 5)
  
  # Extract the confidence intervals as a list for use in the plotting
  FRic.shrub.df.ci <- FRic.shrub.df %>% 
    filter(Covariate %in% c("Icentred_x_variableE2")) %>% 
    dplyr::select("l-95% CI", "u-95% CI") %>% 
    pivot_longer(cols = 1:2, names_to = "Interval", values_to = "Value")
  FRic.forb.df.ci <- FRic.forb.df %>% 
    filter(Covariate %in% c("Icentred_x_variableE2")) %>% 
    dplyr::select("l-95% CI", "u-95% CI") %>% 
    pivot_longer(cols = 1:2, names_to = "Interval", values_to = "Value")
  FRic.gram.df.ci <- FRic.gram.df %>% 
    filter(Covariate %in% c("Icentred_x_variableE2")) %>% 
    dplyr::select("l-95% CI", "u-95% CI") %>% 
    pivot_longer(cols = 1:2, names_to = "Interval", values_to = "Value")
  
  # Save the confidence intervals as a list and remove the intermediate dataframe
  FRic.shrub.ci <- as.list(FRic.shrub.df.ci$Value) # Adding values to the list
  names(FRic.shrub.ci) <- FRic.shrub.df.ci$Interval # Adding names to the values
  
  FRic.forb.ci <- as.list(FRic.forb.df.ci$Value) # Adding values to the list
  names(FRic.forb.ci) <- FRic.forb.df.ci$Interval # Adding names to the values
  
  FRic.gram.ci <- as.list(FRic.gram.df.ci$Value) # Adding values to the list
  names(FRic.gram.ci) <- FRic.gram.df.ci$Interval # Adding names to the values
  
  # Create a ggplot() item combining the three outputs
  (FRic.plot <- ggplot(data = FRic.predictions, aes(x = x_variable, y = .epred)) +
      
      # Ribbons
      stat_lineribbon(data = filter(FRic.predictions, Functional_Group == "Forbs"),
                      aes(y = .epred), .width = 0.95, alpha = 0.22, colour = "#F5D7EA", fill = "#F5D7EA", linewidth = 0.6,
                      linetype = ifelse((FRic.forb.ci$"l-95% CI" < 0 & FRic.forb.ci$"u-95% CI" < 0) | # Automatically assigns significance in subtitle based on CIs
                                          (FRic.forb.ci$"l-95% CI" > 0 & FRic.forb.ci$"u-95% CI" > 0), 1, 2)) +
      
      stat_lineribbon(data = filter(FRic.predictions, Functional_Group == "Graminoids"),
                      aes(y = .epred), .width = 0.95, alpha = 0.22, colour = "#E685C3", fill = "#E685C3", linewidth = 0.6,
                      linetype = ifelse((FRic.shrub.ci$"l-95% CI" < 0 & FRic.shrub.ci$"u-95% CI" < 0) | # Automatically assigns significance in subtitle based on CIs
                                          (FRic.shrub.ci$"l-95% CI" > 0 & FRic.shrub.ci$"u-95% CI" > 0), 1, 2)) +
      
      stat_lineribbon(data = filter(FRic.predictions, Functional_Group == "Shrubs"),
                      aes(y = .epred), .width = 0.95, alpha = 0.22, colour = "#B50675", fill = "#B50675", linewidth = 0.6,
                      linetype = ifelse((FRic.gram.ci$"l-95% CI" < 0 & FRic.gram.ci$"u-95% CI" < 0) | # Automatically assigns significance in subtitle based on CIs
                                          (FRic.gram.ci$"l-95% CI" > 0 & FRic.gram.ci$"u-95% CI" > 0), 1, 2)) +
      
      # Lines
      stat_lineribbon(data = filter(FRic.predictions, Functional_Group == "Forbs"),
                      aes(y = .epred), .width = 0.95, colour = "#F5D7EA", fill = NA, linewidth = 0.6,
                      linetype = ifelse((FRic.forb.ci$"l-95% CI" < 0 & FRic.forb.ci$"u-95% CI" < 0) | # Automatically assigns significance in subtitle based on CIs
                                          (FRic.forb.ci$"l-95% CI" > 0 & FRic.forb.ci$"u-95% CI" > 0), 1, 2)) +
      
      stat_lineribbon(data = filter(FRic.predictions, Functional_Group == "Graminoids"),
                      aes(y = .epred), .width = 0.95, colour = "#E685C3", fill = NA, linewidth = 0.6,
                      linetype = ifelse((FRic.shrub.ci$"l-95% CI" < 0 & FRic.shrub.ci$"u-95% CI" < 0) | # Automatically assigns significance in subtitle based on CIs
                                          (FRic.shrub.ci$"l-95% CI" > 0 & FRic.shrub.ci$"u-95% CI" > 0), 1, 2)) +
      
      stat_lineribbon(data = filter(FRic.predictions, Functional_Group == "Shrubs"),
                      aes(y = .epred), .width = 0.95, colour = "#B50675", fill = NA, linewidth = 0.6,
                      linetype = ifelse((FRic.gram.ci$"l-95% CI" < 0 & FRic.gram.ci$"u-95% CI" < 0) | # Automatically assigns significance in subtitle based on CIs
                                          (FRic.gram.ci$"l-95% CI" > 0 & FRic.gram.ci$"u-95% CI" > 0), 1, 2)) +
      
      labs(y = "Functional Richness \n",
           x = "\nCover (%)",
           fill = "Functional\nGroup",
           colour = "Functional\nGroup") +
      theme_1() +
      theme(legend.position = "bottom",
            legend.key = element_rect(colour = "#000000")))
  
  
  # FEve ----
  
  # Load in model outputs for FEve and process the predictions
  FEve.shrub.mod <- get(load(paste0("data/model_outputs_new/m_FEve_",
                                    "quadratic_ShrubCover", pc.filepath, ".RData")))
  FEve.forb.mod <- get(load(paste0("data/model_outputs_new/m_FEve_",
                                   "quadratic_ForbCover", pc.filepath, ".RData")))
  FEve.gram.mod <- get(load(paste0("data/model_outputs_new/m_FEve_",
                                   "quadratic_GraminoidCover", pc.filepath, ".RData")))
  
  # Determine mean value for centering
  mean.x_variable.shrub <- mean(combo.latest$ShrubCover)
  mean.x_variable.forb <- mean(combo.latest$ForbCover)
  mean.x_variable.gram <- mean(combo.latest$GraminoidCover)
  
  # Create template predictions dataframe
  FEve.shrub.pred.df = data.frame(centred_x_variable = seq(0 - mean.x_variable.shrub, 100 - mean.x_variable.shrub), SurveyedArea = 1)
  FEve.forb.pred.df = data.frame(centred_x_variable = seq(0 - mean.x_variable.forb, 100 - mean.x_variable.forb), SurveyedArea = 1)
  FEve.gram.pred.df = data.frame(centred_x_variable = seq(0 - mean.x_variable.gram, 100 - mean.x_variable.gram), SurveyedArea = 1)
  
  # Add predictions to template datframe
  FEve.shrub.pred = add_epred_draws(FEve.shrub.pred.df, FEve.shrub.mod, re_formula = NA) %>% 
    mutate(x_variable = centred_x_variable + mean.x_variable.shrub) %>% # Uncentre the x variable
    mutate(Functional_Group = "Shrubs")
  FEve.forb.pred = add_epred_draws(FEve.forb.pred.df, FEve.forb.mod, re_formula = NA) %>% 
    mutate(x_variable = centred_x_variable + mean.x_variable.forb) %>% # Uncentre the x variable
    mutate(Functional_Group = "Forbs")
  FEve.gram.pred = add_epred_draws(FEve.gram.pred.df, FEve.gram.mod, re_formula = NA) %>% 
    mutate(x_variable = centred_x_variable + mean.x_variable.gram) %>% # Uncentre the x variable
    mutate(Functional_Group = "Graminoids")
  
  # Combine the predictions outputs into one dataframe and rename variables
  FEve.predictions <- rbind(FEve.shrub.pred, FEve.forb.pred, FEve.gram.pred)
  
  # Convert model output into a dataframe (with 4.d.p.)
  FEve.shrub.df <- brms_SummaryTable(FEve.shrub.mod, formatOptions = list(digits = 5, nsmall = 4), round = 5)
  FEve.forb.df <- brms_SummaryTable(FEve.forb.mod, formatOptions = list(digits = 5, nsmall = 4), round = 5)
  FEve.gram.df <- brms_SummaryTable(FEve.gram.mod, formatOptions = list(digits = 5, nsmall = 4), round = 5)
  
  # Extract the confidence intervals as a list for use in the plotting
  FEve.shrub.df.ci <- FEve.shrub.df %>% 
    filter(Covariate %in% c("Icentred_x_variableE2")) %>% 
    dplyr::select("l-95% CI", "u-95% CI") %>% 
    pivot_longer(cols = 1:2, names_to = "Interval", values_to = "Value")
  FEve.forb.df.ci <- FEve.forb.df %>% 
    filter(Covariate %in% c("Icentred_x_variableE2")) %>% 
    dplyr::select("l-95% CI", "u-95% CI") %>% 
    pivot_longer(cols = 1:2, names_to = "Interval", values_to = "Value")
  FEve.gram.df.ci <- FEve.gram.df %>% 
    filter(Covariate %in% c("Icentred_x_variableE2")) %>% 
    dplyr::select("l-95% CI", "u-95% CI") %>% 
    pivot_longer(cols = 1:2, names_to = "Interval", values_to = "Value")
  
  # Save the confidence intervals as a list and remove the intermediate dataframe
  FEve.shrub.ci <- as.list(FEve.shrub.df.ci$Value) # Adding values to the list
  names(FEve.shrub.ci) <- FEve.shrub.df.ci$Interval # Adding names to the values
  
  FEve.forb.ci <- as.list(FEve.forb.df.ci$Value) # Adding values to the list
  names(FEve.forb.ci) <- FEve.forb.df.ci$Interval # Adding names to the values
  
  FEve.gram.ci <- as.list(FEve.gram.df.ci$Value) # Adding values to the list
  names(FEve.gram.ci) <- FEve.gram.df.ci$Interval # Adding names to the values
  
  # Create a ggplot() item combining the three outputs
  (FEve.plot <- ggplot(data = FEve.predictions, aes(x = x_variable, y = .epred)) +
      
      # Ribbons
      stat_lineribbon(data = filter(FEve.predictions, Functional_Group == "Forbs"),
                      aes(y = .epred), .width = 0.95, alpha = 0.22, colour = "#D1F0D1", fill = "#D1F0D1", linewidth = 0.6,
                      linetype = ifelse((FEve.forb.ci$"l-95% CI" < 0 & FEve.forb.ci$"u-95% CI" < 0) | # Automatically assigns significance in subtitle based on CIs
                                          (FEve.forb.ci$"l-95% CI" > 0 & FEve.forb.ci$"u-95% CI" > 0), 1, 2)) +
      
      stat_lineribbon(data = filter(FEve.predictions, Functional_Group == "Graminoids"),
                      aes(y = .epred), .width = 0.95, alpha = 0.22, colour = "#67B867", fill = "#67B867", linewidth = 0.6,
                      linetype = ifelse((FEve.shrub.ci$"l-95% CI" < 0 & FEve.shrub.ci$"u-95% CI" < 0) | # Automatically assigns significance in subtitle based on CIs
                                          (FEve.shrub.ci$"l-95% CI" > 0 & FEve.shrub.ci$"u-95% CI" > 0), 1, 2)) +
      
      stat_lineribbon(data = filter(FEve.predictions, Functional_Group == "Shrubs"),
                      aes(y = .epred), .width = 0.95, alpha = 0.22, colour = "#0A6E0A", fill = "#0A6E0A", linewidth = 0.6,
                      linetype = ifelse((FEve.gram.ci$"l-95% CI" < 0 & FEve.gram.ci$"u-95% CI" < 0) | # Automatically assigns significance in subtitle based on CIs
                                          (FEve.gram.ci$"l-95% CI" > 0 & FEve.gram.ci$"u-95% CI" > 0), 1, 2)) +
      
      # Lines
      stat_lineribbon(data = filter(FEve.predictions, Functional_Group == "Forbs"),
                      aes(y = .epred), .width = 0.95, colour = "#D1F0D1", fill = NA, linewidth = 0.6,
                      linetype = ifelse((FEve.forb.ci$"l-95% CI" < 0 & FEve.forb.ci$"u-95% CI" < 0) | # Automatically assigns significance in subtitle based on CIs
                                          (FEve.forb.ci$"l-95% CI" > 0 & FEve.forb.ci$"u-95% CI" > 0), 1, 2)) +
      
      stat_lineribbon(data = filter(FEve.predictions, Functional_Group == "Graminoids"),
                      aes(y = .epred), .width = 0.95, colour = "#67B867", fill = NA, linewidth = 0.6,
                      linetype = ifelse((FEve.shrub.ci$"l-95% CI" < 0 & FEve.shrub.ci$"u-95% CI" < 0) | # Automatically assigns significance in subtitle based on CIs
                                          (FEve.shrub.ci$"l-95% CI" > 0 & FEve.shrub.ci$"u-95% CI" > 0), 1, 2)) +
      
      stat_lineribbon(data = filter(FEve.predictions, Functional_Group == "Shrubs"),
                      aes(y = .epred), .width = 0.95, colour = "#0A6E0A", fill = NA, linewidth = 0.6,
                      linetype = ifelse((FEve.gram.ci$"l-95% CI" < 0 & FEve.gram.ci$"u-95% CI" < 0) | # Automatically assigns significance in subtitle based on CIs
                                          (FEve.gram.ci$"l-95% CI" > 0 & FEve.gram.ci$"u-95% CI" > 0), 1, 2)) +
      
      labs(y = "Functional Evenness \n",
           x = "\nCover (%)",
           fill = "Functional\nGroup",
           colour = "Functional\nGroup") +
      theme_1() +
      theme(legend.position = "bottom",
            legend.key = element_rect(colour = "#000000")))
  
  ggsave(FEve.plot, filename = paste0("figures/outputs_new/manuscript_",
                                           "FEve_quadratic_manual", pc.filepath, ".png"), width = 4.5, height = 4.6)
  
  
  # SR ----
  
  # Load in model outputs for SR and process the predictions
  SR.shrub.mod <- get(load(paste0("data/model_outputs_new/m_SR_",
                                  "quadratic_ShrubCover", pc.filepath, ".RData")))
  SR.forb.mod <- get(load(paste0("data/model_outputs_new/m_SR_",
                                 "quadratic_ForbCover", pc.filepath, ".RData")))
  SR.gram.mod <- get(load(paste0("data/model_outputs_new/m_SR_",
                                 "quadratic_GraminoidCover", pc.filepath, ".RData")))
  
  # Determine mean value for centering
  mean.x_variable.shrub <- mean(combo.latest$ShrubCover)
  mean.x_variable.forb <- mean(combo.latest$ForbCover)
  mean.x_variable.gram <- mean(combo.latest$GraminoidCover)
  
  # Create template predictions dataframe
  SR.shrub.pred.df = data.frame(centred_x_variable = seq(0 - mean.x_variable.shrub, 100 - mean.x_variable.shrub), SurveyedArea = 1)
  SR.forb.pred.df = data.frame(centred_x_variable = seq(0 - mean.x_variable.forb, 100 - mean.x_variable.forb), SurveyedArea = 1)
  SR.gram.pred.df = data.frame(centred_x_variable = seq(0 - mean.x_variable.gram, 100 - mean.x_variable.gram), SurveyedArea = 1)
  
  # Add predictions to template datframe
  SR.shrub.pred = add_epred_draws(SR.shrub.pred.df, SR.shrub.mod, re_formula = NA) %>% 
    mutate(x_variable = centred_x_variable + mean.x_variable.shrub) %>% # Uncentre the x variable
    mutate(Functional_Group = "Shrubs")
  SR.forb.pred = add_epred_draws(SR.forb.pred.df, SR.forb.mod, re_formula = NA) %>% 
    mutate(x_variable = centred_x_variable + mean.x_variable.forb) %>% # Uncentre the x variable
    mutate(Functional_Group = "Forbs")
  SR.gram.pred = add_epred_draws(SR.gram.pred.df, SR.gram.mod, re_formula = NA) %>% 
    mutate(x_variable = centred_x_variable + mean.x_variable.gram) %>% # Uncentre the x variable
    mutate(Functional_Group = "Graminoids")
  
  # Combine the predictions outputs into one dataframe and rename variables
  SR.predictions <- rbind(SR.shrub.pred, SR.forb.pred, SR.gram.pred)
  
  # Convert model output into a dataframe (with 4.d.p.)
  SR.shrub.df <- brms_SummaryTable(SR.shrub.mod, formatOptions = list(digits = 5, nsmall = 4), round = 5)
  SR.forb.df <- brms_SummaryTable(SR.forb.mod, formatOptions = list(digits = 5, nsmall = 4), round = 5)
  SR.gram.df <- brms_SummaryTable(SR.gram.mod, formatOptions = list(digits = 5, nsmall = 4), round = 5)
  
  # Extract the confidence intervals as a list for use in the plotting
  SR.shrub.df.ci <- SR.shrub.df %>% 
    filter(Covariate %in% c("Icentred_x_variableE2")) %>% 
    dplyr::select("l-95% CI", "u-95% CI") %>% 
    pivot_longer(cols = 1:2, names_to = "Interval", values_to = "Value")
  SR.forb.df.ci <- SR.forb.df %>% 
    filter(Covariate %in% c("Icentred_x_variableE2")) %>% 
    dplyr::select("l-95% CI", "u-95% CI") %>% 
    pivot_longer(cols = 1:2, names_to = "Interval", values_to = "Value")
  SR.gram.df.ci <- SR.gram.df %>% 
    filter(Covariate %in% c("Icentred_x_variableE2")) %>% 
    dplyr::select("l-95% CI", "u-95% CI") %>% 
    pivot_longer(cols = 1:2, names_to = "Interval", values_to = "Value")
  
  # Save the confidence intervals as a list and remove the intermediate dataframe
  SR.shrub.ci <- as.list(SR.shrub.df.ci$Value) # Adding values to the list
  names(SR.shrub.ci) <- SR.shrub.df.ci$Interval # Adding names to the values
  
  SR.forb.ci <- as.list(SR.forb.df.ci$Value) # Adding values to the list
  names(SR.forb.ci) <- SR.forb.df.ci$Interval # Adding names to the values
  
  SR.gram.ci <- as.list(SR.gram.df.ci$Value) # Adding values to the list
  names(SR.gram.ci) <- SR.gram.df.ci$Interval # Adding names to the values
  
  # Create a ggplot() item combining the three outputs
  (SR.plot <- ggplot(data = SR.predictions, aes(x = x_variable, y = .epred)) +
      
      # Ribbons
      stat_lineribbon(data = filter(SR.predictions, Functional_Group == "Forbs"),
                      aes(y = .epred), .width = 0.95, alpha = 0.22, colour = "#B1CCE8", fill = "#B1CCE8", linewidth = 0.6,
                      linetype = ifelse((SR.forb.ci$"l-95% CI" < 0 & SR.forb.ci$"u-95% CI" < 0) | # Automatically assigns significance in subtitle based on CIs
                                          (SR.forb.ci$"l-95% CI" > 0 & SR.forb.ci$"u-95% CI" > 0), 1, 2)) +
      
      stat_lineribbon(data = filter(SR.predictions, Functional_Group == "Graminoids"),
                      aes(y = .epred), .width = 0.95, alpha = 0.22, colour = "#1C86EE", fill = "#1C86EE", linewidth = 0.6,
                      linetype = ifelse((SR.shrub.ci$"l-95% CI" < 0 & SR.shrub.ci$"u-95% CI" < 0) | # Automatically assigns significance in subtitle based on CIs
                                          (SR.shrub.ci$"l-95% CI" > 0 & SR.shrub.ci$"u-95% CI" > 0), 1, 2)) +
      
      stat_lineribbon(data = filter(SR.predictions, Functional_Group == "Shrubs"),
                      aes(y = .epred), .width = 0.95, alpha = 0.22, colour = "#0D3E6E", fill = "#0D3E6E", linewidth = 0.6,
                      linetype = ifelse((SR.gram.ci$"l-95% CI" < 0 & SR.gram.ci$"u-95% CI" < 0) | # Automatically assigns significance in subtitle based on CIs
                                          (SR.gram.ci$"l-95% CI" > 0 & SR.gram.ci$"u-95% CI" > 0), 1, 2)) +
      
      # Lines
      stat_lineribbon(data = filter(SR.predictions, Functional_Group == "Forbs"),
                      aes(y = .epred), .width = 0.95, colour = "#B1CCE8", fill = NA, linewidth = 0.6,
                      linetype = ifelse((SR.forb.ci$"l-95% CI" < 0 & SR.forb.ci$"u-95% CI" < 0) | # Automatically assigns significance in subtitle based on CIs
                                          (SR.forb.ci$"l-95% CI" > 0 & SR.forb.ci$"u-95% CI" > 0), 1, 2)) +
      
      stat_lineribbon(data = filter(SR.predictions, Functional_Group == "Graminoids"),
                      aes(y = .epred), .width = 0.95, colour = "#1C86EE", fill = NA, linewidth = 0.6,
                      linetype = ifelse((SR.shrub.ci$"l-95% CI" < 0 & SR.shrub.ci$"u-95% CI" < 0) | # Automatically assigns significance in subtitle based on CIs
                                          (SR.shrub.ci$"l-95% CI" > 0 & SR.shrub.ci$"u-95% CI" > 0), 1, 2)) +
      
      stat_lineribbon(data = filter(SR.predictions, Functional_Group == "Shrubs"),
                      aes(y = .epred), .width = 0.95, colour = "#0D3E6E", fill = NA, linewidth = 0.6,
                      linetype = ifelse((SR.gram.ci$"l-95% CI" < 0 & SR.gram.ci$"u-95% CI" < 0) | # Automatically assigns significance in subtitle based on CIs
                                          (SR.gram.ci$"l-95% CI" > 0 & SR.gram.ci$"u-95% CI" > 0), 1, 2)) +
      
      labs(y = "Species Richness \n",
           x = "\nCover (%)",
           fill = "Functional\nGroup",
           colour = "Functional\nGroup") +
      theme_1() +
      theme(legend.position = "none",
            legend.key = element_rect(colour = "#000000")))
  
  
  # FDis ----
  
  # Load in model outputs for FDis and process the predictions
  FDis.shrub.mod <- get(load(paste0("data/model_outputs_new/m_FDis_",
                                    "quadratic_ShrubCover", pc.filepath, ".RData")))
  FDis.forb.mod <- get(load(paste0("data/model_outputs_new/m_FDis_",
                                   "quadratic_ForbCover", pc.filepath, ".RData")))
  FDis.gram.mod <- get(load(paste0("data/model_outputs_new/m_FDis_",
                                   "quadratic_GraminoidCover", pc.filepath, ".RData")))
  
  # Determine mean value for centering
  mean.x_variable.shrub <- mean(combo.latest$ShrubCover)
  mean.x_variable.forb <- mean(combo.latest$ForbCover)
  mean.x_variable.gram <- mean(combo.latest$GraminoidCover)
  
  # Create template predictions dataframe
  FDis.shrub.pred.df = data.frame(centred_x_variable = seq(0 - mean.x_variable.shrub, 100 - mean.x_variable.shrub), SurveyedArea = 1)
  FDis.forb.pred.df = data.frame(centred_x_variable = seq(0 - mean.x_variable.forb, 100 - mean.x_variable.forb), SurveyedArea = 1)
  FDis.gram.pred.df = data.frame(centred_x_variable = seq(0 - mean.x_variable.gram, 100 - mean.x_variable.gram), SurveyedArea = 1)
  
  # Add predictions to template datframe
  FDis.shrub.pred = add_epred_draws(FDis.shrub.pred.df, FDis.shrub.mod, re_formula = NA) %>% 
    mutate(x_variable = centred_x_variable + mean.x_variable.shrub) %>% # Uncentre the x variable
    mutate(Functional_Group = "Shrubs")
  FDis.forb.pred = add_epred_draws(FDis.forb.pred.df, FDis.forb.mod, re_formula = NA) %>% 
    mutate(x_variable = centred_x_variable + mean.x_variable.forb) %>% # Uncentre the x variable
    mutate(Functional_Group = "Forbs")
  FDis.gram.pred = add_epred_draws(FDis.gram.pred.df, FDis.gram.mod, re_formula = NA) %>% 
    mutate(x_variable = centred_x_variable + mean.x_variable.gram) %>% # Uncentre the x variable
    mutate(Functional_Group = "Graminoids")
  
  # Combine the predictions outputs into one dataframe and rename variables
  FDis.predictions <- rbind(FDis.shrub.pred, FDis.forb.pred, FDis.gram.pred)
  
  # Convert model output into a dataframe (with 4.d.p.)
  FDis.shrub.df <- brms_SummaryTable(FDis.shrub.mod, formatOptions = list(digits = 5, nsmall = 4), round = 5)
  FDis.forb.df <- brms_SummaryTable(FDis.forb.mod, formatOptions = list(digits = 5, nsmall = 4), round = 5)
  FDis.gram.df <- brms_SummaryTable(FDis.gram.mod, formatOptions = list(digits = 5, nsmall = 4), round = 5)
  
  # Extract the confidence intervals as a list for use in the plotting
  FDis.shrub.df.ci <- FDis.shrub.df %>% 
    filter(Covariate %in% c("Icentred_x_variableE2")) %>% 
    dplyr::select("l-95% CI", "u-95% CI") %>% 
    pivot_longer(cols = 1:2, names_to = "Interval", values_to = "Value")
  FDis.forb.df.ci <- FDis.forb.df %>% 
    filter(Covariate %in% c("Icentred_x_variableE2")) %>% 
    dplyr::select("l-95% CI", "u-95% CI") %>% 
    pivot_longer(cols = 1:2, names_to = "Interval", values_to = "Value")
  FDis.gram.df.ci <- FDis.gram.df %>% 
    filter(Covariate %in% c("Icentred_x_variableE2")) %>% 
    dplyr::select("l-95% CI", "u-95% CI") %>% 
    pivot_longer(cols = 1:2, names_to = "Interval", values_to = "Value")
  
  # Save the confidence intervals as a list and remove the intermediate dataframe
  FDis.shrub.ci <- as.list(FDis.shrub.df.ci$Value) # Adding values to the list
  names(FDis.shrub.ci) <- FDis.shrub.df.ci$Interval # Adding names to the values
  
  FDis.forb.ci <- as.list(FDis.forb.df.ci$Value) # Adding values to the list
  names(FDis.forb.ci) <- FDis.forb.df.ci$Interval # Adding names to the values
  
  FDis.gram.ci <- as.list(FDis.gram.df.ci$Value) # Adding values to the list
  names(FDis.gram.ci) <- FDis.gram.df.ci$Interval # Adding names to the values
  
  # Create a ggplot() item combining the three outputs
  (FDis.plot <- ggplot(data = FDis.predictions, aes(x = x_variable, y = .epred)) +
      
      # Ribbons
      stat_lineribbon(data = filter(FDis.predictions, Functional_Group == "Forbs"),
                      aes(y = .epred), .width = 0.95, alpha = 0.22, colour = "#F0CDAA", fill = "#F0CDAA", linewidth = 0.6,
                      linetype = ifelse((FDis.forb.ci$"l-95% CI" < 0 & FDis.forb.ci$"u-95% CI" < 0) | # Automatically assigns significance in subtitle based on CIs
                                          (FDis.forb.ci$"l-95% CI" > 0 & FDis.forb.ci$"u-95% CI" > 0), 1, 2)) +
      
      stat_lineribbon(data = filter(FDis.predictions, Functional_Group == "Graminoids"),
                      aes(y = .epred), .width = 0.95, alpha = 0.22, colour = "#FA922A", fill = "#FA922A", linewidth = 0.6,
                      linetype = ifelse((FDis.shrub.ci$"l-95% CI" < 0 & FDis.shrub.ci$"u-95% CI" < 0) | # Automatically assigns significance in subtitle based on CIs
                                          (FDis.shrub.ci$"l-95% CI" > 0 & FDis.shrub.ci$"u-95% CI" > 0), 1, 2)) +
      
      stat_lineribbon(data = filter(FDis.predictions, Functional_Group == "Shrubs"),
                      aes(y = .epred), .width = 0.95, alpha = 0.22, colour = "#E67300", fill = "#E67300", linewidth = 0.6,
                      linetype = ifelse((FDis.gram.ci$"l-95% CI" < 0 & FDis.gram.ci$"u-95% CI" < 0) | # Automatically assigns significance in subtitle based on CIs
                                          (FDis.gram.ci$"l-95% CI" > 0 & FDis.gram.ci$"u-95% CI" > 0), 1, 2)) +
      
      # Lines
      stat_lineribbon(data = filter(FDis.predictions, Functional_Group == "Forbs"),
                      aes(y = .epred), .width = 0.95, colour = "#F0CDAA", fill = NA, linewidth = 0.6,
                      linetype = ifelse((FDis.forb.ci$"l-95% CI" < 0 & FDis.forb.ci$"u-95% CI" < 0) | # Automatically assigns significance in subtitle based on CIs
                                          (FDis.forb.ci$"l-95% CI" > 0 & FDis.forb.ci$"u-95% CI" > 0), 1, 2)) +
      
      stat_lineribbon(data = filter(FDis.predictions, Functional_Group == "Graminoids"),
                      aes(y = .epred), .width = 0.95, colour = "#FA922A", fill = NA, linewidth = 0.6,
                      linetype = ifelse((FDis.shrub.ci$"l-95% CI" < 0 & FDis.shrub.ci$"u-95% CI" < 0) | # Automatically assigns significance in subtitle based on CIs
                                          (FDis.shrub.ci$"l-95% CI" > 0 & FDis.shrub.ci$"u-95% CI" > 0), 1, 2)) +
      
      stat_lineribbon(data = filter(FDis.predictions, Functional_Group == "Shrubs"),
                      aes(y = .epred), .width = 0.95, colour = "#E67300", fill = NA, linewidth = 0.6,
                      linetype = ifelse((FDis.gram.ci$"l-95% CI" < 0 & FDis.gram.ci$"u-95% CI" < 0) | # Automatically assigns significance in subtitle based on CIs
                                          (FDis.gram.ci$"l-95% CI" > 0 & FDis.gram.ci$"u-95% CI" > 0), 1, 2)) +
      
      labs(y = "Functional Dispersion \n",
           x = "\nCover (%)",
           fill = "Functional\nGroup",
           colour = "Functional\nGroup") +
      theme_1() +
      theme(legend.position = "none",
            legend.key = element_rect(colour = "#000000")))
  
  
  # Panel ----
  
  # Create a panel of the four plots
  combined.panel <- grid.arrange(FRic.plot, FEve.plot, FDis.plot, SR.plot, ncol = 4)
  
  # Add in if statements for saving the output
  if (censored == "Yes"){
    
    # Export panel
    ggsave(combined.panel, filename = paste0("figures/outputs_new/manuscript_",
                                             "all_cover_ALL_quadratic_censored", pc.filepath, ".png"), width = 18, height = 4.6)
    
  } else {
    
    # Export panel
    ggsave(combined.panel, filename = paste0("figures/outputs_new/manuscript_",
                                             "all_cover_ALL_quadratic", pc.filepath, ".png"), width = 18, height = 4.6)
    
  }
  
  
} # End of function


# FIGURE 5a: GENERATE THE CHANGE OVER TIME HISTOGRAMS ----

plot.change.histograms <- function(){
  
  # Load in model output
  FRic.mod <- get(load(paste0("data/model_outputs_new/m_FRic_slopes",
                              pc.filepath, ".RData")))
  
  # Save model output as a dataframe
  FRic.df <- brms_SummaryTable(FRic.mod, formatOptions = list(digits = 4, nsmall = 4), round = 4) %>% 
    filter(Covariate == "Intercept")
  
  # Determine bin width for plotting
  FRic_slopes.bin.width <- (max(filter(slopes.input, !is.na(FRic_slopes))$FRic_slopes) - min(filter(slopes.input, !is.na(FRic_slopes))$FRic_slopes))/20
  
  # Determine the estimated change
  FRic.change <- round(as.numeric(FRic.df$Estimate), digits = 3)
  
  # Plot model outputs
  (FRic.plot <- ggplot(slopes.input) +
      geom_histogram(aes(x = FRic_slopes), stat = "bin", binwidth = FRic_slopes.bin.width,
                     fill = "#D02090", colour = "#000000", alpha = 0.65) +
      geom_vline(aes(xintercept = FRic.change), colour = "#000000", linetype = "dashed", size = 0.5) +
      labs(y = "Plots \n",
           x = bquote('Functional Richness Change'~(yr^-1)),
           title = paste0("CIs: ", FRic.df$"l-95% CI", " to ", FRic.df$"u-95% CI"),
           subtitle = paste(ifelse((FRic.df$"l-95% CI" < 0 & FRic.df$"u-95% CI" < 0) | # Automatically asigns significance in subtitle based on CIs
                                     (FRic.df$"l-95% CI" > 0 & FRic.df$"u-95% CI" > 0),
                                   "Significant: CIs do NOT span 0", "Not Significant: CIs span 0"))) +
      theme_1() +
      theme(plot.subtitle = element_text(colour = ifelse((FRic.df$"l-95% CI" < 0 & FRic.df$"u-95% CI" < 0) | # Automatically colours subtitle based on significance from CIs
                                                           (FRic.df$"l-95% CI" > 0 & FRic.df$"u-95% CI" > 0),
                                                         paste("#006400"), paste("#8b0000")))))
  
  # Load in model output
  FEve.mod <- get(load(paste0("data/model_outputs_new/m_FEve_slopes",
                              pc.filepath, ".RData")))
  
  # Save model output as a dataframe
  FEve.df <- brms_SummaryTable(FEve.mod, formatOptions = list(digits = 4, nsmall = 4), round = 4) %>% 
    filter(Covariate == "Intercept")
  
  # Determine bin width for plotting
  FEve_slopes.bin.width <- (max(filter(slopes.input, !is.na(FEve_slopes))$FEve_slopes) - min(filter(slopes.input, !is.na(FEve_slopes))$FEve_slopes))/20
  
  # Determine the estimated change
  FEve.change <- round(as.numeric(FEve.df$Estimate), digits = 3)
  
  # Plot model outputs
  (FEve.plot <- ggplot(slopes.input) +
      geom_histogram(aes(x = FEve_slopes), stat = "bin", binwidth = FEve_slopes.bin.width,
                     fill = "#228B22", colour = "#000000", alpha = 0.65) +
      geom_vline(aes(xintercept = FEve.change), colour = "#000000", linetype = "dashed", size = 0.5) +
      labs(y = "Plots \n",
           x = bquote('Functional Evenness Change'~(yr^-1)),
           title = paste0("CIs: ", FEve.df$"l-95% CI", " to ", FEve.df$"u-95% CI"),
           subtitle = paste(ifelse((FEve.df$"l-95% CI" < 0 & FEve.df$"u-95% CI" < 0) | # Automatically asigns significance in subtitle based on CIs
                                     (FEve.df$"l-95% CI" > 0 & FEve.df$"u-95% CI" > 0),
                                   "Significant: CIs do NOT span 0", "Not Significant: CIs span 0"))) +
      theme_1() +
      theme(plot.subtitle = element_text(colour = ifelse((FEve.df$"l-95% CI" < 0 & FEve.df$"u-95% CI" < 0) | # Automatically colours subtitle based on significance from CIs
                                                           (FEve.df$"l-95% CI" > 0 & FEve.df$"u-95% CI" > 0),
                                                         paste("#006400"), paste("#8b0000")))))
  
  # Load in model output
  FDis.mod <- get(load(paste0("data/model_outputs_new/m_FDis_slopes",
                              pc.filepath, ".RData")))
  
  # Save model output as a dataframe
  FDis.df <- brms_SummaryTable(FDis.mod, formatOptions = list(digits = 4, nsmall = 4), round = 4) %>% 
    filter(Covariate == "Intercept")
  
  # Determine bin width for plotting
  FDis_slopes.bin.width <- (max(filter(slopes.input, !is.na(FDis_slopes))$FDis_slopes) - min(filter(slopes.input, !is.na(FDis_slopes))$FDis_slopes))/20
  
  # Determine the estimated change
  FDis.change <- round(as.numeric(FDis.df$Estimate), digits = 3)
  
  # Plot model outputs
  (FDis.plot <- ggplot(slopes.input) +
      geom_histogram(aes(x = FDis_slopes), stat = "bin", binwidth = FDis_slopes.bin.width,
                     fill = "#EE7600", colour = "#000000", alpha = 0.65) +
      geom_vline(aes(xintercept = FDis.change), colour = "#000000", linetype = "dashed", size = 0.5) +
      labs(y = "Plots \n",
           x = bquote('Functional Dispersion Change'~(yr^-1)),
           title = paste0("CIs: ", FDis.df$"l-95% CI", " to ", FDis.df$"u-95% CI"),
           subtitle = paste(ifelse((FDis.df$"l-95% CI" < 0 & FDis.df$"u-95% CI" < 0) | # Automatically asigns significance in subtitle based on CIs
                                     (FDis.df$"l-95% CI" > 0 & FDis.df$"u-95% CI" > 0),
                                   "Significant: CIs do NOT span 0", "Not Significant: CIs span 0"))) +
      theme_1() +
      theme(plot.subtitle = element_text(colour = ifelse((FDis.df$"l-95% CI" < 0 & FDis.df$"u-95% CI" < 0) | # Automatically colours subtitle based on significance from CIs
                                                           (FDis.df$"l-95% CI" > 0 & FDis.df$"u-95% CI" > 0),
                                                         paste("#006400"), paste("#8b0000")))))
  
  # Load in model output
  SR.mod <- get(load(paste0("data/model_outputs_new/m_SR_slopes",
                            pc.filepath, ".RData")))
  
  # Save model output as a dataframe
  SR.df <- brms_SummaryTable(SR.mod, formatOptions = list(digits = 4, nsmall = 4), round = 4) %>% 
    filter(Covariate == "Intercept")
  
  # Determine bin width for plotting
  SR_slopes.bin.width <- (max(filter(slopes.input, !is.na(SR_slopes))$SR_slopes) - min(filter(slopes.input, !is.na(SR_slopes))$SR_slopes))/20
  
  # Determine the estimated change
  SR.change <- round(as.numeric(SR.df$Estimate), digits = 3)
  
  # Plot model outputs
  (SR.plot <- ggplot(slopes.input) +
      geom_histogram(aes(x = SR_slopes), stat = "bin", binwidth = SR_slopes.bin.width,
                     fill = "#1C86EE", colour = "#000000", alpha = 0.65) +
      geom_vline(aes(xintercept = SR.change), colour = "#000000", linetype = "dashed", size = 0.5) +
      labs(y = "Plots \n",
           x = bquote('Species Richness Change'~(yr^-1)),
           title = paste0("CIs: ", SR.df$"l-95% CI", " to ", SR.df$"u-95% CI"),
           subtitle = paste(ifelse((SR.df$"l-95% CI" < 0 & SR.df$"u-95% CI" < 0) | # Automatically asigns significance in subtitle based on CIs
                                     (SR.df$"l-95% CI" > 0 & SR.df$"u-95% CI" > 0),
                                   "Significant: CIs do NOT span 0", "Not Significant: CIs span 0"))) +
      theme_1() +
      theme(plot.subtitle = element_text(colour = ifelse((SR.df$"l-95% CI" < 0 & SR.df$"u-95% CI" < 0) | # Automatically colours subtitle based on significance from CIs
                                                           (SR.df$"l-95% CI" > 0 & SR.df$"u-95% CI" > 0),
                                                         paste("#006400"), paste("#8b0000")))))
  
  # Load in model output
  ForbCover.mod <- get(load(paste0("data/model_outputs_new/m_ForbCover_slopes",
                                   pc.filepath, ".RData")))
  
  # Save model output as a dataframe
  ForbCover.df <- brms_SummaryTable(ForbCover.mod, formatOptions = list(digits = 4, nsmall = 4), round = 4) %>% 
    filter(Covariate == "Intercept")
  
  # Determine bin width for plotting
  ForbCover_slopes.bin.width <- (max(filter(slopes.input, !is.na(ForbCover_slopes))$ForbCover_slopes) - min(filter(slopes.input, !is.na(ForbCover_slopes))$ForbCover_slopes))/20
  
  # Determine the estimated change
  ForbCover.change <- round(as.numeric(ForbCover.df$Estimate), digits = 3)
  
  # Plot model outputs
  (ForbCover.plot <- ggplot(slopes.input) +
      geom_histogram(aes(x = ForbCover_slopes), stat = "bin", binwidth = ForbCover_slopes.bin.width,
                     fill = "#D02090", colour = "#000000", alpha = 0.35) +
      geom_vline(aes(xintercept = ForbCover.change), colour = "#000000", linetype = "dashed", size = 0.5) +
      labs(y = "Plots \n",
           x = bquote('Forb Cover Change'~(yr^-1)),
           title = paste0("CIs: ", ForbCover.df$"l-95% CI", " to ", ForbCover.df$"u-95% CI"),
           subtitle = paste(ifelse((ForbCover.df$"l-95% CI" < 0 & ForbCover.df$"u-95% CI" < 0) | # Automatically asigns significance in subtitle based on CIs
                                     (ForbCover.df$"l-95% CI" > 0 & ForbCover.df$"u-95% CI" > 0),
                                   "Significant: CIs do NOT span 0", "Not Significant: CIs span 0"))) +
      theme_1() +
      theme(plot.subtitle = element_text(colour = ifelse((ForbCover.df$"l-95% CI" < 0 & ForbCover.df$"u-95% CI" < 0) | # Automatically colours subtitle based on significance from CIs
                                                           (ForbCover.df$"l-95% CI" > 0 & ForbCover.df$"u-95% CI" > 0),
                                                         paste("#006400"), paste("#8b0000")))))
  
  # Load in model output
  GraminoidCover.mod <- get(load(paste0("data/model_outputs_new/m_GraminoidCover_slopes",
                                        pc.filepath, ".RData")))
  
  # Save model output as a dataframe
  GraminoidCover.df <- brms_SummaryTable(GraminoidCover.mod, formatOptions = list(digits = 4, nsmall = 4), round = 4) %>% 
    filter(Covariate == "Intercept")
  
  # Determine bin width for plotting
  GraminoidCover_slopes.bin.width <- (max(filter(slopes.input, !is.na(GraminoidCover_slopes))$GraminoidCover_slopes) - min(filter(slopes.input, !is.na(GraminoidCover_slopes))$GraminoidCover_slopes))/20
  
  # Determine the estimated change
  GraminoidCover.change <- round(as.numeric(GraminoidCover.df$Estimate), digits = 3)
  
  # Plot model outputs
  (GraminoidCover.plot <- ggplot(slopes.input) +
      geom_histogram(aes(x = GraminoidCover_slopes), stat = "bin", binwidth = GraminoidCover_slopes.bin.width,
                     fill = "#228B22", colour = "#000000", alpha = 0.35) +
      geom_vline(aes(xintercept = GraminoidCover.change), colour = "#000000", linetype = "dashed", size = 0.5) +
      labs(y = "Plots \n",
           x = bquote('Graminoid Cover Change'~(yr^-1)),
           title = paste0("CIs: ", GraminoidCover.df$"l-95% CI", " to ", GraminoidCover.df$"u-95% CI"),
           subtitle = paste(ifelse((GraminoidCover.df$"l-95% CI" < 0 & GraminoidCover.df$"u-95% CI" < 0) | # Automatically asigns significance in subtitle based on CIs
                                     (GraminoidCover.df$"l-95% CI" > 0 & GraminoidCover.df$"u-95% CI" > 0),
                                   "Significant: CIs do NOT span 0", "Not Significant: CIs span 0"))) +
      theme_1() +
      theme(plot.subtitle = element_text(colour = ifelse((GraminoidCover.df$"l-95% CI" < 0 & GraminoidCover.df$"u-95% CI" < 0) | # Automatically colours subtitle based on significance from CIs
                                                           (GraminoidCover.df$"l-95% CI" > 0 & GraminoidCover.df$"u-95% CI" > 0),
                                                         paste("#006400"), paste("#8b0000")))))
  
  # Load in model output
  ShrubCover.mod <- get(load(paste0("data/model_outputs_new/m_ShrubCover_slopes",
                                    pc.filepath, ".RData")))
  
  # Save model output as a dataframe
  ShrubCover.df <- brms_SummaryTable(ShrubCover.mod, formatOptions = list(digits = 4, nsmall = 4), round = 4) %>% 
    filter(Covariate == "Intercept")
  
  # Determine bin width for plotting
  ShrubCover_slopes.bin.width <- (max(filter(slopes.input, !is.na(ShrubCover_slopes))$ShrubCover_slopes) - min(filter(slopes.input, !is.na(ShrubCover_slopes))$ShrubCover_slopes))/20
  
  # Determine the estimated change
  ShrubCover.change <- round(as.numeric(ShrubCover.df$Estimate), digits = 3)
  
  # Plot model outputs
  (ShrubCover.plot <- ggplot(slopes.input) +
      geom_histogram(aes(x = ShrubCover_slopes), stat = "bin", binwidth = ShrubCover_slopes.bin.width,
                     fill = "#1C86EE", colour = "#000000", alpha = 0.35) +
      geom_vline(aes(xintercept = ShrubCover.change), colour = "#000000", linetype = "dashed", size = 0.5) +
      labs(y = "Plots \n",
           x = bquote('Shrub Cover Change'~(yr^-1)),
           title = paste0("CIs: ", ShrubCover.df$"l-95% CI", " to ", ShrubCover.df$"u-95% CI"),
           subtitle = paste(ifelse((ShrubCover.df$"l-95% CI" < 0 & ShrubCover.df$"u-95% CI" < 0) | # Automatically asigns significance in subtitle based on CIs
                                     (ShrubCover.df$"l-95% CI" > 0 & ShrubCover.df$"u-95% CI" > 0),
                                   "Significant: CIs do NOT span 0", "Not Significant: CIs span 0"))) +
      theme_1() +
      theme(plot.subtitle = element_text(colour = ifelse((ShrubCover.df$"l-95% CI" < 0 & ShrubCover.df$"u-95% CI" < 0) | # Automatically colours subtitle based on significance from CIs
                                                           (ShrubCover.df$"l-95% CI" > 0 & ShrubCover.df$"u-95% CI" > 0),
                                                         paste("#006400"), paste("#8b0000")))))
  
  # Combine the plots into one output
  combined.metric.panel <- grid.arrange(FRic.plot, FEve.plot, FDis.plot, SR.plot, ncol = 4)
  
  # Export the plot
  ggsave(combined.metric.panel, file = paste0("figures/outputs_new/manuscript_",
                                              "change_over_time_metrics", pc.filepath, ".png"), width = 20, height = 5)
  
  # Combine the plots into one output
  combined.cover.panel <- grid.arrange(ForbCover.plot, GraminoidCover.plot, ShrubCover.plot, ncol = 3)
  
  # Export the plot
  ggsave(combined.cover.panel, file = paste0("figures/outputs_new/manuscript_",
                                             "change_over_time_cover", pc.filepath, ".png"), width = 17.25, height = 5.5)
  
}


# FIGURE 5b: GENERATE COMBINED PLOTS FOR COVER CHANGE ----

plot.combined.cover.change <- function(){
  
    # Load in model output
    FRic.shrub.mod <- get(load(paste0("data/model_outputs_new/m_FRic_slopes_",
                                      "ShrubCover_slopes", pc.filepath, ".RData")))
    FRic.forb.mod <- get(load(paste0("data/model_outputs_new/m_FRic_slopes_",
                                     "ForbCover_slopes", pc.filepath, ".RData")))
    FRic.gram.mod <- get(load(paste0("data/model_outputs_new/m_FRic_slopes_",
                                     "GraminoidCover_slopes", pc.filepath, ".RData")))
    
    # Convert model output into a dataframe (with 4.d.p.)
    FRic.shrub.df <- brms_SummaryTable(FRic.shrub.mod, formatOptions = list(digits = 4, nsmall = 4), round = 4)
    FRic.forb.df <- brms_SummaryTable(FRic.forb.mod, formatOptions = list(digits = 4, nsmall = 4), round = 4)
    FRic.gram.df <- brms_SummaryTable(FRic.gram.mod, formatOptions = list(digits = 4, nsmall = 4), round = 4)
    
    # Extract the confidence intervals as a list for use in the plotting
    FRic.shrub.df.ci <- FRic.shrub.df %>% 
      filter(Covariate %in% c("x_variable")) %>% 
      dplyr::select("l-95% CI", "u-95% CI") %>% 
      pivot_longer(cols = 1:2, names_to = "Interval", values_to = "Value")
    FRic.forb.df.ci <- FRic.forb.df %>% 
      filter(Covariate %in% c("x_variable")) %>% 
      dplyr::select("l-95% CI", "u-95% CI") %>% 
      pivot_longer(cols = 1:2, names_to = "Interval", values_to = "Value")
    FRic.gram.df.ci <- FRic.gram.df %>% 
      filter(Covariate %in% c("x_variable")) %>% 
      dplyr::select("l-95% CI", "u-95% CI") %>% 
      pivot_longer(cols = 1:2, names_to = "Interval", values_to = "Value")
    
    # Save the confidence intervals as a list and remove the intermediate dataframe
    FRic.shrub.ci <- as.list(FRic.shrub.df.ci$Value) # Adding values to the list
    names(FRic.shrub.ci) <- FRic.shrub.df.ci$Interval # Adding names to the values

    FRic.forb.ci <- as.list(FRic.forb.df.ci$Value) # Adding values to the list
    names(FRic.forb.ci) <- FRic.forb.df.ci$Interval # Adding names to the values
    
    FRic.gram.ci <- as.list(FRic.gram.df.ci$Value) # Adding values to the list
    names(FRic.gram.ci) <- FRic.gram.df.ci$Interval # Adding names to the values
    
    # Extract the prediction data frame
    FRic.shrub.pred <- ggpredict(FRic.shrub.mod, terms = "x_variable [-5:5, sample = 50]", back.transform = TRUE) %>% 
      mutate(Functional_Group = "Shrubs")
    FRic.forb.pred <- ggpredict(FRic.forb.mod, terms = "x_variable [-5:5, sample = 50]", back.transform = TRUE) %>% 
      mutate(Functional_Group = "Forbs")
    FRic.gram.pred <- ggpredict(FRic.gram.mod, terms = "x_variable [-5:5, sample = 50]", back.transform = TRUE) %>% 
      mutate(Functional_Group = "Graminoids")
    
    # Combine the predictions into one dataframe
    FRic.predictions <- rbind(FRic.shrub.pred, FRic.forb.pred, FRic.gram.pred) %>% 
      rename(Cover_Change = x) %>% 
      dplyr::select(-group)
    
      # NOTE: plotting only between -5 and 5, but all values in model
    
    # Create a ggplot() item combining the three outputs
    (FRic.plot <- ggplot(data = FRic.predictions) +
        geom_ribbon(aes(x = Cover_Change, ymin = conf.low, ymax = conf.high, fill = Functional_Group),
                    alpha = 0.35, linetype = 2) +
        scale_fill_manual(values = c("#F5D7EA", "#E685C3", "#B50675")) +
        
        geom_line(data = filter(FRic.predictions, Functional_Group == "Forbs"),
                  aes(x = Cover_Change, y = predicted), colour = "#B50675",
                  linetype = ifelse((FRic.forb.ci$"l-95% CI" < 0 & FRic.forb.ci$"u-95% CI" < 0) | # Automatically asigns significance in subtitle based on CIs
                                      (FRic.forb.ci$"l-95% CI" > 0 & FRic.forb.ci$"u-95% CI" > 0), 1, 2)) +
        
        geom_line(data = filter(FRic.predictions, Functional_Group == "Shrubs"),
                  aes(x = Cover_Change, y = predicted), colour = "#B50675",
                  linetype = ifelse((FRic.shrub.ci$"l-95% CI" < 0 & FRic.shrub.ci$"u-95% CI" < 0) | # Automatically asigns significance in subtitle based on CIs
                                      (FRic.shrub.ci$"l-95% CI" > 0 & FRic.shrub.ci$"u-95% CI" > 0), 1, 2)) +
        
        geom_line(data = filter(FRic.predictions, Functional_Group == "Graminoids"),
                  aes(x = Cover_Change, y = predicted), colour = "#B50675",
                  linetype = ifelse((FRic.gram.ci$"l-95% CI" < 0 & FRic.gram.ci$"u-95% CI" < 0) | # Automatically asigns significance in subtitle based on CIs
                                      (FRic.gram.ci$"l-95% CI" > 0 & FRic.gram.ci$"u-95% CI" > 0), 1, 2)) +
        
        labs(y = bquote('Functional Richness Change'~(yr^-1)),
             x = bquote('Cover Change'~(yr^-1)),
             fill = "Functional\nGroup") +
        theme_1() +
        theme(legend.position = "bottom",
              legend.title = element_blank(),
              legend.key = element_rect(colour = "#000000")))
    
  
    # Load in model output
    FEve.shrub.mod <- get(load(paste0("data/model_outputs_new/m_FEve_slopes_",
                                      "ShrubCover_slopes", pc.filepath, ".RData")))
    FEve.forb.mod <- get(load(paste0("data/model_outputs_new/m_FEve_slopes_",
                                     "ForbCover_slopes", pc.filepath, ".RData")))
    FEve.gram.mod <- get(load(paste0("data/model_outputs_new/m_FEve_slopes_",
                                     "GraminoidCover_slopes", pc.filepath, ".RData")))
    
    # Convert model output into a dataframe (with 4.d.p.)
    FEve.shrub.df <- brms_SummaryTable(FEve.shrub.mod, formatOptions = list(digits = 4, nsmall = 4), round = 4)
    FEve.forb.df <- brms_SummaryTable(FEve.forb.mod, formatOptions = list(digits = 4, nsmall = 4), round = 4)
    FEve.gram.df <- brms_SummaryTable(FEve.gram.mod, formatOptions = list(digits = 4, nsmall = 4), round = 4)
    
    # Extract the confidence intervals as a list for use in the plotting
    FEve.shrub.df.ci <- FEve.shrub.df %>% 
      filter(Covariate %in% c("x_variable")) %>% 
      dplyr::select("l-95% CI", "u-95% CI") %>% 
      pivot_longer(cols = 1:2, names_to = "Interval", values_to = "Value")
    FEve.forb.df.ci <- FEve.forb.df %>% 
      filter(Covariate %in% c("x_variable")) %>% 
      dplyr::select("l-95% CI", "u-95% CI") %>% 
      pivot_longer(cols = 1:2, names_to = "Interval", values_to = "Value")
    FEve.gram.df.ci <- FEve.gram.df %>% 
      filter(Covariate %in% c("x_variable")) %>% 
      dplyr::select("l-95% CI", "u-95% CI") %>% 
      pivot_longer(cols = 1:2, names_to = "Interval", values_to = "Value")
    
    # Save the confidence intervals as a list and remove the intermediate dataframe
    FEve.shrub.ci <- as.list(FEve.shrub.df.ci$Value) # Adding values to the list
    names(FEve.shrub.ci) <- FEve.shrub.df.ci$Interval # Adding names to the values
    
    FEve.forb.ci <- as.list(FEve.forb.df.ci$Value) # Adding values to the list
    names(FEve.forb.ci) <- FEve.forb.df.ci$Interval # Adding names to the values
    
    FEve.gram.ci <- as.list(FEve.gram.df.ci$Value) # Adding values to the list
    names(FEve.gram.ci) <- FEve.gram.df.ci$Interval # Adding names to the values
    
    # Extract the prediction data frame
    FEve.shrub.pred <- ggpredict(FEve.shrub.mod, terms = "x_variable [-5:5, sample = 50]", back.transform = TRUE) %>% 
      mutate(Functional_Group = "Shrubs")
    FEve.forb.pred <- ggpredict(FEve.forb.mod, terms = "x_variable [-5:5, sample = 50]", back.transform = TRUE) %>% 
      mutate(Functional_Group = "Forbs")
    FEve.gram.pred <- ggpredict(FEve.gram.mod, terms = "x_variable [-5:5, sample = 50]", back.transform = TRUE) %>% 
      mutate(Functional_Group = "Graminoids")
    
    # Combine the predictions into one dataframe
    FEve.predictions <- rbind(FEve.shrub.pred, FEve.forb.pred, FEve.gram.pred) %>% 
      rename(Cover_Change = x) %>% 
      dplyr::select(-group)
    
      # NOTE: plotting only between -5 and 5, but all values in model
    
    # Create a ggplot() item combining the three outputs
    (FEve.plot <- ggplot(data = FEve.predictions) +
        geom_ribbon(aes(x = Cover_Change, ymin = conf.low, ymax = conf.high, fill = Functional_Group),
                    alpha = 0.35, linetype = 2) +
        scale_fill_manual(values = c("#D1F0D1", "#67B867", "#0A6E0A")) +
        
        geom_line(data = filter(FEve.predictions, Functional_Group == "Forbs"),
                  aes(x = Cover_Change, y = predicted), colour = "#0A6E0A",
                  linetype = ifelse((FEve.forb.ci$"l-95% CI" < 0 & FEve.forb.ci$"u-95% CI" < 0) | # Automatically asigns significance in subtitle based on CIs
                                      (FEve.forb.ci$"l-95% CI" > 0 & FEve.forb.ci$"u-95% CI" > 0), 1, 2)) +
        
        geom_line(data = filter(FEve.predictions, Functional_Group == "Shrubs"),
                  aes(x = Cover_Change, y = predicted), colour = "#0A6E0A",
                  linetype = ifelse((FEve.shrub.ci$"l-95% CI" < 0 & FEve.shrub.ci$"u-95% CI" < 0) | # Automatically asigns significance in subtitle based on CIs
                                      (FEve.shrub.ci$"l-95% CI" > 0 & FEve.shrub.ci$"u-95% CI" > 0), 1, 2)) +
        
        geom_line(data = filter(FEve.predictions, Functional_Group == "Graminoids"),
                  aes(x = Cover_Change, y = predicted), colour = "#0A6E0A",
                  linetype = ifelse((FEve.gram.ci$"l-95% CI" < 0 & FEve.gram.ci$"u-95% CI" < 0) | # Automatically asigns significance in subtitle based on CIs
                                      (FEve.gram.ci$"l-95% CI" > 0 & FEve.gram.ci$"u-95% CI" > 0), 1, 2)) +
        
        labs(y = bquote('Functional Evenness Change'~(yr^-1)),
             x = bquote('Cover Change'~(yr^-1)),
             fill = "Functional\nGroup") +
        theme_1() +
        theme(legend.position = "bottom",
              legend.title = element_blank(),
              legend.key = element_rect(colour = "#000000")))
    
    
    # Load in model output
    FDis.shrub.mod <- get(load(paste0("data/model_outputs_new/m_FDis_slopes_",
                                      "ShrubCover_slopes", pc.filepath, ".RData")))
    FDis.forb.mod <- get(load(paste0("data/model_outputs_new/m_FDis_slopes_",
                                     "ForbCover_slopes", pc.filepath, ".RData")))
    FDis.gram.mod <- get(load(paste0("data/model_outputs_new/m_FDis_slopes_",
                                     "GraminoidCover_slopes", pc.filepath, ".RData")))
    
    # Convert model output into a dataframe (with 4.d.p.)
    FDis.shrub.df <- brms_SummaryTable(FDis.shrub.mod, formatOptions = list(digits = 4, nsmall = 4), round = 4)
    FDis.forb.df <- brms_SummaryTable(FDis.forb.mod, formatOptions = list(digits = 4, nsmall = 4), round = 4)
    FDis.gram.df <- brms_SummaryTable(FDis.gram.mod, formatOptions = list(digits = 4, nsmall = 4), round = 4)
    
    # Extract the confidence intervals as a list for use in the plotting
    FDis.shrub.df.ci <- FDis.shrub.df %>% 
      filter(Covariate %in% c("x_variable")) %>% 
      dplyr::select("l-95% CI", "u-95% CI") %>% 
      pivot_longer(cols = 1:2, names_to = "Interval", values_to = "Value")
    FDis.forb.df.ci <- FDis.forb.df %>% 
      filter(Covariate %in% c("x_variable")) %>% 
      dplyr::select("l-95% CI", "u-95% CI") %>% 
      pivot_longer(cols = 1:2, names_to = "Interval", values_to = "Value")
    FDis.gram.df.ci <- FDis.gram.df %>% 
      filter(Covariate %in% c("x_variable")) %>% 
      dplyr::select("l-95% CI", "u-95% CI") %>% 
      pivot_longer(cols = 1:2, names_to = "Interval", values_to = "Value")
    
    # Save the confidence intervals as a list and remove the intermediate dataframe
    FDis.shrub.ci <- as.list(FDis.shrub.df.ci$Value) # Adding values to the list
    names(FDis.shrub.ci) <- FDis.shrub.df.ci$Interval # Adding names to the values
    
    FDis.forb.ci <- as.list(FDis.forb.df.ci$Value) # Adding values to the list
    names(FDis.forb.ci) <- FDis.forb.df.ci$Interval # Adding names to the values
    
    FDis.gram.ci <- as.list(FDis.gram.df.ci$Value) # Adding values to the list
    names(FDis.gram.ci) <- FDis.gram.df.ci$Interval # Adding names to the values
    
    # Extract the prediction data frame
    FDis.shrub.pred <- ggpredict(FDis.shrub.mod, terms = "x_variable [-5:5, sample = 50]", back.transform = TRUE) %>% 
      mutate(Functional_Group = "Shrubs")
    FDis.forb.pred <- ggpredict(FDis.forb.mod, terms = "x_variable [-5:5, sample = 50]", back.transform = TRUE) %>% 
      mutate(Functional_Group = "Forbs")
    FDis.gram.pred <- ggpredict(FDis.gram.mod, terms = "x_variable [-5:5, sample = 50]", back.transform = TRUE) %>% 
      mutate(Functional_Group = "Graminoids")
    
    # Combine the predictions into one dataframe
    FDis.predictions <- rbind(FDis.shrub.pred, FDis.forb.pred, FDis.gram.pred) %>% 
      rename(Cover_Change = x) %>% 
      dplyr::select(-group)
    
    # NOTE: plotting only between -5 and 5, but all values in model
    
    # Create a ggplot() item combining the three outputs
    (FDis.plot <- ggplot(data = FDis.predictions) +
        geom_ribbon(aes(x = Cover_Change, ymin = conf.low, ymax = conf.high, fill = Functional_Group),
                    alpha = 0.35, linetype = 2) +
        scale_fill_manual(values = c("#F0CDAA", "#FA922A", "#E67300")) +
        
        geom_line(data = filter(FDis.predictions, Functional_Group == "Forbs"),
                  aes(x = Cover_Change, y = predicted), colour = "#E67300",
                  linetype = ifelse((FDis.forb.ci$"l-95% CI" < 0 & FDis.forb.ci$"u-95% CI" < 0) | # Automatically asigns significance in subtitle based on CIs
                                      (FDis.forb.ci$"l-95% CI" > 0 & FDis.forb.ci$"u-95% CI" > 0), 1, 2)) +
        
        geom_line(data = filter(FDis.predictions, Functional_Group == "Shrubs"),
                  aes(x = Cover_Change, y = predicted), colour = "#E67300",
                  linetype = ifelse((FDis.shrub.ci$"l-95% CI" < 0 & FDis.shrub.ci$"u-95% CI" < 0) | # Automatically asigns significance in subtitle based on CIs
                                      (FDis.shrub.ci$"l-95% CI" > 0 & FDis.shrub.ci$"u-95% CI" > 0), 1, 2)) +
        
        geom_line(data = filter(FDis.predictions, Functional_Group == "Graminoids"),
                  aes(x = Cover_Change, y = predicted), colour = "#E67300",
                  linetype = ifelse((FDis.gram.ci$"l-95% CI" < 0 & FDis.gram.ci$"u-95% CI" < 0) | # Automatically asigns significance in subtitle based on CIs
                                      (FDis.gram.ci$"l-95% CI" > 0 & FDis.gram.ci$"u-95% CI" > 0), 1, 2)) +
        
        labs(y = bquote('Functional Dispersion Change'~(yr^-1)),
             x = bquote('Cover Change'~(yr^-1)),
             fill = "Functional\nGroup") +
        theme_1() +
        theme(legend.position = "bottom",
              legend.title = element_blank(),
              legend.key = element_rect(colour = "#000000")))
    
    
    # Load in model output
    SR.shrub.mod <- get(load(paste0("data/model_outputs_new/m_SR_slopes_",
                                    "ShrubCover_slopes", pc.filepath, ".RData")))
    SR.forb.mod <- get(load(paste0("data/model_outputs_new/m_SR_slopes_",
                                   "ForbCover_slopes", pc.filepath, ".RData")))
    SR.gram.mod <- get(load(paste0("data/model_outputs_new/m_SR_slopes_",
                                   "GraminoidCover_slopes", pc.filepath, ".RData")))
    
    # Convert model output into a dataframe (with 4.d.p.)
    SR.shrub.df <- brms_SummaryTable(SR.shrub.mod, formatOptions = list(digits = 4, nsmall = 4), round = 4)
    SR.forb.df <- brms_SummaryTable(SR.forb.mod, formatOptions = list(digits = 4, nsmall = 4), round = 4)
    SR.gram.df <- brms_SummaryTable(SR.gram.mod, formatOptions = list(digits = 4, nsmall = 4), round = 4)
    
    # Extract the confidence intervals as a list for use in the plotting
    SR.shrub.df.ci <- SR.shrub.df %>% 
      filter(Covariate %in% c("x_variable")) %>% 
      dplyr::select("l-95% CI", "u-95% CI") %>% 
      pivot_longer(cols = 1:2, names_to = "Interval", values_to = "Value")
    SR.forb.df.ci <- SR.forb.df %>% 
      filter(Covariate %in% c("x_variable")) %>% 
      dplyr::select("l-95% CI", "u-95% CI") %>% 
      pivot_longer(cols = 1:2, names_to = "Interval", values_to = "Value")
    SR.gram.df.ci <- SR.gram.df %>% 
      filter(Covariate %in% c("x_variable")) %>% 
      dplyr::select("l-95% CI", "u-95% CI") %>% 
      pivot_longer(cols = 1:2, names_to = "Interval", values_to = "Value")
    
    # Save the confidence intervals as a list and remove the intermediate dataframe
    SR.shrub.ci <- as.list(SR.shrub.df.ci$Value) # Adding values to the list
    names(SR.shrub.ci) <- SR.shrub.df.ci$Interval # Adding names to the values
    
    SR.forb.ci <- as.list(SR.forb.df.ci$Value) # Adding values to the list
    names(SR.forb.ci) <- SR.forb.df.ci$Interval # Adding names to the values
    
    SR.gram.ci <- as.list(SR.gram.df.ci$Value) # Adding values to the list
    names(SR.gram.ci) <- SR.gram.df.ci$Interval # Adding names to the values
    
    # Extract the prediction data frame
    SR.shrub.pred <- ggpredict(SR.shrub.mod, terms = "x_variable [-5:5, sample = 50]", back.transform = TRUE) %>% 
      mutate(Functional_Group = "Shrubs")
    SR.forb.pred <- ggpredict(SR.forb.mod, terms = "x_variable [-5:5, sample = 50]", back.transform = TRUE) %>% 
      mutate(Functional_Group = "Forbs")
    SR.gram.pred <- ggpredict(SR.gram.mod, terms = "x_variable [-5:5, sample = 50]", back.transform = TRUE) %>% 
      mutate(Functional_Group = "Graminoids")
    
    # Combine the predictions into one dataframe
    SR.predictions <- rbind(SR.shrub.pred, SR.forb.pred, SR.gram.pred) %>% 
      rename(Cover_Change = x) %>% 
      dplyr::select(-group)
    
      # NOTE: plotting only between -5 and 5, but all values in model
    
    # Create a ggplot() item combining the three outputs
    (SR.plot <- ggplot(data = SR.predictions) +
        geom_ribbon(aes(x = Cover_Change, ymin = conf.low, ymax = conf.high, fill = Functional_Group),
                    alpha = 0.35, linetype = 2) +
        scale_fill_manual(values = c("#B1CCE8", "#1C86EE", "#0D3E6E")) +
        
        geom_line(data = filter(SR.predictions, Functional_Group == "Forbs"),
                  aes(x = Cover_Change, y = predicted), colour = "#0D3E6E",
                  linetype = ifelse((SR.forb.ci$"l-95% CI" < 0 & SR.forb.ci$"u-95% CI" < 0) | # Automatically asigns significance in subtitle based on CIs
                                      (SR.forb.ci$"l-95% CI" > 0 & SR.forb.ci$"u-95% CI" > 0), 1, 2)) +
        
        geom_line(data = filter(SR.predictions, Functional_Group == "Shrubs"),
                  aes(x = Cover_Change, y = predicted), colour = "#0D3E6E",
                  linetype = ifelse((SR.shrub.ci$"l-95% CI" < 0 & SR.shrub.ci$"u-95% CI" < 0) | # Automatically asigns significance in subtitle based on CIs
                                      (SR.shrub.ci$"l-95% CI" > 0 & SR.shrub.ci$"u-95% CI" > 0), 1, 2)) +
        
        geom_line(data = filter(SR.predictions, Functional_Group == "Graminoids"),
                  aes(x = Cover_Change, y = predicted), colour = "#0D3E6E",
                  linetype = ifelse((SR.gram.ci$"l-95% CI" < 0 & SR.gram.ci$"u-95% CI" < 0) | # Automatically asigns significance in subtitle based on CIs
                                      (SR.gram.ci$"l-95% CI" > 0 & SR.gram.ci$"u-95% CI" > 0), 1, 2)) +
        
        labs(y = bquote('Species Richness Change'~(yr^-1)),
             x = bquote('Cover Change'~(yr^-1)),
             fill = "Functional\nGroup") +
        theme_1() +
        theme(legend.position = "bottom",
              legend.title = element_blank(),
              legend.key = element_rect(colour = "#000000")))
    
    # Create a panel of all three plots
    combined.panel <- grid.arrange(FRic.plot, FEve.plot, FDis.plot, SR.plot, ncol = 4, widths = c(1,1,1,1))
    
    # Output the saved panel
    ggsave(combined.panel, filename = paste0("figures/outputs_new/manuscript_",
                                             "all_cover_slopes", pc.filepath, ".png"), width = 20, height = 5.5)
    
} # End of function


# Figure 5a?: GENERATE COMBINED PLOTS FOR TEMP CHANGS VS METRICS ----

plot.temp.change <- function(){
  
  # Load in model output
  FRic.mod <- get(load(paste0("data/model_outputs_new/m_FRic_slopes_",
                              "WarmQSlope", pc.filepath, ".RData")))
  
  # Range to predict variable
  range.to.predict <- paste0("x_variable [", min(slopes.input$WarmQSlope, na.rm = TRUE), ":", max(slopes.input$WarmQSlope, na.rm = TRUE), 
                             " by =", (max(slopes.input$WarmQSlope, na.rm = TRUE) - min(slopes.input$WarmQSlope, na.rm = TRUE))/50, "]")
  
  # Extract the prediction data frame
  FRic.pred <- ggpredict(FRic.mod, terms = range.to.predict)
  
  # Convert model output into a dataframe (with 4.d.p.)
  FRic.df <- brms_SummaryTable(FRic.mod, formatOptions = list(digits = 4, nsmall = 4), round = 4)
  
  # Extract the confidence intervals as a list for use in the plotting
  FRic.df.ci <- FRic.df %>% 
    filter(Covariate %in% c("x_variable")) %>% 
    dplyr::select("l-95% CI", "u-95% CI") %>% 
    pivot_longer(cols = 1:2, names_to = "Interval", values_to = "Value")
  
  # Save the confidence intervals as a list and remove the intermediate dataframe
  FRic.ci <- as.list(FRic.df.ci$Value) # Adding values to the list
  names(FRic.ci) <- FRic.df.ci$Interval # Adding names to the values
  rm(FRic.df, FRic.df.ci) # Remove unnecessary dataframes of summary and CIs
  
  # Plot model outputs
  (FRic.plot <- ggplot(FRic.pred) +
      
      # Hexplot
      geom_hex(data = slopes.input, aes(x = WarmQSlope, y = FRic_slopes, bins = 18)) +
      scale_fill_gradient(low = "#FFE1FF", high = "#D02090") +
      geom_ribbon(aes(x = x, ymin = conf.low, ymax = conf.high),
                  fill = "#D02090", alpha = 0.2, colour = "#D02090", linewdith = 0.05) + # Adds c.intervals for predictions as ribbon
      geom_line(aes(x = x, y = predicted), colour = "#000000", linewidth = 1.1,
                linetype = ifelse((FRic.ci$"l-95% CI" < 0 & FRic.ci$"u-95% CI" < 0) | # Automatically asigns significance in subtitle based on CIs
                                    (FRic.ci$"l-95% CI" > 0 & FRic.ci$"u-95% CI" > 0), 1, 2)) + # Adds line for predicted values of FRic vs LAT
      
      labs(y = bquote('Functional Richness Change'~(yr^-1)),
           x = paste0("\nTemperature Change (", "\u00B0", "C yr", "\U207B", "\U00B9", ")"),
           title = paste0("CIs: ", FRic.ci$"l-95% CI", " to ", FRic.ci$"u-95% CI"),
           subtitle = paste(ifelse((FRic.ci$"l-95% CI" < 0 & FRic.ci$"u-95% CI" < 0) | # Automatically asigns significance in subtitle based on CIs
                                     (FRic.ci$"l-95% CI" > 0 & FRic.ci$"u-95% CI" > 0),
                                   "Significant: CIs do NOT span 0", "Not Significant: CIs span 0"))) +
      theme_1() +
      theme(legend.position = "none",
            plot.subtitle = element_text(colour = ifelse((FRic.ci$"l-95% CI" < 0 & FRic.ci$"u-95% CI" < 0) | # Automatically colours subtitle based on significance from CIs
                                                           (FRic.ci$"l-95% CI" > 0 & FRic.ci$"u-95% CI" > 0),
                                                         paste("#006400"), paste("#8b0000")))))
  
  # Load in model output
  FEve.mod <- get(load(paste0("data/model_outputs_new/m_FEve_slopes_",
                              "WarmQSlope", pc.filepath, ".RData")))
  
  # Range to predict variable
  range.to.predict <- paste0("x_variable [", min(slopes.input$WarmQSlope, na.rm = TRUE), ":", max(slopes.input$WarmQSlope, na.rm = TRUE), 
                             " by =", (max(slopes.input$WarmQSlope, na.rm = TRUE) - min(slopes.input$WarmQSlope, na.rm = TRUE))/50, "]")
  
  # Extract the prediction data frame
  FEve.pred <- ggpredict(FEve.mod, terms = range.to.predict)
  
  # Convert model output into a dataframe (with 4.d.p.)
  FEve.df <- brms_SummaryTable(FEve.mod, formatOptions = list(digits = 4, nsmall = 4), round = 4)
  
  # Extract the confidence intervals as a list for use in the plotting
  FEve.df.ci <- FEve.df %>% 
    filter(Covariate %in% c("x_variable")) %>% 
    dplyr::select("l-95% CI", "u-95% CI") %>% 
    pivot_longer(cols = 1:2, names_to = "Interval", values_to = "Value")
  
  # Save the confidence intervals as a list and remove the intermediate dataframe
  FEve.ci <- as.list(FEve.df.ci$Value) # Adding values to the list
  names(FEve.ci) <- FEve.df.ci$Interval # Adding names to the values
  rm(FEve.df, FEve.df.ci) # Remove unnecessary dataframes of summary and CIs
  
  # Plot model outputs
  (FEve.plot <- ggplot(FEve.pred) +
      
      # Hexplot
      geom_hex(data = slopes.input, aes(x = WarmQSlope, y = FEve_slopes, bins = 18)) +
      scale_fill_gradient(low = "#CFF099", high = "#228B22") +
      geom_ribbon(aes(x = x, ymin = conf.low, ymax = conf.high),
                  fill = "#228B22", alpha = 0.2, colour = "#228B22", linewdith = 0.05) + # Adds c.intervals for predictions as ribbon
      geom_line(aes(x = x, y = predicted), colour = "#000000", linewidth = 1.1,
                linetype = ifelse((FEve.ci$"l-95% CI" < 0 & FEve.ci$"u-95% CI" < 0) | # Automatically asigns significance in subtitle based on CIs
                                    (FEve.ci$"l-95% CI" > 0 & FEve.ci$"u-95% CI" > 0), 1, 2)) + # Adds line for predicted values of FEve vs LAT
      
      labs(y = bquote('Functional Evenness Change'~(yr^-1)),
           x = paste0("\nTemperature Change (", "\u00B0", "C yr", "\U207B", "\U00B9", ")"),
           title = paste0("CIs: ", FEve.ci$"l-95% CI", " to ", FEve.ci$"u-95% CI"),
           subtitle = paste(ifelse((FEve.ci$"l-95% CI" < 0 & FEve.ci$"u-95% CI" < 0) | # Automatically asigns significance in subtitle based on CIs
                                     (FEve.ci$"l-95% CI" > 0 & FEve.ci$"u-95% CI" > 0),
                                   "Significant: CIs do NOT span 0", "Not Significant: CIs span 0"))) +
      theme_1() +
      theme(legend.position = "none",
            plot.subtitle = element_text(colour = ifelse((FEve.ci$"l-95% CI" < 0 & FEve.ci$"u-95% CI" < 0) | # Automatically colours subtitle based on significance from CIs
                                                           (FEve.ci$"l-95% CI" > 0 & FEve.ci$"u-95% CI" > 0),
                                                         paste("#006400"), paste("#8b0000")))))
  
  
  # Load in model output
  SR.mod <- get(load(paste0("data/model_outputs_new/m_SR_slopes_",
                            "WarmQSlope", pc.filepath, ".RData")))
  
  # Range to predict variable
  range.to.predict <- paste0("x_variable [", min(slopes.input$WarmQSlope, na.rm = TRUE), ":", max(slopes.input$WarmQSlope, na.rm = TRUE), 
                             " by =", (max(slopes.input$WarmQSlope, na.rm = TRUE) - min(slopes.input$WarmQSlope, na.rm = TRUE))/50, "]")
  
  # Extract the prediction data frame
  SR.pred <- ggpredict(SR.mod, terms = range.to.predict)
  
  # Convert model output into a dataframe (with 4.d.p.)
  SR.df <- brms_SummaryTable(SR.mod, formatOptions = list(digits = 4, nsmall = 4), round = 4)
  
  # Extract the confidence intervals as a list for use in the plotting
  SR.df.ci <- SR.df %>% 
    filter(Covariate %in% c("x_variable")) %>% 
    dplyr::select("l-95% CI", "u-95% CI") %>% 
    pivot_longer(cols = 1:2, names_to = "Interval", values_to = "Value")
  
  # Save the confidence intervals as a list and remove the intermediate dataframe
  SR.ci <- as.list(SR.df.ci$Value) # Adding values to the list
  names(SR.ci) <- SR.df.ci$Interval # Adding names to the values
  rm(SR.df, SR.df.ci) # Remove unnecessary dataframes of summary and CIs
  
  # Plot model outputs
  (SR.plot <- ggplot(SR.pred) +
      
      # Hexplot
      geom_hex(data = slopes.input, aes(x = WarmQSlope, y = SR_slopes, bins = 18)) +
      scale_fill_gradient(low = "#ABE0EB", high = "#1C86EE") +
      geom_ribbon(aes(x = x, ymin = conf.low, ymax = conf.high),
                  fill = "#1C86EE", alpha = 0.2, colour = "#1C86EE", linewdith = 0.05) + # Adds c.intervals for predictions as ribbon
      geom_line(aes(x = x, y = predicted), colour = "#000000", linewidth = 1.1,
                linetype = ifelse((SR.ci$"l-95% CI" < 0 & SR.ci$"u-95% CI" < 0) | # Automatically asigns significance in subtitle based on CIs
                                    (SR.ci$"l-95% CI" > 0 & SR.ci$"u-95% CI" > 0), 1, 2)) + # Adds line for predicted values of SR vs LAT
      
      labs(y = bquote('Species Richness Change'~(yr^-1)),
           x = paste0("\nTemperature Change (", "\u00B0", "C yr", "\U207B", "\U00B9", ")"),
           title = paste0("CIs: ", SR.ci$"l-95% CI", " to ", SR.ci$"u-95% CI"),
           subtitle = paste(ifelse((SR.ci$"l-95% CI" < 0 & SR.ci$"u-95% CI" < 0) | # Automatically asigns significance in subtitle based on CIs
                                     (SR.ci$"l-95% CI" > 0 & SR.ci$"u-95% CI" > 0),
                                   "Significant: CIs do NOT span 0", "Not Significant: CIs span 0"))) +
      theme_1() +
      theme(legend.position = "none",
            plot.subtitle = element_text(colour = ifelse((SR.ci$"l-95% CI" < 0 & SR.ci$"u-95% CI" < 0) | # Automatically colours subtitle based on significance from CIs
                                                           (SR.ci$"l-95% CI" > 0 & SR.ci$"u-95% CI" > 0),
                                                         paste("#006400"), paste("#8b0000")))))
  
  
  # Load in model output
  FDis.mod <- get(load(paste0("data/model_outputs_new/m_FDis_slopes_",
                              "WarmQSlope", pc.filepath, ".RData")))
  
  # Range to predict variable
  range.to.predict <- paste0("x_variable [", min(slopes.input$WarmQSlope, na.rm = TRUE), ":", max(slopes.input$WarmQSlope, na.rm = TRUE), 
                             " by =", (max(slopes.input$WarmQSlope, na.rm = TRUE) - min(slopes.input$WarmQSlope, na.rm = TRUE))/50, "]")
  
  # Extract the prediction data frame
  FDis.pred <- ggpredict(FDis.mod, terms = range.to.predict)
  
  # Convert model output into a dataframe (with 4.d.p.)
  FDis.df <- brms_SummaryTable(FDis.mod, formatOptions = list(digits = 4, nsmall = 4), round = 4)
  
  # Extract the confidence intervals as a list for use in the plotting
  FDis.df.ci <- FDis.df %>% 
    filter(Covariate %in% c("x_variable")) %>% 
    dplyr::select("l-95% CI", "u-95% CI") %>% 
    pivot_longer(cols = 1:2, names_to = "Interval", values_to = "Value")
  
  # Save the confidence intervals as a list and remove the intermediate dataframe
  FDis.ci <- as.list(FDis.df.ci$Value) # Adding values to the list
  names(FDis.ci) <- FDis.df.ci$Interval # Adding names to the values
  rm(FDis.df, FDis.df.ci) # Remove unnecessary dataframes of summary and CIs
  
  # Plot model outputs
  (FDis.plot <- ggplot(FDis.pred) +
      
      # Hexplot
      geom_hex(data = slopes.input, aes(x = WarmQSlope, y = FDis_slopes, bins = 18)) +
      scale_fill_gradient(low = "#FFCF91", high = "#EE7600") +
      geom_ribbon(aes(x = x, ymin = conf.low, ymax = conf.high),
                  fill = "#EE7600", alpha = 0.2, colour = "#EE7600", linewdith = 0.05) + # Adds c.intervals for predictions as ribbon
      geom_line(aes(x = x, y = predicted), colour = "#000000", linewidth = 1.1,
                linetype = ifelse((FDis.ci$"l-95% CI" < 0 & FDis.ci$"u-95% CI" < 0) | # Automatically asigns significance in subtitle based on CIs
                                    (FDis.ci$"l-95% CI" > 0 & FDis.ci$"u-95% CI" > 0), 1, 2)) + # Adds line for predicted values of FDis vs LAT
      
      labs(y = bquote('Functional Dispersion Change'~(yr^-1)),
           x = paste0("\nTemperature Change (", "\u00B0", "C yr", "\U207B", "\U00B9", ")"),
           title = paste0("CIs: ", FDis.ci$"l-95% CI", " to ", FDis.ci$"u-95% CI"),
           subtitle = paste(ifelse((FDis.ci$"l-95% CI" < 0 & FDis.ci$"u-95% CI" < 0) | # Automatically asigns significance in subtitle based on CIs
                                     (FDis.ci$"l-95% CI" > 0 & FDis.ci$"u-95% CI" > 0),
                                   "Significant: CIs do NOT span 0", "Not Significant: CIs span 0"))) +
      theme_1() +
      theme(legend.position = "none",
            plot.subtitle = element_text(colour = ifelse((FDis.ci$"l-95% CI" < 0 & FDis.ci$"u-95% CI" < 0) | # Automatically colours subtitle based on significance from CIs
                                                           (FDis.ci$"l-95% CI" > 0 & FDis.ci$"u-95% CI" > 0),
                                                         paste("#006400"), paste("#8b0000")))))
  
  
  # Create a panel of all three plots
  combined.panel <- grid.arrange(FRic.plot, FEve.plot, FDis.plot, SR.plot, ncol = 4)
  
  # Export panel
  ggsave(combined.panel, filename = paste0("figures/outputs_new/manuscript_",
                                           "TempChange", pc.filepath, ".png"), width = 20, height = 5.5)
  
} # End of function


# EXTRA 1: PRECIP ----

# FIGURE 3b: GENERATE PrecipAnnITUDE PLOTS WITH COLOUR SCHEME ----

plot.precip <- function(censored){
  
  # Determine which version of FRic to import (censored or not)
  
  if (censored == "No"){
    
    # Load in model output
    FRic.mod <- get(load(paste0("data/model_outputs_new/m_FRic_",
                                FRic.distribution, "_PrecipAnn", pc.filepath, ".RData")))
    
    # Extract the prediction data frame
    FRic.pred <- ggpredict(FRic.mod, terms = "x_variable [0:2000, sample = 30]", back.transform = TRUE)
    
    # Convert model output into a dataframe (with 4.d.p.)
    FRic.df <- brms_SummaryTable(FRic.mod, formatOptions = list(digits = 4, nsmall = 4), round = 4)
    
    # Extract the confidence intervals as a list for use in the plotting
    FRic.df.ci <- FRic.df %>% 
      filter(Covariate %in% c("x_variable")) %>% 
      dplyr::select("l-95% CI", "u-95% CI") %>% 
      pivot_longer(cols = 1:2, names_to = "Interval", values_to = "Value")
    
  } # End of import and not censored
  
  if (censored == "Yes"){
    
    # Load in model output
    FRic.mod <- get(load(paste0("data/model_outputs_new/m_FRic_",
                                "censored_PrecipAnn", pc.filepath, ".RData")))
    
    # Extract the prediction data frame and exponentiate the outputs
    FRic.pred <- ggpredict(FRic.mod, terms = "x_variable [0:2000, sample = 30]", back.transform = FALSE) %>% 
      mutate(conf.low = exp(conf.low),
             conf.high = exp(conf.high),
             predicted = exp(predicted))
    
    # Convert model output into a dataframe (with 4.d.p.)
    FRic.df <- brms_SummaryTable(FRic.mod, formatOptions = list(digits = 4, nsmall = 4), round = 4)
    
    # Extract the confidence intervals as a list for use in the plotting
    FRic.df.ci <- FRic.df %>% 
      filter(Covariate %in% c("x_variable")) %>%
      dplyr::select("l-95% CI", "u-95% CI") %>% 
      pivot_longer(cols = 1:2, names_to = "Interval", values_to = "Value")
    
  } # End of import and censored
  
  # Save the confidence intervals as a list and remove the intermediate dataframe
  FRic.ci <- as.list(FRic.df.ci$Value) # Adding values to the list
  names(FRic.ci) <- FRic.df.ci$Interval # Adding names to the values
  rm(FRic.df, FRic.df.ci) # Remove unnecessary dataframes of summary and CIs
  
  # Plot model outputs
  (FRic.plot <- ggplot(FRic.pred) +
      
      # Scatterplot
      # geom_point(data = combo.latest, aes(x = PrecipAnn, y = FRic), # Adds original FRic vs PrecipAnn data points and colours by region
      #            fill = "#D02090", colour = c("#000000"), alpha = 0.2, shape = 21, size = 3) +
      # geom_ribbon(aes(x = x, ymin = conf.low, ymax = conf.high),
      #             fill = "#D02090", alpha = 0.27, colour = "#D02090", linewdith = 0.1) + # Adds c.intervals for predictions as ribbon
      # geom_line(aes(x = x, y = predicted), colour = "#D02090", linewidth = 1.1,
      #           linetype = ifelse((FRic.ci$"l-95% CI" < 0 & FRic.ci$"u-95% CI" < 0) | # Automatically asigns significance in subtitle based on CIs
      # (FRic.ci$"l-95% CI" > 0 & FRic.ci$"u-95% CI" > 0), 1, 2)) + # Adds line for predicted values of FRic vs PrecipAnn
      
      # Hexplot
    geom_hex(data = combo.latest, aes(x = PrecipAnn, y = FRic), bins = 18) +
      scale_fill_gradient(low = "#FFE1FF", high = "#D02090") +
      geom_ribbon(aes(x = x, ymin = conf.low, ymax = conf.high),
                  fill = "#D02090", alpha = 0.2, colour = "#D02090", linewdith = 0.05) + # Adds c.intervals for predictions as ribbon
      geom_line(aes(x = x, y = predicted), colour = "#000000", linewidth = 1.1,
                linetype = ifelse((FRic.ci$"l-95% CI" < 0 & FRic.ci$"u-95% CI" < 0) | # Automatically asigns significance in subtitle based on CIs
                                    (FRic.ci$"l-95% CI" > 0 & FRic.ci$"u-95% CI" > 0), 1, 2)) + # Adds line for predicted values of FRic vs PrecipAnn
      
      labs(y = "Functional Richness\n",
           x = paste0("\nPrecipitation (mm)"),
           title = paste0("CIs: ", FRic.ci$"l-95% CI", " to ", FRic.ci$"u-95% CI"),
           subtitle = paste(ifelse((FRic.ci$"l-95% CI" < 0 & FRic.ci$"u-95% CI" < 0) | # Automatically asigns significance in subtitle based on CIs
                                     (FRic.ci$"l-95% CI" > 0 & FRic.ci$"u-95% CI" > 0),
                                   "Significant: CIs do NOT span 0", "Not Significant: CIs span 0"))) +
      theme_1() +
      theme(legend.position = "none",
            plot.subtitle = element_text(colour = ifelse((FRic.ci$"l-95% CI" < 0 & FRic.ci$"u-95% CI" < 0) | # Automatically colours subtitle based on significance from CIs
                                                           (FRic.ci$"l-95% CI" > 0 & FRic.ci$"u-95% CI" > 0),
                                                         paste("#006400"), paste("#8b0000")))))
  
  
  # Load in model output
  FEve.mod <- get(load(paste0("data/model_outputs_new/m_FEve_",
                              "PrecipAnn", pc.filepath, ".RData")))
  
  # Convert model output into a dataframe (with 4.d.p.)
  FEve.df <- brms_SummaryTable(FEve.mod, formatOptions = list(digits = 4, nsmall = 4), round = 4)
  
  # Extract the confidence intervals as a list for use in the plotting
  FEve.df.ci <- FEve.df %>% 
    filter(Covariate %in% c("x_variable")) %>% 
    dplyr::select("l-95% CI", "u-95% CI") %>% 
    pivot_longer(cols = 1:2, names_to = "Interval", values_to = "Value")
  
  # Save the confidence intervals as a list and remove the intermediate dataframe
  FEve.ci <- as.list(FEve.df.ci$Value) # Adding values to the list
  names(FEve.ci) <- FEve.df.ci$Interval # Adding names to the values
  rm(FEve.df, FEve.df.ci) # Remove unnecessary dataframes of summary and CIs  
  
  # Extract the prediction data frame
  FEve.pred <- ggpredict(FEve.mod, terms = "x_variable [0:2000, sample = 30]", back.transform = FALSE)
  
  # Generate plot
  (FEve.plot <- ggplot(FEve.pred) +
      
      
      # # SCatterplot
      # geom_point(data = combo.latest, aes(x = PrecipAnn, y = FEve), # Adds original FEve vs PrecipAnn data points and colours by region
      #            fill = "#228B22", colour = c("#000000"), alpha = 0.2, shape = 21, size = 3) +
      # geom_ribbon(aes(x = x, ymin = conf.low, ymax = conf.high),
      #             fill = "#228B22", alpha = 0.27, colour = "#228B22", linewdith = 0.1) + # Adds c.intervals for predictions as ribbon
      # geom_line(aes(x = x, y = predicted), colour = "#228B22", linewidth = 1.1,
      #           linetype = ifelse((FEve.ci$"l-95% CI" < 0 & FEve.ci$"u-95% CI" < 0) | # Automatically asigns significance in subtitle based on CIs
      #                               (FEve.ci$"l-95% CI" > 0 & FEve.ci$"u-95% CI" > 0), 1, 2)) + # Adds line for predicted values of FEve vs PrecipAnn
      
    # Hexplot
    geom_hex(data = combo.latest, aes(x = PrecipAnn, y = FEve), bins = 18) +
      scale_fill_gradient(low = "#CFF099", high = "#228B22") +
      geom_ribbon(aes(x = x, ymin = conf.low, ymax = conf.high),
                  fill = "#228B22", alpha = 0.25, colour = "#228B22", linewdith = 0.05) + # Adds c.intervals for predictions as ribbon
      geom_line(aes(x = x, y = predicted), colour = "#000000", linewidth = 1.1,
                linetype = ifelse((FEve.ci$"l-95% CI" < 0 & FEve.ci$"u-95% CI" < 0) | # Automatically asigns significance in subtitle based on CIs
                                    (FEve.ci$"l-95% CI" > 0 & FEve.ci$"u-95% CI" > 0), 1, 2)) + # Adds line for predicted values of FRic vs PrecipAnn
      
      labs(y = "Functional Evenness\n",
           x = paste0("\nPrecipitation (mm)"),
           title = paste0("CIs: ", FEve.ci$"l-95% CI", " to ", FEve.ci$"u-95% CI"),
           subtitle = paste(ifelse((FEve.ci$"l-95% CI" < 0 & FEve.ci$"u-95% CI" < 0) | # Automatically asigns significance in subtitle based on CIs
                                     (FEve.ci$"l-95% CI" > 0 & FEve.ci$"u-95% CI" > 0),
                                   "Significant: CIs do NOT span 0", "Not Significant: CIs span 0"))) +
      theme_1() +
      theme(legend.position = "none",
            plot.subtitle = element_text(colour = ifelse((FEve.ci$"l-95% CI" < 0 & FEve.ci$"u-95% CI" < 0) | # Automatically colours subtitle based on significance from CIs
                                                           (FEve.ci$"l-95% CI" > 0 & FEve.ci$"u-95% CI" > 0),
                                                         paste("#006400"), paste("#8b0000")))))
  
  
  # Load in model output
  SR.mod <- get(load(paste0("data/model_outputs_new/m_SR_",
                            "PrecipAnn", pc.filepath, ".RData")))
  
  # Convert model output into a dataframe (with 4.d.p.)
  SR.df <- brms_SummaryTable(SR.mod, formatOptions = list(digits = 4, nsmall = 4), round = 4)
  
  # Extract the confidence intervals as a list for use in the plotting
  SR.df.ci <- SR.df %>% 
    filter(Covariate %in% c("x_variable")) %>% 
    dplyr::select("l-95% CI", "u-95% CI") %>% 
    pivot_longer(cols = 1:2, names_to = "Interval", values_to = "Value")
  
  # Save the confidence intervals as a list and remove the intermediate dataframe
  SR.ci <- as.list(SR.df.ci$Value) # Adding values to the list
  names(SR.ci) <- SR.df.ci$Interval # Adding names to the values
  rm(SR.df, SR.df.ci) # Remove unnecessary dataframes of summary and CIs  
  
  # Extract the prediction data frame
  SR.pred <- ggpredict(SR.mod, terms = "x_variable [0:2000, sample = 30]", back.transform = FALSE)
  
  # Generate plot
  (SR.plot <- ggplot(SR.pred) +
      
      # Scatterplot
      # geom_point(data = combo.latest, aes(x = PrecipAnn, y = SR), # Adds original SR vs PrecipAnn data points and colours by region
      #            fill = "#1C86EE", colour = c("#000000"), alpha = 0.2, shape = 21, size = 3) +
      # geom_ribbon(aes(x = x, ymin = conf.low, ymax = conf.high),
      #             fill = "#1C86EE", alpha = 0.27, colour = "#1C86EE", linewdith = 0.1) + # Adds c.intervals for predictions as ribbon
      # geom_line(aes(x = x, y = predicted), colour = "#1C86EE", linewidth = 1.1,
      #           linetype = ifelse((SR.ci$"l-95% CI" < 0 & SR.ci$"u-95% CI" < 0) | # Automatically asigns significance in subtitle based on CIs
      #                               (SR.ci$"l-95% CI" > 0 & SR.ci$"u-95% CI" > 0), 1, 2)) + # Adds line for predicted values of SR vs PrecipAnn
      
      # Hexplot
    geom_hex(data = combo.latest, aes(x = PrecipAnn, y = SR), bins = 18) +
      scale_fill_gradient(low = "#ABE0EB", high = "#1C86EE") +
      geom_ribbon(aes(x = x, ymin = conf.low, ymax = conf.high),
                  fill = "#1C86EE", alpha = 0.25, colour = "#1C86EE", linewdith = 0.05) + # Adds c.intervals for predictions as ribbon
      geom_line(aes(x = x, y = predicted), colour = "#000000", linewidth = 1.1,
                linetype = ifelse((SR.ci$"l-95% CI" < 0 & SR.ci$"u-95% CI" < 0) | # Automatically asigns significance in subtitle based on CIs
                                    (SR.ci$"l-95% CI" > 0 & SR.ci$"u-95% CI" > 0), 1, 2)) + # Adds line for predicted values of SR vs PrecipAnn
      
      labs(y = "Species Richness\n",
           x = paste0("\nPrecipitation (mm)"),
           title = paste0("CIs: ", SR.ci$"l-95% CI", " to ", SR.ci$"u-95% CI"),
           subtitle = paste(ifelse((SR.ci$"l-95% CI" < 0 & SR.ci$"u-95% CI" < 0) | # Automatically asigns significance in subtitle based on CIs
                                     (SR.ci$"l-95% CI" > 0 & SR.ci$"u-95% CI" > 0),
                                   "Significant: CIs do NOT span 0", "Not Significant: CIs span 0"))) +
      theme_1() +
      theme(legend.position = "none",
            plot.subtitle = element_text(colour = ifelse((SR.ci$"l-95% CI" < 0 & SR.ci$"u-95% CI" < 0) | # Automatically colours subtitle based on significance from CIs
                                                           (SR.ci$"l-95% CI" > 0 & SR.ci$"u-95% CI" > 0),
                                                         paste("#006400"), paste("#8b0000")))))
  
  
  # Load in model output
  FDis.mod <- get(load(paste0("data/model_outputs_new/m_FDis_",
                              "PrecipAnn", pc.filepath, ".RData")))
  
  # Convert model output into a dataframe (with 4.d.p.)
  FDis.df <- brms_SummaryTable(FDis.mod, formatOptions = list(digits = 4, nsmall = 4), round = 4)
  
  # Extract the confidence intervals as a list for use in the plotting
  FDis.df.ci <- FDis.df %>% 
    filter(Covariate %in% c("x_variable")) %>% 
    dplyr::select("l-95% CI", "u-95% CI") %>% 
    pivot_longer(cols = 1:2, names_to = "Interval", values_to = "Value")
  
  # Save the confidence intervals as a list and remove the intermediate dataframe
  FDis.ci <- as.list(FDis.df.ci$Value) # Adding values to the list
  names(FDis.ci) <- FDis.df.ci$Interval # Adding names to the values
  rm(FDis.df, FDis.df.ci) # Remove unnecessary dataframes of summary and CIs  
  
  # Extract the prediction data frame
  FDis.pred <- ggpredict(FDis.mod, terms = "x_variable [0:2000, sample = 30]", back.transform = FALSE)
  
  # Generate plot
  (FDis.plot <- ggplot(FDis.pred) +
      
      
      # # SCatterplot
      # geom_point(data = combo.latest, aes(x = PrecipAnn, y = FDis), # Adds original FDis vs PrecipAnn data points and colours by region
      #            fill = "#228B22", colour = c("#000000"), alpha = 0.2, shape = 21, size = 3) +
      # geom_ribbon(aes(x = x, ymin = conf.low, ymax = conf.high),
      #             fill = "#228B22", alpha = 0.27, colour = "#228B22", linewdith = 0.1) + # Adds c.intervals for predictions as ribbon
      # geom_line(aes(x = x, y = predicted), colour = "#228B22", linewidth = 1.1,
      #           linetype = ifelse((FDis.ci$"l-95% CI" < 0 & FDis.ci$"u-95% CI" < 0) | # Automatically asigns significance in subtitle based on CIs
      #                               (FDis.ci$"l-95% CI" > 0 & FDis.ci$"u-95% CI" > 0), 1, 2)) + # Adds line for predicted values of FDis vs PrecipAnn
      
    # Hexplot
    geom_hex(data = combo.latest, aes(x = PrecipAnn, y = FDis), bins = 18) +
      scale_fill_gradient(low = "#FFCF91", high = "#EE7600") +
      geom_ribbon(aes(x = x, ymin = conf.low, ymax = conf.high),
                  fill = "#EE7600", alpha = 0.25, colour = "#EE7600", linewdith = 0.05) + # Adds c.intervals for predictions as ribbon
      geom_line(aes(x = x, y = predicted), colour = "#000000", linewidth = 1.1,
                linetype = ifelse((FDis.ci$"l-95% CI" < 0 & FDis.ci$"u-95% CI" < 0) | # Automatically asigns significance in subtitle based on CIs
                                    (FDis.ci$"l-95% CI" > 0 & FDis.ci$"u-95% CI" > 0), 1, 2)) + # Adds line for predicted values of FRic vs PrecipAnn
      
      labs(y = "Functional Dispersion\n",
           x = paste0("\nPrecipitation (mm)"),
           title = paste0("CIs: ", FDis.ci$"l-95% CI", " to ", FDis.ci$"u-95% CI"),
           subtitle = paste(ifelse((FDis.ci$"l-95% CI" < 0 & FDis.ci$"u-95% CI" < 0) | # Automatically asigns significance in subtitle based on CIs
                                     (FDis.ci$"l-95% CI" > 0 & FDis.ci$"u-95% CI" > 0),
                                   "Significant: CIs do NOT span 0", "Not Significant: CIs span 0"))) +
      theme_1() +
      theme(legend.position = "none",
            plot.subtitle = element_text(colour = ifelse((FDis.ci$"l-95% CI" < 0 & FDis.ci$"u-95% CI" < 0) | # Automatically colours subtitle based on significance from CIs
                                                           (FDis.ci$"l-95% CI" > 0 & FDis.ci$"u-95% CI" > 0),
                                                         paste("#006400"), paste("#8b0000")))))
  
  
  # Create a panel of all three plots
  combined.panel <- grid.arrange(FRic.plot, FEve.plot, FDis.plot, SR.plot, ncol = 2)
  
  # Add in if statements for saving the output
  if (censored == "No"){
    
    # Export panel
    ggsave(combined.panel, filename = paste0("figures/outputs_new/manuscript_",
                                             "PrecipAnn", pc.filepath, ".png"), width = 12, height = 12)
    
  }
  
  # Add in if statements for saving the output
  if (censored == "Yes"){
    
    # Export panel
    ggsave(combined.panel, filename = paste0("figures/outputs_new/manuscript_",
                                             "PrecipAnn_censored", pc.filepath, ".png"), width = 12, height = 12)
    
  }
  
}


# EXTRA 2: SOIL MOISTURE ----

plot.moisture <- function(){
  
  # FRic ----
  
  # Load in model output
  FRic.mod <- get(load(paste0("data/model_outputs_new/m_FRic_censored_MOISTURE",
                              pc.filepath, ".RData")))
  
  summary(FRic.mod)
  
  # Extract the prediction data frame
  FRic.pred <- ggpredict(FRic.mod, terms = "x_variable", back.transform = FALSE)
  
  # Plot the model outputs
  (FRic.plot <- ggplot(FRic.pred) +
      geom_errorbar(aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high), width = 0.2) +
      geom_point(aes(x = x, y = predicted), stat = "identity", fill = "#D02090",
                 colour = c("#000000"), alpha = 0.9, shape = 21, size = 5) +
      labs(y = "FRic \n",
           x = "\n Soil Moisture") +
      coord_flip() +
      theme_1() +
      theme(legend.position = "none",
            plot.title = element_text(size = 12),
            plot.subtitle = element_text(size = 6)))
  
  
  # FEve ----
  
  # Load in model output
  FEve.mod <- get(load(paste0("data/model_outputs_new/m_FEve_MOISTURE",
                              pc.filepath, ".RData")))
  
  summary(FEve.mod)
  
  # Extract the prediction data frame
  FEve.pred <- ggpredict(FEve.mod, terms = "x_variable", back.transform = FALSE)
  
  # Plot the model outputs
  (FEve.plot <- ggplot(FEve.pred) +
      geom_errorbar(aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high), width = 0.2) +
      geom_point(aes(x = x, y = predicted), stat = "identity", fill = "#228B22",
                 colour = c("#000000"), alpha = 0.9, shape = 21, size = 5) +
      labs(y = "FEve \n",
           x = "\n Soil Moisture") +
      coord_flip() +
      theme_1() +
      theme(legend.position = "none",
            plot.title = element_text(size = 12),
            plot.subtitle = element_text(size = 6)))
  
  # SR ----
  
  # Load in model output
  SR.mod <- get(load(paste0("data/model_outputs_new/m_SR_MOISTURE",
                            pc.filepath, ".RData")))
  
  summary(SR.mod)
  
  # Extract the prediction data frame
  SR.pred <- ggpredict(SR.mod, terms = "x_variable", back.transform = FALSE)
  
  # Plot the model outputs
  (SR.plot <- ggplot(SR.pred) +
      geom_errorbar(aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high), width = 0.2) +
      geom_point(aes(x = x, y = predicted), stat = "identity", fill = "#1C86EE",
                 colour = c("#000000"), alpha = 0.9, shape = 21, size = 5) +
      labs(y = "SR \n",
           x = "\n Soil Moisture") +
      coord_flip() +
      theme_1() +
      theme(legend.position = "none",
            plot.title = element_text(size = 12),
            plot.subtitle = element_text(size = 6)))
  
  # FDis ----
  
  # Load in model output
  FDis.mod <- get(load(paste0("data/model_outputs_new/m_FDis_MOISTURE",
                              pc.filepath, ".RData")))
  
  summary(FDis.mod)
  
  # Extract the prediction data frame
  FDis.pred <- ggpredict(FDis.mod, terms = "x_variable", back.transform = FALSE)
  
  # Plot the model outputs
  (FDis.plot <- ggplot(FDis.pred) +
      geom_errorbar(aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high), width = 0.2) +
      geom_point(aes(x = x, y = predicted), stat = "identity", fill = "#EE7600",
                 colour = c("#000000"), alpha = 0.9, shape = 21, size = 5) +
      labs(y = "FDis \n",
           x = "\n Soil Moisture") +
      coord_flip() +
      theme_1() +
      theme(legend.position = "none",
            plot.title = element_text(size = 12),
            plot.subtitle = element_text(size = 6)))
  
  # Create a panel of all three plots
  combined.panel <- grid.arrange(FRic.plot, FEve.plot, FDis.plot, SR.plot, ncol = 2)
  
  # Export panel
  ggsave(combined.panel, filename = paste0("figures/outputs_new/combined_",
                                           "MOISTURE_censored", pc.filepath, ".png"), width = 10, height = 10)
  
  
}


# EXTRA 3: REGION ----

plot.region <- function(){
  
  # FRic ----
  
  # Load in model output
  FRic.mod <- get(load(paste0("data/model_outputs_new/m_FRic_censored_Region",
                              pc.filepath, ".RData")))
  
  summary(FRic.mod)
  
  # Extract the prediction data frame
  FRic.pred <- ggpredict(FRic.mod, terms = "x_variable", back.transform = FALSE)
  
  # Plot the model outputs
  (FRic.plot <- ggplot(FRic.pred) +
      geom_errorbar(aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high), width = 0.2) +
      geom_point(aes(x = x, y = predicted), stat = "identity", fill = "#D02090",
                 colour = c("#000000"), alpha = 0.9, shape = 21, size = 5) +
      labs(y = "FRic \n",
           x = "\n Region") +
      coord_flip() +
      theme_1() +
      theme(legend.position = "none",
            plot.title = element_text(size = 12),
            plot.subtitle = element_text(size = 6)))
  
  
  # FEve ----
  
  # Load in model output
  FEve.mod <- get(load(paste0("data/model_outputs_new/m_FEve_Region",
                              pc.filepath, ".RData")))
  
  summary(FEve.mod)
  
  # Extract the prediction data frame
  FEve.pred <- ggpredict(FEve.mod, terms = "x_variable", back.transform = FALSE)
  
  # Plot the model outputs
  (FEve.plot <- ggplot(FEve.pred) +
      geom_errorbar(aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high), width = 0.2) +
      geom_point(aes(x = x, y = predicted), stat = "identity", fill = "#228B22",
                 colour = c("#000000"), alpha = 0.9, shape = 21, size = 5) +
      labs(y = "FEve \n",
           x = "\n Region") +
      coord_flip() +
      theme_1() +
      theme(legend.position = "none",
            plot.title = element_text(size = 12),
            plot.subtitle = element_text(size = 6)))
  
  # SR ----
  
  # Load in model output
  SR.mod <- get(load(paste0("data/model_outputs_new/m_SR_Region",
                            pc.filepath, ".RData")))
  
  summary(SR.mod)
  
  # Extract the prediction data frame
  SR.pred <- ggpredict(SR.mod, terms = "x_variable", back.transform = FALSE)
  
  # Plot the model outputs
  (SR.plot <- ggplot(SR.pred) +
      geom_errorbar(aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high), width = 0.2) +
      geom_point(aes(x = x, y = predicted), stat = "identity", fill = "#1C86EE",
                 colour = c("#000000"), alpha = 0.9, shape = 21, size = 5) +
      labs(y = "SR \n",
           x = "\n Region") +
      coord_flip() +
      theme_1() +
      theme(legend.position = "none",
            plot.title = element_text(size = 12),
            plot.subtitle = element_text(size = 6)))
  
  # FDis ----
  
  # Load in model output
  FDis.mod <- get(load(paste0("data/model_outputs_new/m_FDis_Region",
                              pc.filepath, ".RData")))
  
  summary(FDis.mod)
  
  # Extract the prediction data frame
  FDis.pred <- ggpredict(FDis.mod, terms = "x_variable", back.transform = FALSE)
  
  # Plot the model outputs
  (FDis.plot <- ggplot(FDis.pred) +
      geom_errorbar(aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high), width = 0.2) +
      geom_point(aes(x = x, y = predicted), stat = "identity", fill = "#EE7600",
                 colour = c("#000000"), alpha = 0.9, shape = 21, size = 5) +
      labs(y = "FDis \n",
           x = "\n Region") +
      coord_flip() +
      theme_1() +
      theme(legend.position = "none",
            plot.title = element_text(size = 12),
            plot.subtitle = element_text(size = 6)))
  
  # Create a panel of all three plots
  combined.panel <- grid.arrange(FRic.plot, FEve.plot, FDis.plot, SR.plot, ncol = 2)
  
  # Export panel
  ggsave(combined.panel, filename = paste0("figures/outputs_new/combined_",
                                           "Region_censored", pc.filepath, ".png"), width = 15, height = 10)
  
  
}


# EXTRA 4: LATITUDE TEMPORAL ----

plot.latitude.time <- function(){
  
  # Produce dataframe for running models on
  input.data <- slopes.input %>% 
    dplyr::select(LAT, FRic_slopes, FEve_slopes, SR_slopes, FDis_slopes,
                  SurveyedArea, gridcell, SUBSITE, PlotDominatingFG)
  
  # FRic ----
  
  # Load in model output
  FRic.mod <- get(load(paste0("data/model_outputs_new/m_FRic_slopes_LAT",
                              pc.filepath, ".RData")))
  
  # Convert model output into a dataframe (with 4.d.p.)
  FRic.df <- brms_SummaryTable(FRic.mod, formatOptions = list(digits = 5, nsmall = 5), round = 5)
  
  # Extract the confidence intervals as a list for use in the plotting
  FRic.df.ci <- FRic.df %>% 
    filter(Covariate %in% c("x_variable")) %>% 
    dplyr::select("l-95% CI", "u-95% CI") %>% 
    pivot_longer(cols = 1:2, names_to = "Interval", values_to = "Value")
  
  # Save the confidence intervals as a list and remove the intermediate dataframe
  FRic.ci <- as.list(FRic.df.ci$Value) # Adding values to the list
  names(FRic.ci) <- FRic.df.ci$Interval # Adding names to the values
  rm(FRic.df, FRic.df.ci) # Remove unnecessary dataframes of summary and CIs  
  
  # Range to predict variable
  range.to.predict <- paste0("x_variable [", min(input.data$LAT, na.rm = TRUE), ":", max(input.data$LAT, na.rm = TRUE), 
                             " by =", (max(input.data$LAT, na.rm = TRUE) - min(input.data$LAT, na.rm = TRUE))/50, "]")
  
  # Extract the prediction data frame
  FRic.pred <- ggpredict(FRic.mod, terms = range.to.predict)
  
  # Extract dataframe with only the x-variable for plotting
  data.x.var.only <- dplyr::select(input.data, LAT)
  
  # Plot model outputs
  (FRic.plot <- ggplot(FRic.pred) +
      geom_line(aes(x = x, y = predicted), colour = "grey10") + # Adds line for predicted values of FRic vs LAT
      geom_ribbon(aes(x = x, ymin = conf.low, ymax = conf.high), fill = "grey10", alpha = 0.2) + # Adds c.intervals for predictions as ribbon
      geom_point(data = input.data, aes(x = data.x.var.only[,1], y = FRic_slopes), # Adds original FRic vs LAT data points and colours by region
                 fill = "#D02090", colour = c("#000000"), alpha = 0.5, shape = 21, size = 3) +
      labs(y = "Change in FRic \n",
           x = paste0("\n Latitude"),
           title = paste0("CIs: ", FRic.ci$"l-95% CI", " to ", FRic.ci$"u-95% CI"),
           subtitle = paste(ifelse((FRic.ci$"l-95% CI" < 0 & FRic.ci$"u-95% CI" < 0) | # Automatically asigns significance in subtitle based on CIs
                                     (FRic.ci$"l-95% CI" > 0 & FRic.ci$"u-95% CI" > 0),
                                   "Significant: CIs do NOT span 0", "Not Significant: CIs span 0"))) +
      theme_1() +
      theme(legend.position = "bottom",
            plot.subtitle = element_text(colour = ifelse((FRic.ci$"l-95% CI" < 0 & FRic.ci$"u-95% CI" < 0) | # Automatically colours subtitle based on significance from CIs
                                                           (FRic.ci$"l-95% CI" > 0 & FRic.ci$"u-95% CI" > 0),
                                                         paste("#006400"), paste("#8b0000")))))
  
  # FEve ----
  
  # Load in model output
  FEve.mod <- get(load(paste0("data/model_outputs_new/m_FEve_slopes_LAT",
                              pc.filepath, ".RData")))
  
  # Convert model output into a dataframe (with 4.d.p.)
  FEve.df <- brms_SummaryTable(FEve.mod, formatOptions = list(digits = 5, nsmall = 5), round = 5)
  
  # Extract the confidence intervals as a list for use in the plotting
  FEve.df.ci <- FEve.df %>% 
    filter(Covariate %in% c("x_variable")) %>% 
    dplyr::select("l-95% CI", "u-95% CI") %>% 
    pivot_longer(cols = 1:2, names_to = "Interval", values_to = "Value")
  
  # Save the confidence intervals as a list and remove the intermediate dataframe
  FEve.ci <- as.list(FEve.df.ci$Value) # Adding values to the list
  names(FEve.ci) <- FEve.df.ci$Interval # Adding names to the values
  rm(FEve.df, FEve.df.ci) # Remove unnecessary dataframes of summary and CIs  
  
  # Range to predict variable
  range.to.predict <- paste0("x_variable [", min(input.data$LAT, na.rm = TRUE), ":", max(input.data$LAT, na.rm = TRUE), 
                             " by =", (max(input.data$LAT, na.rm = TRUE) - min(input.data$LAT, na.rm = TRUE))/50, "]")
  
  # Extract the prediction data frame
  FEve.pred <- ggpredict(FEve.mod, terms = range.to.predict)
  
  # Extract dataframe with only the x-variable for plotting
  data.x.var.only <- dplyr::select(input.data, LAT)
  
  # Plot model outputs
  (FEve.plot <- ggplot(FEve.pred) +
      geom_line(aes(x = x, y = predicted), colour = "grey10") + # Adds line for predicted values of FEve vs LAT
      geom_ribbon(aes(x = x, ymin = conf.low, ymax = conf.high), fill = "grey10", alpha = 0.2) + # Adds c.intervals for predictions as ribbon
      geom_point(data = input.data, aes(x = data.x.var.only[,1], y = FEve_slopes), # Adds original FEve vs LAT data points and colours by region
                 fill = "#228B22", colour = c("#000000"), alpha = 0.5, shape = 21, size = 3) +
      labs(y = "Change in FEve \n",
           x = paste0("\n Latitude"),
           title = paste0("CIs: ", FEve.ci$"l-95% CI", " to ", FEve.ci$"u-95% CI"),
           subtitle = paste(ifelse((FEve.ci$"l-95% CI" < 0 & FEve.ci$"u-95% CI" < 0) | # Automatically asigns significance in subtitle based on CIs
                                     (FEve.ci$"l-95% CI" > 0 & FEve.ci$"u-95% CI" > 0),
                                   "Significant: CIs do NOT span 0", "Not Significant: CIs span 0"))) +
      theme_1() +
      theme(legend.position = "bottom",
            plot.subtitle = element_text(colour = ifelse((FEve.ci$"l-95% CI" < 0 & FEve.ci$"u-95% CI" < 0) | # Automatically colours subtitle based on significance from CIs
                                                           (FEve.ci$"l-95% CI" > 0 & FEve.ci$"u-95% CI" > 0),
                                                         paste("#006400"), paste("#8b0000")))))
  
  # SR ----
  
  # Load in model output
  SR.mod <- get(load(paste0("data/model_outputs_new/m_SR_slopes_LAT",
                            pc.filepath, ".RData")))
  
  # Convert model output into a dataframe (with 4.d.p.)
  SR.df <- brms_SummaryTable(SR.mod, formatOptions = list(digits = 5, nsmall = 5), round = 5)
  
  # Extract the confidence intervals as a list for use in the plotting
  SR.df.ci <- SR.df %>% 
    filter(Covariate %in% c("x_variable")) %>% 
    dplyr::select("l-95% CI", "u-95% CI") %>% 
    pivot_longer(cols = 1:2, names_to = "Interval", values_to = "Value")
  
  # Save the confidence intervals as a list and remove the intermediate dataframe
  SR.ci <- as.list(SR.df.ci$Value) # Adding values to the list
  names(SR.ci) <- SR.df.ci$Interval # Adding names to the values
  rm(SR.df, SR.df.ci) # Remove unnecessary dataframes of summary and CIs  
  
  # Range to predict variable
  range.to.predict <- paste0("x_variable [", min(input.data$LAT, na.rm = TRUE), ":", max(input.data$LAT, na.rm = TRUE), 
                             " by =", (max(input.data$LAT, na.rm = TRUE) - min(input.data$LAT, na.rm = TRUE))/50, "]")
  
  # Extract the prediction data frame
  SR.pred <- ggpredict(SR.mod, terms = range.to.predict)
  
  # Extract dataframe with only the x-variable for plotting
  data.x.var.only <- dplyr::select(input.data, LAT)
  
  # Plot model outputs
  (SR.plot <- ggplot(SR.pred) +
      geom_line(aes(x = x, y = predicted), colour = "grey10") + # Adds line for predicted values of SR vs LAT
      geom_ribbon(aes(x = x, ymin = conf.low, ymax = conf.high), fill = "grey10", alpha = 0.2) + # Adds c.intervals for predictions as ribbon
      geom_point(data = input.data, aes(x = data.x.var.only[,1], y = SR_slopes), # Adds original SR vs LAT data points and colours by region
                 fill = "#1C86EE", colour = c("#000000"), alpha = 0.5, shape = 21, size = 3) +
      labs(y = "Change in SR \n",
           x = paste0("\n Latitude"),
           title = paste0("CIs: ", SR.ci$"l-95% CI", " to ", SR.ci$"u-95% CI"),
           subtitle = paste(ifelse((SR.ci$"l-95% CI" < 0 & SR.ci$"u-95% CI" < 0) | # Automatically asigns significance in subtitle based on CIs
                                     (SR.ci$"l-95% CI" > 0 & SR.ci$"u-95% CI" > 0),
                                   "Significant: CIs do NOT span 0", "Not Significant: CIs span 0"))) +
      theme_1() +
      theme(legend.position = "bottom",
            plot.subtitle = element_text(colour = ifelse((SR.ci$"l-95% CI" < 0 & SR.ci$"u-95% CI" < 0) | # Automatically colours subtitle based on significance from CIs
                                                           (SR.ci$"l-95% CI" > 0 & SR.ci$"u-95% CI" > 0),
                                                         paste("#006400"), paste("#8b0000")))))
  
  # FDis ----
  
  # Load in model output
  FDis.mod <- get(load(paste0("data/model_outputs_new/m_FDis_slopes_LAT",
                              pc.filepath, ".RData")))
  
  # Convert model output into a dataframe (with 4.d.p.)
  FDis.df <- brms_SummaryTable(FDis.mod, formatOptions = list(digits = 5, nsmall = 5), round = 5)
  
  # Extract the confidence intervals as a list for use in the plotting
  FDis.df.ci <- FDis.df %>% 
    filter(Covariate %in% c("x_variable")) %>% 
    dplyr::select("l-95% CI", "u-95% CI") %>% 
    pivot_longer(cols = 1:2, names_to = "Interval", values_to = "Value")
  
  # Save the confidence intervals as a list and remove the intermediate dataframe
  FDis.ci <- as.list(FDis.df.ci$Value) # Adding values to the list
  names(FDis.ci) <- FDis.df.ci$Interval # Adding names to the values
  rm(FDis.df, FDis.df.ci) # Remove unnecessary dataframes of summary and CIs  
  
  # Range to predict variable
  range.to.predict <- paste0("x_variable [", min(input.data$LAT, na.rm = TRUE), ":", max(input.data$LAT, na.rm = TRUE), 
                             " by =", (max(input.data$LAT, na.rm = TRUE) - min(input.data$LAT, na.rm = TRUE))/50, "]")
  
  # Extract the prediction data frame
  FDis.pred <- ggpredict(FDis.mod, terms = range.to.predict)
  
  # Extract dataframe with only the x-variable for plotting
  data.x.var.only <- dplyr::select(input.data, LAT)
  
  # Plot model outputs
  (FDis.plot <- ggplot(FDis.pred) +
      geom_line(aes(x = x, y = predicted), colour = "grey10") + # Adds line for predicted values of FDis vs LAT
      geom_ribbon(aes(x = x, ymin = conf.low, ymax = conf.high), fill = "grey10", alpha = 0.2) + # Adds c.intervals for predictions as ribbon
      geom_point(data = input.data, aes(x = data.x.var.only[,1], y = FDis_slopes), # Adds original FDis vs LAT data points and colours by region
                 fill = "#EE7600", colour = c("#000000"), alpha = 0.5, shape = 21, size = 3) +
      labs(y = "Change in FDis \n",
           x = paste0("\n Latitude"),
           title = paste0("CIs: ", FDis.ci$"l-95% CI", " to ", FDis.ci$"u-95% CI"),
           subtitle = paste(ifelse((FDis.ci$"l-95% CI" < 0 & FDis.ci$"u-95% CI" < 0) | # Automatically asigns significance in subtitle based on CIs
                                     (FDis.ci$"l-95% CI" > 0 & FDis.ci$"u-95% CI" > 0),
                                   "Significant: CIs do NOT span 0", "Not Significant: CIs span 0"))) +
      theme_1() +
      theme(legend.position = "bottom",
            plot.subtitle = element_text(colour = ifelse((FDis.ci$"l-95% CI" < 0 & FDis.ci$"u-95% CI" < 0) | # Automatically colours subtitle based on significance from CIs
                                                           (FDis.ci$"l-95% CI" > 0 & FDis.ci$"u-95% CI" > 0),
                                                         paste("#006400"), paste("#8b0000")))))
  
  # Create a panel of all three plots
  combined.panel <- grid.arrange(FRic.plot, FEve.plot, FDis.plot, SR.plot, ncol = 2)
  
  # Export panel
  ggsave(combined.panel, filename = paste0("figures/outputs_new/combined_slopes_LAT",
                                           pc.filepath, ".png"), width = 15, height = 14)
  
}


# EXTRA 5: TEMPERATURE TEMPORAL ----

plot.temperature.time <- function(){
  
  # Produce dataframe for running models on
  input.data <- slopes.input %>% 
    dplyr::select(TempAvSum, FRic_slopes, FEve_slopes, SR_slopes, FDis_slopes,
                  SurveyedArea, gridcell, SUBSITE, PlotDominatingFG)
  
  # FRic ----
  
  # Load in model output
  FRic.mod <- get(load(paste0("data/model_outputs_new/m_FRic_slopes_TempAvSum",
                              pc.filepath, ".RData")))
  
  # Convert model output into a dataframe (with 4.d.p.)
  FRic.df <- brms_SummaryTable(FRic.mod, formatOptions = list(digits = 5, nsmall = 5), round = 5)
  
  # Extract the confidence intervals as a list for use in the plotting
  FRic.df.ci <- FRic.df %>% 
    filter(Covariate %in% c("x_variable")) %>% 
    dplyr::select("l-95% CI", "u-95% CI") %>% 
    pivot_longer(cols = 1:2, names_to = "Interval", values_to = "Value")
  
  # Save the confidence intervals as a list and remove the intermediate dataframe
  FRic.ci <- as.list(FRic.df.ci$Value) # Adding values to the list
  names(FRic.ci) <- FRic.df.ci$Interval # Adding names to the values
  rm(FRic.df, FRic.df.ci) # Remove unnecessary dataframes of summary and CIs  
  
  # Range to predict variable
  range.to.predict <- paste0("x_variable [", min(input.data$TempAvSum, na.rm = TRUE), ":", max(input.data$TempAvSum, na.rm = TRUE), 
                             " by =", (max(input.data$TempAvSum, na.rm = TRUE) - min(input.data$TempAvSum, na.rm = TRUE))/50, "]")
  
  # Extract the prediction data frame
  FRic.pred <- ggpredict(FRic.mod, terms = range.to.predict)
  
  # Extract dataframe with only the x-variable for plotting
  data.x.var.only <- dplyr::select(input.data, TempAvSum)
  
  # Plot model outputs
  (FRic.plot <- ggplot(FRic.pred) +
      geom_line(aes(x = x, y = predicted), colour = "grey10") + # Adds line for predicted values of FRic vs TempAvSum
      geom_ribbon(aes(x = x, ymin = conf.low, ymax = conf.high), fill = "grey10", alpha = 0.2) + # Adds c.intervals for predictions as ribbon
      geom_point(data = input.data, aes(x = data.x.var.only[,1], y = FRic_slopes), # Adds original FRic vs TempAvSum data points and colours by region
                 fill = "#D02090", colour = c("#000000"), alpha = 0.5, shape = 21, size = 3) +
      labs(y = "Change in FRic \n",
           x = paste0("\n Temperature"),
           title = paste0("CIs: ", FRic.ci$"l-95% CI", " to ", FRic.ci$"u-95% CI"),
           subtitle = paste(ifelse((FRic.ci$"l-95% CI" < 0 & FRic.ci$"u-95% CI" < 0) | # Automatically asigns significance in subtitle based on CIs
                                     (FRic.ci$"l-95% CI" > 0 & FRic.ci$"u-95% CI" > 0),
                                   "Significant: CIs do NOT span 0", "Not Significant: CIs span 0"))) +
      theme_1() +
      theme(legend.position = "bottom",
            plot.subtitle = element_text(colour = ifelse((FRic.ci$"l-95% CI" < 0 & FRic.ci$"u-95% CI" < 0) | # Automatically colours subtitle based on significance from CIs
                                                           (FRic.ci$"l-95% CI" > 0 & FRic.ci$"u-95% CI" > 0),
                                                         paste("#006400"), paste("#8b0000")))))
  
  # FEve ----
  
  # Load in model output
  FEve.mod <- get(load(paste0("data/model_outputs_new/m_FEve_slopes_TempAvSum",
                              pc.filepath, ".RData")))
  
  # Convert model output into a dataframe (with 4.d.p.)
  FEve.df <- brms_SummaryTable(FEve.mod, formatOptions = list(digits = 5, nsmall = 5), round = 5)
  
  # Extract the confidence intervals as a list for use in the plotting
  FEve.df.ci <- FEve.df %>% 
    filter(Covariate %in% c("x_variable")) %>% 
    dplyr::select("l-95% CI", "u-95% CI") %>% 
    pivot_longer(cols = 1:2, names_to = "Interval", values_to = "Value")
  
  # Save the confidence intervals as a list and remove the intermediate dataframe
  FEve.ci <- as.list(FEve.df.ci$Value) # Adding values to the list
  names(FEve.ci) <- FEve.df.ci$Interval # Adding names to the values
  rm(FEve.df, FEve.df.ci) # Remove unnecessary dataframes of summary and CIs  
  
  # Range to predict variable
  range.to.predict <- paste0("x_variable [", min(input.data$TempAvSum, na.rm = TRUE), ":", max(input.data$TempAvSum, na.rm = TRUE), 
                             " by =", (max(input.data$TempAvSum, na.rm = TRUE) - min(input.data$TempAvSum, na.rm = TRUE))/50, "]")
  
  # Extract the prediction data frame
  FEve.pred <- ggpredict(FEve.mod, terms = range.to.predict)
  
  # Extract dataframe with only the x-variable for plotting
  data.x.var.only <- dplyr::select(input.data, TempAvSum)
  
  # Plot model outputs
  (FEve.plot <- ggplot(FEve.pred) +
      geom_line(aes(x = x, y = predicted), colour = "grey10") + # Adds line for predicted values of FEve vs TempAvSum
      geom_ribbon(aes(x = x, ymin = conf.low, ymax = conf.high), fill = "grey10", alpha = 0.2) + # Adds c.intervals for predictions as ribbon
      geom_point(data = input.data, aes(x = data.x.var.only[,1], y = FEve_slopes), # Adds original FEve vs TempAvSum data points and colours by region
                 fill = "#228B22", colour = c("#000000"), alpha = 0.5, shape = 21, size = 3) +
      labs(y = "Change in FEve \n",
           x = paste0("\n Temperature"),
           title = paste0("CIs: ", FEve.ci$"l-95% CI", " to ", FEve.ci$"u-95% CI"),
           subtitle = paste(ifelse((FEve.ci$"l-95% CI" < 0 & FEve.ci$"u-95% CI" < 0) | # Automatically asigns significance in subtitle based on CIs
                                     (FEve.ci$"l-95% CI" > 0 & FEve.ci$"u-95% CI" > 0),
                                   "Significant: CIs do NOT span 0", "Not Significant: CIs span 0"))) +
      theme_1() +
      theme(legend.position = "bottom",
            plot.subtitle = element_text(colour = ifelse((FEve.ci$"l-95% CI" < 0 & FEve.ci$"u-95% CI" < 0) | # Automatically colours subtitle based on significance from CIs
                                                           (FEve.ci$"l-95% CI" > 0 & FEve.ci$"u-95% CI" > 0),
                                                         paste("#006400"), paste("#8b0000")))))
  
  # SR ----
  
  # Load in model output
  SR.mod <- get(load(paste0("data/model_outputs_new/m_SR_slopes_TempAvSum",
                            pc.filepath, ".RData")))
  
  # Convert model output into a dataframe (with 4.d.p.)
  SR.df <- brms_SummaryTable(SR.mod, formatOptions = list(digits = 5, nsmall = 5), round = 5)
  
  # Extract the confidence intervals as a list for use in the plotting
  SR.df.ci <- SR.df %>% 
    filter(Covariate %in% c("x_variable")) %>% 
    dplyr::select("l-95% CI", "u-95% CI") %>% 
    pivot_longer(cols = 1:2, names_to = "Interval", values_to = "Value")
  
  # Save the confidence intervals as a list and remove the intermediate dataframe
  SR.ci <- as.list(SR.df.ci$Value) # Adding values to the list
  names(SR.ci) <- SR.df.ci$Interval # Adding names to the values
  rm(SR.df, SR.df.ci) # Remove unnecessary dataframes of summary and CIs  
  
  # Range to predict variable
  range.to.predict <- paste0("x_variable [", min(input.data$TempAvSum, na.rm = TRUE), ":", max(input.data$TempAvSum, na.rm = TRUE), 
                             " by =", (max(input.data$TempAvSum, na.rm = TRUE) - min(input.data$TempAvSum, na.rm = TRUE))/50, "]")
  
  # Extract the prediction data frame
  SR.pred <- ggpredict(SR.mod, terms = range.to.predict)
  
  # Extract dataframe with only the x-variable for plotting
  data.x.var.only <- dplyr::select(input.data, TempAvSum)
  
  # Plot model outputs
  (SR.plot <- ggplot(SR.pred) +
      geom_line(aes(x = x, y = predicted), colour = "grey10") + # Adds line for predicted values of SR vs TempAvSum
      geom_ribbon(aes(x = x, ymin = conf.low, ymax = conf.high), fill = "grey10", alpha = 0.2) + # Adds c.intervals for predictions as ribbon
      geom_point(data = input.data, aes(x = data.x.var.only[,1], y = SR_slopes), # Adds original SR vs TempAvSum data points and colours by region
                 fill = "#1C86EE", colour = c("#000000"), alpha = 0.5, shape = 21, size = 3) +
      labs(y = "Change in SR \n",
           x = paste0("\n Temperature"),
           title = paste0("CIs: ", SR.ci$"l-95% CI", " to ", SR.ci$"u-95% CI"),
           subtitle = paste(ifelse((SR.ci$"l-95% CI" < 0 & SR.ci$"u-95% CI" < 0) | # Automatically asigns significance in subtitle based on CIs
                                     (SR.ci$"l-95% CI" > 0 & SR.ci$"u-95% CI" > 0),
                                   "Significant: CIs do NOT span 0", "Not Significant: CIs span 0"))) +
      theme_1() +
      theme(legend.position = "bottom",
            plot.subtitle = element_text(colour = ifelse((SR.ci$"l-95% CI" < 0 & SR.ci$"u-95% CI" < 0) | # Automatically colours subtitle based on significance from CIs
                                                           (SR.ci$"l-95% CI" > 0 & SR.ci$"u-95% CI" > 0),
                                                         paste("#006400"), paste("#8b0000")))))
  
  # FDis ----
  
  # Load in model output
  FDis.mod <- get(load(paste0("data/model_outputs_new/m_FDis_slopes_TempAvSum",
                              pc.filepath, ".RData")))
  
  # Convert model output into a dataframe (with 4.d.p.)
  FDis.df <- brms_SummaryTable(FDis.mod, formatOptions = list(digits = 5, nsmall = 5), round = 5)
  
  # Extract the confidence intervals as a list for use in the plotting
  FDis.df.ci <- FDis.df %>% 
    filter(Covariate %in% c("x_variable")) %>% 
    dplyr::select("l-95% CI", "u-95% CI") %>% 
    pivot_longer(cols = 1:2, names_to = "Interval", values_to = "Value")
  
  # Save the confidence intervals as a list and remove the intermediate dataframe
  FDis.ci <- as.list(FDis.df.ci$Value) # Adding values to the list
  names(FDis.ci) <- FDis.df.ci$Interval # Adding names to the values
  rm(FDis.df, FDis.df.ci) # Remove unnecessary dataframes of summary and CIs  
  
  # Range to predict variable
  range.to.predict <- paste0("x_variable [", min(input.data$TempAvSum, na.rm = TRUE), ":", max(input.data$TempAvSum, na.rm = TRUE), 
                             " by =", (max(input.data$TempAvSum, na.rm = TRUE) - min(input.data$TempAvSum, na.rm = TRUE))/50, "]")
  
  # Extract the prediction data frame
  FDis.pred <- ggpredict(FDis.mod, terms = range.to.predict)
  
  # Extract dataframe with only the x-variable for plotting
  data.x.var.only <- dplyr::select(input.data, TempAvSum)
  
  # Plot model outputs
  (FDis.plot <- ggplot(FDis.pred) +
      geom_line(aes(x = x, y = predicted), colour = "grey10") + # Adds line for predicted values of FDis vs TempAvSum
      geom_ribbon(aes(x = x, ymin = conf.low, ymax = conf.high), fill = "grey10", alpha = 0.2) + # Adds c.intervals for predictions as ribbon
      geom_point(data = input.data, aes(x = data.x.var.only[,1], y = FDis_slopes), # Adds original FDis vs TempAvSum data points and colours by region
                 fill = "#EE7600", colour = c("#000000"), alpha = 0.5, shape = 21, size = 3) +
      labs(y = "Change in FDis \n",
           x = paste0("\n Temperature"),
           title = paste0("CIs: ", FDis.ci$"l-95% CI", " to ", FDis.ci$"u-95% CI"),
           subtitle = paste(ifelse((FDis.ci$"l-95% CI" < 0 & FDis.ci$"u-95% CI" < 0) | # Automatically asigns significance in subtitle based on CIs
                                     (FDis.ci$"l-95% CI" > 0 & FDis.ci$"u-95% CI" > 0),
                                   "Significant: CIs do NOT span 0", "Not Significant: CIs span 0"))) +
      theme_1() +
      theme(legend.position = "bottom",
            plot.subtitle = element_text(colour = ifelse((FDis.ci$"l-95% CI" < 0 & FDis.ci$"u-95% CI" < 0) | # Automatically colours subtitle based on significance from CIs
                                                           (FDis.ci$"l-95% CI" > 0 & FDis.ci$"u-95% CI" > 0),
                                                         paste("#006400"), paste("#8b0000")))))
  
  # Create a panel of all three plots
  combined.panel <- grid.arrange(FRic.plot, FEve.plot, FDis.plot, SR.plot, ncol = 2)
  
  # Export panel
  ggsave(combined.panel, filename = paste0("figures/outputs_new/combined_slopes_TempAvSum",
                                           pc.filepath, ".png"), width = 15, height = 14)
  
}


# EXTRA 6: PRECIPITATION CHANGE TEMPORAL ----

plot.precip.change <- function(){
  
  # Produce dataframe for running models on
  input.data <- slopes.input %>% 
    dplyr::select(PrecSlope, FRic_slopes, FEve_slopes, SR_slopes, FDis_slopes,
                  SurveyedArea, gridcell, SUBSITE, PlotDominatingFG)
  
  # FRic ----
  
  # Load in model output
  FRic.mod <- get(load(paste0("data/model_outputs_new/m_FRic_slopes_PrecSlope",
                              pc.filepath, ".RData")))
  
  # Convert model output into a dataframe (with 4.d.p.)
  FRic.df <- brms_SummaryTable(FRic.mod, formatOptions = list(digits = 5, nsmall = 5), round = 5)
  
  # Extract the confidence intervals as a list for use in the plotting
  FRic.df.ci <- FRic.df %>% 
    filter(Covariate %in% c("x_variable")) %>% 
    dplyr::select("l-95% CI", "u-95% CI") %>% 
    pivot_longer(cols = 1:2, names_to = "Interval", values_to = "Value")
  
  # Save the confidence intervals as a list and remove the intermediate dataframe
  FRic.ci <- as.list(FRic.df.ci$Value) # Adding values to the list
  names(FRic.ci) <- FRic.df.ci$Interval # Adding names to the values
  rm(FRic.df, FRic.df.ci) # Remove unnecessary dataframes of summary and CIs  
  
  # Range to predict variable
  range.to.predict <- paste0("x_variable [", min(input.data$PrecSlope, na.rm = TRUE), ":", max(input.data$PrecSlope, na.rm = TRUE), 
                             " by =", (max(input.data$PrecSlope, na.rm = TRUE) - min(input.data$PrecSlope, na.rm = TRUE))/50, "]")
  
  # Extract the prediction data frame
  FRic.pred <- ggpredict(FRic.mod, terms = range.to.predict)
  
  # Extract dataframe with only the x-variable for plotting
  data.x.var.only <- dplyr::select(input.data, PrecSlope)
  
  # Plot model outputs
  (FRic.plot <- ggplot(FRic.pred) +
      geom_line(aes(x = x, y = predicted), colour = "grey10") + # Adds line for predicted values of FRic vs PrecSlope
      geom_ribbon(aes(x = x, ymin = conf.low, ymax = conf.high), fill = "grey10", alpha = 0.2) + # Adds c.intervals for predictions as ribbon
      geom_point(data = input.data, aes(x = data.x.var.only[,1], y = FRic_slopes), # Adds original FRic vs PrecSlope data points and colours by region
                 fill = "#D02090", colour = c("#000000"), alpha = 0.5, shape = 21, size = 3) +
      labs(y = "Change in FRic \n",
           x = paste0("\n Change in Precipitation"),
           title = paste0("CIs: ", FRic.ci$"l-95% CI", " to ", FRic.ci$"u-95% CI"),
           subtitle = paste(ifelse((FRic.ci$"l-95% CI" < 0 & FRic.ci$"u-95% CI" < 0) | # Automatically asigns significance in subtitle based on CIs
                                     (FRic.ci$"l-95% CI" > 0 & FRic.ci$"u-95% CI" > 0),
                                   "Significant: CIs do NOT span 0", "Not Significant: CIs span 0"))) +
      theme_1() +
      theme(legend.position = "bottom",
            plot.subtitle = element_text(colour = ifelse((FRic.ci$"l-95% CI" < 0 & FRic.ci$"u-95% CI" < 0) | # Automatically colours subtitle based on significance from CIs
                                                           (FRic.ci$"l-95% CI" > 0 & FRic.ci$"u-95% CI" > 0),
                                                         paste("#006400"), paste("#8b0000")))))
  
  # FEve ----
  
  # Load in model output
  FEve.mod <- get(load(paste0("data/model_outputs_new/m_FEve_slopes_PrecSlope",
                              pc.filepath, ".RData")))
  
  # Convert model output into a dataframe (with 4.d.p.)
  FEve.df <- brms_SummaryTable(FEve.mod, formatOptions = list(digits = 5, nsmall = 5), round = 5)
  
  # Extract the confidence intervals as a list for use in the plotting
  FEve.df.ci <- FEve.df %>% 
    filter(Covariate %in% c("x_variable")) %>% 
    dplyr::select("l-95% CI", "u-95% CI") %>% 
    pivot_longer(cols = 1:2, names_to = "Interval", values_to = "Value")
  
  # Save the confidence intervals as a list and remove the intermediate dataframe
  FEve.ci <- as.list(FEve.df.ci$Value) # Adding values to the list
  names(FEve.ci) <- FEve.df.ci$Interval # Adding names to the values
  rm(FEve.df, FEve.df.ci) # Remove unnecessary dataframes of summary and CIs  
  
  # Range to predict variable
  range.to.predict <- paste0("x_variable [", min(input.data$PrecSlope, na.rm = TRUE), ":", max(input.data$PrecSlope, na.rm = TRUE), 
                             " by =", (max(input.data$PrecSlope, na.rm = TRUE) - min(input.data$PrecSlope, na.rm = TRUE))/50, "]")
  
  # Extract the prediction data frame
  FEve.pred <- ggpredict(FEve.mod, terms = range.to.predict)
  
  # Extract dataframe with only the x-variable for plotting
  data.x.var.only <- dplyr::select(input.data, PrecSlope)
  
  # Plot model outputs
  (FEve.plot <- ggplot(FEve.pred) +
      geom_line(aes(x = x, y = predicted), colour = "grey10") + # Adds line for predicted values of FEve vs PrecSlope
      geom_ribbon(aes(x = x, ymin = conf.low, ymax = conf.high), fill = "grey10", alpha = 0.2) + # Adds c.intervals for predictions as ribbon
      geom_point(data = input.data, aes(x = data.x.var.only[,1], y = FEve_slopes), # Adds original FEve vs PrecSlope data points and colours by region
                 fill = "#228B22", colour = c("#000000"), alpha = 0.5, shape = 21, size = 3) +
      labs(y = "Change in FEve \n",
           x = paste0("\n Change in Precipitation"),
           title = paste0("CIs: ", FEve.ci$"l-95% CI", " to ", FEve.ci$"u-95% CI"),
           subtitle = paste(ifelse((FEve.ci$"l-95% CI" < 0 & FEve.ci$"u-95% CI" < 0) | # Automatically asigns significance in subtitle based on CIs
                                     (FEve.ci$"l-95% CI" > 0 & FEve.ci$"u-95% CI" > 0),
                                   "Significant: CIs do NOT span 0", "Not Significant: CIs span 0"))) +
      theme_1() +
      theme(legend.position = "bottom",
            plot.subtitle = element_text(colour = ifelse((FEve.ci$"l-95% CI" < 0 & FEve.ci$"u-95% CI" < 0) | # Automatically colours subtitle based on significance from CIs
                                                           (FEve.ci$"l-95% CI" > 0 & FEve.ci$"u-95% CI" > 0),
                                                         paste("#006400"), paste("#8b0000")))))
  
  # SR ----
  
  # Load in model output
  SR.mod <- get(load(paste0("data/model_outputs_new/m_SR_slopes_PrecSlope",
                            pc.filepath, ".RData")))
  
  # Convert model output into a dataframe (with 4.d.p.)
  SR.df <- brms_SummaryTable(SR.mod, formatOptions = list(digits = 5, nsmall = 5), round = 5)
  
  # Extract the confidence intervals as a list for use in the plotting
  SR.df.ci <- SR.df %>% 
    filter(Covariate %in% c("x_variable")) %>% 
    dplyr::select("l-95% CI", "u-95% CI") %>% 
    pivot_longer(cols = 1:2, names_to = "Interval", values_to = "Value")
  
  # Save the confidence intervals as a list and remove the intermediate dataframe
  SR.ci <- as.list(SR.df.ci$Value) # Adding values to the list
  names(SR.ci) <- SR.df.ci$Interval # Adding names to the values
  rm(SR.df, SR.df.ci) # Remove unnecessary dataframes of summary and CIs  
  
  # Range to predict variable
  range.to.predict <- paste0("x_variable [", min(input.data$PrecSlope, na.rm = TRUE), ":", max(input.data$PrecSlope, na.rm = TRUE), 
                             " by =", (max(input.data$PrecSlope, na.rm = TRUE) - min(input.data$PrecSlope, na.rm = TRUE))/50, "]")
  
  # Extract the prediction data frame
  SR.pred <- ggpredict(SR.mod, terms = range.to.predict)
  
  # Extract dataframe with only the x-variable for plotting
  data.x.var.only <- dplyr::select(input.data, PrecSlope)
  
  # Plot model outputs
  (SR.plot <- ggplot(SR.pred) +
      geom_line(aes(x = x, y = predicted), colour = "grey10") + # Adds line for predicted values of SR vs PrecSlope
      geom_ribbon(aes(x = x, ymin = conf.low, ymax = conf.high), fill = "grey10", alpha = 0.2) + # Adds c.intervals for predictions as ribbon
      geom_point(data = input.data, aes(x = data.x.var.only[,1], y = SR_slopes), # Adds original SR vs PrecSlope data points and colours by region
                 fill = "#1C86EE", colour = c("#000000"), alpha = 0.5, shape = 21, size = 3) +
      labs(y = "Change in SR \n",
           x = paste0("\n Change in Precipitation"),
           title = paste0("CIs: ", SR.ci$"l-95% CI", " to ", SR.ci$"u-95% CI"),
           subtitle = paste(ifelse((SR.ci$"l-95% CI" < 0 & SR.ci$"u-95% CI" < 0) | # Automatically asigns significance in subtitle based on CIs
                                     (SR.ci$"l-95% CI" > 0 & SR.ci$"u-95% CI" > 0),
                                   "Significant: CIs do NOT span 0", "Not Significant: CIs span 0"))) +
      theme_1() +
      theme(legend.position = "bottom",
            plot.subtitle = element_text(colour = ifelse((SR.ci$"l-95% CI" < 0 & SR.ci$"u-95% CI" < 0) | # Automatically colours subtitle based on significance from CIs
                                                           (SR.ci$"l-95% CI" > 0 & SR.ci$"u-95% CI" > 0),
                                                         paste("#006400"), paste("#8b0000")))))
  
  # FDis ----
  
  # Load in model output
  FDis.mod <- get(load(paste0("data/model_outputs_new/m_FDis_slopes_PrecSlope",
                              pc.filepath, ".RData")))
  
  # Convert model output into a dataframe (with 4.d.p.)
  FDis.df <- brms_SummaryTable(FDis.mod, formatOptions = list(digits = 5, nsmall = 5), round = 5)
  
  # Extract the confidence intervals as a list for use in the plotting
  FDis.df.ci <- FDis.df %>% 
    filter(Covariate %in% c("x_variable")) %>% 
    dplyr::select("l-95% CI", "u-95% CI") %>% 
    pivot_longer(cols = 1:2, names_to = "Interval", values_to = "Value")
  
  # Save the confidence intervals as a list and remove the intermediate dataframe
  FDis.ci <- as.list(FDis.df.ci$Value) # Adding values to the list
  names(FDis.ci) <- FDis.df.ci$Interval # Adding names to the values
  rm(FDis.df, FDis.df.ci) # Remove unnecessary dataframes of summary and CIs  
  
  # Range to predict variable
  range.to.predict <- paste0("x_variable [", min(input.data$PrecSlope, na.rm = TRUE), ":", max(input.data$PrecSlope, na.rm = TRUE), 
                             " by =", (max(input.data$PrecSlope, na.rm = TRUE) - min(input.data$PrecSlope, na.rm = TRUE))/50, "]")
  
  # Extract the prediction data frame
  FDis.pred <- ggpredict(FDis.mod, terms = range.to.predict)
  
  # Extract dataframe with only the x-variable for plotting
  data.x.var.only <- dplyr::select(input.data, PrecSlope)
  
  # Plot model outputs
  (FDis.plot <- ggplot(FDis.pred) +
      geom_line(aes(x = x, y = predicted), colour = "grey10") + # Adds line for predicted values of FDis vs PrecSlope
      geom_ribbon(aes(x = x, ymin = conf.low, ymax = conf.high), fill = "grey10", alpha = 0.2) + # Adds c.intervals for predictions as ribbon
      geom_point(data = input.data, aes(x = data.x.var.only[,1], y = FDis_slopes), # Adds original FDis vs PrecSlope data points and colours by region
                 fill = "#EE7600", colour = c("#000000"), alpha = 0.5, shape = 21, size = 3) +
      labs(y = "Change in FDis \n",
           x = paste0("\n Change in Precipitation"),
           title = paste0("CIs: ", FDis.ci$"l-95% CI", " to ", FDis.ci$"u-95% CI"),
           subtitle = paste(ifelse((FDis.ci$"l-95% CI" < 0 & FDis.ci$"u-95% CI" < 0) | # Automatically asigns significance in subtitle based on CIs
                                     (FDis.ci$"l-95% CI" > 0 & FDis.ci$"u-95% CI" > 0),
                                   "Significant: CIs do NOT span 0", "Not Significant: CIs span 0"))) +
      theme_1() +
      theme(legend.position = "bottom",
            plot.subtitle = element_text(colour = ifelse((FDis.ci$"l-95% CI" < 0 & FDis.ci$"u-95% CI" < 0) | # Automatically colours subtitle based on significance from CIs
                                                           (FDis.ci$"l-95% CI" > 0 & FDis.ci$"u-95% CI" > 0),
                                                         paste("#006400"), paste("#8b0000")))))
  
  # Create a panel of all three plots
  combined.panel <- grid.arrange(FRic.plot, FEve.plot, FDis.plot, SR.plot, ncol = 2)
  
  # Export panel
  ggsave(combined.panel, filename = paste0("figures/outputs_new/combined_slopes_PrecSlope",
                                           pc.filepath, ".png"), width = 15, height = 14)
  
}


# END ----


# OLD VERSIONS OF FUNCTIONS ----

# # FUNCTION: SR VS FRic ----
# 
# bayesian.metric.comparison.FRic <- function(run.model, quadratic, x.var, data, colour.by, colour.discrete,
#                                             chains, cores, warmup, iterations, delta, treedepth){
#   
#   # Modify combo.latest for colouring points
#   combo.latest <- combo.latest %>% 
#     mutate(PlotDominatingFG = ifelse(PlotDominatingFG == "Shrub-Dominated",
#                                      "Shrub", PlotDominatingFG),
#            PlotDominatingFG = ifelse(PlotDominatingFG == "Graminoid-Dominated",
#                                      "Graminoid", PlotDominatingFG),
#            PlotDominatingFG = ifelse(PlotDominatingFG == "Forb-Dominated",
#                                      "Forb", PlotDominatingFG)) %>% 
#     rename(colour.variable = colour.by)
#   
#   # Run if loop for if linear
#   if (quadratic == "No"){
#     
#     # Generate input dataframe
#     input.data <- data %>% 
#       dplyr::select(x.var, FRic, FEve, SurveyedArea, gridcell, SUBSITE) %>% 
#       rename(x_variable = x.var) # Rename first column to x_variable
#     
#     # Run model in if loop
#     if (run.model == TRUE){
#       
#       # Function for running the Bayesian model
#       SR.FRic.mod <- brm(FRic ~ x_variable + log(SurveyedArea) + (1 | gridcell / SUBSITE),
#                          data = input.data, family = FRic.distribution, chains = chains,
#                          warmup = warmup, iter = iterations, cores = cores,
#                          control = list(adapt_delta = delta,
#                                         max_treedepth = treedepth))
#       
#       # Export model output
#       save(SR.FRic.mod, file = paste0("data/model_outputs_new/m_FRic_",
#                                       FRic.distribution, "_SR", pc.filepath, ".RData"))
#       
#     } else {
#       
#       # Load in model output
#       SR.FRic.mod <- get(load(paste0("data/model_outputs_new/m_FRic_",
#                                      FRic.distribution, "_SR", pc.filepath, ".RData")))
#       
#     } # End of if else loop
#     
#     # Convert model output into a dataframe (with 4.d.p.)
#     SR.FRic.df <- brms_SummaryTable(SR.FRic.mod, formatOptions = list(digits = 4, nsmall = 4), round = 4)
#     
#     # Extract the confidence intervals as a list for use in the plotting
#     SR.FRic.df.ci <- SR.FRic.df %>% 
#       filter(Covariate %in% c("x_variable")) %>% 
#       dplyr::select("l-95% CI", "u-95% CI") %>% 
#       pivot_longer(cols = 1:2, names_to = "Interval", values_to = "Value")
#     
#     # Save the confidence intervals as a list and remove the intermediate dataframe
#     SR.FRic.ci <- as.list(SR.FRic.df.ci$Value) # Adding values to the list
#     names(SR.FRic.ci) <- SR.FRic.df.ci$Interval # Adding names to the values
#     rm(SR.FRic.df, SR.FRic.df.ci) # Remove unnecessary dataframes of summary and CIs  
#     
#     # Extract the prediction data frame
#     SR.FRic.pred <- ggpredict(SR.FRic.mod, terms = "x_variable [4:20]",  back.transform = TRUE)
#     
#     # Extract dataframe with only the x-variable for plotting
#     data.x.var.only <- dplyr::select(data, SR)
#     
#     # Plot model outputs
#     (SR.FRic.plot <- ggplot(SR.FRic.pred) +
#         geom_line(aes(x = x, y = predicted), colour = "grey10") + # Adds line for predicted values of FRic vs LAT
#         geom_ribbon(aes(x = x, ymin = conf.low, ymax = conf.high), fill = "grey10", alpha = 0.2) + # Adds c.intervals for predictions as ribbon
#         geom_point(data = combo.latest, aes(x = data.x.var.only[,1], y = FRic, fill = colour.variable), # Adds original FRic vs LAT data points and colours by region
#                    colour = c("#000000"), alpha = 0.5, shape = 21, size = 3) +
#         # scale_fill_viridis(option = FRic.colour, begin = 0.2, end = 1, direction = -1, discrete = colour.discrete,
#         #                    guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")) +
#         scale_fill_gradient2(low = "#D02090",
#                              # mid = "#228B22",
#                              mid = "#FFFFFF",
#                              high = "#1E90FF", midpoint = 0.5,
#                              guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")) +
#         scale_y_continuous(limits = c(0, max(combo.latest$FRic)*1.05)) +
#         labs(y = "FRic \n",
#              x = "\n SR",
#              fill = paste0(colour.by),
#              title = paste0("CIs: ", SR.FRic.ci$"l-95% CI", " to ", SR.FRic.ci$"u-95% CI"),
#              subtitle = paste(ifelse((SR.FRic.ci$"l-95% CI" < 0 & SR.FRic.ci$"u-95% CI" < 0) | # Automatically asigns significance in subtitle based on CIs
#                                        (SR.FRic.ci$"l-95% CI" > 0 & SR.FRic.ci$"u-95% CI" > 0),
#                                      "Significant: CIs do NOT span 0", "Not Significant: CIs span 0"))) +
#         theme_1() +
#         theme(legend.position = "right",
#               plot.subtitle = element_text(colour = ifelse((SR.FRic.ci$"l-95% CI" < 0 & SR.FRic.ci$"u-95% CI" < 0) | # Automatically colours subtitle based on significance from CIs
#                                                              (SR.FRic.ci$"l-95% CI" > 0 & SR.FRic.ci$"u-95% CI" > 0),
#                                                            paste("#006400"), paste("#8b0000")))))
#     
#     # Export plot
#     ggsave(SR.FRic.plot, filename = paste0("figures/outputs_new/SR_vs_FRic",
#                                            pc.filepath, ".png"), width = 8, height = 7)
#     
#   } # End of not quadratic if loop
#   
#   # Run if loop for is quadratic
#   if (quadratic == "Yes"){
#     
#     # Generate input dataframe
#     input.data <- data %>% 
#       dplyr::select(x.var, FRic, FEve, SurveyedArea, gridcell, SUBSITE) %>% 
#       rename(SR = x.var) # Rename first column to x_variable
#     
#     # Determine mean value for centering
#     mean.SR <- mean(input.data$SR, na.rm = TRUE)
#     
#     # Center the x_variable
#     input.data <- input.data %>% 
#       mutate(centred_SR = SR - mean.SR)
#     
#     # Run model in if loop
#     if (run.model == TRUE){
#       
#       # Function for running the Bayesian model
#       SR.FRic.mod <- brm(FRic ~ centred_SR + I(centred_SR^2) + (1 | gridcell / SUBSITE),
#                          data = input.data, family = FRic.distribution, chains = chains,
#                          warmup = warmup, iter = iterations, cores = cores,
#                          control = list(adapt_delta = delta,
#                                         max_treedepth = treedepth))
#       
#       # Export model output
#       save(SR.FRic.mod, file = paste0("data/model_outputs_new/m_FRic_",
#                                       FRic.distribution, "_quadratic_SR", pc.filepath, ".RData"))
#       
#     } else {
#       
#       # Load in model output
#       SR.FRic.mod <- get(load(paste0("data/model_outputs_new/m_FRic_",
#                                      FRic.distribution, "_quadratic_SR", pc.filepath, ".RData")))
#       
#     } # End of if else loop
#     
#     # Check summary of model
#     summary(SR.FRic.mod)
#     
#     # Convert model output into a dataframe (with 4.d.p.)
#     SR.FRic.df <- brms_SummaryTable(SR.FRic.mod, formatOptions = list(digits = 4, nsmall = 4), round = 4)
#     
#     # Extract the confidence intervals as a list for use in the plotting
#     SR.FRic.df.ci <- SR.FRic.df %>% 
#       filter(Covariate %in% c("centred_SR")) %>% 
#       dplyr::select("l-95% CI", "u-95% CI") %>% 
#       pivot_longer(cols = 1:2, names_to = "Interval", values_to = "Value")
#     
#     # Save the confidence intervals as a list and remove the intermediate dataframe
#     SR.FRic.ci <- as.list(SR.FRic.df.ci$Value) # Adding values to the list
#     names(SR.FRic.ci) <- SR.FRic.df.ci$Interval # Adding names to the values
#     rm(SR.FRic.df, SR.FRic.df.ci) # Remove unnecessary dataframes of summary and CIs  
#     
#     # Create template predictions dataframe
#     SR.FRic.pred.df = data.frame(centred_SR = seq(0 - mean.SR, 100 - mean.SR))
#     
#     # Add predictions to template datframe
#     SR.FRic.pred = add_epred_draws(SR.FRic.pred.df, SR.FRic.mod, re_formula = NA) %>% 
#       mutate(SR = centred_SR + mean.SR) # Uncentre the x variable
#     
#     # Plot model outputs
#     (SR.FRic.plot <- ggplot(data = SR.FRic.pred, aes(x = SR, y = .epred)) +
#         geom_point(data = input.data, aes(x = SR, y = FRic), alpha = 0) +
#         stat_lineribbon(aes(y = .epred), .width = 0.95, fill = "grey10",
#                         alpha = 0.2, colour = "#000000", linetype = 2, linewidth = 0.1) +
#         # scale_fill_viridis(option = FEve.colour, begin = 0.2, end = 1, direction = -1, discrete = colour.discrete,
#         #                    guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")) +
#         scale_x_continuous(limits = c(0,max(combo.latest$SR, na.rm = TRUE))) +
#         scale_y_continuous(limits = c(0, 61)) +
#         labs(y = "FRic \n",
#              x = "SR \n",
#              fill = "FEve \n",
#              title = paste0("CIs: ", SR.FRic.ci$"l-95% CI", " to ", SR.FRic.ci$"u-95% CI"),
#              subtitle = paste(ifelse((SR.FRic.ci$"l-95% CI" < 0 & SR.FRic.ci$"u-95% CI" < 0) | # Automatically asigns significance in subtitle based on CIs
#                                        (SR.FRic.ci$"l-95% CI" > 0 & SR.FRic.ci$"u-95% CI" > 0),
#                                      "Significant: CIs do NOT span 0", "Not Significant: CIs span 0"))) +
#         theme_1() +
#         theme(legend.position = "right",
#               plot.subtitle = element_text(colour = ifelse((SR.FRic.ci$"l-95% CI" < 0 & SR.FRic.ci$"u-95% CI" < 0) | # Automatically colours subtitle based on significance from CIs
#                                                              (SR.FRic.ci$"l-95% CI" > 0 & SR.FRic.ci$"u-95% CI" > 0),
#                                                            paste("#006400"), paste("#8b0000")))))
#     
#     # Export plot
#     ggsave(SR.FRic.plot, filename = paste0("figures/outputs_new/SR_vs_FRic_quadratic",
#                                            pc.filepath, ".png"), width = 8, height = 7)
#     
#     # Plot another one without line just for making final plot with points
#     (SR.FRic.plot.points <- ggplot() +
#         geom_point(data = input.data, aes(x = SR, y = FRic, fill = FEve),
#                    colour = c("#000000"), alpha = 0.5, shape = 21, size = 3) +
#         scale_fill_gradient2(low = "#D02090",
#                              # mid = "#228B22",
#                              mid = "#FFFFFF",
#                              high = "#1E90FF", midpoint = 0.5,
#                              guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")) +
#         scale_x_continuous(limits = c(0, max(combo.latest$SR, na.rm = TRUE))) +
#         scale_y_continuous(limits = c(0, 61)) +
#         labs(y = "FRic \n",
#              x = "SR \n",
#              fill = "FEve \n",
#              title = paste0("CIs: ", SR.FRic.ci$"l-95% CI", " to ", SR.FRic.ci$"u-95% CI"),
#              subtitle = paste(ifelse((SR.FRic.ci$"l-95% CI" < 0 & SR.FRic.ci$"u-95% CI" < 0) | # Automatically asigns significance in subtitle based on CIs
#                                        (SR.FRic.ci$"l-95% CI" > 0 & SR.FRic.ci$"u-95% CI" > 0),
#                                      "Significant: CIs do NOT span 0", "Not Significant: CIs span 0"))) +
#         theme_1() +
#         theme(legend.position = "none",
#               plot.subtitle = element_text(colour = ifelse((SR.FRic.ci$"l-95% CI" < 0 & SR.FRic.ci$"u-95% CI" < 0) | # Automatically colours subtitle based on significance from CIs
#                                                              (SR.FRic.ci$"l-95% CI" > 0 & SR.FRic.ci$"u-95% CI" > 0),
#                                                            paste("#006400"), paste("#8b0000")))))
#     
#     # Export plot
#     ggsave(SR.FRic.plot.points, filename = paste0("figures/outputs_new/SR_vs_FRic_quadratic_points",
#                                                   pc.filepath, ".png"), width = 8, height = 7)
#     
#   } # End of quadratic if loop
#   
# } # End of function
# 
# # FUNCTION: SR VS FEve ----
# bayesian.metric.comparison.FEve <- function(run.model, quadratic, x.var, data, colour.by, colour.discrete,
#                                             chains, cores, warmup, iterations, delta, treedepth){
#   
#   # Modify combo.latest for colouring points
#   combo.latest <- combo.latest %>% 
#     mutate(PlotDominatingFG = ifelse(PlotDominatingFG == "Shrub-Dominated",
#                                      "Shrub", PlotDominatingFG),
#            PlotDominatingFG = ifelse(PlotDominatingFG == "Graminoid-Dominated",
#                                      "Graminoid", PlotDominatingFG),
#            PlotDominatingFG = ifelse(PlotDominatingFG == "Forb-Dominated",
#                                      "Forb", PlotDominatingFG)) %>% 
#     rename(colour.variable = colour.by)
#   
#   
#   # Generate input dataframe
#   input.data <- data %>% 
#     dplyr::select(x.var, FEve, FEve, SurveyedArea, gridcell, SUBSITE) %>% 
#     rename(x_variable = x.var) # Rename first column to x_variable
#   
#   # Run model in if loop
#   if (run.model == TRUE){
#     
#     # Function for running the Bayesian model
#     SR.FEve.mod <- brm(FEve ~ x_variable + log(SurveyedArea) + (1 | gridcell / SUBSITE),
#                        data = input.data, family = FEve.distribution, chains = chains,
#                        warmup = warmup, iter = iterations, cores = cores,
#                        control = list(adapt_delta = delta,
#                                       max_treedepth = treedepth))
#     
#     # Export model output
#     save(SR.FEve.mod, file = paste0("data/model_outputs_new/m_FEve_",
#                                     FEve.distribution, "_SR", pc.filepath, ".RData"))
#     
#   } else {
#     
#     # Load in model output
#     SR.FEve.mod <- get(load(paste0("data/model_outputs_new/m_FEve_",
#                                    FEve.distribution, "_SR", pc.filepath, ".RData")))
#     
#   } # End of if else loop
#   
#   # Convert model output into a dataframe (with 4.d.p.)
#   SR.FEve.df <- brms_SummaryTable(SR.FEve.mod, formatOptions = list(digits = 4, nsmall = 4), round = 4)
#   
#   # Extract the confidence intervals as a list for use in the plotting
#   SR.FEve.df.ci <- SR.FEve.df %>% 
#     filter(Covariate %in% c("x_variable")) %>% 
#     dplyr::select("l-95% CI", "u-95% CI") %>% 
#     pivot_longer(cols = 1:2, names_to = "Interval", values_to = "Value")
#   
#   # Save the confidence intervals as a list and remove the intermediate dataframe
#   SR.FEve.ci <- as.list(SR.FEve.df.ci$Value) # Adding values to the list
#   names(SR.FEve.ci) <- SR.FEve.df.ci$Interval # Adding names to the values
#   rm(SR.FEve.df, SR.FEve.df.ci) # Remove unnecessary dataframes of summary and CIs  
#   
#   # Extract the prediction data frame
#   SR.FEve.pred <- ggpredict(SR.FEve.mod, terms = "x_variable [4:20]",  back.transform = TRUE)
#   
#   # Extract dataframe with only the x-variable for plotting
#   data.x.var.only <- dplyr::select(data, SR)
#   
#   # Plot model outputs
#   (SR.FEve.plot <- ggplot(SR.FEve.pred) +
#       geom_line(aes(x = x, y = predicted), colour = "grey10") + # Adds line for predicted values of FEve vs LAT
#       geom_ribbon(aes(x = x, ymin = conf.low, ymax = conf.high), fill = "grey10", alpha = 0.2) + # Adds c.intervals for predictions as ribbon
#       geom_point(data = combo.latest, aes(x = SR, y = FEve, fill = colour.variable), # Adds original FEve vs LAT data points and colours by region
#                  colour = c("#000000"), alpha = 0.5, shape = 21, size = 3) +
#       # scale_fill_viridis(option = FEve.colour, begin = 0.2, end = 1, direction = -1, discrete = colour.discrete,
#       #                    guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")) +
#       scale_fill_gradient2(low = "#D02090",
#                            # mid = "#228B22",
#                            mid = "#FFFFFF",
#                            high = "#1E90FF", midpoint = 0.5,
#                            guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")) +
#       scale_y_continuous(limits = c(0, max(combo.latest$FEve)*1.05)) +
#       labs(y = "FEve \n",
#            x = "\n SR",
#            fill = paste0(colour.by),
#            title = paste0("CIs: ", SR.FEve.ci$"l-95% CI", " to ", SR.FEve.ci$"u-95% CI"),
#            subtitle = paste(ifelse((SR.FEve.ci$"l-95% CI" < 0 & SR.FEve.ci$"u-95% CI" < 0) | # Automatically asigns significance in subtitle based on CIs
#                                      (SR.FEve.ci$"l-95% CI" > 0 & SR.FEve.ci$"u-95% CI" > 0),
#                                    "Significant: CIs do NOT span 0", "Not Significant: CIs span 0"))) +
#       theme_1() +
#       theme(legend.position = "right",
#             plot.subtitle = element_text(colour = ifelse((SR.FEve.ci$"l-95% CI" < 0 & SR.FEve.ci$"u-95% CI" < 0) | # Automatically colours subtitle based on significance from CIs
#                                                            (SR.FEve.ci$"l-95% CI" > 0 & SR.FEve.ci$"u-95% CI" > 0),
#                                                          paste("#006400"), paste("#8b0000")))))
#   
#   # Export plot
#   ggsave(SR.FEve.plot, filename = paste0("figures/outputs_new/SR_vs_FEve",
#                                          pc.filepath, ".png"), width = 8, height = 7)
#   
#   
# } # End of function
