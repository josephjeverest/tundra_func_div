# 08 - FUNCTIONS: FD Calculations
# Joseph Everest
# December 2021, adapted March 2022, July 2023


# FUNCTION: Create a function for running the a FD calculation loop on multiple subset dataframes ----

FD.calc <- function(x){
  
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
  
  # Create output dataframe
  fd.output <- data.frame()
  
  # Get actual name of vector x
  vector.name = str_remove((deparse(substitute(x))), pattern = "itex_rownames_")

  # Run loop to calulate FD indices
  for (i in x){
    
    # Split up ITEX input into smaller chunks (full itex.input too large to run at once)
    itex.input.cut <- itex.input %>%
      filter(row.names(.) %in% i) %>%
      select_if(colSums(.) != 0)
    
    # Determine number of species in plot ("nbsp" and "sing.sp")
    nbsp <- ncol(itex.input.cut)
    
    if (nbsp %in% c(1)){ # If 1 species in the plot, breaks calculations
      
      # Create dataframe to append to FD output
      single.species.append <- c(nbsp, nbsp, NA, NA, NA, NA, NA, NA, NA) # FRic = NA and FEve = NA when only one species
      
      single.species.append <- data.frame(t(single.species.append))
      
      colnames(single.species.append) <- c("nbsp", "sing.sp", "FRic", "qual.FRic", "FEve", "FDis", "RaoQ", "Warning", "Error")
      
      rownames(single.species.append) <- c(i)
      
      # Join dataframe to FD output
      fd.output <- rbind(fd.output, single.species.append)
      
    } else { # End of if nbsp == 1 if statement
     
      # Generate list of species (41,461) to retain
      retain.col <- colnames(itex.input.cut)
      
      # Retain only the rows for these species in trait input
      trait.input.cut <- trait.input %>%
        filter(row.names(trait.input) %in% retain.col)
      
      # Set trait input row names back to those filtered by (e.g. not starting from 1)
      rownames(trait.input.cut) <- retain.col
      
      # Calculate FRic and FEve from FD package inside tryCatch function
      fd.raw <- myTryCatch(dbFD(trait.input.cut, itex.input.cut,
                                w.abun = TRUE, stand.x = TRUE, calc.FRic = TRUE,
                                calc.CWM = FALSE, calc.FDiv = FALSE, ord = c("podani"),
                                corr = c("cailliez"), m = m.input, print.pco = TRUE)) # Lower m if R keeps crashing during processing
      
      # Add if statement for correcting the qual.FRic score for species == 3
      if (nbsp >= 4){
        
        # Append output values including warnings and errors to dataframe
        fd.raw.df <- data.frame("nbsp" = ifelse(is.null(fd.raw$value$nbsp[1]), NA, as.numeric(fd.raw$value$nbsp[1])),
                                "sing.sp" = ifelse(is.null(fd.raw$value$sing.sp[1]), NA, as.numeric(fd.raw$value$sing.sp[1])),
                                "FRic" = ifelse(is.null(fd.raw$value$FRic[1]), NA, as.numeric(fd.raw$value$FRic[1])),
                                "qual.FRic" = ifelse(is.null(fd.raw$value$qual.FRic[1]), NA, as.numeric(fd.raw$value$qual.FRic[1])),
                                "FEve" = ifelse(is.null(fd.raw$value$FEve[1]), NA, as.numeric(fd.raw$value$FEve[1])),
                                "FDis" = ifelse(is.null(fd.raw$value$FDis[1]), NA, as.numeric(fd.raw$value$FDis[1])),
                                "RaoQ" = ifelse(is.null(fd.raw$value$RaoQ[1]), NA, as.numeric(fd.raw$value$RaoQ[1])),
                                "Warning" = ifelse(is.null(fd.raw$warning$message[1]), NA, as.character(fd.raw$warning$message[1])),
                                "Error" = ifelse(is.null(fd.raw$error$message[1]), NA, as.character(fd.raw$error$message[1])))
        
      } else { # End of if loop for qual.FRic corrections
        
        # Append output values including warnings and errors to dataframe for when nbsp == 3
        fd.raw.df <- data.frame("nbsp" = ifelse(is.null(fd.raw$value$nbsp[1]), NA, as.numeric(fd.raw$value$nbsp[1])),
                                "sing.sp" = ifelse(is.null(fd.raw$value$sing.sp[1]), NA, as.numeric(fd.raw$value$sing.sp[1])),
                                "FRic" = ifelse(is.null(fd.raw$value$FRic[1]), NA, as.numeric(fd.raw$value$FRic[1])),
                                "qual.FRic" = NA,
                                # "qual.FRic" = ifelse(is.null(fd.raw$value$qual.FRic[1]), NA, as.numeric(fd.raw$value$qual.FRic[1])),
                                "FEve" = ifelse(is.null(fd.raw$value$FEve[1]), NA, as.numeric(fd.raw$value$FEve[1])),
                                "FDis" = ifelse(is.null(fd.raw$value$FDis[1]), NA, as.numeric(fd.raw$value$FDis[1])),
                                "RaoQ" = ifelse(is.null(fd.raw$value$RaoQ[1]), NA, as.numeric(fd.raw$value$RaoQ[1])),
                                "Warning" = ifelse(is.null(fd.raw$warning$message[1]), NA, as.character(fd.raw$warning$message[1])),
                                "Error" = ifelse(is.null(fd.raw$error$message[1]), NA, as.character(fd.raw$error$message[1])))
        
      } # End of else loop for qual.FRic corrections
      
      # Append ID number (i) (of 6346) as row numbers
      rownames(fd.raw.df) <- i
      
      # Append to output dataframe
      fd.output <- rbind(fd.output, fd.raw.df)
       
    } # End of else statement
    
  } # End of for loop
  
  # Add further corrections to output
  fd.output.clean <- fd.output %>% 
    mutate(FRic = ifelse(nbsp == 2, NA, FRic))
  
  # Create output filepath
  path.output <- paste0("data/output_08_fd_output_", vector.name, pc.filepath, ".csv")
  
  # Save output dataframe
  write.csv(fd.output.clean, path.output)
  
  # Remove unwanted variables
  rm(itex.input.cut, retain.col, trait.input.cut, fd.raw, fd.raw.df, fd.output, path.output, path.output.csv,
     vector.name, duplicates.trait, duplicates.true, duplicates.false, duplicates.append, nbsp, single.species.append)
  
} # End of function











# 
# # OLD FUNCTION: Function for running FD calculations but has multiple stipulations involved ----
# 
# FD.calc.old <- function(x){
#   
#   # Load bespoke tryCatch function into the local environment
#     # Found at: https://stackoverflow.com/questions/4948361/how-do-i-save-warnings-and-errors-as-output-from-a-function
#   myTryCatch <- function(expr) {
#     warn <- err <- NULL
#     value <- withCallingHandlers(
#       tryCatch(expr, error=function(e) {
#         err <<- e
#         NULL
#       }), warning=function(w) {
#         warn <<- w
#         invokeRestart("muffleWarning")
#       })
#     list(value=value, warning=warn, error=err)
#   }
#   
#   # Create output dataframe
#   fd.output <- data.frame()
#   
#   # Get actual name of vector x
#   vector.name = str_remove((deparse(substitute(x))), pattern = "itex_rownames_")
#   
#   # Run loop to calulate FD indices
#   for (i in x){
#     
#     # Split up ITEX input into smaller chunks (full itex.input too large to run at once)
#     itex.input.cut <- itex.input %>%
#       filter(row.names(.) %in% i) %>%
#       select_if(colSums(.) != 0)
#     
#     # Determine number of species in plot ("nbsp" and "sing.sp")
#     nbsp <- ncol(itex.input.cut)
#     
#     if (nbsp == 1 | i == "2428"){ # Row 2428 breaks convex hull creation with a super confusing error!
#       
#       # Create dataframe to append to FD output
#       single.species.append <- c(nbsp, nbsp, NA, NA, NA, NA, NA, NA, NA) # FRic = NA and FEve = NA when only one species
#       
#       single.species.append <- data.frame(t(single.species.append))
#       
#       colnames(single.species.append) <- c("nbsp", "sing.sp", "FRic", "qual.FRic", "FEve", "FDis", "RaoQ", "Warning", "Error")
#       
#       rownames(single.species.append) <- c(i)
#       
#       # Join dataframe to FD output
#       fd.output <- rbind(fd.output, single.species.append)
#       
#     } else {
#       
#       # Generate list of species (41,461) to retain
#       retain.col <- colnames(itex.input.cut)
#       
#       # Retain only the rows for these species in trait input
#       trait.input.cut <- trait.input %>%
#         filter(row.names(trait.input) %in% retain.col)
#       
#       # Set trait input row names back to those filtered by (e.g. not starting from 1)
#       rownames(trait.input.cut) <- retain.col
#       
#       # Check for duplicates in traits dataframe
#       duplicates.trait <- duplicated(trait.input.cut)
#       
#       # Remove first record as always says FALSE as no previous records to be duplicated with
#       duplicates.trait <- duplicates.trait[-1]
#       
#       # Test to see whether there are any duplicates in the dataframe (e.g. includes TRUE)
#       duplicates.true <- TRUE %in% duplicates.trait
#       
#       # See whether contains any non-duplicates
#       duplicates.false <- FALSE %in% duplicates.trait
#       
#       # IF ONLY CONTAINS DUPLICATES
#       if (duplicates.true == TRUE & duplicates.false == FALSE){
#         
#         # IF LESS THAN THREE SPECIES
#         if (nbsp < 3){
#           
#           # Create dataframe to append to FD output
#           duplicates.append <- c(nbsp, nbsp, "0", NA, NA, NA, NA, NA, NA) # FRic = 0 and FEve = NA when traits identical but < 3 species
#           
#           duplicates.append <- data.frame(t(duplicates.append))
#           
#           colnames(duplicates.append) <- c("nbsp", "sing.sp", "FRic", "qual.FRic", "FEve", "FDis", "RaoQ", "Warning", "Error")
#           
#           rownames(duplicates.append) <- c(i)
#           
#           # Join dataframe to FD output
#           fd.output <- rbind(fd.output, duplicates.append)
#           
#         } else {
#           
#           # Create dataframe to append to FD output
#           duplicates.append <- c(nbsp, nbsp, "0", NA, "1", NA, NA, NA, NA) # FRic = 0 and FEve = 1 when traits identical and 3 or more species
#           
#           duplicates.append <- data.frame(t(duplicates.append))
#           
#           colnames(duplicates.append) <- c("nbsp", "sing.sp", "FRic", "qual.FRic", "FEve", "FDis", "RaoQ", "Warning", "Error")
#           
#           rownames(duplicates.append) <- c(i)
#           
#           # Join dataframe to FD output
#           fd.output <- rbind(fd.output, duplicates.append)
#           
#         }
#         
#       } else {
#         
#         # Calculate FRic and FEve from FD package inside tryCatch function
#         fd.raw <- myTryCatch(dbFD(trait.input.cut, itex.input.cut,
#                                   w.abun = TRUE, stand.x = TRUE, calc.FRic = TRUE,
#                                   calc.CWM = FALSE, calc.FDiv = FALSE, ord = c("podani"),
#                                   corr = c("cailliez"), m = m.input, print.pco = TRUE)) # Lower m if R keeps crashing during processing
#         
#         # Append output values including warnings and errors to dataframe
#         fd.raw.df <- data.frame("nbsp" = ifelse(is.null(fd.raw$value$nbsp[1]), NA, as.numeric(fd.raw$value$nbsp[1])),
#                                 "sing.sp" = ifelse(is.null(fd.raw$value$sing.sp[1]), NA, as.numeric(fd.raw$value$sing.sp[1])),
#                                 "FRic" = ifelse(is.null(fd.raw$value$FRic[1]), NA, as.numeric(fd.raw$value$FRic[1])),
#                                 "qual.FRic" = ifelse(is.null(fd.raw$value$qual.FRic[1]), NA, as.numeric(fd.raw$value$qual.FRic[1])),
#                                 "FEve" = ifelse(is.null(fd.raw$value$FEve[1]), NA, as.numeric(fd.raw$value$FEve[1])),
#                                 "FDis" = ifelse(is.null(fd.raw$value$FDis[1]), NA, as.numeric(fd.raw$value$FDis[1])),
#                                 "RaoQ" = ifelse(is.null(fd.raw$value$RaoQ[1]), NA, as.numeric(fd.raw$value$RaoQ[1])),
#                                 "Warning" = ifelse(is.null(fd.raw$warning$message[1]), NA, as.character(fd.raw$warning$message[1])),
#                                 "Error" = ifelse(is.null(fd.raw$error$message[1]), NA, as.character(fd.raw$error$message[1])))
# 
#         # Append ID number (i) (of 6346) as row numbers
#         rownames(fd.raw.df) <- i
#         
#         # Append to output dataframe
#         fd.output <- rbind(fd.output, fd.raw.df)
#         
#       }
#       
#     }
#     
#   }
#   
#   # Create output filepath
#   path.output <- paste0("data/output_08_fd_output_", vector.name, pc.filepath, ".csv")
#   
#   # Save output dataframe
#   write.csv(fd.output, path.output)
#   
#   # Remove unwanted variables
#   rm(itex.input.cut, retain.col, trait.input.cut, fd.raw, fd.raw.df, fd.output, path.output, path.output.csv,
#      vector.name, duplicates.trait, duplicates.true, duplicates.false, duplicates.append, nbsp, single.species.append)
#   
# }
