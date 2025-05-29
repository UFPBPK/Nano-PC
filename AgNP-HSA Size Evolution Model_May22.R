### AgNP-HSA Size Evolution Model ###
### Author: Xinyue Chen a,b,c, Zhoumeng Lin a,b,c,*
##### a Department of Environmental and Global Health, College of Public Health and Health Professions, University of Florida, Gainesville, FL 32611, United States 
##### b Center for Environmental and Human Toxicology, University of Florida, Gainesville, FL 32611, United States 
##### c Center for Pharmacometrics and Systems Pharmacology, University of Florida, Orlando, FL 32827, United States 
##### * Corresponding author at: Department of Environmental and Global Health, College of Public Health and Health Professions, University of Florida, 2187 Mowry Rd, Gainesville, FL 32611, United States. Phone: +1-352-273-6160. 
##### E-mail address: linzhoumeng@ufl.edu (Z. Lin).


# Load required libraries
library(deSolve) # Purpose: Solving differential equations (ODE system for HSA corona formation kinetics on silver nanoparticles (AgNPs))
library(dplyr)   # Purpose: Data manipulation and transformation (processing simulation results, filtering datasets, and combining multiple nanoparticle data)
library(R6)      # Purpose: Object-oriented programming framework (creating the AgNP_HSASizeSimulator class with encapsulated methods and parameters)


# Define the AgNP-HSA Size Evolution Simulator class
AgNP_HSASizeSimulator <- R6::R6Class(
  "AgNP_HSASizeSimulator",
  public = list(
    # Calibration data (AgNP20) - Zhao et al., 2021
    Agnp20_data = list(
      period1 = list(initial = 23.0, final = 28.5, time = 0.3),
      period2 = list(initial = 25.1, final = 36.9, time = 1.0),
      period3 = list(initial = 26.7, final = 40.3, time = 3.0)
    ),
    # Validation data (AgNP40) - Zhao et al., 2021
    Agnp40_data = list(
      period1 = list(initial = 45.2, final = 58.5, time = 0.3),
      period2 = list(initial = 40.5, final = 74.9, time = 1.0),
      period3 = list(initial = 46.6, final = 72.3, time = 3.0)
    ),
    # Validation data (AgNP80) - Zhao et al., 2021
    Agnp80_data = list(
      period1 = list(initial = 88.9, final = 95.0, time = 0.3),
      period2 = list(initial = 95.1, final = 124.4, time = 1.0),
      period3 = list(initial = 93.7, final = 133.5, time = 3.0)
    ),
    
    # Validation data (AgNP19) - Shannahan et al., 2015
    shannahan_data = list(
      initial = 19.13,
      final = 69.99,
      time = 8.0  # 8 hour equilibrium point
    ),
    
    # Basic parameters
    protein_sizes = list(HSA = 8),
    k_on_base = 2.4e3,                    # Base association rate constant for HSA (M⁻¹·s⁻¹ ),  cited from Dell’Orco et al. (2010)
    k_off_base = 2.0e-3,                  # Base dissociation rate constant for HSA (s⁻¹ ), cited from Dell’Orco et al. (2010)
    c_init = list(HSA = 6.0e-4),            # Initial HSA concentration (M), cited from Dell’Orco et al. (2010)
    
    

    # Store protein layer thickness for each period
    L = list(),
    
    # Improved size-dependent scaling factors
    get_size_factors = function(diameter, time_h) {
      # Enhanced size-dependent effects
      relative_size <- diameter / 20  # normalized to SNP20
      
      # Progressive size effect over time
      time_factor <- 1 + 0.2 * log(time_h + 0.01)
      
      # Stronger size dependence for larger particles
      k_on_factor <- relative_size^0.8 * time_factor  # Increased size dependence
      k_off_factor <- relative_size^(-0.5)  # Stronger size effect on dissociation
      
      # Additional size-dependent multiplier for larger particles
      size_multiplier <- ifelse(diameter > 40, 
                                1 + 0.2 * log(diameter/40), 
                                1)
      
      return(list(
        k_on = self$k_on_base * k_on_factor * size_multiplier,
        k_off = self$k_off_base * k_off_factor
      ))
    },
    
    initialize = function() {
      self$calibrate_model()
    },
    
    calibrate_model = function() {
      # Calculate protein layer thickness using SNP20 data
      self$L$period1 <- 
        (self$Agnp20_data$period1$final - self$Agnp20_data$period1$initial) / 2
      
      self$L$period2 <- 
        (self$Agnp20_data$period2$final - self$Agnp20_data$period2$initial) / 2
      
      self$L$period3 <- 
        (self$Agnp20_data$period3$final - self$Agnp20_data$period3$initial) / 2
      
      self$L$long_term <- mean(c(
        (self$shannahan_data$final - self$shannahan_data$initial) / 2
      ))
      
      cat("Calibrated protein layer thickness (from SNP20):\n")
      cat("Period 1 (0-0.3h):", round(self$L$period1, 2), "nm\n")
      cat("Period 2 (0-1.0h):", round(self$L$period2, 2), "nm\n")
      cat("Period 3 (0-3.0h):", round(self$L$period3, 2), "nm\n")
      cat("Long-term (8.0h):", round(self$L$long_term, 2), "nm\n")
    },
    
    simulate_period = function(initial_size, period_end, protein_layer) {
      # Define the model with enhanced mechanisms
      model <- function(t, state, parameters) {
        with(as.list(c(state, parameters)), {
          time_h <- t/3600  # Convert time to hours for factor calculation
          
          # Get updated kinetic parameters based on current time
          size_factors <- self$get_size_factors(initial_size, time_h)
          
          # Enhanced surface availability model
          available_surface <- (1 - theta)^2
          
          # Improved crowding effects
          crowding_factor <- (1 - 0.1 * theta) * (1 + 0.7442 * theta^2)  # Non-linear crowding
          
          # Size-dependent cooperative binding
          cooperative_factor <- 1 + 0.5 * theta * (initial_size/20)^0.3
          
          dtheta <- size_factors$k_on * self$c_init$HSA * available_surface * 
            crowding_factor * cooperative_factor - 
            size_factors$k_off * theta^1.8  # Adjusted exponent
          
          list(c(dtheta))
        })
      }
      
      times <- seq(0, period_end * 3600, by = 60)  # One minute time steps
      state <- c(theta = 0.0)
      out <- ode(y = state, times = times, func = model, parms = NULL)
      
      data <- as.data.frame(out)
      names(data) <- c("time", "coverage")
      
      # Enhanced size-dependent layer thickness calculation
      size_factor <- (initial_size / 20)^0.3  # Increased size dependence
      base_layer <- protein_layer * size_factor #base_layer = 
      
      data <- data %>%
        mutate(
          time_h = time / 3600,
          # Diameter calculation
          diameter = initial_size + 2 * coverage * base_layer * 
            (1 + 0.15 * log(time_h + 0.001)) * 
            (1 + 0.1 * (initial_size/20)^0.5)     # L = base_layer * (1 + 0.15 * log(time_h + 0.001)) * (1 + 0.1 * (initial_size/20)^0.5)
        )
      
      return(data)
    },
    
    # Run all simulations with consistent approach
    simulate_all_data = function() {
      # Create empty list to store all results
      all_results <- list()
      
      # Define simulation parameters for each dataset
      sim_params <- list(
        # AgNP20 data - run separate simulations for each time period
        "Zhao et al., 2021 (20 nm, Calibration)" = list(
          initial_sizes = c(
            self$Agnp20_data$period1$initial,
            self$Agnp20_data$period2$initial,
            self$Agnp20_data$period3$initial
          ),
          end_times = c(
            self$Agnp20_data$period1$time,
            self$Agnp20_data$period2$time,
            self$Agnp20_data$period3$time
          ),
          protein_layers = c(
            self$L$period1,
            self$L$period2,
            self$L$period3
          )
        ),
        
        # AgNP40 data
        "Zhao et al., 2021 (40 nm)" = list(
          initial_sizes = c(
            self$Agnp40_data$period1$initial,
            self$Agnp40_data$period2$initial,
            self$Agnp40_data$period3$initial
          ),
          end_times = c(
            self$Agnp40_data$period1$time,
            self$Agnp40_data$period2$time,
            self$Agnp40_data$period3$time
          ),
          protein_layers = c(
            self$L$period1,
            self$L$period2,
            self$L$period3
          )
        ),
        
        # AgNP80 data
        "Zhao et al., 2021 (80 nm)" = list(
          initial_sizes = c(
            self$Agnp80_data$period1$initial,
            self$Agnp80_data$period2$initial,
            self$Agnp80_data$period3$initial
          ),
          end_times = c(
            self$Agnp80_data$period1$time,
            self$Agnp80_data$period2$time,
            self$Agnp80_data$period3$time
          ),
          protein_layers = c(
            self$L$period1,
            self$L$period2,
            self$L$period3
          )
        ),
        
        # Shannahan data - only one time point
        "Shannahan et al., 2015 (19 nm)" = list(
          initial_sizes = self$shannahan_data$initial,
          end_times = self$shannahan_data$time,
          protein_layers = self$L$long_term
        )
      )
      
      # Run simulations for each dataset and each time period
      for (dataset_name in names(sim_params)) {
        params <- sim_params[[dataset_name]]
        dataset_results <- list()
        
        # Loop through all time periods for this dataset
        for (i in 1:length(params$end_times)) {
          # Run simulation
          result <- self$simulate_period(
            params$initial_sizes[i],
            params$end_times[i],
            params$protein_layers[i]
          )
          
          # Add dataset info
          result$group <- dataset_name
          result$dataset_type <- ifelse(grepl("Calibration", dataset_name), "Calibration", "Validation")
          result$period <- paste0("0-", params$end_times[i], "h")
          result$initial_size <- params$initial_sizes[i]  # Store initial size for validation
          result$combo_id <- paste(dataset_name, paste0("0-", params$end_times[i], "h"), sep = "_")
          
          # Add to results list
          dataset_results[[i]] <- result
        }
        
        # Combine all periods for this dataset
        all_results[[dataset_name]] <- bind_rows(dataset_results)
      }
      
      # Combine all datasets
      all_data <- bind_rows(all_results)
      return(all_data)
    },
    
    # Create experimental data points
    create_experimental_data = function() {
      # AgNP20 data (Calibration)
      df_agnp20 <- data.frame(
        time_h = c(
          0, self$Agnp20_data$period1$time,
          0, self$Agnp20_data$period2$time,
          0, self$Agnp20_data$period3$time
        ),
        diameter = c(
          self$Agnp20_data$period1$initial, self$Agnp20_data$period1$final,
          self$Agnp20_data$period2$initial, self$Agnp20_data$period2$final,
          self$Agnp20_data$period3$initial, self$Agnp20_data$period3$final
        ),
        group = rep("Zhao et al., 2021 (20 nm, Calibration)", 6),
        point_type = rep(c("Initial", "Final"), 3),
        period = rep(c("0-0.3h", "0-1.0h", "0-3.0h"), each = 2)
      )
      
      # AgNP40 data (Validation)
      df_agnp40 <- data.frame(
        time_h = c(
          0, self$Agnp40_data$period1$time,
          0, self$Agnp40_data$period2$time,
          0, self$Agnp40_data$period3$time
        ),
        diameter = c(
          self$Agnp40_data$period1$initial, self$Agnp40_data$period1$final,
          self$Agnp40_data$period2$initial, self$Agnp40_data$period2$final,
          self$Agnp40_data$period3$initial, self$Agnp40_data$period3$final
        ),
        group = rep("Zhao et al., 2021 (40 nm)", 6),
        point_type = rep(c("Initial", "Final"), 3),
        period = rep(c("0-0.3h", "0-1.0h", "0-3.0h"), each = 2)
      )
      
      # AgNP80 data (Validation)
      df_agnp80 <- data.frame(
        time_h = c(
          0, self$Agnp80_data$period1$time,
          0, self$Agnp80_data$period2$time,
          0, self$Agnp80_data$period3$time
        ),
        diameter = c(
          self$Agnp80_data$period1$initial, self$Agnp80_data$period1$final,
          self$Agnp80_data$period2$initial, self$Agnp80_data$period2$final,
          self$Agnp80_data$period3$initial, self$Agnp80_data$period3$final
        ),
        group = rep("Zhao et al., 2021 (80 nm)", 6),
        point_type = rep(c("Initial", "Final"), 3),
        period = rep(c("0-0.3h", "0-1.0h", "0-3.0h"), each = 2)
      )
      
      # Shannahan data (Validation)
      df_shannahan <- data.frame(
        time_h = c(0, self$shannahan_data$time),
        diameter = c(
          self$shannahan_data$initial, 
          self$shannahan_data$final
        ),
        group = rep("Shannahan et al., 2015 (19 nm)", 2),
        point_type = c("Initial", "Final"),
        period = rep("0-8.0h", 2)
      )
      
      # Combine all data
      exp_data <- bind_rows(df_agnp20, df_agnp40, df_agnp80, df_shannahan)
      return(exp_data)
    },
    
    # Function to create a summary table of validation results
    create_summary_table = function() {
      # Define dataset information - using the correct initial and final values
      summary_df <- data.frame(
        Dataset = c(
          "Zhao et al., 2021 (20 nm, Calibration) - 0.3h",
          "Zhao et al., 2021 (20 nm, Calibration) - 1.0h", 
          "Zhao et al., 2021 (20 nm, Calibration) - 3.0h",
          "Zhao et al., 2021 (40 nm) - 0.3h",
          "Zhao et al., 2021 (40 nm) - 1.0h", 
          "Zhao et al., 2021 (40 nm) - 3.0h",
          "Zhao et al., 2021 (80 nm) - 0.3h",
          "Zhao et al., 2021 (80 nm) - 1.0h", 
          "Zhao et al., 2021 (80 nm) - 3.0h", 
          "Shannahan et al., 2015 (19 nm) - 8.0h"
        ),
        Type = c(
          rep("Calibration", 3),
          rep("Validation", 7)
        ),
        Initial_Size_nm = c(
          self$Agnp20_data$period1$initial,
          self$Agnp20_data$period2$initial, 
          self$Agnp20_data$period3$initial,
          self$Agnp40_data$period1$initial,
          self$Agnp40_data$period2$initial, 
          self$Agnp40_data$period3$initial,
          self$Agnp80_data$period1$initial,
          self$Agnp80_data$period2$initial, 
          self$Agnp80_data$period3$initial,
          self$shannahan_data$initial
        ),
        Final_Size_nm = c(
          self$Agnp20_data$period1$final,
          self$Agnp20_data$period2$final, 
          self$Agnp20_data$period3$final,
          self$Agnp40_data$period1$final,
          self$Agnp40_data$period2$final, 
          self$Agnp40_data$period3$final,
          self$Agnp80_data$period1$final,
          self$Agnp80_data$period2$final, 
          self$Agnp80_data$period3$final,
          self$shannahan_data$final
        ),
        Measurement_Time_h = c(
          0.3, 1.0, 3.0,
          0.3, 1.0, 3.0,
          0.3, 1.0, 3.0,
          8.0
        )
      )
      
      # Calculate size increase and percentage for each dataset
      summary_df$Size_Increase_nm <- summary_df$Final_Size_nm - summary_df$Initial_Size_nm
      summary_df$Size_Increase_Percent <- round((summary_df$Size_Increase_nm / summary_df$Initial_Size_nm) * 100, 1)
      
      return(summary_df)
    },
    
    # Helper function to print validation statistics in a clean format
    print_validation_stats = function(results) {
      # Get the summary table
      summary_table <- self$create_summary_table()
      
      # Run simulations to get the predicted values at each time point
      simulated_values <- list()
      
      # For each dataset and time, get the simulated value
      for (i in 1:nrow(summary_table)) {
        dataset_info <- summary_table[i, ]
        
        # Extract parts from dataset name
        dataset_parts <- strsplit(as.character(dataset_info$Dataset), " - ")[[1]]
        dataset_name <- dataset_parts[1]
        time_str <- gsub("h", "", dataset_parts[2])
        time_val <- as.numeric(time_str)
        
        # Find the matching results
        dataset_results <- results %>%
          filter(group == dataset_name & abs(time_h - time_val) < 0.01)
        
        if (nrow(dataset_results) > 0) {
          # Take the closest time point
          closest_result <- dataset_results %>%
            arrange(abs(time_h - time_val)) %>%
            slice(1)
          
          simulated_values[[i]] <- closest_result$diameter
        } else {
          simulated_values[[i]] <- NA
        }
      }
      
      # Add simulated values to the data frame
      summary_table$Final_Size_Simulated_nm <- round(unlist(simulated_values), 2)
      
      # Calculate error percentage
      summary_table$Error_Percent <- round(abs(summary_table$Final_Size_Simulated_nm - summary_table$Final_Size_nm) / 
                                             summary_table$Final_Size_nm * 100, 1)
      
      # Reorganize columns for better readability
      summary_table <- summary_table %>%
        select(Dataset, Type, Initial_Size_nm, Final_Size_nm, Final_Size_Simulated_nm, 
               Measurement_Time_h, Size_Increase_nm, Size_Increase_Percent, Error_Percent)
      
      # Print the summary table
      cat("\n===================================================\n")
      cat("MODEL VALIDATION RESULTS\n")
      cat("===================================================\n\n")
      
      print(summary_table)
      
      cat("\n===================================================\n")
      
      # Return the summary table for potential further use
      return(summary_table)
    }
  )
)

# Main execution function
run_corona_simulation <- function() {
  # Create simulator instance
  simulator <- AgNP_HSASizeSimulator$new()
  
  # Run simulations
  results <- simulator$simulate_all_data()
  
  # Print validation statistics
  summary_table <- simulator$print_validation_stats(results)
  
  # Return results for potential further use
  return(list(
    simulator = simulator,
    results = results,
    summary_table = summary_table
  ))
}

# Use this function to run the simulation:
simulation_output <- run_corona_simulation()