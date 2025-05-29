### AgNP-BSA Size Evolution Model ###
### Author: Xinyue Chen a,b,c, Zhoumeng Lin a,b,c,*
##### a Department of Environmental and Global Health, College of Public Health and Health Professions, University of Florida, Gainesville, FL 32611, United States 
##### b Center for Environmental and Human Toxicology, University of Florida, Gainesville, FL 32611, United States 
##### c Center for Pharmacometrics and Systems Pharmacology, University of Florida, Orlando, FL 32827, United States 
##### * Corresponding author at: Department of Environmental and Global Health, College of Public Health and Health Professions, University of Florida, 2187 Mowry Rd, Gainesville, FL 32611, United States. Phone: +1-352-273-6160. 
##### E-mail address: linzhoumeng@ufl.edu (Z. Lin).


# Load required libraries
library(deSolve) # Purpose: Solving differential equations (ODE system for BSA corona formation kinetics on silver nanoparticles (AgNPs))
library(dplyr)   # Purpose: Data manipulation and transformation (processing simulation results, filtering datasets, and combining multiple nanoparticle data)
library(R6)      # Purpose: Object-oriented programming framework (creating the AgNP_BSASizeSimulator class with encapsulated methods and parameters)


# Define the AgNP-BSA Size Evolution Simulator class
AgNP_BSASizeSimulator <- R6::R6Class(
  "AgNP_BSASizeSimulator",
  public = list(
    # Basic Parameters
    # Calibration data
    D_np_1 = 38.2,                     # AgNP-citrate (Calibration) (nm)
    D_np_2 = 40.5,                     # AgNP-PVP (Calibration) (nm)
    
    # Validation data
    D_np_3 = 144.2,                    # AgNP-citrate (Large) (nm)
    D_np_4 = 165.8,                    # AgNP-PVP (Large) (nm)
    D_np_5 = 23.6,                     # Shannahan et al., 2015 (validation) (nm)
    D_np_6 = 89.3,                     # Zhang et al., 2023 (nm)
    
    # BSA protein parameters
    protein_sizes = list(BSA = 7),      # Unit (nm) for BSA
    k_on_base = 2.0e3,                  # Base association rate constant for BSA (M⁻¹·s⁻¹ ), estimated value based on Hanigan-Diebel et al. (2024)
    k_off_base = 1.5e-3,                # Base dissociation rate constant for BSA (s⁻¹ ), estimated value based on Hanigan-Diebel et al. (2024)
    c_init = 6.0e-4,                      # Initial BSA concentration (M), based on Dell’Orco et al. (2010)
    
    # Surface coating factors
    coating_factors = list(
      "citrate" = 0.922,     # F_coat=(d_(final,1)-d_(initial,1))/2÷L_base=((91.5-38.2)/2)÷28.9=0.922
      "PVP" = 1.078         # F_coat=(d_(final,2)-d_(initial,2))/2÷L_base=((102.8-40.5)/2)÷28.9=1.078
    ),
    
    # Size scaling factors
    size_scaling = list(
      "small" = 1.0,      # Base scaling for particles < 50 nm
      "medium" = 0.85,    # Reduced scaling for 50-100 nm
      "large" = 0.5       # Further reduced scaling for > 100 nm
    ),
    
    # Calibration parameters
    L_base = NULL, # Will be calculated from calibration data
    L = list(),  # For each set of data
    
    initialize = function(initial_conditions = list()) {
      if (length(initial_conditions) > 0) {
        for (param in names(initial_conditions)) {
          if (param %in% names(self)) {
            self[[param]] <- initial_conditions[[param]]
          }
        }
      }
      # Calculate protein layer thickness from calibration data
      self$calibrate_model()
    },
    
    calibrate_model = function() {
      # Using AgNP-citrate and AgNP-PVP data for calibration
      final_size_1 = 91.5   # AgNP-citrate + BSA (1h)
      final_size_2 = 102.8  # AgNP-PVP + BSA (1h)
      
      # Calculate protein layer thickness for full coverage (θ = 1)
      # Then we'll adjust for actual coverage in the simulation
      thickness_1 = (final_size_1 - self$D_np_1) / 2
      thickness_2 = (final_size_2 - self$D_np_2) / 2
      
      # Use average thickness as base
      self$L_base = mean(c(thickness_1, thickness_2))
      
      # Store individual thicknesses for comparison
      self$L = list(
        "AgNP-citrate" = thickness_1,
        "AgNP-PVP" = thickness_2
      )
      
      cat("Calibrated protein layer thickness (base):", self$L_base, "nm\n")
      cat("Individual thicknesses for different AgNP coatings:\n")
      for (type in names(self$L)) {
        cat(type, ":", self$L[[type]], "nm\n")
      }
    },
    
    # Determine size category for scaling
    get_size_category = function(d_initial) {
      if (d_initial < 50) {
        return("small")
      } else if (d_initial < 100) {
        return("medium")
      } else {
        return("large")
      }
    },
    
    # Determine coating type based on dataset name
    get_coating_type = function(dataset_name) {
      if (grepl("citrate", dataset_name, ignore.case = TRUE)) {
        return("citrate")
      } else if (grepl("PVP", dataset_name, ignore.case = TRUE)) {
        return("PVP")
      } else {
        return("citrate")  # Zhang et al.,2023 & Shannahan et al.,2015 (validation) both used citrate-coated AgNPs
      }
    },
    
    # Calculate the size-adjusted protein layer thickness
    calculate_layer_thickness = function(d_initial, coating_type) {
      # Base thickness adjusted by coating type and size
      size_cat <- self$get_size_category(d_initial)
      size_factor <- self$size_scaling[[size_cat]]
      coating_factor <- self$coating_factors[[coating_type]]
      
      # Adjustment formula with coating factor influence and size-dependent scaling
      L = self$L_base * coating_factor * size_factor * (1 + 0.15 * log(d_initial/40 + 1))
      return(L)
    },
    
    # Calculate size effect on k_on
    calculate_size_effect = function(d_initial) {
      # Modified size effect - larger particles have relatively smaller surface area to volume ratio
      # More significant reduction for larger particles
      if (d_initial < 50) {
        size_effect = (d_initial/40)^0.25
      } else if (d_initial < 100) {
        size_effect = (d_initial/40)^0.15
      } else {
        size_effect = (d_initial/40)^0.1
      }
      return(size_effect)
    },
    
    # Calculate time-dependent k_on scaling
    calculate_time_effect = function(t) {
      # Modified time effect - slower initial binding rate for silver
      time_effect = (1 + 0.6*t)/(0.3 + t)
      return(time_effect)
    },
    
    simulate = function(sim_time = 8.0, time_step = 0.1) {
      if (is.null(self$L_base)) {
        stop("Model must be calibrated before simulation")
      }
      
      # Run for all diameters and combine results
      results_list <- list()
      
      # Calibration data
      results_list[[1]] <- self$simulate_single(self$D_np_1, "AgNP-citrate (Calibration)", sim_time, time_step)
      results_list[[2]] <- self$simulate_single(self$D_np_2, "AgNP-PVP (Calibration)", sim_time, time_step)
      
      # Validation data
      results_list[[3]] <- self$simulate_single(self$D_np_3, "AgNP-citrate (Large)", sim_time, time_step)
      results_list[[4]] <- self$simulate_single(self$D_np_4, "AgNP-PVP (Large)", sim_time, time_step)
      results_list[[5]] <- self$simulate_single(self$D_np_5, "Shannahan et al., 2015 (validation)", sim_time, time_step)
      results_list[[6]] <- self$simulate_single(self$D_np_6, "Zhang et al., 2023", sim_time, time_step)
      
      # Combine all results
      results <- bind_rows(results_list)
      return(results)
    },
    
    simulate_single = function(d_initial, dataset_name, sim_time = 8.0, time_step = 0.1) {
      # Determine coating type from dataset name
      coating_type <- self$get_coating_type(dataset_name)
      
      # Calculate specific layer thickness for this diameter and coating
      L = self$calculate_layer_thickness(d_initial, coating_type)
      
      # Calculate size effect on k_on
      size_effect = self$calculate_size_effect(d_initial)
      
      # Define the differential equations
      model <- function(t, state, parameters) {
        with(as.list(c(state, parameters)), {
          # Core model with modifying factors
          # Surface Availability Term (S)
          S = (1 - theta)^1.5
          
          # Cooperative Binding Term (F_coop) - modified for silver
          F_coop = 1 + theta * (d_initial/40)^0.4
          
          # Crowding Effect Term (F_crowd) - modified for silver
          F_crowd = exp(-0.25 * theta)
          
          # Surface coating effect on binding (F_coat)
          F_coat = parameters$coating_factor
          
          # Time-dependent k_on
          k_on = self$k_on_base * size_effect * self$calculate_time_effect(t) * F_coat
          
          # Size-dependent k_off
          k_off = self$k_off_base / (1 + 0.15 * log(d_initial/40 + 1))
          
          # Core equation: dθ/dt = k_on * C_BSA * S * F_coop * F_crowd * F_coat - k_off * θ
          dtheta <- k_on * self$c_init * S * F_coop * F_crowd - k_off * theta
          
          list(c(dtheta))
        })
      }
      
      # Convert simulation time to hours
      times <- seq(0, sim_time, by = time_step)
      state <- c(theta = 0.0)  # Start with no coverage
      
      # Add coating factor to parameters
      model_params <- list(
        d_initial = d_initial,
        coating_factor = self$coating_factors[[coating_type]]
      )
      
      out <- ode(y = state, times = times, func = model, parms = model_params)
      
      data <- as.data.frame(out)
      names(data) <- c("time", "coverage")
      
      # Calculate diameter over time
      data <- data %>%
        mutate(
          time_h = time,
          diameter = d_initial + 2 * L * coverage,
          dataset = dataset_name,
          coating = coating_type
        )
      
      return(data)
    },
    
 
    
    report_validation = function(results) {
      # Define validation times and datasets
      validation_times <- c(1.0, 1.0, 1.0, 1.0, 8.0, 1.0)
      dataset_names <- c(
        "AgNP-citrate (Calibration)", 
        "AgNP-PVP (Calibration)",
        "AgNP-citrate (Large)",
        "AgNP-PVP (Large)",
        "Shannahan et al., 2015 (validation)",
        "Zhang et al., 2023"
      )
      
      # Expected final sizes
      expected_sizes <- list(
        "AgNP-citrate (Calibration)" = 91.5,
        "AgNP-PVP (Calibration)" = 102.8,
        "AgNP-citrate (Large)" = 190.0,
        "AgNP-PVP (Large)" = 225.0,
        "Shannahan et al., 2015 (validation)" = 38.0,
        "Zhang et al., 2023" = 145.0
      )
      
      # Initial sizes
      initial_sizes <- list(
        "AgNP-citrate (Calibration)" = self$D_np_1,
        "AgNP-PVP (Calibration)" = self$D_np_2,
        "AgNP-citrate (Large)" = self$D_np_3,
        "AgNP-PVP (Large)" = self$D_np_4,
        "Shannahan et al., 2015 (validation)" = self$D_np_5,
        "Zhang et al., 2023" = self$D_np_6
      )
      
      cat("\nModel Validation Results:\n")
      cat("=======================\n")
      
      # Create a data frame to store validation results
      val_results <- data.frame(
        Dataset = character(),
        Initial_Size_nm = numeric(),
        Final_Time_h = numeric(),
        Simulated_Size_nm = numeric(),
        Experimental_Size_nm = numeric(),
        Error_percent = numeric(),
        Size_Increase_nm = numeric(),
        stringsAsFactors = FALSE
      )
      
      # Calculate results at specific validation times
      for (i in 1:length(dataset_names)) {
        dataset <- dataset_names[i]
        time_val <- validation_times[i]
        
        if (!is.na(time_val)) {
          # Find the closest time point in the simulation
          dataset_results <- results %>% filter(dataset == !!dataset)
          
          if (nrow(dataset_results) > 0) {
            closest_time <- dataset_results %>%
              mutate(time_diff = abs(time_h - time_val)) %>%
              filter(time_diff == min(time_diff)) %>%
              slice(1)
            
            # Calculate error and size increase
            error_percent <- abs(closest_time$diameter - expected_sizes[[dataset]]) / 
              expected_sizes[[dataset]] * 100
            size_increase <- closest_time$diameter - initial_sizes[[dataset]]
            
            # Add to results data frame
            val_results <- rbind(val_results, data.frame(
              Dataset = dataset,
              Initial_Size_nm = initial_sizes[[dataset]],
              Final_Time_h = time_val,
              Simulated_Size_nm = round(closest_time$diameter, 2),
              Experimental_Size_nm = expected_sizes[[dataset]],
              Error_percent = round(error_percent, 2),
              Size_Increase_nm = round(size_increase, 2)
            ))
            
            # Print results
            cat("\nDataset:", dataset, "\n")
            cat("Coating type:", closest_time$coating, "\n")
            cat("Initial Size:", initial_sizes[[dataset]], "nm\n")
            cat("Simulated Final Size at t =", time_val, "h:", 
                round(closest_time$diameter, 2), "nm\n")
            cat("Experimental Final Size:", expected_sizes[[dataset]], "nm\n")
            cat("Error:", round(error_percent, 2), "%\n")
            cat("Size Increase:", round(size_increase, 2), "nm\n")
          }
        }
      }
      
      # Return the validation results table
      return(val_results)
    },
    
    # Create a function to show mathematical model details
    show_model_details = function() {
      cat("\nAgNP-BSA Size Evolution Mathematical Model\n")
      cat("============================================\n\n")
      
      cat("Core Equation:\n")
      cat("dθ/dt = k_on * C_BSA * S * F_coop * F_crowd * F_coat - k_off * θ\n\n")
      
      cat("Where:\n")
      cat("θ = BSA surface coverage (0 to 1)\n")
      cat("C_BSA = Initial BSA concentration (", self$c_init, ")\n")
      cat("k_on = Association rate constant (varies with size, coating and time)\n")
      cat("k_off = Dissociation rate constant (varies with size)\n\n")
      
      cat("Modifying Factors:\n")
      cat("- Surface Availability Term: S = (1 - θ)^1.5\n")
      cat("- Cooperative Binding Term: F_coop = 1 + θ * (d_initial/40)^0.4\n")
      cat("- Crowding Effect Term: F_crowd = exp(-0.25 * θ)\n")
      cat("- Coating Effect Term: F_coat = coating_factor (citrate: ", self$coating_factors$citrate, 
          ", PVP: ", self$coating_factors$PVP, ")\n")
      cat("- Time-dependent k_on: k_on = k_on_base * size_effect * time_effect * F_coat\n")
      cat("  where size_effect = (d_initial/40)^0.25\n")
      cat("  and time_effect = (1 + 0.6t)/(0.3 + t)\n")
      cat("- Size-dependent k_off: k_off = k_off_base / (1 + 0.15 * ln(d_initial/40 + 1))\n\n")
      
      cat("Size Scaling Factors:\n")
      cat("- Small particles (<50nm): ", self$size_scaling$small, "\n")
      cat("- Medium particles (50-100nm): ", self$size_scaling$medium, "\n")  
      cat("- Large particles (>100nm): ", self$size_scaling$large, "\n\n")
      
      cat("Final Diameter:\n")
      cat("d(t) = d_initial + 2L * θ(t)\n")
      cat("where L = L_base * coating_factor * size_scaling * (1 + 0.15 * ln(d_initial/40 + 1))\n")
      cat("L_base =", self$L_base, "nm (calibrated value)\n\n")
      
      cat("Model Parameters:\n")
      cat("- k_on_base =", self$k_on_base, "\n")
      cat("- k_off_base =", self$k_off_base, "\n")
      cat("- BSA concentration =", self$c_init, "M\n")
    }
  )
)

# Function to create a summary table of all datasets
create_summary_table <- function(simulator) {
  # Create data frame with all datasets
  summary_df <- data.frame(
    Dataset = c(
      "AgNP-citrate (Calibration)", 
      "AgNP-PVP (Calibration)",
      "AgNP-citrate (Large)",
      "AgNP-PVP (Large)",
      "Shannahan et al., 2015 (validation)",
      "Zhang et al., 2023"
    ),
    Coating = c("citrate", "PVP", "citrate", "PVP", "citrate", "PVP"),
    Type = c(rep("Calibration", 2), rep("Validation", 4)),
    Initial_Size_nm = c(
      simulator$D_np_1, simulator$D_np_2, 
      simulator$D_np_3, simulator$D_np_4, 
      simulator$D_np_5, simulator$D_np_6
    ),
    Final_Size_nm = c(91.5, 102.8, 190.0, 225.0, 38.0, 145.0),
    Measurement_Time_h = c(1.0, 1.0, 1.0, 1.0, 8.0, 1.0)
  )
  
  # Calculate size increase for each dataset
  summary_df$Size_Increase_nm <- summary_df$Final_Size_nm - summary_df$Initial_Size_nm
  summary_df$Size_Increase_Percent <- round((summary_df$Size_Increase_nm / summary_df$Initial_Size_nm) * 100, 1)
  
  return(summary_df)
}

# Function to run the simulation and create plots
run_agnp_bsa_simulation <- function() {
  # Create simulator instance
  simulator <- AgNP_BSASizeSimulator$new()
  
  # Show model details
  simulator$show_model_details()
  
  # Run simulation for all datasets
  results <- simulator$simulate(sim_time = 8.0, time_step = 0.1)

  # Generate validation report
  validation_results <- simulator$report_validation(results)
  
  # Return simulation results and plots
  return(list(
    simulator = simulator,
    results = results,
    validation = validation_results
  ))
}

# Run the simulation and store results
simulation_output <- run_agnp_bsa_simulation()

# Create and print summary table
summary_table <- create_summary_table(simulation_output$simulator)
print(summary_table)



