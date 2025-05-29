### AuNP-BSA Size Evolution Model ###
### Author: Xinyue Chen a,b,c, Zhoumeng Lin a,b,c,*
##### a Department of Environmental and Global Health, College of Public Health and Health Professions, University of Florida, Gainesville, FL 32611, United States 
##### b Center for Environmental and Human Toxicology, University of Florida, Gainesville, FL 32611, United States 
##### c Center for Pharmacometrics and Systems Pharmacology, University of Florida, Orlando, FL 32827, United States 
##### * Corresponding author at: Department of Environmental and Global Health, College of Public Health and Health Professions, University of Florida, 2187 Mowry Rd, Gainesville, FL 32611, United States. Phone: +1-352-273-6160. 
##### E-mail address: linzhoumeng@ufl.edu (Z. Lin).


# Load required libraries
library(deSolve) # Purpose: Solving differential equations (ODE system for BSA corona formation kinetics on gold nanoparticles (AuNPs))
library(dplyr)   # Purpose: Data manipulation and transformation (processing simulation results, filtering datasets, and combining multiple nanoparticle data)
library(R6)      # Purpose: Object-oriented programming framework (creating the AuNP_BSASizeSimulator class with encapsulated methods and parameters)


# Define the AuNP-BSA Size Evolution Simulator class
AuNP_BSASizeSimulator <- R6::R6Class(
  "AuNP_BSASizeSimulator",
  public = list(
    # Basic Parameters
    # Calibration data from Hanigan-Diebel et al., 2024
    D_np_1 = 9.4,                      # AuNP-MPTMA (nm)
    D_np_2 = 9.8,                      # AuNP-MEEE (nm)
    D_np_3 = 10.3,                     # AuNP-MHA (nm)
    
    # Validation data
    D_np_4 = 24,                       # Kumar et al., 2024 (nm)
    D_np_5 = 16,                       # Mishra & Das, 2022 (nm)
    D_np_6 = 28,                       # Mishra & Das, 2022 (nm)
    D_np_7 = 41,                       # Mishra & Das, 2022 (nm)
    D_np_8 = 69,                       # Mishra & Das, 2022 (nm)
    
    # BSA protein parameters
    protein_sizes = list(BSA = 7),      # Unit (nm) for BSA
    k_on_base = 2.0e3,                  # Base association rate constant for BSA (M⁻¹·s⁻¹ ), estimated value based on Hanigan-Diebel et al. (2024)
    k_off_base = 1.5e-3,                # Base dissociation rate constant for BSA (s⁻¹ ), estimated value based on Hanigan-Diebel et al. (2024)
    c_init = 6.0e-4,                      # Initial BSA concentration (M), based on Dell’Orco et al. (2010)
    
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
      # Using Hanigan-Diebel data for calibration
      final_size_1 = 22.3  # AuNP-MPTMA + BSA
      final_size_2 = 13.7  # AuNP-MEEE + BSA
      final_size_3 = 13.3  # AuNP-MHA + BSA
      
      # Calculate protein layer thickness for full coverage (θ = 1)
      # Then we'll adjust for actual coverage in the simulation
      thickness_1 = (final_size_1 - self$D_np_1) / 2
      thickness_2 = (final_size_2 - self$D_np_2) / 2
      thickness_3 = (final_size_3 - self$D_np_3) / 2
      
      # Use average thickness as base
      self$L_base = mean(c(thickness_1, thickness_2, thickness_3))
      
      # Store individual thicknesses for comparison
      self$L = list(
        "AuNP-MPTMA" = thickness_1,
        "AuNP-MEEE" = thickness_2,
        "AuNP-MHA" = thickness_3
      )
      
      cat("Calibrated protein layer thickness (base):", self$L_base, "nm\n")
      cat("Individual thicknesses for different AuNP types:\n")
      for (type in names(self$L)) {
        cat(type, ":", self$L[[type]], "nm\n")
      }
    },
    
    # Calculate the size-adjusted protein layer thickness
    calculate_layer_thickness = function(d_initial) {
      L = self$L_base * (1 + 0.2 * log(d_initial/10 + 1))
      return(L)
    },
    
    # Calculate size effect on k_on
    calculate_size_effect = function(d_initial) {
      size_effect = (d_initial/10)^0.3
      return(size_effect)
    },
    
    # Calculate time-dependent k_on scaling
    calculate_time_effect = function(t) {
      time_effect = (1 + 0.8*t)/(0.2 + t)
      return(time_effect)
    },
    
    simulate = function(d_initial = NULL, sim_time = 1.5, time_step = 0.01) {
      if (is.null(self$L_base)) {
        stop("Model must be calibrated before simulation")
      }
      
      # If no specific diameter is provided, run simulations for all datasets
      if (is.null(d_initial)) {
        # Run for all diameters and combine results
        results_list <- list()
        
        # Calibration data
        results_list[[1]] <- self$simulate_single(self$D_np_1, "AuNP-MPTMA (Calibration)", sim_time, time_step)
        results_list[[2]] <- self$simulate_single(self$D_np_2, "AuNP-MEEE (Calibration)", sim_time, time_step)
        results_list[[3]] <- self$simulate_single(self$D_np_3, "AuNP-MHA (Calibration)", sim_time, time_step)
        
        # Validation data
        results_list[[4]] <- self$simulate_single(self$D_np_4, "Kumar et al., 2024 (24 nm)", sim_time, time_step)
        results_list[[5]] <- self$simulate_single(self$D_np_5, "Mishra & Das, 2022 (16 nm)", sim_time, time_step)
        results_list[[6]] <- self$simulate_single(self$D_np_6, "Mishra & Das, 2022 (28 nm)", sim_time, time_step)
        results_list[[7]] <- self$simulate_single(self$D_np_7, "Mishra & Das, 2022 (41 nm)", sim_time, time_step)
        results_list[[8]] <- self$simulate_single(self$D_np_8, "Mishra & Das, 2022 (69 nm)", sim_time, time_step)
        
        # Combine all results
        results <- bind_rows(results_list)
        return(results)
      } else {
        # Run for specific diameter
        return(self$simulate_single(d_initial, paste0("Custom (", d_initial, "nm)"), sim_time, time_step))
      }
    },
    
    simulate_single = function(d_initial, dataset_name, sim_time = 1.5, time_step = 0.01) {
      # Calculate specific layer thickness for this diameter
      L = self$calculate_layer_thickness(d_initial)
      
      # Calculate size effect on k_on
      size_effect = self$calculate_size_effect(d_initial)
      
      # Define the differential equations
      model <- function(t, state, parameters) {
        with(as.list(c(state, parameters)), {
          # Core model with modifying factors
          # Surface Availability Term (S)
          S = (1 - theta)^1.5
          
          # Cooperative Binding Term (F_coop)
          F_coop = 1 + theta * (d_initial/10)^0.5
          
          # Crowding Effect Term (F_crowd)
          F_crowd = exp(-0.2 * theta)
          
          # Time-dependent k_on
          k_on = self$k_on_base * size_effect * self$calculate_time_effect(t)
          
          # Size-dependent k_off
          k_off = self$k_off_base / (1 + 0.2 * log(d_initial/10 + 1))
          
          # Core equation: dθ/dt = k_on * C_BSA * S * F_coop * F_crowd - k_off * θ
          dtheta <- k_on * self$c_init * S * F_coop * F_crowd - k_off * theta
          
          list(c(dtheta))
        })
      }
      
      # Convert simulation time to hours
      times <- seq(0, sim_time, by = time_step)
      state <- c(theta = 0.0)  # Start with no coverage
      out <- ode(y = state, times = times, func = model, parms = list(d_initial = d_initial))
      
      data <- as.data.frame(out)
      names(data) <- c("time", "coverage")
      
      # Calculate diameter over time
      data <- data %>%
        mutate(
          time_h = time,
          diameter = d_initial + 2 * L * coverage,
          dataset = dataset_name
        )
      
      return(data)
    },
    
 
    report_validation = function(results) {
      # Define validation times
      validation_times <- c(1.5, 1.5, 1.5, 1.5, 1/6, 1/6, 1/6, 1/6)
      dataset_names <- c(
        "AuNP-MPTMA (Calibration)", 
        "AuNP-MEEE (Calibration)", 
        "AuNP-MHA (Calibration)",
        "Kumar et al., 2024 (24 nm)", 
        "Mishra & Das, 2022 (16 nm)", 
        "Mishra & Das, 2022 (28 nm)", 
        "Mishra & Das, 2022 (41 nm)", 
        "Mishra & Das, 2022 (69 nm)"
      )
      
      # Expected final sizes
      expected_sizes <- list(
        "AuNP-MPTMA (Calibration)" = 22.3,
        "AuNP-MEEE (Calibration)" = 13.7,
        "AuNP-MHA (Calibration)" = 13.3,
        "Kumar et al., 2024 (24 nm)" = 33,
        "Mishra & Das, 2022 (16 nm)" = 21.2,
        "Mishra & Das, 2022 (28 nm)" = 37.8,
        "Mishra & Das, 2022 (41 nm)" = 53.4,
        "Mishra & Das, 2022 (69 nm)" = 82.0
      )
      
      # Initial sizes
      initial_sizes <- list(
        "AuNP-MPTMA (Calibration)" = self$D_np_1,
        "AuNP-MEEE (Calibration)" = self$D_np_2,
        "AuNP-MHA (Calibration)" = self$D_np_3,
        "Kumar et al., 2024 (24 nm)" = self$D_np_4,
        "Mishra & Das, 2022 (16 nm)" = self$D_np_5,
        "Mishra & Das, 2022 (28 nm)" = self$D_np_6,
        "Mishra & Das, 2022 (41 nm)" = self$D_np_7,
        "Mishra & Das, 2022 (69 nm)" = self$D_np_8
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
      cat("\nAuNP-BSA Size Evolution Mathematical Model\n")
      cat("============================================\n\n")
      
      cat("Core Equation:\n")
      cat("dθ/dt = k_on * C_BSA * S * F_coop * F_crowd - k_off * θ\n\n")
      
      cat("Where:\n")
      cat("θ = BSA surface coverage (0 to 1)\n")
      cat("C_BSA = Initial BSA concentration (", self$c_init, ")\n")
      cat("k_on = Association rate constant (varies with size and time)\n")
      cat("k_off = Dissociation rate constant (varies with size)\n\n")
      
      cat("Modifying Factors:\n")
      cat("- Surface Availability Term: S = (1 - θ)^1.5\n")
      cat("- Cooperative Binding Term: F_coop = 1 + θ * (d_initial/10)^0.5\n")
      cat("- Crowding Effect Term: F_crowd = exp(-0.2 * θ)\n")
      cat("- Time-dependent k_on: k_on = k_on_base * size_effect * ((1 + 0.8t)/(0.2 + t))\n")
      cat("  where size_effect = (d_initial/10)^0.3\n")
      cat("- Size-dependent k_off: k_off = k_off_base / (1 + 0.2 * ln(d/10 + 1))\n\n")
      
      cat("Final Diameter:\n")
      cat("d(t) = d_initial + 2L * θ(t)\n")
      cat("where L = L_base * (1 + 0.2 * ln(d_initial/10 + 1))\n")
      cat("L_base =", self$L_base, "nm (calibrated value)\n\n")
      
      cat("Model Parameters:\n")
      cat("- k_on_base =", self$k_on_base, "\n")
      cat("- k_off_base =", self$k_off_base, "\n")
      cat("- BSA concentration =", self$c_init, "M\n")
    }
  )
)

# Function to run the simulation and create plots
run_aunp_bsa_simulation <- function() {
  # Create simulator instance
  simulator <- AuNP_BSASizeSimulator$new()
  
  # Show model details
  simulator$show_model_details()
  
  # Run simulation for all datasets
  results <- simulator$simulate(sim_time = 1.5, time_step = 0.01)
  

  # Generate validation report
  validation_results <- simulator$report_validation(results)
  
  # Return simulation results and plots
  return(list(
    simulator = simulator,
    results = results,
    validation = validation_results
  ))
}

# Function to create a summary table of all datasets
create_summary_table <- function(simulator) {
  # Create data frame with all datasets
  summary_df <- data.frame(
    Dataset = c(
      "AuNP-MPTMA (Calibration)", 
      "AuNP-MEEE (Calibration)",
      "AuNP-MHA (Calibration)",
      "Kumar et al., 2024 (24 nm)",
      "Mishra & Das, 2022 (16 nm)",
      "Mishra & Das, 2022 (28 nm)",
      "Mishra & Das, 2022 (41 nm)",
      "Mishra & Das, 2022 (69 nm)"
    ),
    Type = c(rep("Calibration", 3), rep("Validation", 5)),
    Initial_Size_nm = c(
      simulator$D_np_1, simulator$D_np_2, simulator$D_np_3,
      simulator$D_np_4, simulator$D_np_5, simulator$D_np_6, 
      simulator$D_np_7, simulator$D_np_8
    ),
    Final_Size_nm = c(22.3, 13.7, 13.3, 33.0, 21.2, 37.8, 53.4, 82.0),
    Measurement_Time_h = c(1.5, 1.5, 1.5, 1.5, 1/6, 1/6, 1/6, 1/6)
  )
  
  # Calculate size increase for each dataset
  summary_df$Size_Increase_nm <- summary_df$Final_Size_nm - summary_df$Initial_Size_nm
  summary_df$Size_Increase_Percent <- round((summary_df$Size_Increase_nm / summary_df$Initial_Size_nm) * 100, 1)
  
  return(summary_df)
}

# Run the simulation and store results
simulation_output <- run_aunp_bsa_simulation()

# Create and print summary table
summary_table <- create_summary_table(simulation_output$simulator)
print(summary_table)



