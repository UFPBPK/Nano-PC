### AuNP-HSA Size Evolution Model ###
### Author: Xinyue Chen a,b,c, Zhoumeng Lin a,b,c,*
##### a Department of Environmental and Global Health, College of Public Health and Health Professions, University of Florida, Gainesville, FL 32611, United States 
##### b Center for Environmental and Human Toxicology, University of Florida, Gainesville, FL 32611, United States 
##### c Center for Pharmacometrics and Systems Pharmacology, University of Florida, Orlando, FL 32827, United States 
##### * Corresponding author at: Department of Environmental and Global Health, College of Public Health and Health Professions, University of Florida, 2187 Mowry Rd, Gainesville, FL 32611, United States. Phone: +1-352-273-6160. 
##### E-mail address: linzhoumeng@ufl.edu (Z. Lin).


# Load required libraries
library(deSolve) # Purpose: Solving differential equations (ODE system for HSA corona formation kinetics on gold nanoparticles (AuNPs) across multiple size ranges)
library(dplyr)   # Purpose: Data manipulation and transformation (processing multi-dataset simulation results, filtering calibration vs validation data, and combining nanoparticle size series)
library(R6)      # Purpose: Object-oriented programming framework (creating the AuNP_HSASizeSimulator class with comprehensive dataset management, calibration methods, and validation reporting)


# Define the Corona Simulator class for AuNP-HSA
AuNP_HSASizeSimulator <- R6::R6Class(
  "AuNP_HSASizeSimulator",
  public = list(
    # Protein-specific parameters
    protein_sizes = list(HSA = 8),      # Unit (nm); Experimentally determined data
    k_on = list(HSA = 2.4e3),           # Base association rate constant for HSA (M⁻¹·s⁻¹ ),  cited from Dell’Orco et al. (2010)
    k_off = list(HSA = 2.0e-3),         # Base dissociation rate constant for HSA (s⁻¹ ), cited from Dell’Orco et al. (2010)
    c_init = list(HSA = 6e-4),          # Initial HSA concentration (M), cited from Dell’Orco et al. (2010)
    
    # Nanoparticle datasets and calibration parameters
    nanoparticle_datasets = NULL,
    L_base = NULL,     # Will be calculated from calibration data
    
    # Initialization method
    initialize = function() {
      # Define comprehensive nanoparticle datasets
      self$nanoparticle_datasets <- list(
        # Original calibration and validation datasets
        list(
          name = "Guglielmelli et al., 2023 (Calibration)",
          dataset_type = "Calibration",
          initial_size = 10.0,
          final_size = 21.6,
          measurement_time = 0.1,
          incubation_time = 1 # assume 1 hour here
        ),
        list(
          name = "Capomaccio et al., 2015 (18.5 nm)",
          dataset_type = "Validation",
          initial_size = 18.5,
          final_size = 24.4,
          measurement_time = 0.1,
          incubation_time = 5/60  # 5 minutes converted to hours
        ),
        # New Goy-Lopez et al., 2012 datasets with multiple nanoparticle sizes
        list(
          name = "Goy-Lopez et al., 2012 (5 nm)",
          dataset_type = "Validation",
          initial_size = 5.2,
          final_size = 10.2,
          measurement_time = 0.1,
          incubation_time = 1 # assume 1 hour here
        ),
        list(
          name = "Goy-Lopez et al., 2012 (10 nm)",
          dataset_type = "Validation",
          initial_size = 13.6,
          final_size = 26.3,
          measurement_time = 0.1,
          incubation_time = 1 # assume 1 hour here
        ),
        list(
          name = "Goy-Lopez et al., 2012 (20 nm)",
          dataset_type = "Validation",
          initial_size = 20.7,
          final_size = 33.4,
          measurement_time = 0.1,
          incubation_time = 1 # assume 1 hour here
        ),
        list(
          name = "Goy-Lopez et al., 2012 (40 nm)",
          dataset_type = "Validation",
          initial_size = 44.7,
          final_size = 58.4,
          measurement_time = 0.1,
          incubation_time = 1 # assume 1 hour here
        ),
        list(
          name = "Goy-Lopez et al., 2012 (60 nm)",
          dataset_type = "Validation",
          initial_size = 64.7,
          final_size = 79.6,
          measurement_time = 0.1,
          incubation_time = 1 # assume 1 hour here
        ),
        list(
          name = "Goy-Lopez et al., 2012 (80 nm)",
          dataset_type = "Validation",
          initial_size = 80.1,
          final_size = 95.3,
          measurement_time = 0.1,
          incubation_time = 1 # assume 1 hour here
        ),
        list(
          name = "Goy-Lopez et al., 2012 (100 nm)",
          dataset_type = "Validation",
          initial_size = 102.9,
          final_size = 118.6,
          measurement_time = 0.1,
          incubation_time = 1 # assume 1 hour here
        )
      )
      
      # Calibrate using the first dataset
      self$calibrate_model(self$nanoparticle_datasets[[1]])
    },
    
    # Calibration method
    calibrate_model = function(calibration_dataset) {
      # Calculate protein layer thickness from the calibration dataset
      final_size = calibration_dataset$final_size
      initial_size = calibration_dataset$initial_size
      
      self$L_base = (final_size - initial_size) / 2
      cat("Calibrated protein layer thickness:", round(self$L_base, 2), "nm\n")
    },
    
    # Differential equation model for single dataset simulation
    simulate_single = function(dataset_param, sim_time = 1, time_step = 0.05) {
      # Check if the model is calibrated
      if (is.null(self$L_base)) {
        stop("Model must be calibrated before simulation")
      }
      
      # Define the differential equations
      model <- function(t, state, parameters) {
        with(as.list(c(state, parameters)), {
          # Cooperative binding effect
          F_coop = 1 + theta * (d_initial/10)^0.5
          
          # Surface availability
          S = (1 - theta)^1.5
          
          # Crowding effect
          F_crowd = exp(-0.2 * theta)
          
          # Size effect on k_on
          size_effect = (d_initial/10)^0.3
          
          # Time effect on k_on
          time_effect = (1 + 0.8*t)/(0.2 + t)
          
          # Adjusted binding rate
          k_on = self$k_on$HSA * size_effect * time_effect
          
          # Adjusted dissociation rate
          k_off = self$k_off$HSA / (1 + 0.2 * log(d_initial/10 + 1))
          
          # Final equation: dtheta/dt = k_on * c_init * S * F_coop * F_crowd - k_off * theta
          dtheta = k_on * self$c_init$HSA * S * F_coop * F_crowd - k_off * theta
          
          # Make sure to return a list with names matching the state variables
          return(list(c(dtheta)))
        })
      }
      
      # Prepare simulation parameters
      times <- seq(0, sim_time, by = time_step)
      state <- c(theta = 0.0)  # Start with no coverage
      d_initial <- dataset_param$initial_size
      
      # Solve the ODE
      out <- ode(y = state, times = times, func = model, parms = list(d_initial = d_initial))
      
      # Convert to dataframe and add metadata
      data <- as.data.frame(out)
      names(data) <- c("time", "coverage")
      
      # Calculate adjusted protein layer thickness based on nanoparticle size
      L = self$L_base * (1 + 0.2 * log(d_initial/10 + 1))
      
      # Calculate diameter over time
      data <- data %>%
        mutate(
          time_h = time,
          diameter = d_initial + 2 * L * coverage,
          dataset = dataset_param$name,
          initial_size = d_initial,
          dataset_type = dataset_param$dataset_type
        )
      
      return(data)
    },
    
    # Simulate all datasets
    simulate = function(sim_time = 1, time_step = 0.05) {
      # Verify model is calibrated
      if (is.null(self$L_base)) {
        stop("Model must be calibrated before simulation")
      }
      
      # Run simulation for all datasets
      results_list <- lapply(self$nanoparticle_datasets, function(dataset) {
        self$simulate_single(dataset, sim_time, time_step)
      })
      
      # Combine all results
      results <- bind_rows(results_list)
      return(results)
    },
    

# Validation reporting method
report_validation = function(results) {
  # Create validation report
  val_results <- do.call(rbind, lapply(self$nanoparticle_datasets, function(ds) {
    # Find the closest time point to the measurement time
    dataset_results <- results %>% filter(dataset == ds$name)
    
    if (nrow(dataset_results) > 0) {
      # Find closest time point to incubation time
      closest_time_idx <- which.min(abs(dataset_results$time_h - ds$incubation_time))
      closest_time <- dataset_results[closest_time_idx, ]
      
      # Calculate error and size increase
      error_percent <- abs(closest_time$diameter - ds$final_size) / ds$final_size * 100
      size_increase <- closest_time$diameter - ds$initial_size
      
      data.frame(
        Dataset = ds$name,
        Type = ds$dataset_type,
        Initial_Size_nm = ds$initial_size,
        Final_Time_h = ds$incubation_time,
        Simulated_Size_nm = round(closest_time$diameter, 2),
        Experimental_Size_nm = ds$final_size,
        Error_percent = round(error_percent, 2),
        Size_Increase_nm = round(size_increase, 2)
      )
    }
  }))
  
  # Print detailed validation results
  cat("\nModel Validation Results:\n")
  cat("=======================\n")
  print(val_results)
  
  return(val_results)
},

# Model details method
show_model_details = function() {
  cat("\nAuNP-HSA Size Evolution Mathematical Model\n")
  cat("============================================\n\n")
  
  cat("Core Equation:\n")
  cat("dθ/dt = k_on * C_HSA * S * F_coop * F_crowd - k_off * θ\n\n")
  
  cat("Where:\n")
  cat("θ = HSA surface coverage (0 to 1)\n")
  cat("C_HSA = Initial HSA concentration (", self$c_init$HSA, ")\n")
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
  cat("L_base =", round(self$L_base, 2), "nm (calibrated value)\n\n")
  
  cat("Model Parameters:\n")
  cat("- k_on_base =", self$k_on$HSA, "\n")
  cat("- k_off_base =", self$k_off$HSA, "\n")
  cat("- HSA concentration =", self$c_init$HSA, "M\n")
}
  )
)

# Function to run the simulation and create plots
run_aunp_hsa_simulation <- function() {
  # Create simulator instance
  simulator <- AuNP_HSASizeSimulator$new()
  
  # Show model details
  simulator$show_model_details()
  
  # Run simulation for all datasets
  results <- simulator$simulate(sim_time = 1, time_step = 0.05)

  
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
  # Convert nanoparticle datasets to a data frame
  summary_df <- do.call(rbind, lapply(simulator$nanoparticle_datasets, function(dataset) {
    data.frame(
      Dataset = dataset$name,
      Type = dataset$dataset_type,
      Initial_Size_nm = dataset$initial_size,
      Final_Size_nm = dataset$final_size,
      Measurement_Time_h = dataset$incubation_time,
      Size_Increase_nm = dataset$final_size - dataset$initial_size,
      Size_Increase_Percent = round(((dataset$final_size - dataset$initial_size) / dataset$initial_size) * 100, 1)
    )
  }))
  
  return(summary_df)
}

# Run the simulation and store results
simulation_output <- run_aunp_hsa_simulation()

# Create and print summary table
summary_table <- create_summary_table(simulation_output$simulator)
print(summary_table)






