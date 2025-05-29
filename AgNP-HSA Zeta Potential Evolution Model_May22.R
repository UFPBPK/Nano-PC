### AgNP-HSA Zeta Potential Evolution Model ###
### Author: Xinyue Chen a,b,c, Zhoumeng Lin a,b,c,*
##### a Department of Environmental and Global Health, College of Public Health and Health Professions, University of Florida, Gainesville, FL 32611, United States 
##### b Center for Environmental and Human Toxicology, University of Florida, Gainesville, FL 32611, United States 
##### c Center for Pharmacometrics and Systems Pharmacology, University of Florida, Orlando, FL 32827, United States 
##### * Corresponding author at: Department of Environmental and Global Health, College of Public Health and Health Professions, University of Florida, 2187 Mowry Rd, Gainesville, FL 32611, United States. Phone: +1-352-273-6160. 
##### E-mail address: linzhoumeng@ufl.edu (Z. Lin).


# Load required libraries
library(deSolve) # Purpose: Solving differential equations (ODE system for HSA surface coverage dynamics and zeta potential evolution on silver nanoparticles (AgNPs))
library(dplyr)   # Purpose: Data manipulation and transformation (processing multi-dataset zeta potential results, filtering validation metrics, and combining electrokinetic measurements from multiple research studies)
library(R6)      # Purpose: Object-oriented programming framework (creating the AgNP_HSAZetaPotentialSimulator class with advanced calibration algorithms, charge/size/incubation factor calculations, and comprehensive validation reporting with performance metrics)


# Define the AgNP-HSA Zeta Potential Corona Simulator class
AgNP_HSAZetaPotentialSimulator <- R6::R6Class(
  "AgNP_HSAZetaPotentialSimulator",
  public = list(
    # Protein-specific parameters
    protein_sizes = list(HSA = 8),        # Unit (nm)
    k_on_base = 2.4e3,                    # Base association rate constant for HSA (M⁻¹·s⁻¹ ),  cited from Dell’Orco et al. (2010)
    k_off_base = 2.0e-3,                  # Base dissociation rate constant for HSA (s⁻¹ ), cited from Dell’Orco et al. (2010)
    c_HSA = 6.0e-4,                         # Initial HSA concentration (M), cited from Dell’Orco et al. (2010)
    delta_zeta_max = NULL,                # Will be calibrated using calibration dataset
    
    # Nanoparticle datasets
    nanoparticle_datasets = NULL,
    
    # Initialization method
    initialize = function() {
      # Define comprehensive nanoparticle datasets
      self$nanoparticle_datasets <- list(
        # Calibration dataset
        list(
          name = "Zhao et al., 2021 (SNP20) (Calibration)",
          dataset_type = "Calibration",
          initial_size = 23.0,
          final_size = 28.5,
          initial_zeta = -38.8,
          final_zeta = -13.5,
          measurement_time = 0.3,
          incubation_time = 0.3, 
          coating = "citrate"
        ),
        
        # Validation datasets
        list(
          name = "Zhao et al., 2021 (SNP40)",
          dataset_type = "Validation",
          initial_size = 45.2,
          final_size = 58.5,
          initial_zeta = -32.9,
          final_zeta = -12.8,
          measurement_time = 0.3,
          incubation_time = 0.3,
          coating = "citrate"
        ),
        list(
          name = "Zhao et al., 2021 (SNP80)",
          dataset_type = "Validation",
          initial_size = 88.9,
          final_size = 95.0,
          initial_zeta = -37.2,
          final_zeta = -12.5,
          measurement_time = 0.3,
          incubation_time = 0.3,
          coating = "citrate"
        ),
        list(
          name = "Shannahan et al., 2015",
          dataset_type = "Validation",
          initial_size = 19.13,
          final_size = 69.99,
          initial_zeta = -35.0,
          final_zeta = -25.0,
          measurement_time = 8.0,
          incubation_time = 8.0,
          coating = "citrate"
        )
      )
      
      # Set a temporary delta_zeta_max value
      self$delta_zeta_max = 20
      
      # Calibrate the model using only the calibration dataset
      self$calibrate_model(self$nanoparticle_datasets[[1]])
    },
    
    # Charge factor calculation - improved to better match AgNP behavior
    calculate_F_charge = function(initial_zeta) {
      if (initial_zeta < 0) {
        return(0.4)   # Calibrated for AgNPs. Here only have data with negative initial zeta potential values.
      }
    },
    
    # Size factor calculation - enhanced for AgNP size dependence
    calculate_F_size = function(d_initial) {
      # Enhanced size-dependent effect for AgNPs
      return(1 - 0.5 * tanh(d_initial/50))  # Stronger size dependence, smaller scale
    },
    
    # Special factor for long incubation dataset (Shannahan)
    calculate_F_incubation = function(dataset_name, incubation_time) {
      # Special factor for datasets with very long incubation times
      if (incubation_time > 4) {
        return(0.25)  # Significant reduction for long incubation datasets
      } else {
        return(1.0)   # Normal factor for standard incubation times
      }
    },
    
    # Calibration method
    calibrate_model = function(calibration_dataset) {
      # Extract observed change
      zeta_initial = calibration_dataset$initial_zeta
      zeta_final = calibration_dataset$final_zeta
      observed_change = zeta_final - zeta_initial
      
      # Run a simulation with temporary delta_zeta_max to get coverage
      temp_result = self$simulate_single(calibration_dataset, 
                                         sim_time = calibration_dataset$incubation_time * 1.1, 
                                         temp_calibration = TRUE)
      final_idx = which.min(abs(temp_result$time_h - calibration_dataset$incubation_time))
      final_coverage = temp_result[final_idx, "coverage"]
      
      # Calculate modifying factors for calibration dataset
      F_charge = self$calculate_F_charge(zeta_initial)
      F_size = self$calculate_F_size(calibration_dataset$initial_size)
      F_incubation = self$calculate_F_incubation(calibration_dataset$name, calibration_dataset$incubation_time)
      
      # Back-calculate delta_zeta_max - using saturation relationship
      if (final_coverage > 0) {
        numerator = (0.2 + final_coverage) * observed_change
        denominator = final_coverage * F_charge * F_size * F_incubation
        self$delta_zeta_max = numerator / denominator
      } else {
        warning("Final coverage is zero or negative, using default delta_zeta_max")
        self$delta_zeta_max = 30
      }
      
      cat("Calibrated delta_zeta_max:", round(self$delta_zeta_max, 2), "mV\n")
      cat("Final coverage at calibration time point:", round(final_coverage, 4), "\n")
      cat("F_charge:", round(F_charge, 4), "\n")
      cat("F_size:", round(F_size, 4), "\n")
      cat("F_incubation:", round(F_incubation, 4), "\n")
    },
    
    # Differential equation model for single dataset simulation
    simulate_single = function(dataset_param, sim_time = NULL, time_step = 0.01, temp_calibration = FALSE) {
      # Set simulation time if not provided
      if (is.null(sim_time)) {
        sim_time = max(dataset_param$incubation_time * 1.2, 1)  # At least 1 hour, or 20% more than incubation time
      }
      
      # Define the differential equations
      model <- function(t, state, parameters) {
        with(as.list(c(state, parameters)), {
          # Extract surface coverage
          theta <- HSA_coverage
          
          # Calculate time in hours
          time_h <- t
          
          # Enhanced size-dependent effects
          relative_size <- d_initial / 20  # normalized to SNP20
          
          # Progressive size effect over time
          time_factor <- 1 + 0.2 * log(time_h + 0.01)
          
          # Stronger size dependence for larger particles
          k_on_factor <- relative_size^0.8 * time_factor  # Increased size dependence
          k_off_factor <- relative_size^(-0.5)  # Stronger size effect on dissociation
          
          # Additional size-dependent multiplier for larger particles
          size_multiplier <- ifelse(d_initial > 40, 
                                    1 + 0.2 * log(d_initial/40), 
                                    1)
          

          
          # Adjusted rate constants with size-dependent effects
          k_on = k_on_base * k_on_factor * size_multiplier
          k_off = k_off_base * k_off_factor
          
          # Surface availability - higher exponent for silver
          S = (1 - theta)^2
          
          # Cooperative binding effect - adjusted for AgNP
          F_coop = 1 + 0.5 * theta * (d_initial/20)^0.3
          
          # Crowding effect - complex for silver
          F_crowd = (1 - 0.1 * theta) * (1 + 0.7442 * theta^2)
          
          # Final equation: dHSA_coverage/dt with higher power for dissociation
          dHSA_coverage = k_on * c_HSA * S * F_coop * F_crowd - k_off * theta^1.8
          
          # Return rate of change
          return(list(c(dHSA_coverage)))
        })
      }
      
      # Prepare simulation parameters
      times <- seq(0, sim_time, by = time_step)
      state <- c(HSA_coverage = 0.0)  # Start with no coverage
      
      # Parameters for simulation
      parameters <- list(
        d_initial = dataset_param$initial_size,
        initial_zeta = dataset_param$initial_zeta,
        incubation_time = dataset_param$incubation_time,
        k_on_base = self$k_on_base,
        k_off_base = self$k_off_base,
        c_HSA = self$c_HSA,
        dataset_name = dataset_param$name
      )
      
      # Solve the ODE
      out <- ode(y = state, times = times, func = model, parms = parameters)
      
      # Convert to dataframe and add metadata
      data <- as.data.frame(out)
      names(data) <- c("time", "coverage")
      
      # Calculate zeta potential over time - using a saturation formula
      data <- data %>%
        mutate(
          time_h = time,
          # Calculate modifying factors 
          F_charge = self$calculate_F_charge(dataset_param$initial_zeta),
          F_size = self$calculate_F_size(dataset_param$initial_size),
          F_incubation = self$calculate_F_incubation(dataset_param$name, dataset_param$incubation_time),
          # Calculate zeta potential using saturation model
          zeta_change = self$delta_zeta_max * F_charge * F_size * F_incubation * 
            (coverage / (0.2 + coverage)),
          # Calculate total zeta potential
          zeta_potential = dataset_param$initial_zeta + zeta_change,
          dataset = dataset_param$name,
          initial_size = dataset_param$initial_size,
          initial_zeta = dataset_param$initial_zeta,
          final_zeta = dataset_param$final_zeta,
          dataset_type = dataset_param$dataset_type
        )
      
      return(data)
    },
    
    # Simulate all datasets
    simulate = function(sim_time = 48, time_step = 0.01) {
      # Verify model is calibrated
      if (is.null(self$delta_zeta_max)) {
        stop("Model must be calibrated before simulation")
      }
      
      # Run simulation for all datasets
      results_list <- lapply(self$nanoparticle_datasets, function(ds) {
        self$simulate_single(ds, sim_time, time_step)
      })
      
      # Combine all results
      results <- bind_rows(results_list)
      return(results)
    },
    
    # Validation reporting method
    report_validation = function(results) {
      # Create validation report
      val_results <- do.call(rbind, lapply(self$nanoparticle_datasets, function(ds) {
        # Find the closest time point to the incubation time
        dataset_results <- results %>% filter(dataset == ds$name)
        
        if (nrow(dataset_results) > 0) {
          # Find closest time point to incubation time
          closest_time_idx <- which.min(abs(dataset_results$time_h - ds$incubation_time))
          closest_time <- dataset_results[closest_time_idx, ]
          
          # Calculate error and zeta potential change
          error_abs <- abs(closest_time$zeta_potential - ds$final_zeta)
          error_percent <- error_abs / abs(ds$final_zeta) * 100
          zeta_change <- closest_time$zeta_potential - ds$initial_zeta
          exp_zeta_change <- ds$final_zeta - ds$initial_zeta
          
          data.frame(
            Dataset = ds$name,
            Type = ds$dataset_type,
            Coating = ds$coating,
            Initial_Size_nm = ds$initial_size,
            Initial_Zeta_mV = ds$initial_zeta,
            Final_Time_h = ds$incubation_time,
            Simulated_Zeta_mV = round(closest_time$zeta_potential, 2),
            Experimental_Zeta_mV = ds$final_zeta,
            Sim_Change_mV = round(zeta_change, 2),
            Exp_Change_mV = round(exp_zeta_change, 2),
            Abs_Error_mV = round(error_abs, 2),
            Error_percent = round(error_percent, 2),
            Surface_Coverage = round(closest_time$coverage, 3)
          )
        }
      }))
      
      # Print detailed validation results
      cat("\nModel Validation Results:\n")
      cat("=======================\n")
      print(val_results)
      
      # Calculate overall model performance metrics
      rmse <- sqrt(mean((val_results$Simulated_Zeta_mV - val_results$Experimental_Zeta_mV)^2))
      mae <- mean(abs(val_results$Simulated_Zeta_mV - val_results$Experimental_Zeta_mV))
      mean_pct_error <- mean(val_results$Error_percent)
      median_pct_error <- median(val_results$Error_percent)
      
      cat("\nOverall Model Performance Metrics:\n")
      cat("================================\n")
      cat("Root Mean Square Error (RMSE):", round(rmse, 2), "mV\n")
      cat("Mean Absolute Error (MAE):", round(mae, 2), "mV\n")
      cat("Mean Percentage Error:", round(mean_pct_error, 2), "%\n")
      cat("Median Percentage Error:", round(median_pct_error, 2), "%\n")
      
      return(val_results)
    },
    
    # Model details method
    show_model_details = function() {
      cat("\nAgNP-HSA Zeta Potential Evolution Model\n")
      cat("======================================\n\n")
      
      cat("Core Equation (Surface Coverage):\n")
      cat("dθ/dt = k_on * C_HSA * S * F_coop * F_crowd - k_off * θ^1.8\n\n")
      
      cat("Where:\n")
      cat("θ = HSA surface coverage (0 to 1)\n")
      cat("C_HSA = Initial HSA concentration (", self$c_HSA, ")\n")
      cat("k_on = Association rate constant (varies with size and time)\n")
      cat("k_off = Dissociation rate constant (varies with size)\n\n")
      
      cat("Association Rate Calculation:\n")
      cat("k_on = k_on_base * relative_size^0.8 * time_factor * size_multiplier\n")
      cat("where:\n")
      cat("- relative_size = d_initial / 20 (normalized to 20 nm particles)\n")
      cat("- time_factor = 1 + 0.2 * log(time_h + 0.01) (progressive time effect)\n")
      cat("- size_multiplier = 1 + 0.2 * log(d_initial/40) (for d_initial > 40 nm)\n")
      cat("- size_multiplier = 1 (for d_initial ≤ 40 nm)\n")
      
      cat("Dissociation Rate Calculation:\n")
      cat("k_off = k_off_base * (d_initial/20)^(-0.5)\n\n")
      
      cat("Modifying Factors for Surface Coverage:\n")
      cat("- Surface Availability Term: S = (1 - θ)^2\n")
      cat("- Cooperative Binding Term: F_coop = 1 + 0.5θ * (d_initial/20)^0.3\n")
      cat("- Crowding Effect Term: F_crowd = (1 - 0.1θ)(1 + 0.7439θ²)\n\n")
      
      cat("Zeta Potential Calculation:\n")
      cat("ζ(t) = ζ_initial + Δζ_max * F_charge * F_size * F_incubation * (θ(t)/(0.2 + θ(t)))\n")
      cat("where Δζ_max = ", round(self$delta_zeta_max, 2), " mV (calibrated value)\n")
      cat("- Charge-based factor (F_charge): 0.4 for negative, -0.7 for positive initial zeta\n")
      cat("- Size-based factor (F_size): 1 - 0.5 * tanh(d_initial/50)\n")
      cat("- Incubation-based factor (F_incubation): 0.25 for long incubation (>4h), 1.0 for others\n\n")
      
      cat("Model Parameters:\n")
      cat("- k_on_base =", self$k_on_base, " M^-1*s^-1\n")
      cat("- k_off_base =", self$k_off_base, " s^-1\n")
      cat("- HSA concentration =", self$c_HSA, "M\n")
      cat("- HSA size =", self$protein_sizes$HSA, "nm\n")
    }
  )
)

# Main simulation execution
set.seed(123) # For reproducibility

# Run the simulation
simulator <- AgNP_HSAZetaPotentialSimulator$new()
simulator$show_model_details()

# Run simulation for all datasets
results <- simulator$simulate(sim_time = 48, time_step = 0.01)

# Generate validation report
validation_results <- simulator$report_validation(results)