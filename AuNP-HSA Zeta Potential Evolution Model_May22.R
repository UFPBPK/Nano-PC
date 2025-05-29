### AuNP-HSA Zeta Potential Evolution Model ###
### Author: Xinyue Chen a,b,c, Zhoumeng Lin a,b,c,*
##### a Department of Environmental and Global Health, College of Public Health and Health Professions, University of Florida, Gainesville, FL 32611, United States 
##### b Center for Environmental and Human Toxicology, University of Florida, Gainesville, FL 32611, United States 
##### c Center for Pharmacometrics and Systems Pharmacology, University of Florida, Orlando, FL 32827, United States 
##### * Corresponding author at: Department of Environmental and Global Health, College of Public Health and Health Professions, University of Florida, 2187 Mowry Rd, Gainesville, FL 32611, United States. Phone: +1-352-273-6160. 
##### E-mail address: linzhoumeng@ufl.edu (Z. Lin).


# Load required libraries
library(deSolve) # Purpose: Solving differential equations (ODE system for HSA surface coverage dynamics and comprehensive zeta potential evolution modeling on gold nanoparticles (AuNPs) with multiple coating types)
library(dplyr)   # Purpose: Data manipulation and transformation (processing multi-study zeta potential datasets, performing coating-specific analysis, charge-type grouping, size-range categorization, and generating comprehensive validation statistics across diverse experimental conditions)
library(R6)      # Purpose: Object-oriented programming framework (creating the ZetaPotentialSimulator class with advanced multi-coating calibration, charge-dependent factor calculations, coating-specific effects modeling, and sophisticated validation reporting with performance metrics stratified by coating type, charge polarity, and nanoparticle size ranges)


# Define the Zeta Potential Corona Simulator class for AuNP-HSA
ZetaPotentialSimulator <- R6::R6Class(
  "ZetaPotentialSimulator",
  public = list(
    # Protein-specific parameters
    protein_sizes = list(HSA = 8),      # Unit (nm)
    k_on = list(HSA = 2.4e3),           # Base association rate constant for HSA (M⁻¹·s⁻¹ ),  cited from Dell’Orco et al. (2010)
    k_off = list(HSA = 2.0e-3),           # Base dissociation rate constant for HSA (s⁻¹ ), cited from Dell’Orco et al. (2010)
    c_init = list(HSA = 6.0e-4),          # Initial HSA concentration (M), cited from Dell’Orco et al. (2010)
    
    # Nanoparticle datasets and calibration parameters
    nanoparticle_datasets = NULL,
    delta_zeta_max = NULL,              # Maximum zeta potential change (to be calibrated)
    
    # Initialization method
    initialize = function() {
      # Define comprehensive nanoparticle datasets
      self$nanoparticle_datasets <- list(
        # Calibration dataset - Dai et al. 2023 is the calibration dataset.
        list(
          name = "Dai et al., 2023 (Calibration)",
          dataset_type = "Calibration",
          initial_size = 28.20,
          final_size = 164.20,
          initial_zeta = 60.60,
          final_zeta = 26.27,
          measurement_time = 0.1,
          incubation_time = 24,
          coating = "CTAB"  
        ),
        
        # Validation datasets
        list(
          name = "Dell'Aglio et al., 2021 (10.5 nmol/L)",
          dataset_type = "Validation",
          initial_size = 10.00,
          final_size = 21.70,
          initial_zeta = -48.93,
          final_zeta = -28.25,
          measurement_time = 0.1,
          incubation_time = 24,
          coating = "citrate"
        ),
        list(
          name = "Dell'Aglio et al., 2021 (28.5 nmol/L)",
          dataset_type = "Validation",
          initial_size = 10.00,
          final_size = 10.00,  # No size data provided
          initial_zeta = -41.81,
          final_zeta = -22.37,
          measurement_time = 0.1,
          incubation_time = 24,
          coating = "citrate"
        ),
        
        # Goy-Lopez datasets
        list(
          name = "Goy-Lopez et al., 2012 (5 nm)",
          dataset_type = "Validation",
          initial_size = 5.20,
          final_size = 10.20,
          initial_zeta = -50.80,
          final_zeta = -37.40,
          measurement_time = 0.1,
          incubation_time = 48,
          coating = "citrate"
        ),
        list(
          name = "Goy-Lopez et al., 2012 (10 nm)",
          dataset_type = "Validation",
          initial_size = 13.60,
          final_size = 26.30,
          initial_zeta = -57.30,
          final_zeta = -29.80,
          measurement_time = 0.1,
          incubation_time = 48,
          coating = "citrate"
        ),
        list(
          name = "Goy-Lopez et al., 2012 (20 nm)",
          dataset_type = "Validation",
          initial_size = 20.70,
          final_size = 33.40,
          initial_zeta = -50.10,
          final_zeta = -32.20,
          measurement_time = 0.1,
          incubation_time = 48,
          coating = "citrate"
        ),
        list(
          name = "Goy-Lopez et al., 2012 (40 nm)",
          dataset_type = "Validation",
          initial_size = 44.70,
          final_size = 58.40,
          initial_zeta = -56.80,
          final_zeta = -31.10,
          measurement_time = 0.1,
          incubation_time = 48,
          coating = "citrate"
        ),
        list(
          name = "Goy-Lopez et al., 2012 (60 nm)",
          dataset_type = "Validation",
          initial_size = 64.70,
          final_size = 79.60,
          initial_zeta = -55.10,
          final_zeta = -30.90,
          measurement_time = 0.1,
          incubation_time = 48,
          coating = "citrate"
        ),
        list(
          name = "Goy-Lopez et al., 2012 (80 nm)",
          dataset_type = "Validation",
          initial_size = 80.10,
          final_size = 95.30,
          initial_zeta = -54.80,
          final_zeta = -29.10,
          measurement_time = 0.1,
          incubation_time = 48,
          coating = "citrate"
        ),
        list(
          name = "Goy-Lopez et al., 2012 (100 nm)",
          dataset_type = "Validation",
          initial_size = 102.90,
          final_size = 118.60,
          initial_zeta = -52.70,
          final_zeta = -26.20,
          measurement_time = 0.1,
          incubation_time = 48,
          coating = "citrate"
        )
      )
      
      # Set a temporary delta_zeta_max value
      self$delta_zeta_max = 20
      
      # Calibrate using the first dataset
      self$calibrate_model(self$nanoparticle_datasets[[1]])
    },
    
    # Charge factor calculation
    calculate_F_charge = function(initial_zeta) {
      if (initial_zeta > 0) {
        # For positive initial charge (protein binding decreases the charge)
        return(-0.8)  # Negative factor to reduce positive charge
      } else {
        # For negative initial charge (protein binding makes it less negative)
        return(0.6)   # Positive factor to make it less negative
      }
    },
    
    # Size factor calculation
    calculate_F_size = function(d_initial) {
      # Larger particles show less relative zeta potential change
      return(1 - 0.3 * tanh(d_initial/100))  # Range: 0.7-1.0
    },
    
    # Coating factor calculation
    calculate_F_coat = function(coating) {
      if (coating == "CTAB") {
        return(1.2)  # Enhanced zeta change for positive CTAB
      } else {
        return(1.0)  # Standard for citrate/default coatings
      }
    },
    
    # Calibration method
    calibrate_model = function(calibration_dataset) {
      # Extract observed change
      zeta_initial = calibration_dataset$initial_zeta
      zeta_final = calibration_dataset$final_zeta
      observed_change = zeta_final - zeta_initial
      
      # Run a simulation with temporary delta_zeta_max to get coverage
      temp_result = self$simulate_single(calibration_dataset, sim_time = calibration_dataset$incubation_time * 1.1, temp_calibration = TRUE)
      final_idx = which.min(abs(temp_result$time_h - calibration_dataset$incubation_time))
      final_coverage = temp_result[final_idx, "coverage"]
      
      # Calculate F_charge for calibration dataset
      F_charge = self$calculate_F_charge(zeta_initial)
      
      # Calculate F_size for calibration dataset
      F_size = self$calculate_F_size(calibration_dataset$initial_size)
      
      # Calculate F_coating for calibration dataset
      F_coat = self$calculate_F_coat(calibration_dataset$coating)
      
      # Back-calculate delta_zeta_max
      if (final_coverage > 0) {
        self$delta_zeta_max = observed_change / (F_charge * F_size * F_coat * final_coverage)
      } else {
        warning("Final coverage is zero or negative, using default delta_zeta_max")
        self$delta_zeta_max = 30
      }
      
      cat("Calibrated delta_zeta_max:", round(self$delta_zeta_max, 2), "mV\n")
      cat("Final coverage at calibration time point:", round(final_coverage, 4), "\n")
      cat("F_charge:", round(F_charge, 4), "\n")
      cat("F_size:", round(F_size, 4), "\n")
      cat("F_coat:", round(F_coat, 4), "\n")
    },
    
    # Differential equation model for single dataset simulation
    simulate_single = function(dataset_param, sim_time = NULL, time_step = 0.1, temp_calibration = FALSE) {
      # Set simulation time if not provided
      if (is.null(sim_time)) {
        sim_time = max(dataset_param$incubation_time * 1.2, 1)  # At least 1 hour, or 20% more than incubation time
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
          
          # Coating effect on k_on
          if (coating == "CTAB") {
            coat_effect = 1.2  # Enhance binding for positive CTAB
          } else {
            coat_effect = 1.0  # Default for citrate coating
          }
          
          # Adjusted binding rate
          k_on = self$k_on$HSA * size_effect * time_effect * coat_effect
          
          # Adjusted dissociation rate - size and charge dependent
          k_off = self$k_off$HSA / (1 + 0.2 * log(d_initial/10 + 1))
          
          # Modification for charge effects on dissociation rate
          if (abs(initial_zeta) > 40) {
            # High absolute charge decreases dissociation (stronger electrostatic binding)
            k_off = k_off * 0.8
          } else if (abs(initial_zeta) < 10) {
            # Low absolute charge increases dissociation (weaker binding)
            k_off = k_off * 1.2
          }
          
          # Final equation: dtheta/dt
          dtheta = k_on * self$c_init$HSA * S * F_coop * F_crowd - k_off * theta
          
          # Return rate of change
          return(list(c(dtheta)))
        })
      }
      
      # Prepare simulation parameters
      times <- seq(0, sim_time, by = time_step)
      state <- c(theta = 0.0)  # Start with no coverage
      
      # Parameters for simulation
      parameters <- list(
        d_initial = dataset_param$initial_size,
        initial_zeta = dataset_param$initial_zeta,
        coating = dataset_param$coating
      )
      
      # Solve the ODE
      out <- ode(y = state, times = times, func = model, parms = parameters)
      
      # Convert to dataframe and add metadata
      data <- as.data.frame(out)
      names(data) <- c("time", "coverage")
      
      # Calculate zeta potential over time
      data <- data %>%
        mutate(
          time_h = time,
          F_charge = self$calculate_F_charge(dataset_param$initial_zeta),
          F_size = self$calculate_F_size(dataset_param$initial_size),
          F_coat = self$calculate_F_coat(dataset_param$coat),
          # Calculate zeta potential change
          zeta_change = self$delta_zeta_max * F_charge * F_size * F_coat * coverage,
          # Calculate total zeta potential
          zeta_potential = dataset_param$initial_zeta + zeta_change,
          dataset = dataset_param$name,
          initial_size = dataset_param$initial_size,
          initial_zeta = dataset_param$initial_zeta,
          dataset_type = dataset_param$dataset_type,
          coating = dataset_param$coating
        )
      
      return(data)
    },
    
    # Simulate all datasets
    simulate = function(sim_time = 48, time_step = 0.1) {
      # Verify model is calibrated
      if (is.null(self$delta_zeta_max)) {
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
      
      # Separate performance metrics for positive and negative initial zeta
      pos_results <- val_results[val_results$Initial_Zeta_mV > 0, ]
      neg_results <- val_results[val_results$Initial_Zeta_mV < 0, ]
      
      if (nrow(pos_results) > 0) {
        pos_rmse <- sqrt(mean((pos_results$Simulated_Zeta_mV - pos_results$Experimental_Zeta_mV)^2))
        pos_mean_pct_error <- mean(pos_results$Error_percent)
        cat("\nPositive Initial Zeta Potential (n=", nrow(pos_results), "):\n", sep="")
        cat("RMSE:", round(pos_rmse, 2), "mV, Mean % Error:", round(pos_mean_pct_error, 2), "%\n")
      }
      
      if (nrow(neg_results) > 0) {
        neg_rmse <- sqrt(mean((neg_results$Simulated_Zeta_mV - neg_results$Experimental_Zeta_mV)^2))
        neg_mean_pct_error <- mean(neg_results$Error_percent)
        cat("Negative Initial Zeta Potential (n=", nrow(neg_results), "):\n", sep="")
        cat("RMSE:", round(neg_rmse, 2), "mV, Mean % Error:", round(neg_mean_pct_error, 2), "%\n")
      }
      
      return(val_results)
    },
    
    # Model details method
    show_model_details = function() {
      cat("\nAuNP-HSA Zeta Potential Evolution Mathematical Model\n")
      cat("==================================================\n\n")
      
      cat("Core Equation (Surface Coverage):\n")
      cat("dθ/dt = k_on * C_HSA * S * F_coop * F_crowd - k_off * θ\n\n")
      
      cat("Where:\n")
      cat("θ = HSA surface coverage (0 to 1)\n")
      cat(paste0("C_HSA = Initial HSA concentration (", self$c_init$HSA, ")\n"))
      cat("k_on = Association rate constant (varies with size, time, and coating)\n")
      cat("k_off = Dissociation rate constant (varies with size and initial charge)\n\n")
      
      cat("Modifying Factors for Surface Coverage:\n")
      cat("- Surface Availability Term: S = (1 - θ)^1.5\n")
      cat("- Cooperative Binding Term: F_coop = 1 + θ * (d_initial/10)^0.5\n")
      cat("- Crowding Effect Term: F_crowd = exp(-0.2 * θ)\n")
      cat("- Time-dependent k_on: k_on = k_on_base * size_effect * time_effect * coat_effect\n")
      cat("  where size_effect = (d_initial/10)^0.3\n")
      cat("  where time_effect = (1 + 0.8t)/(0.2 + t)\n")
      cat("  where coat_effect = 1.2 for CTAB or 1.0 for citrate\n")
      cat("- Size and charge-dependent k_off: k_off = k_off_base / (1 + 0.2 * log(d/10 + 1))\n")
      cat("  with adjustments based on initial zeta potential\n\n")
      
      cat("Zeta Potential Calculation:\n")
      cat("ζ(t) = ζ_initial + Δζ_max * F_charge * F_size * F_coat * θ(t)\n")
      cat(paste0("where Δζ_max = ", round(self$delta_zeta_max, 2), " mV (calibrated value)\n"))
      cat("Charge-based factor: F_charge = 0.6 for negative, -0.8 for positive initial zeta\n")
      cat("Size-based factor: F_size = 1 - 0.3 * tanh(d_initial/100)\n")
      cat("Coating-specific factor: F_coat = 1.2 for CTAB or 1.0 for citrate\n\n")
      
      cat("Model Parameters:\n")
      cat(paste0("- k_on_base = ", self$k_on$HSA, " M^-1*s^-1\n"))
      cat(paste0("- k_off_base = ", self$k_off$HSA, " s^-1\n"))
      cat(paste0("- HSA concentration = ", self$c_init$HSA, " M\n"))
    }
  )
)

# Function to run the simulation and create plots
run_aunp_hsa_zeta_simulation <- function() {
  # Create simulator instance
  simulator <- ZetaPotentialSimulator$new()
  
  # Show model details
  simulator$show_model_details()
  
  # Run simulation for all datasets
  results <- simulator$simulate(sim_time = 48, time_step = 0.1)
  
  # Generate validation report
  validation_results <- simulator$report_validation(results)
  
  # Return simulation results and plots
  return(list(
    simulator = simulator,
    results = results,
    validation = validation_results
  ))
}

# Create function to summarize validation results
summarize_validation <- function(validation_results) {
  # Group by coating type
  coating_summary <- validation_results %>%
    group_by(Coating) %>%
    summarize(
      Count = n(),
      Mean_Error_mV = mean(Abs_Error_mV),
      Median_Error_mV = median(Abs_Error_mV),
      Mean_Error_Pct = mean(Error_percent),
      Median_Error_Pct = median(Error_percent)
    )
  
  # Group by initial charge type (positive/negative)
  charge_summary <- validation_results %>%
    mutate(Charge_Type = ifelse(Initial_Zeta_mV > 0, "Positive", "Negative")) %>%
    group_by(Charge_Type) %>%
    summarize(
      Count = n(),
      Mean_Error_mV = mean(Abs_Error_mV),
      Median_Error_mV = median(Abs_Error_mV),
      Mean_Error_Pct = mean(Error_percent),
      Median_Error_Pct = median(Error_percent)
    )
  
  # Group by size ranges
  size_summary <- validation_results %>%
    mutate(Size_Range = case_when(
      Initial_Size_nm < 20 ~ "<20 nm",
      Initial_Size_nm < 50 ~ "20-50 nm",
      Initial_Size_nm < 80 ~ "50-80 nm",
      TRUE ~ ">80 nm"
    )) %>%
    group_by(Size_Range) %>%
    summarize(
      Count = n(),
      Mean_Error_mV = mean(Abs_Error_mV),
      Median_Error_mV = median(Abs_Error_mV),
      Mean_Error_Pct = mean(Error_percent),
      Median_Error_Pct = median(Error_percent)
    )
  
  return(list(
    by_coating = coating_summary,
    by_charge = charge_summary,
    by_size = size_summary
  ))
}

# Run the simulation
simulation_output <- run_aunp_hsa_zeta_simulation()

# Create summary statistics
validation_summary <- summarize_validation(simulation_output$validation)
print(validation_summary$by_coating)
print(validation_summary$by_charge)
print(validation_summary$by_size)