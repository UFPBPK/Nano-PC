### Parameter Optimization Script for AgNP-HSA Corona Formation Model ###
### Author: Xinyue Chen a,b,c, Zhoumeng Lin a,b,c,*
##### a Department of Environmental and Global Health, College of Public Health and Health Professions, University of Florida, Gainesville, FL 32611, United States 
##### b Center for Environmental and Human Toxicology, University of Florida, Gainesville, FL 32611, United States 
##### c Center for Pharmacometrics and Systems Pharmacology, University of Florida, Orlando, FL 32827, United States 
##### * Corresponding author at: Department of Environmental and Global Health, College of Public Health and Health Professions, University of Florida, 2187 Mowry Rd, Gainesville, FL 32611, United States. Phone: +1-352-273-6160. 
##### E-mail address: linzhoumeng@ufl.edu (Z. Lin).


# Load required libraries
library(deSolve)   # Purpose: Solving differential equations (ODE system for corona formation kinetics)
library(ggplot2)   # Purpose: Creating high-quality data visualizations and plots
library(tidyr)     # Purpose: Data reshaping and tidying (converting between wide/long formats)
library(dplyr)     # Purpose: Data manipulation and transformation (filtering, selecting, summarizing)
library(gridExtra) # Purpose: Arranging and combining multiple ggplot objects into layouts
library(grid)      # Purpose: Low-level graphics operations and plot customization
library(stats)     # Purpose: Statistical functions and optimization routines


# Define calibration data for AgNP-HSA system from Zhao et al., 2021 (20 nm only)
calibration_data <- data.frame(
  time_h = c(0, 0.3, 0, 1.0, 0, 3.0),  # AgNP20 timepoints
  diameter = c(23.0, 28.5, 25.1, 36.9, 26.7, 40.3),  # AgNP20 diameters
  particle_size = rep("AgNP20", 6),
  timepoint = c("Initial", "Final", "Initial", "Final", "Initial", "Final"),
  experiment = c("Period1", "Period1", "Period2", "Period2", "Period3", "Period3")
)

# Base parameters for AgNP-HSA model
base_params <- list(
  k_on_base = 2.4e3,                    # Base association rate constant for HSA (M⁻¹·s⁻¹ ),  cited from Dell’Orco et al. (2010)
  k_off_base = 2.0e-3,                  # Base dissociation rate constant for HSA (s⁻¹ ), cited from Dell’Orco et al. (2010)
  c_init = 6.0e-4                         # Initial HSA concentration (M), cited from Dell’Orco et al. (2010)
)

# Define the differential equation model for corona formation
corona_ode_model <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Extract parameters
    d_initial <- initial_diameter
    
    # Calculate time-dependent and size-dependent factors
    size_effect <- (d_initial/20)^0.8
    time_effect <- (1 + 0.2 * log(time + 0.01))
    
    # Calculate adjusted rate constants
    k_on <- k_on_base * size_effect * time_effect
    k_off <- k_off_base * (d_initial/20)^(-0.5)
    
    # Surface availability
    S <- (1 - theta)^2
    
    # Cooperative binding factor
    F_coop <- 1 + 0.5 * theta * (d_initial/20)^0.3
    
    # Crowding effect factor - this is what we're optimizing
    # Using parameters a and b for the quadratic form
    F_crowd <- (1 - a * theta) * (1 + b * theta^2)
    
    # Core differential equation
    dtheta <- k_on * c_init * S * F_coop * F_crowd - k_off * theta^1.8
    
    # Return the rate of change
    return(list(c(dtheta)))
  })
}

# Function to calculate the simulated diameter given theta and parameters
calculate_diameter <- function(theta, d_initial) {
  # Protein layer thickness calculations
  L_base <- 5.8  # Base thickness calibrated from experiments
  
  # Size-dependent thickness
  L <- L_base * (d_initial/20)^0.3
  
  # Time and size dependent factors (matching full model)
  T <- 1 + 0.15 * log(1 + 0.001)  # Time factor at t=1h
  S <- 1 + 0.1 * (d_initial/20)^0.5  # Size factor
  
  # Calculate diameter
  d <- d_initial + 2 * L * theta * T * S
  
  return(d)
}

# Function to simulate the model and get diameters at specified timepoints
simulate_model <- function(params, experimental_data) {
  # Get unique experimental setups
  experiments <- unique(experimental_data$experiment)
  
  # Initialize results dataframe
  results <- data.frame()
  
  # Process each experiment separately
  for (exp in experiments) {
    # Get experiment data
    exp_data <- experimental_data[experimental_data$experiment == exp, ]
    
    # Get initial and final timepoints
    initial_data <- exp_data[exp_data$timepoint == "Initial", ]
    final_data <- exp_data[exp_data$timepoint == "Final", ]
    
    # Set parameters for simulation
    parameters <- c(
      k_on_base = base_params$k_on_base,
      k_off_base = base_params$k_off_base,
      c_init = base_params$c_init,
      a = params[1],  # Crowding factor parameter a
      b = params[2],  # Crowding factor parameter b
      initial_diameter = initial_data$diameter
    )
    
    # Initial state (no coverage)
    initial_state <- c(theta = 0.0)
    
    # Simulation timepoints
    times <- seq(0, final_data$time_h + 0.1, by = 0.01)
    
    # Solve the ODE system
    out <- ode(y = initial_state, times = times, func = corona_ode_model, parms = parameters)
    
    # Convert to data frame
    sim_data <- as.data.frame(out)
    colnames(sim_data) <- c("time", "coverage")
    
    # Calculate diameters
    sim_data$diameter <- calculate_diameter(sim_data$coverage, initial_data$diameter)
    
    # Find model value at experimental final timepoint
    idx <- which.min(abs(sim_data$time - final_data$time_h))
    sim_diameter <- sim_data$diameter[idx]
    
    # Add to results
    results <- rbind(results, data.frame(
      experiment = exp,
      time_h = final_data$time_h,
      exp_diameter = final_data$diameter,
      sim_diameter = sim_diameter,
      initial_diameter = initial_data$diameter
    ))
  }
  
  return(results)
}

# Function to calculate error metrics
calculate_error <- function(params, experimental_data) {
  # Get simulation results
  sim_results <- simulate_model(params, experimental_data)
  
  # Calculate sum of squared errors
  sse <- sum((sim_results$sim_diameter - sim_results$exp_diameter)^2)
  
  # Calculate RMSE
  rmse <- sqrt(sse / nrow(sim_results))
  
  return(rmse)
}

# Function to perform parameter sweep and visualize results (sequential version)
parameter_sweep <- function(a_range, b_range, steps, experimental_data) {
  # Create parameter grid
  a_values <- seq(a_range[1], a_range[2], length.out = steps)
  b_values <- seq(b_range[1], b_range[2], length.out = steps)
  
  # Initialize results matrix
  results <- matrix(NA, nrow = steps, ncol = steps)
  
  # Perform parameter sweep
  cat("Starting parameter sweep...\n")
  total_steps <- steps * steps
  step_count <- 0
  
  for (i in 1:steps) {
    for (j in 1:steps) {
      params <- c(a_values[i], b_values[j])
      results[i, j] <- calculate_error(params, experimental_data)
      
      # Progress update
      step_count <- step_count + 1
      if (step_count %% 10 == 0 || step_count == total_steps) {
        cat("\rProgress: ", round(step_count/total_steps*100), "% complete", sep="")
        flush.console()
      }
    }
  }
  cat("\n")
  
  # Convert to data frame for plotting
  results_df <- expand.grid(a = a_values, b = b_values)
  results_df$rmse <- as.vector(results)
  
  # Create heatmap plot
  p <- ggplot(results_df, aes(x = a, y = b, fill = rmse)) +
    geom_tile() +
    scale_fill_gradient2(low = "blue", mid = "yellow", high = "red", midpoint = median(results_df$rmse)) +
    labs(
      title = "Parameter Sweep Error Heatmap",
      subtitle = "Root Mean Square Error for different parameter combinations",
      x = "Parameter a (Crowding Linear Term)",
      y = "Parameter b (Crowding Quadratic Term)",
      fill = "RMSE"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5),
      axis.title = element_text(size = 12, face = "bold"),
      legend.position = "right",
      legend.title = element_text(size = 10, face = "bold")
    )
  
  # Find best parameters
  best_idx <- which.min(results_df$rmse)
  best_params <- c(results_df$a[best_idx], results_df$b[best_idx])
  best_rmse <- results_df$rmse[best_idx]
  
  cat("Best parameters from sweep:\n")
  cat("a =", best_params[1], "\n")
  cat("b =", best_params[2], "\n")
  cat("RMSE =", best_rmse, "\n")
  
  return(list(
    params = best_params,
    rmse = best_rmse,
    results = results_df,
    plot = p
  ))
}

# Function to perform gradient descent optimization
optimize_parameters <- function(initial_params, experimental_data, max_iterations = 100, learning_rate = 0.05, tolerance = 1e-4) {
  # Initialize parameters and error
  current_params <- initial_params
  current_error <- calculate_error(current_params, experimental_data)
  
  # Store history
  history <- data.frame(
    iteration = 0,
    a = current_params[1],
    b = current_params[2],
    rmse = current_error
  )
  
  # Run gradient descent
  for (iter in 1:max_iterations) {
    # Calculate gradient using finite differences
    h <- 0.01  # Step size for finite difference
    
    # Gradient for parameter a
    params_plus_a <- current_params
    params_plus_a[1] <- current_params[1] + h
    error_plus_a <- calculate_error(params_plus_a, experimental_data)
    
    # Gradient for parameter b
    params_plus_b <- current_params
    params_plus_b[2] <- current_params[2] + h
    error_plus_b <- calculate_error(params_plus_b, experimental_data)
    
    # Calculate gradient
    gradient <- c(
      (error_plus_a - current_error) / h,
      (error_plus_b - current_error) / h
    )
    
    # Update parameters
    new_params <- current_params - learning_rate * gradient
    
    # Ensure parameters remain in reasonable range
    new_params[1] <- max(0.1, min(0.9, new_params[1]))  # a between 0.1 and 0.9
    new_params[2] <- max(0.0, min(1.0, new_params[2]))  # b between 0 and 1
    
    # Calculate new error
    new_error <- calculate_error(new_params, experimental_data)
    
    # Update history
    history <- rbind(history, data.frame(
      iteration = iter,
      a = new_params[1],
      b = new_params[2],
      rmse = new_error
    ))
    
    # Check for convergence
    if (abs(new_error - current_error) < tolerance) {
      cat("Converged after", iter, "iterations\n")
      break
    }
    
    # Update current values
    current_params <- new_params
    current_error <- new_error
    
    # Print progress every 10 iterations
    if (iter %% 10 == 0) {
      cat("Iteration", iter, "- RMSE:", round(current_error, 6), 
          "- params:", round(current_params[1], 4), round(current_params[2], 4), "\n")
    }
  }
  
  # Plot convergence
  p <- ggplot(history, aes(x = iteration, y = rmse)) +
    geom_line() +
    geom_point() +
    labs(
      title = "Gradient Descent Optimization",
      subtitle = "RMSE vs. Iteration",
      x = "Iteration",
      y = "Root Mean Square Error"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5),
      axis.title = element_text(size = 12, face = "bold")
    )
  
  return(list(
    optimal_params = current_params,
    final_rmse = current_error,
    history = history,
    plot = p
  ))
}

# Function to validate model with optimized parameters
validate_model <- function(params, experimental_data) {
  # Get simulation results
  sim_results <- simulate_model(params, experimental_data)
  
  # Generate simulation curves for plotting
  plot_data <- data.frame()
  
  # For each experimental setup
  for (i in 1:nrow(sim_results)) {
    # Get experimental setup
    exp <- sim_results$experiment[i]
    exp_data <- experimental_data[experimental_data$experiment == exp, ]
    
    # Get initial and final timepoints
    initial_data <- exp_data[exp_data$timepoint == "Initial", ]
    final_data <- exp_data[exp_data$timepoint == "Final", ]
    
    # Set parameters for simulation
    parameters <- c(
      k_on_base = base_params$k_on_base,
      k_off_base = base_params$k_off_base,
      c_init = base_params$c_init,
      a = params[1],
      b = params[2],
      initial_diameter = initial_data$diameter
    )
    
    # Initial state (no coverage)
    initial_state <- c(theta = 0.0)
    
    # Simulation timepoints
    times <- seq(0, final_data$time_h + 0.1, by = 0.01)
    
    # Solve the ODE system
    out <- ode(y = initial_state, times = times, func = corona_ode_model, parms = parameters)
    
    # Convert to data frame
    sim_data <- as.data.frame(out)
    colnames(sim_data) <- c("time", "coverage")
    
    # Calculate diameters
    sim_data$diameter <- calculate_diameter(sim_data$coverage, initial_data$diameter)
    sim_data$experiment <- exp
    
    # Add to plot data
    plot_data <- rbind(plot_data, sim_data)
  }
  
  # Create plot
  p <- ggplot() +
    # Add simulation curves
    geom_line(data = plot_data, aes(x = time, y = diameter, color = experiment)) +
    # Add experimental points
    geom_point(data = experimental_data[experimental_data$timepoint == "Final", ], 
               aes(x = time_h, y = diameter, color = experiment), size = 3, shape = 18) +
    geom_point(data = experimental_data[experimental_data$timepoint == "Initial", ], 
               aes(x = time_h, y = diameter, color = experiment), size = 3, shape = 15) +
    labs(
      title = "AgNP-HSA Corona Formation Model Validation",
      subtitle = paste("Crowding parameters: a =", round(params[1], 3), "b =", round(params[2], 3)),
      x = "Time (h)",
      y = "Diameter (nm)",
      color = "Experiment"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5),
      axis.title = element_text(size = 12, face = "bold"),
      legend.position = "right",
      legend.title = element_text(size = 10, face = "bold")
    )
  
  # Calculate error metrics
  sim_results$error_nm <- sim_results$sim_diameter - sim_results$exp_diameter
  sim_results$error_pct <- abs(sim_results$error_nm) / sim_results$exp_diameter * 100
  
  # Print validation results
  cat("\nValidation Results:\n")
  print(sim_results)
  
  # Calculate overall RMSE
  overall_rmse <- sqrt(mean(sim_results$error_nm^2))
  cat("\nOverall RMSE:", overall_rmse, "nm\n")
  
  # Calculate average percent error
  avg_error <- mean(sim_results$error_pct)
  cat("Average percent error:", round(avg_error, 2), "%\n")
  
  return(list(
    validation_results = sim_results,
    plot = p,
    overall_rmse = overall_rmse,
    avg_error = avg_error
  ))
}

# Main script execution
cat("Starting AgNP-HSA Parameter Optimization\n")
cat("----------------------------------------\n")

# For simplicity, let's just try a small grid for a proof of concept
# STAGE 1: Broad parameter sweep
cat("\nSTAGE 1: Broad Parameter Sweep (4x4 grid for quick demonstration)\n")
sweep_results <- parameter_sweep(
  a_range = c(0.4, 0.8),   # Linear term parameter range
  b_range = c(0.1, 0.5),   # Quadratic term parameter range
  steps = 4,               # Small grid for quick demonstration
  experimental_data = calibration_data
)
print(sweep_results$plot)

# STAGE 2: Focused parameter sweep around promising region
cat("\nSTAGE 2: Focused Parameter Sweep\n")
# Use the results from stage 1 to define a narrower range
a_center <- sweep_results$params[1]
b_center <- sweep_results$params[2]
a_range <- c(max(0.3, a_center - 0.1), min(0.9, a_center + 0.1))
b_range <- c(max(0.0, b_center - 0.1), min(0.6, b_center + 0.1))

focused_sweep_results <- parameter_sweep(
  a_range = a_range,
  b_range = b_range,
  steps = 5,            # Small grid for quick demonstration
  experimental_data = calibration_data
)
print(focused_sweep_results$plot)

# STAGE 3: Gradient descent optimization
cat("\nSTAGE 3: Gradient Descent Optimization\n")
# Use the results from stage 2 as starting point
initial_params <- focused_sweep_results$params
optimization_results <- optimize_parameters(
  initial_params = initial_params,
  experimental_data = calibration_data,
  max_iterations = 20,     # Reduced for demonstration
  learning_rate = 0.02,
  tolerance = 1e-5
)
print(optimization_results$plot)

# Final validation
cat("\nFinal Model Validation\n")
validation_results <- validate_model(
  params = optimization_results$optimal_params,
  experimental_data = calibration_data
)
print(validation_results$plot)

# Display final optimized parameters
cat("\nFinal Optimized Parameters:\n")
cat("a (linear crowding term) =", round(optimization_results$optimal_params[1], 4), "\n")
cat("b (quadratic crowding term) =", round(optimization_results$optimal_params[2], 4), "\n")
cat("Final RMSE =", round(optimization_results$final_rmse, 4), "nm\n")

