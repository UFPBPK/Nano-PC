### Complete Sensitivity Analysis Framework for Nanoparticle Protein Corona Evolution Models ###
### Purpose: Comprehensive parameter sensitivity analysis for all NP-PC systems:
###          - AgNP-HSA (size and zeta potential evolution)
###          - AuNP-HSA (size and zeta potential evolution) 
###          - AuNP-BSA (size evolution)
###          - AgNP-BSA (size evolution)
### Author: Xinyue Chen a,b,c, Zhoumeng Lin a,b,c,*
##### a Department of Environmental and Global Health, College of Public Health and Health Professions, University of Florida, Gainesville, FL 32611, United States 
##### b Center for Environmental and Human Toxicology, University of Florida, Gainesville, FL 32611, United States 
##### c Center for Pharmacometrics and Systems Pharmacology, University of Florida, Orlando, FL 32827, United States 
##### * Corresponding author at: Department of Environmental and Global Health, College of Public Health and Health Professions, University of Florida, 2187 Mowry Rd, Gainesville, FL 32611, United States. Phone: +1-352-273-6160. 
##### E-mail address: linzhoumeng@ufl.edu (Z. Lin).

# Load required libraries
library(deSolve)    # Purpose: Solving differential equations for protein corona formation kinetics
# Required for running modified ODE models with varied parameters
# Provides ode() function for numerical integration of rate equations

library(dplyr)      # Purpose: Data manipulation and transformation for sensitivity analysis results
# Used for: grouping parameters, calculating sensitivity indices, filtering datasets, summarizing results across parameter variations

library(ggplot2)    # Purpose: Creating publication-quality visualizations for sensitivity results
# Used for: tornado plots (parameter ranking), response curves, parameter correlation plots, and sensitivity heatmaps

library(gridExtra)  # Purpose: Arranging multiple ggplot objects into combined layouts
# Used for: creating multi-panel figures showing sensitivity results for both size and zeta potential in single publication figure

library(reshape2)   # Purpose: Data reshaping between wide and long formats for visualization
# Used for: converting sensitivity matrices for heatmap plotting, transforming parameter sweep results for faceted plotting

library(R6)         # Purpose: Object-oriented programming framework

# =============================================================================
# 1. AgNP-HSA SENSITIVITY ANALYSIS CLASS
# =============================================================================

AgNP_HSA_SensitivityAnalyzer <- R6::R6Class(
  "AgNP_HSA_SensitivityAnalyzer",
  public = list(
    base_params = NULL,
    param_ranges = NULL,
    sensitivity_results = NULL,
    
    initialize = function() {
      # Define base parameters (from AgNP-HSA models)
      self$base_params <- list(
        # Kinetic parameters
        k_on_base = 2.4e3,
        k_off_base = 2.0e-3,
        c_init = 6.0e-4,
        
        # Physical parameters
        L_base = 5.8,               # Base layer thickness
        
        # Size evolution specific parameters (from original model)
        a = 0.1,                    # Crowding parameter
        b = 0.74,                   # Crowding parameter  
        cooperative_factor = 0.5,   # Cooperative binding factor
        size_power = 0.8,           # Size dependence power (relative_size^0.8)
        size_power_coop = 0.3,      # Size power for cooperative binding
        size_power_layer = 0.3,     # Size power for layer thickness
        k_off_size_power = -0.5,    # k_off size dependence power
        time_log_factor = 0.2,      # Time logarithmic factor
        layer_time_factor = 0.15,   # Layer time adjustment
        layer_size_factor = 0.1,    # Layer size adjustment
        theta_power = 1.8,          # Theta power in dissociation
        
        # Zeta potential model parameters
        delta_zeta_max = 149.8,     # From zeta potential model calibration
        F_charge = 0.4,             # Charge factor for negative initial zeta
        F_size_factor = 0.5,        # Size factor coefficient
        F_incubation_long = 0.25    # Factor for long incubation times
      )
      
      # Define parameter variation ranges (±50% around base values)
      self$param_ranges <- list(
        k_on_base = c(1.2e3, 3.6e3),
        k_off_base = c(1.0e-3, 3.0e-3),
        c_init = c(3.0e-4, 9.0e-4),
        L_base = c(2.9, 8.7),
        a = c(0.05, 0.15),
        b = c(0.37, 1.11),
        cooperative_factor = c(0.25, 0.75),
        size_power = c(0.4, 1.2),
        size_power_coop = c(0.15, 0.45),
        size_power_layer = c(0.15, 0.45),
        k_off_size_power = c(-0.75, -0.25),
        time_log_factor = c(0.1, 0.3),
        layer_time_factor = c(0.075, 0.225),
        layer_size_factor = c(0.05, 0.15),
        theta_power = c(0.9, 2.7),
        delta_zeta_max = c(74.9, 224.7),
        F_charge = c(0.2, 0.6),
        F_size_factor = c(0.25, 0.75),
        F_incubation_long = c(0.125, 0.375)
      )
    },
    
    # One-at-a-time sensitivity analysis (TORNADO PLOT METHOD)
    run_oat_sensitivity = function(target_variable = "size", time_point = 1.0, 
                                   initial_size = 23.0, n_points = 21) {
      
      results_list <- list()
      
      for (param_name in names(self$param_ranges)) {
        cat("Analyzing AgNP-HSA parameter:", param_name, "\n")
        
        # Get baseline output
        base_output <- if (target_variable == "size") {
          self$simulate_size_single(self$base_params, time_point, initial_size)
        } else if (target_variable == "zeta") {
          self$simulate_zeta_single(self$base_params, time_point, initial_size)
        }
        
        # Calculate outputs at +50% and -50% of parameter value
        param_range <- self$param_ranges[[param_name]]
        param_low <- param_range[1]   # -50%
        param_high <- param_range[2]  # +50%
        
        # Low value (-50%)
        modified_params_low <- self$base_params
        modified_params_low[[param_name]] <- param_low
        output_low <- if (target_variable == "size") {
          self$simulate_size_single(modified_params_low, time_point, initial_size)
        } else if (target_variable == "zeta") {
          self$simulate_zeta_single(modified_params_low, time_point, initial_size)
        }
        
        # High value (+50%)
        modified_params_high <- self$base_params
        modified_params_high[[param_name]] <- param_high
        output_high <- if (target_variable == "size") {
          self$simulate_size_single(modified_params_high, time_point, initial_size)
        } else if (target_variable == "zeta") {
          self$simulate_zeta_single(modified_params_high, time_point, initial_size)
        }
        
        # Calculate tornado plot values
        low_effect <- (output_low - base_output) / base_output * 100    # % change from baseline
        high_effect <- (output_high - base_output) / base_output * 100  # % change from baseline
        
        # Store results for tornado plot
        results_list[[param_name]] <- data.frame(
          parameter = param_name,
          baseline_output = base_output,
          low_output = output_low,
          high_output = output_high,
          low_effect = low_effect,      # Negative side of tornado
          high_effect = high_effect,    # Positive side of tornado
          total_swing = abs(high_effect - low_effect),  # Total sensitivity
          stringsAsFactors = FALSE
        )
      }
      
      # Combine all results
      sensitivity_df <- bind_rows(results_list)
      self$sensitivity_results <- sensitivity_df
      
      return(sensitivity_df)
    },
    
    # Simulate size with modified parameters
    simulate_size_single = function(params, time_point, initial_size) {
      # AgNP-HSA size evolution model (matching original model)
      model <- function(t, state, parameters) {
        with(as.list(c(state, parameters)), {
          theta <- state[1]
          time_h <- t/3600  # Convert to hours
          
          # Enhanced size-dependent effects (from original model)
          relative_size <- d_initial / 20  # normalized to SNP20
          
          # Progressive size effect over time
          time_factor_calc <- 1 + time_log_factor * log(time_h + 0.01)
          
          # Stronger size dependence for larger particles
          k_on_factor <- relative_size^size_power * time_factor_calc  # From original model
          k_off_factor <- relative_size^k_off_size_power  # From original model
          
          # Additional size-dependent multiplier for larger particles
          size_multiplier <- ifelse(d_initial > 40, 
                                    1 + 0.2 * log(d_initial/40), 
                                    1)
          
          # Adjusted rate constants with size-dependent effects
          k_on <- k_on_base * k_on_factor * size_multiplier
          k_off <- k_off_base * k_off_factor
          
          # Enhanced surface availability model
          available_surface <- (1 - theta)^2
          
          # Improved crowding effects (using optimized parameters)
          crowding_factor_calc <- (1 - a * theta) * (1 + b * theta^2)
          
          # Size-dependent cooperative binding
          cooperative_factor_calc <- 1 + cooperative_factor * theta * (d_initial/20)^size_power_coop
          
          # Final equation
          dtheta <- k_on * c_init * available_surface * 
            crowding_factor_calc * cooperative_factor_calc - 
            k_off * theta^theta_power  # Adjustable exponent
          
          return(list(c(dtheta)))
        })
      }
      
      # Solve ODE with increased tolerances to reduce warnings
      times <- seq(0, time_point * 3600, by = 60)  # Convert to seconds, minute steps
      out <- ode(y = c(0), times = times, func = model, 
                 parms = c(params, d_initial = initial_size),
                 rtol = 1e-4, atol = 1e-6)  # Increased tolerances
      
      # Calculate final size with size-dependent layer thickness
      final_coverage <- out[nrow(out), 2]
      size_factor <- (initial_size / 20)^params$size_power_layer  # From original model
      base_layer <- params$L_base * size_factor
      
      # Enhanced diameter calculation (from original model)
      layer_thickness <- base_layer * 
        (1 + params$layer_time_factor * log(time_point + 0.001)) * 
        (1 + params$layer_size_factor * (initial_size/20)^0.5)
      
      final_size <- initial_size + 2 * final_coverage * layer_thickness
      
      return(final_size)
    },
    
    # Simulate zeta potential with modified parameters
    simulate_zeta_single = function(params, time_point, initial_size, initial_zeta = -38.8) {
      # Get coverage from size model
      coverage <- (self$simulate_size_single(params, time_point, initial_size) - initial_size) / 
        (2 * params$L_base * (initial_size / 20)^params$size_power_layer * 
           (1 + params$layer_time_factor * log(time_point + 0.001)) * 
           (1 + params$layer_size_factor * (initial_size/20)^0.5))
      
      # Calculate zeta potential change using saturation formula
      F_size <- 1 - params$F_size_factor * tanh(initial_size/50)
      F_incubation <- if (time_point > 4) params$F_incubation_long else 1.0
      zeta_change <- params$delta_zeta_max * params$F_charge * F_size * F_incubation * (coverage / (0.2 + coverage))
      final_zeta <- initial_zeta + zeta_change
      
      return(final_zeta)
    },
    
    # Calculate sensitivity indices for tornado plot
    calculate_sensitivity_indices = function() {
      if (is.null(self$sensitivity_results)) {
        stop("Run sensitivity analysis first")
      }
      
      # Calculate tornado plot metrics
      indices <- self$sensitivity_results %>%
        mutate(
          max_effect = pmax(abs(low_effect), abs(high_effect)),
          total_range = abs(high_effect - low_effect)
        ) %>%
        arrange(desc(total_range)) %>%
        mutate(importance_rank = row_number())
      
      return(indices)
    },
    
    # Create proper tornado plot
    plot_tornado = function(target_variable = "Size") {
      indices <- self$calculate_sensitivity_indices()
      
      # Prepare data for tornado plot
      tornado_data <- indices %>%
        select(parameter, low_effect, high_effect, total_range) %>%
        arrange(total_range) %>%  # Order by sensitivity (smallest to largest)
        mutate(
          parameter = factor(parameter, levels = parameter),  # Preserve order
          low_effect = ifelse(low_effect > 0, 0, low_effect),  # Ensure low effects are ≤ 0
          high_effect = ifelse(high_effect < 0, 0, high_effect)  # Ensure high effects are ≥ 0
        )
      
      # Create tornado plot
      p <- ggplot(tornado_data) +
        # Left side (negative effects)
        geom_col(aes(x = parameter, y = low_effect), 
                 fill = "coral", width = 0.7, alpha = 0.8) +
        # Right side (positive effects)  
        geom_col(aes(x = parameter, y = high_effect), 
                 fill = "steelblue", width = 0.7, alpha = 0.8) +
        coord_flip() +
        geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
        labs(
          title = paste("AgNP-HSA Parameter Sensitivity Analysis -", target_variable),
          subtitle = "Tornado Plot: Impact of ±50% parameter variation on model output",
          x = "Model Parameters",
          y = "% Change in Output from Baseline",
          caption = "Orange: When parameter decreases by 50%  | Blue: When parameter increases by 50%"
        ) +
        theme_minimal() +
        theme(
          plot.title = element_text(size = 14, face = "bold"),
          plot.subtitle = element_text(size = 11),
          axis.title = element_text(size = 12),
          plot.caption = element_text(size = 9, color = "gray50")
        ) +
        scale_y_continuous(labels = function(x) paste0(x, "%"))
      
      return(p)
    },
    
    # Create parameter response curves
    plot_response_curves = function(top_n = 6, target_variable = "size") {
      if (is.null(self$sensitivity_results)) {
        stop("Run sensitivity analysis first")
      }
      
      # Get top N most sensitive parameters based on total range
      indices <- self$calculate_sensitivity_indices()
      top_params <- indices$parameter[1:min(top_n, nrow(indices))]
      
      # For response curves, we need to generate detailed parameter sweeps
      response_data <- list()
      
      for (param_name in top_params) {
        # Create detailed parameter sweep
        param_range <- self$param_ranges[[param_name]]
        param_values <- seq(param_range[1], param_range[2], length.out = 21)
        
        outputs <- numeric(length(param_values))
        base_output <- indices$baseline_output[indices$parameter == param_name]
        
        for (i in 1:length(param_values)) {
          modified_params <- self$base_params
          modified_params[[param_name]] <- param_values[i]
          
          outputs[i] <- self$simulate_size_single(modified_params, 1.0, 23.0)
        }
        
        response_data[[param_name]] <- data.frame(
          parameter = param_name,
          param_normalized = (param_values - self$base_params[[param_name]]) / self$base_params[[param_name]] * 100,
          output_normalized = (outputs - base_output) / base_output * 100
        )
      }
      
      # Combine and plot
      plot_data <- bind_rows(response_data) %>%
        mutate(parameter = factor(parameter, levels = top_params))
      target_label <- if(target_variable == "size") "Size Evolution" else "Zeta Potential Evolution"
      y_label <- if(target_variable == "size") "Size Change (%)" else "Zeta Potential Change (%)"
      
      p <- ggplot(plot_data, aes(x = param_normalized, y = output_normalized)) +
        geom_line(color = "steelblue", linewidth = 1) +
        geom_point(color = "darkblue", size = 1.5) +
        geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
        geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
        facet_wrap(~parameter, scales = "free", ncol = 3) +
        labs(
          title = paste("AgNP-HSA Parameter Response Curves -", target_label),
          subtitle = paste("Impact of parameter changes on", tolower(target_label)),
          x = "Parameter Change (%)",
          y = y_label
        ) +
        theme_minimal() +
        theme(
          plot.title = element_text(size = 14, face = "bold"),
          strip.text = element_text(size = 10, face = "bold")
        )
      
      return(p)
    },
    
    # Generate comprehensive sensitivity report
    generate_report = function(target_variable = "size") {
      cat("\n=== AgNP-HSA SENSITIVITY ANALYSIS REPORT ===\n")
      cat("Target Variable:", toupper(target_variable), "\n")
      cat("Analysis Date:", Sys.Date(), "\n\n")
      
      indices <- self$calculate_sensitivity_indices()
      
      cat("TOP 5 MOST INFLUENTIAL PARAMETERS:\n")
      cat("===================================\n")
      for (i in 1:min(5, nrow(indices))) {
        param <- indices[i, ]
        cat(sprintf("%d. %s\n", i, param$parameter))
        cat(sprintf("   Total Sensitivity Range: %.1f%%\n", param$total_range))
        cat(sprintf("   Low Effect (-50%%): %.1f%%\n", param$low_effect))
        cat(sprintf("   High Effect (+50%%): %.1f%%\n", param$high_effect))
        cat("\n")
      }
      
      cat("PARAMETER CATEGORIES:\n")
      cat("====================\n")
      kinetic_params <- c("k_on_base", "k_off_base", "c_init")
      physical_params <- c("L_base", "F_size_factor")
      size_model_params <- c("a", "b", "cooperative_factor", "size_power", "size_power_coop", "size_power_layer", 
                             "k_off_size_power", "time_log_factor", "layer_time_factor", "layer_size_factor", "theta_power")
      zeta_params <- c("delta_zeta_max", "F_charge", "F_incubation_long")
      
      for (category in list(
        list(name = "Kinetic Parameters", params = kinetic_params),
        list(name = "Physical Parameters", params = physical_params),
        list(name = "Size Model Parameters", params = size_model_params),
        list(name = "Zeta Potential Parameters", params = zeta_params)
      )) {
        cat(sprintf("%s:\n", category$name))
        cat_params <- indices[indices$parameter %in% category$params, ]
        if (nrow(cat_params) > 0) {
          for (j in 1:nrow(cat_params)) {
            cat(sprintf("  - %s (Range: %.1f%%, Rank: %d)\n", 
                        cat_params$parameter[j], 
                        cat_params$total_range[j],
                        cat_params$importance_rank[j]))
          }
        }
        cat("\n")
      }
      
      return(indices)
    }
  )
)

# =============================================================================
# 2. AuNP-HSA SENSITIVITY ANALYSIS CLASS
# =============================================================================

AuNP_HSA_SensitivityAnalyzer <- R6::R6Class(
  "AuNP_HSA_SensitivityAnalyzer",
  public = list(
    base_params = NULL,
    param_ranges = NULL,
    sensitivity_results = NULL,
    
    initialize = function() {
      # Define base parameters (from AuNP-HSA models)
      self$base_params <- list(
        # Kinetic parameters
        k_on_base = 2.4e3,
        k_off_base = 2.0e-3,
        c_init = 6.0e-4,
        
        # Physical parameters
        L_base = 5.8,           # From calibration
        
        # Size evolution model parameters (from original model)
        cooperative_exp = 0.5,      # Exponent in (d_initial/10)^0.5
        surface_exp = 1.5,          # Surface availability exponent
        crowding_factor = 0.2,      # Crowding effect parameter
        size_effect_exp = 0.3,      # Size effect exponent (d_initial/10)^0.3
        time_factor = 0.8,          # Time-dependent factor
        time_offset = 0.2,          # Time offset parameter
        layer_adjust_factor = 0.2,  # Layer thickness adjustment (log factor)
        k_off_adjust = 0.2,         # k_off size adjustment factor
        
        # Zeta potential model parameters
        delta_zeta_max = 39.14,     # From calibration
        F_charge_pos = -0.8,        # For positive initial zeta
        F_charge_neg = 0.6,         # For negative initial zeta
        F_size_factor = 0.3,        # Size factor coefficient
        F_size_scale = 100,         # Size scaling parameter
        F_coat_CTAB = 1.2,          # CTAB coating factor
        F_coat_citrate = 1.0        # Citrate coating factor
      )
      
      # Define parameter variation ranges (±50% around base values)
      self$param_ranges <- list(
        k_on_base = c(1.2e3, 3.6e3),
        k_off_base = c(1.0e-3, 3.0e-3),
        c_init = c(3.0e-4, 9.0e-4),
        L_base = c(2.9, 8.7),
        cooperative_exp = c(0.25, 0.75),
        surface_exp = c(0.75, 2.25),
        crowding_factor = c(0.1, 0.3),
        size_effect_exp = c(0.15, 0.45),
        time_factor = c(0.4, 1.2),
        time_offset = c(0.1, 0.3),
        layer_adjust_factor = c(0.1, 0.3),
        k_off_adjust = c(0.1, 0.3),
        delta_zeta_max = c(19.57, 58.71),
        F_charge_pos = c(-1.2, -0.4),
        F_charge_neg = c(0.3, 0.9),
        F_size_factor = c(0.15, 0.45),
        F_size_scale = c(50, 150),
        F_coat_CTAB = c(0.6, 1.8),
        F_coat_citrate = c(0.5, 1.5)
      )
    },
    
    # One-at-a-time sensitivity analysis (TORNADO PLOT METHOD)
    run_oat_sensitivity = function(target_variable = "size", time_point = 1.0, 
                                   initial_size = 10.0, initial_zeta = -48.93, n_points = 21) {
      
      results_list <- list()
      
      for (param_name in names(self$param_ranges)) {
        cat("Analyzing AuNP-HSA parameter:", param_name, "\n")
        
        # Get baseline output
        base_output <- if (target_variable == "size") {
          self$simulate_size_single(self$base_params, time_point, initial_size)
        } else if (target_variable == "zeta") {
          self$simulate_zeta_single(self$base_params, time_point, initial_size, initial_zeta)
        }
        
        # Calculate outputs at +50% and -50% of parameter value
        param_range <- self$param_ranges[[param_name]]
        param_low <- param_range[1]   # -50%
        param_high <- param_range[2]  # +50%
        
        # Low value (-50%)
        modified_params_low <- self$base_params
        modified_params_low[[param_name]] <- param_low
        output_low <- if (target_variable == "size") {
          self$simulate_size_single(modified_params_low, time_point, initial_size)
        } else if (target_variable == "zeta") {
          self$simulate_zeta_single(modified_params_low, time_point, initial_size, initial_zeta)
        }
        
        # High value (+50%)
        modified_params_high <- self$base_params
        modified_params_high[[param_name]] <- param_high
        output_high <- if (target_variable == "size") {
          self$simulate_size_single(modified_params_high, time_point, initial_size)
        } else if (target_variable == "zeta") {
          self$simulate_zeta_single(modified_params_high, time_point, initial_size, initial_zeta)
        }
        
        # Calculate tornado plot values
        low_effect <- (output_low - base_output) / base_output * 100    # % change from baseline
        high_effect <- (output_high - base_output) / base_output * 100  # % change from baseline
        
        # Store results for tornado plot
        results_list[[param_name]] <- data.frame(
          parameter = param_name,
          baseline_output = base_output,
          low_output = output_low,
          high_output = output_high,
          low_effect = low_effect,      # Negative side of tornado
          high_effect = high_effect,    # Positive side of tornado
          total_swing = abs(high_effect - low_effect),  # Total sensitivity
          stringsAsFactors = FALSE
        )
      }
      
      # Combine all results
      sensitivity_df <- bind_rows(results_list)
      self$sensitivity_results <- sensitivity_df
      
      return(sensitivity_df)
    },
    
    # Simulate size with modified parameters
    simulate_size_single = function(params, time_point, initial_size) {
      # AuNP-HSA size evolution model (matching original model)
      model <- function(t, state, parameters) {
        with(as.list(c(state, parameters)), {
          theta <- state[1]
          
          # Cooperative binding effect (from original model)
          F_coop <- 1 + theta * (d_initial/10)^cooperative_exp
          
          # Surface availability (from original model)
          S <- (1 - theta)^surface_exp
          
          # Crowding effect (from original model)
          F_crowd <- exp(-crowding_factor * theta)
          
          # Size effect on k_on (from original model)
          size_effect <- (d_initial/10)^size_effect_exp
          
          # Time effect on k_on (from original model)
          time_effect <- (1 + time_factor*t)/(time_offset + t)
          
          # Adjusted binding rate
          k_on <- k_on_base * size_effect * time_effect
          
          # Adjusted dissociation rate (from original model)
          k_off <- k_off_base / (1 + k_off_adjust * log(d_initial/10 + 1))
          
          # Final equation from original model
          dtheta <- k_on * c_init * S * F_coop * F_crowd - k_off * theta
          
          return(list(c(dtheta)))
        })
      }
      
      # Solve ODE with increased tolerances to reduce warnings
      times <- seq(0, time_point, by = 0.01)
      out <- ode(y = c(0), times = times, func = model, 
                 parms = c(params, d_initial = initial_size),
                 rtol = 1e-4, atol = 1e-6)  # Increased tolerances
      
      # Calculate final size with size-adjusted layer thickness (from original model)
      final_coverage <- out[nrow(out), 2]
      # Adjusted protein layer thickness based on nanoparticle size
      L <- params$L_base * (1 + params$layer_adjust_factor * log(initial_size/10 + 1))
      final_size <- initial_size + 2 * L * final_coverage
      
      return(final_size)
    },
    
    # Simulate zeta potential with modified parameters
    simulate_zeta_single = function(params, time_point, initial_size, initial_zeta) {
      # Get coverage from size model
      coverage <- (self$simulate_size_single(params, time_point, initial_size) - initial_size) / 
        (2 * params$L_base * (1 + params$layer_adjust_factor * log(initial_size/10 + 1)))
      
      # Calculate modifying factors
      F_charge <- if (initial_zeta > 0) params$F_charge_pos else params$F_charge_neg
      F_size <- 1 - params$F_size_factor * tanh(initial_size/params$F_size_scale)
      F_coat <- if (initial_zeta > 0) params$F_coat_CTAB else params$F_coat_citrate
      
      # Calculate zeta potential change
      zeta_change <- params$delta_zeta_max * F_charge * F_size * F_coat * coverage
      final_zeta <- initial_zeta + zeta_change
      
      return(final_zeta)
    },
    
    # Calculate sensitivity indices for tornado plot
    calculate_sensitivity_indices = function() {
      if (is.null(self$sensitivity_results)) {
        stop("Run sensitivity analysis first")
      }
      
      # Calculate tornado plot metrics
      indices <- self$sensitivity_results %>%
        mutate(
          max_effect = pmax(abs(low_effect), abs(high_effect)),
          total_range = abs(high_effect - low_effect)
        ) %>%
        arrange(desc(total_range)) %>%
        mutate(importance_rank = row_number())
      
      return(indices)
    },
    
    # Create proper tornado plot
    plot_tornado = function(target_variable = "Size") {
      indices <- self$calculate_sensitivity_indices()
      
      # Prepare data for tornado plot
      tornado_data <- indices %>%
        select(parameter, low_effect, high_effect, total_range) %>%
        arrange(total_range) %>%  # Order by sensitivity (smallest to largest)
        mutate(
          parameter = factor(parameter, levels = parameter),  # Preserve order
          low_effect = ifelse(low_effect > 0, 0, low_effect),  # Ensure low effects are ≤ 0
          high_effect = ifelse(high_effect < 0, 0, high_effect)  # Ensure high effects are ≥ 0
        )
      
      # Create tornado plot
      p <- ggplot(tornado_data) +
        # Left side (negative effects)
        geom_col(aes(x = parameter, y = low_effect), 
                 fill = "coral", width = 0.7, alpha = 0.8) +
        # Right side (positive effects)  
        geom_col(aes(x = parameter, y = high_effect), 
                 fill = "steelblue", width = 0.7, alpha = 0.8) +
        coord_flip() +
        geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
        labs(
          title = paste("AuNP-HSA Parameter Sensitivity Analysis -", target_variable),
          subtitle = "Tornado Plot: Impact of ±50% parameter variation on model output",
          x = "Model Parameters",
          y = "% Change in Output from Baseline",
          caption = "Orange: When parameter decreases by 50%  | Blue: When parameter increases by 50%"
        ) +
        theme_minimal() +
        theme(
          plot.title = element_text(size = 14, face = "bold"),
          plot.subtitle = element_text(size = 11),
          axis.title = element_text(size = 12),
          plot.caption = element_text(size = 9, color = "gray50")
        ) +
        scale_y_continuous(labels = function(x) paste0(x, "%"))
      
      return(p)
    },
    
    # Create parameter response curves
    plot_response_curves = function(top_n = 6, target_variable = "size") {
      if (is.null(self$sensitivity_results)) {
        stop("Run sensitivity analysis first")
      }
      
      # Get top N most sensitive parameters based on total range
      indices <- self$calculate_sensitivity_indices()
      top_params <- indices$parameter[1:min(top_n, nrow(indices))]
      
      # For response curves, we need to generate detailed parameter sweeps
      response_data <- list()
      
      for (param_name in top_params) {
        # Create detailed parameter sweep
        param_range <- self$param_ranges[[param_name]]
        param_values <- seq(param_range[1], param_range[2], length.out = 21)
        
        outputs <- numeric(length(param_values))
        base_output <- indices$baseline_output[indices$parameter == param_name]
        
        for (i in 1:length(param_values)) {
          modified_params <- self$base_params
          modified_params[[param_name]] <- param_values[i]
          
          outputs[i] <- self$simulate_size_single(modified_params, 1.0, 10.0)
        }
        
        response_data[[param_name]] <- data.frame(
          parameter = param_name,
          param_normalized = (param_values - self$base_params[[param_name]]) / self$base_params[[param_name]] * 100,
          output_normalized = (outputs - base_output) / base_output * 100
        )
      }
      
      # Combine and plot
      plot_data <- bind_rows(response_data) %>%
        mutate(parameter = factor(parameter, levels = top_params))
      
      target_label <- if(target_variable == "size") "Size Evolution" else "Zeta Potential Evolution"
      y_label <- if(target_variable == "size") "Size Change (%)" else "Zeta Potential Change (%)"
      
      
      p <- ggplot(plot_data, aes(x = param_normalized, y = output_normalized)) +
        geom_line(color = "steelblue", linewidth = 1) +
        geom_point(color = "darkblue", size = 1.5) +
        geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
        geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
        facet_wrap(~parameter, scales = "free", ncol = 3) +
        labs(
          title = paste("AuNP-HSA Parameter Response Curves -", target_label),
          subtitle = paste("Impact of parameter changes on", tolower(target_label)),
          x = "Parameter Change (%)",
          y = y_label
        ) +
        theme_minimal() +
        theme(
          plot.title = element_text(size = 14, face = "bold"),
          strip.text = element_text(size = 10, face = "bold")
        )
      
      return(p)
    },
    
    # Generate comprehensive sensitivity report
    generate_report = function(target_variable = "size") {
      cat("\n=== AuNP-HSA SENSITIVITY ANALYSIS REPORT ===\n")
      cat("Target Variable:", toupper(target_variable), "\n")
      cat("Analysis Date:", Sys.Date(), "\n\n")
      
      indices <- self$calculate_sensitivity_indices()
      
      cat("TOP 5 MOST INFLUENTIAL PARAMETERS:\n")
      cat("===================================\n")
      for (i in 1:min(5, nrow(indices))) {
        param <- indices[i, ]
        cat(sprintf("%d. %s\n", i, param$parameter))
        cat(sprintf("   Total Sensitivity Range: %.1f%%\n", param$total_range))
        cat(sprintf("   Low Effect (-50%%): %.1f%%\n", param$low_effect))
        cat(sprintf("   High Effect (+50%%): %.1f%%\n", param$high_effect))
        cat("\n")
      }
      
      cat("PARAMETER CATEGORIES:\n")
      cat("====================\n")
      kinetic_params <- c("k_on_base", "k_off_base", "c_init")
      physical_params <- c("L_base", "layer_adjust_factor", "k_off_adjust")
      size_model_params <- c("cooperative_exp", "surface_exp", "crowding_factor", "size_effect_exp", "time_factor", "time_offset")
      zeta_model_params <- c("delta_zeta_max", "F_charge_pos", "F_charge_neg", "F_size_factor", "F_size_scale", "F_coat_CTAB", "F_coat_citrate")
      
      for (category in list(
        list(name = "Kinetic Parameters", params = kinetic_params),
        list(name = "Physical Parameters", params = physical_params),
        list(name = "Size Model Parameters", params = size_model_params),
        list(name = "Zeta Potential Parameters", params = zeta_model_params)
      )) {
        cat(sprintf("%s:\n", category$name))
        cat_params <- indices[indices$parameter %in% category$params, ]
        if (nrow(cat_params) > 0) {
          for (j in 1:nrow(cat_params)) {
            cat(sprintf("  - %s (Range: %.1f%%, Rank: %d)\n", 
                        cat_params$parameter[j], 
                        cat_params$total_range[j],
                        cat_params$importance_rank[j]))
          }
        }
        cat("\n")
      }
      
      return(indices)
    }
  )
)

# =============================================================================
# 3. AuNP-BSA SENSITIVITY ANALYSIS CLASS
# =============================================================================

AuNP_BSA_SensitivityAnalyzer <- R6::R6Class(
  "AuNP_BSA_SensitivityAnalyzer",
  public = list(
    base_params = NULL,
    param_ranges = NULL,
    sensitivity_results = NULL,
    
    initialize = function() {
      # Define base parameters (from AuNP-BSA size evolution model)
      self$base_params <- list(
        # Kinetic parameters
        k_on_base = 2.0e3,
        k_off_base = 1.5e-3,
        c_init = 6.0e-4,
        
        # Physical parameters
        L_base = 3.3,               # From calibration
        
        # Model parameters
        cooperative_factor = 0.5,   # From F_coop calculation (d_initial/10)^0.5
        surface_exp = 1.5,          # Surface availability exponent
        crowding_factor = 0.2,      # Crowding effect parameter
        size_effect_exp = 0.3,      # Size effect exponent (d_initial/10)^0.3
        time_factor = 0.8,          # Time-dependent factor
        time_offset = 0.2,          # Time offset parameter
        layer_adjust_factor = 0.2,  # Layer thickness adjustment
        k_off_adjust = 0.2          # k_off size adjustment factor
      )
      
      # Define parameter variation ranges (±50% around base values)
      self$param_ranges <- list(
        k_on_base = c(1.0e3, 3.0e3),
        k_off_base = c(0.75e-3, 2.25e-3),
        c_init = c(3.0e-4, 9.0e-4),
        L_base = c(1.65, 4.95),
        cooperative_factor = c(0.25, 0.75),
        surface_exp = c(0.75, 2.25),
        crowding_factor = c(0.1, 0.3),
        size_effect_exp = c(0.15, 0.45),
        time_factor = c(0.4, 1.2),
        time_offset = c(0.1, 0.3),
        layer_adjust_factor = c(0.1, 0.3),
        k_off_adjust = c(0.1, 0.3)
      )
    },
    
    # One-at-a-time sensitivity analysis (TORNADO PLOT METHOD)
    run_oat_sensitivity = function(target_variable = "size", time_point = 1.5, 
                                   initial_size = 9.4, n_points = 21) {
      
      results_list <- list()
      
      for (param_name in names(self$param_ranges)) {
        cat("Analyzing AuNP-BSA parameter:", param_name, "\n")
        
        # Get baseline output
        base_output <- self$simulate_size_single(self$base_params, time_point, initial_size)
        
        # Calculate outputs at +50% and -50% of parameter value
        param_range <- self$param_ranges[[param_name]]
        param_low <- param_range[1]   # -50%
        param_high <- param_range[2]  # +50%
        
        # Low value (-50%)
        modified_params_low <- self$base_params
        modified_params_low[[param_name]] <- param_low
        output_low <- self$simulate_size_single(modified_params_low, time_point, initial_size)
        
        # High value (+50%)
        modified_params_high <- self$base_params
        modified_params_high[[param_name]] <- param_high
        output_high <- self$simulate_size_single(modified_params_high, time_point, initial_size)
        
        # Calculate tornado plot values
        low_effect <- (output_low - base_output) / base_output * 100    # % change from baseline
        high_effect <- (output_high - base_output) / base_output * 100  # % change from baseline
        
        # Store results for tornado plot
        results_list[[param_name]] <- data.frame(
          parameter = param_name,
          baseline_output = base_output,
          low_output = output_low,
          high_output = output_high,
          low_effect = low_effect,      # Negative side of tornado
          high_effect = high_effect,    # Positive side of tornado
          total_swing = abs(high_effect - low_effect),  # Total sensitivity
          stringsAsFactors = FALSE
        )
      }
      
      # Combine all results
      sensitivity_df <- bind_rows(results_list)
      self$sensitivity_results <- sensitivity_df
      
      return(sensitivity_df)
    },
    
    # Simulate size with modified parameters
    simulate_size_single = function(params, time_point, initial_size) {
      # AuNP-BSA size evolution model
      model <- function(t, state, parameters) {
        with(as.list(c(state, parameters)), {
          theta <- state[1]
          
          # Core model with modifying factors
          # Surface Availability Term (S)
          S <- (1 - theta)^surface_exp
          
          # Cooperative Binding Term (F_coop)
          F_coop <- 1 + theta * (d_initial/10)^cooperative_factor
          
          # Crowding Effect Term (F_crowd)
          F_crowd <- exp(-crowding_factor * theta)
          
          # Time-dependent k_on
          k_on <- k_on_base * size_effect * time_effect
          
          # Size-dependent k_off
          k_off <- k_off_base / (1 + k_off_adjust * log(d_initial/10 + 1))
          
          # Core equation
          dtheta <- k_on * c_init * S * F_coop * F_crowd - k_off * theta
          
          return(list(c(dtheta)))
        })
      }
      
      # Calculate size and time effects
      size_effect <- (initial_size/10)^params$size_effect_exp
      time_effect <- (1 + params$time_factor*time_point)/(params$time_offset + time_point)
      
      # Solve ODE with increased tolerances to reduce warnings
      times <- seq(0, time_point, by = 0.01)
      out <- ode(y = c(0), times = times, func = model, 
                 parms = c(params, d_initial = initial_size, size_effect = size_effect, time_effect = time_effect),
                 rtol = 1e-4, atol = 1e-6)  # Increased tolerances
      
      # Calculate final size
      final_coverage <- out[nrow(out), 2]
      # Calculate adjusted protein layer thickness
      L <- params$L_base * (1 + params$layer_adjust_factor * log(initial_size/10 + 1))
      final_size <- initial_size + 2 * L * final_coverage
      
      return(final_size)
    },
    
    # Calculate sensitivity indices for tornado plot
    calculate_sensitivity_indices = function() {
      if (is.null(self$sensitivity_results)) {
        stop("Run sensitivity analysis first")
      }
      
      # Calculate tornado plot metrics
      indices <- self$sensitivity_results %>%
        mutate(
          max_effect = pmax(abs(low_effect), abs(high_effect)),
          total_range = abs(high_effect - low_effect)
        ) %>%
        arrange(desc(total_range)) %>%
        mutate(importance_rank = row_number())
      
      return(indices)
    },
    
    # Create proper tornado plot
    plot_tornado = function(target_variable = "Size") {
      indices <- self$calculate_sensitivity_indices()
      
      # Prepare data for tornado plot
      tornado_data <- indices %>%
        select(parameter, low_effect, high_effect, total_range) %>%
        arrange(total_range) %>%  # Order by sensitivity (smallest to largest)
        mutate(
          parameter = factor(parameter, levels = parameter),  # Preserve order
          low_effect = ifelse(low_effect > 0, 0, low_effect),  # Ensure low effects are ≤ 0
          high_effect = ifelse(high_effect < 0, 0, high_effect)  # Ensure high effects are ≥ 0
        )
      
      # Create tornado plot
      p <- ggplot(tornado_data) +
        # Left side (negative effects)
        geom_col(aes(x = parameter, y = low_effect), 
                 fill = "coral", width = 0.7, alpha = 0.8) +
        # Right side (positive effects)  
        geom_col(aes(x = parameter, y = high_effect), 
                 fill = "steelblue", width = 0.7, alpha = 0.8) +
        coord_flip() +
        geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
        labs(
          title = paste("AuNP-BSA Parameter Sensitivity Analysis -", target_variable),
          subtitle = "Tornado Plot: Impact of ±50% parameter variation on model output",
          x = "Model Parameters",
          y = "% Change in Output from Baseline",
          caption = "Orange: When parameter decreases by 50%  | Blue: When parameter increases by 50%"
        ) +
        theme_minimal() +
        theme(
          plot.title = element_text(size = 14, face = "bold"),
          plot.subtitle = element_text(size = 11),
          axis.title = element_text(size = 12),
          plot.caption = element_text(size = 9, color = "gray50")
        ) +
        scale_y_continuous(labels = function(x) paste0(x, "%"))
      
      return(p)
    },
    
    # Create parameter response curves
    plot_response_curves = function(top_n = 6, target_variable = "size") {
      if (is.null(self$sensitivity_results)) {
        stop("Run sensitivity analysis first")
      }
      
      # Get top N most sensitive parameters based on total range
      indices <- self$calculate_sensitivity_indices()
      top_params <- indices$parameter[1:min(top_n, nrow(indices))]
      
      # For response curves, we need to generate detailed parameter sweeps
      response_data <- list()
      
      for (param_name in top_params) {
        # Create detailed parameter sweep
        param_range <- self$param_ranges[[param_name]]
        param_values <- seq(param_range[1], param_range[2], length.out = 21)
        
        outputs <- numeric(length(param_values))
        base_output <- indices$baseline_output[indices$parameter == param_name]
        
        for (i in 1:length(param_values)) {
          modified_params <- self$base_params
          modified_params[[param_name]] <- param_values[i]
          
          outputs[i] <- self$simulate_size_single(modified_params, 1.5, 9.4)
        }
        
        response_data[[param_name]] <- data.frame(
          parameter = param_name,
          param_normalized = (param_values - self$base_params[[param_name]]) / self$base_params[[param_name]] * 100,
          output_normalized = (outputs - base_output) / base_output * 100
        )
      }
      
      # Combine and plot
      plot_data <- bind_rows(response_data) %>%
        mutate(parameter = factor(parameter, levels = top_params))
      
      target_label <- "Size Evolution"
      y_label <- "Size Change (%)"
      
      p <- ggplot(plot_data, aes(x = param_normalized, y = output_normalized)) +
        geom_line(color = "steelblue", linewidth = 1) +
        geom_point(color = "darkblue", size = 1.5) +
        geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
        geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
        facet_wrap(~parameter, scales = "free", ncol = 3) +
        labs(
          title = paste("AuNP-BSA Parameter Response Curves -", target_label),
          subtitle = "Impact of parameter changes on size evolution",
          x = "Parameter Change (%)",
          y = y_label
        ) +
        theme_minimal() +
        theme(
          plot.title = element_text(size = 14, face = "bold"),
          strip.text = element_text(size = 10, face = "bold")
        )
      
      return(p)
    },
    
    # Generate comprehensive sensitivity report
    generate_report = function(target_variable = "size") {
      cat("\n=== AuNP-BSA SENSITIVITY ANALYSIS REPORT ===\n")
      cat("Target Variable:", toupper(target_variable), "\n")
      cat("Analysis Date:", Sys.Date(), "\n\n")
      
      indices <- self$calculate_sensitivity_indices()
      
      cat("TOP 5 MOST INFLUENTIAL PARAMETERS:\n")
      cat("===================================\n")
      for (i in 1:min(5, nrow(indices))) {
        param <- indices[i, ]
        cat(sprintf("%d. %s\n", i, param$parameter))
        cat(sprintf("   Total Sensitivity Range: %.1f%%\n", param$total_range))
        cat(sprintf("   Low Effect (-50%%): %.1f%%\n", param$low_effect))
        cat(sprintf("   High Effect (+50%%): %.1f%%\n", param$high_effect))
        cat("\n")
      }
      
      cat("PARAMETER CATEGORIES:\n")
      cat("====================\n")
      kinetic_params <- c("k_on_base", "k_off_base", "c_init")
      physical_params <- c("L_base", "layer_adjust_factor", "k_off_adjust")
      model_params <- c("cooperative_factor", "surface_exp", "crowding_factor", "size_effect_exp", "time_factor", "time_offset")
      
      for (category in list(
        list(name = "Kinetic Parameters", params = kinetic_params),
        list(name = "Physical Parameters", params = physical_params),
        list(name = "Model Parameters", params = model_params)
      )) {
        cat(sprintf("%s:\n", category$name))
        cat_params <- indices[indices$parameter %in% category$params, ]
        if (nrow(cat_params) > 0) {
          for (j in 1:nrow(cat_params)) {
            cat(sprintf("  - %s (Range: %.1f%%, Rank: %d)\n", 
                        cat_params$parameter[j], 
                        cat_params$total_range[j],
                        cat_params$importance_rank[j]))
          }
        }
        cat("\n")
      }
      
      return(indices)
    }
  )
)

# =============================================================================
# 4. AgNP-BSA SENSITIVITY ANALYSIS CLASS
# =============================================================================

AgNP_BSA_SensitivityAnalyzer <- R6::R6Class(
  "AgNP_BSA_SensitivityAnalyzer",
  public = list(
    base_params = NULL,
    param_ranges = NULL,
    sensitivity_results = NULL,
    
    initialize = function() {
      # Define base parameters (from AgNP-BSA size evolution model)
      self$base_params <- list(
        # Kinetic parameters
        k_on_base = 2.0e3,
        k_off_base = 1.5e-3,
        c_init = 6.0e-4,
        
        # Physical parameters
        L_base = 28.9,              # From calibration
        
        # Surface coating factors
        coating_factor_citrate = 0.922, #F_coat
        coating_factor_PVP = 1.078,   #F_coat
        
        # Size scaling factors
        size_scaling_small = 1.0,   # Base scaling for particles < 50 nm
        size_scaling_medium = 0.85, # Reduced scaling for 50-100 nm
        size_scaling_large = 0.5,   # Further reduced scaling for > 100 nm
        
        # Model parameters
        cooperative_factor = 0.4,   # Exponent for (d_initial/40)^0.4
        surface_exp = 1.5,          # Surface availability exponent
        crowding_factor = 0.25,     # Crowding effect parameter
        size_effect_exp_small = 0.25, # Size effect exponent for small particles
        size_effect_exp_medium = 0.15, # Size effect exponent for medium particles
        size_effect_exp_large = 0.1,   # Size effect exponent for large particles
        time_factor = 0.6,          # Time-dependent factor
        time_offset = 0.3,          # Time offset parameter
        layer_adjust_factor = 0.15, # Layer thickness adjustment
        k_off_adjust = 0.15         # k_off size adjustment factor
      )
      
      # Define parameter variation ranges (±50% around base values)
      self$param_ranges <- list(
        k_on_base = c(1.0e3, 3.0e3),
        k_off_base = c(0.75e-3, 2.25e-3),
        c_init = c(3.0e-4, 9.0e-4),
        L_base = c(14.45, 43.35),
        coating_factor_citrate = c(0.461, 1.383),
        coating_factor_PVP = c(0.539, 1.617),
        size_scaling_small = c(0.5, 1.5),
        size_scaling_medium = c(0.425, 1.275),
        size_scaling_large = c(0.25, 0.75),
        cooperative_factor = c(0.2, 0.6),
        surface_exp = c(0.75, 2.25),
        crowding_factor = c(0.125, 0.375),
        size_effect_exp_small = c(0.125, 0.375),
        size_effect_exp_medium = c(0.075, 0.225),
        size_effect_exp_large = c(0.05, 0.15),
        time_factor = c(0.3, 0.9),
        time_offset = c(0.15, 0.45),
        layer_adjust_factor = c(0.075, 0.225),
        k_off_adjust = c(0.075, 0.225)
      )
    },
    
    # One-at-a-time sensitivity analysis (TORNADO PLOT METHOD)
    run_oat_sensitivity = function(target_variable = "size", time_point = 1.0, 
                                   initial_size = 38.2, coating = "citrate", n_points = 21) {
      
      results_list <- list()
      
      for (param_name in names(self$param_ranges)) {
        cat("Analyzing AgNP-BSA parameter:", param_name, "\n")
        
        # Get baseline output
        base_output <- self$simulate_size_single(self$base_params, time_point, initial_size, coating)
        
        # Calculate outputs at +50% and -50% of parameter value
        param_range <- self$param_ranges[[param_name]]
        param_low <- param_range[1]   # -50%
        param_high <- param_range[2]  # +50%
        
        # Low value (-50%)
        modified_params_low <- self$base_params
        modified_params_low[[param_name]] <- param_low
        output_low <- self$simulate_size_single(modified_params_low, time_point, initial_size, coating)
        
        # High value (+50%)
        modified_params_high <- self$base_params
        modified_params_high[[param_name]] <- param_high
        output_high <- self$simulate_size_single(modified_params_high, time_point, initial_size, coating)
        
        # Calculate tornado plot values
        low_effect <- (output_low - base_output) / base_output * 100    # % change from baseline
        high_effect <- (output_high - base_output) / base_output * 100  # % change from baseline
        
        # Store results for tornado plot
        results_list[[param_name]] <- data.frame(
          parameter = param_name,
          baseline_output = base_output,
          low_output = output_low,
          high_output = output_high,
          low_effect = low_effect,      # Negative side of tornado
          high_effect = high_effect,    # Positive side of tornado
          total_swing = abs(high_effect - low_effect),  # Total sensitivity
          stringsAsFactors = FALSE
        )
      }
      
      # Combine all results
      sensitivity_df <- bind_rows(results_list)
      self$sensitivity_results <- sensitivity_df
      
      return(sensitivity_df)
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
    
    # Simulate size with modified parameters
    simulate_size_single = function(params, time_point, initial_size, coating) {
      # Determine size category and get appropriate parameters
      size_cat <- self$get_size_category(initial_size)
      
      # Get size-specific parameters
      if (size_cat == "small") {
        size_scaling <- params$size_scaling_small
        size_effect_exp <- params$size_effect_exp_small
      } else if (size_cat == "medium") {
        size_scaling <- params$size_scaling_medium
        size_effect_exp <- params$size_effect_exp_medium
      } else {
        size_scaling <- params$size_scaling_large
        size_effect_exp <- params$size_effect_exp_large
      }
      
      # Get coating factor
      coating_factor <- if (coating == "citrate") params$coating_factor_citrate else params$coating_factor_PVP
      
      # AgNP-BSA size evolution model
      model <- function(t, state, parameters) {
        with(as.list(c(state, parameters)), {
          theta <- state[1]
          
          # Core model with modifying factors
          # Surface Availability Term (S)
          S <- (1 - theta)^surface_exp
          
          # Cooperative Binding Term (F_coop) - modified for silver
          F_coop <- 1 + theta * (d_initial/40)^cooperative_factor
          
          # Crowding Effect Term (F_crowd) - modified for silver
          F_crowd <- exp(-crowding_factor * theta)
          
          # Surface coating effect on binding (F_coat)
          F_coat <- coating_factor
          
          # Time-dependent k_on
          k_on <- k_on_base * size_effect * time_effect * F_coat
          
          # Size-dependent k_off
          k_off <- k_off_base / (1 + k_off_adjust * log(d_initial/40 + 1))
          
          # Core equation
          dtheta <- k_on * c_init * S * F_coop * F_crowd - k_off * theta
          
          return(list(c(dtheta)))
        })
      }
      
      # Calculate size and time effects
      size_effect <- (initial_size/40)^size_effect_exp
      time_effect <- (1 + params$time_factor*time_point)/(params$time_offset + time_point)
      
      # Solve ODE with increased tolerances to reduce warnings
      times <- seq(0, time_point, by = 0.01)
      out <- ode(y = c(0), times = times, func = model, 
                 parms = c(params, d_initial = initial_size, size_effect = size_effect, 
                           time_effect = time_effect, coating_factor = coating_factor),
                 rtol = 1e-4, atol = 1e-6)  # Increased tolerances
      
      # Calculate final size
      final_coverage <- out[nrow(out), 2]
      
      # Calculate adjusted protein layer thickness
      L <- params$L_base * coating_factor * size_scaling * 
        (1 + params$layer_adjust_factor * log(initial_size/40 + 1))
      
      final_size <- initial_size + 2 * L * final_coverage
      
      return(final_size)
    },
    
    # Calculate sensitivity indices for tornado plot
    calculate_sensitivity_indices = function() {
      if (is.null(self$sensitivity_results)) {
        stop("Run sensitivity analysis first")
      }
      
      # Calculate tornado plot metrics
      indices <- self$sensitivity_results %>%
        mutate(
          max_effect = pmax(abs(low_effect), abs(high_effect)),
          total_range = abs(high_effect - low_effect)
        ) %>%
        arrange(desc(total_range)) %>%
        mutate(importance_rank = row_number())
      
      return(indices)
    },
    
    # Create proper tornado plot
    plot_tornado = function(target_variable = "Size") {
      indices <- self$calculate_sensitivity_indices()
      
      # Prepare data for tornado plot
      tornado_data <- indices %>%
        select(parameter, low_effect, high_effect, total_range) %>%
        arrange(total_range) %>%  # Order by sensitivity (smallest to largest)
        mutate(
          parameter = factor(parameter, levels = parameter),  # Preserve order
          low_effect = ifelse(low_effect > 0, 0, low_effect),  # Ensure low effects are ≤ 0
          high_effect = ifelse(high_effect < 0, 0, high_effect)  # Ensure high effects are ≥ 0
        )
      
      # Create tornado plot
      p <- ggplot(tornado_data) +
        geom_col(aes(x = parameter, y = low_effect), 
                 fill = "coral", width = 0.7, alpha = 0.8) +
        # Right side (positive effects)  
        geom_col(aes(x = parameter, y = high_effect), 
                 fill = "steelblue", width = 0.7, alpha = 0.8) +
        coord_flip() +
        geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
        labs(
          title = paste("AgNP-BSA Parameter Sensitivity Analysis -", target_variable),
          subtitle = "Tornado Plot: Impact of ±50% parameter variation on model output",
          x = "Model Parameters",
          y = "% Change in Output from Baseline",
          caption = "Orange: When parameter decreases by 50%  | Blue: When parameter increases by 50%"
        ) +
        theme_minimal() +
        theme(
          plot.title = element_text(size = 14, face = "bold"),
          plot.subtitle = element_text(size = 11),
          axis.title = element_text(size = 12),
          plot.caption = element_text(size = 9, color = "gray50")
        ) +
        scale_y_continuous(labels = function(x) paste0(x, "%"))
      
      return(p)
    },
    
    # Create parameter response curves
    plot_response_curves = function(top_n = 6, target_variable = "size", coating = "citrate") {
      if (is.null(self$sensitivity_results)) {
        stop("Run sensitivity analysis first")
      }
      
      # Get top N most sensitive parameters based on total range
      indices <- self$calculate_sensitivity_indices()
      top_params <- indices$parameter[1:min(top_n, nrow(indices))]
      
      # For response curves, we need to generate detailed parameter sweeps
      response_data <- list()
      
      for (param_name in top_params) {
        # Create detailed parameter sweep
        param_range <- self$param_ranges[[param_name]]
        param_values <- seq(param_range[1], param_range[2], length.out = 21)
        
        outputs <- numeric(length(param_values))
        base_output <- indices$baseline_output[indices$parameter == param_name]
        
        for (i in 1:length(param_values)) {
          modified_params <- self$base_params
          modified_params[[param_name]] <- param_values[i]
          
          outputs[i] <- self$simulate_size_single(modified_params, 1.0, 38.2, "citrate")
        }
        
        response_data[[param_name]] <- data.frame(
          parameter = param_name,
          param_normalized = (param_values - self$base_params[[param_name]]) / self$base_params[[param_name]] * 100,
          output_normalized = (outputs - base_output) / base_output * 100
        )
      }
      
      # Combine and plot
      plot_data <- bind_rows(response_data) %>%
        mutate(parameter = factor(parameter, levels = top_params))
      
      target_label <- "Size Evolution"
      coating_label <- toupper(coating)
      y_label <- "Size Change (%)"
      
      p <- ggplot(plot_data, aes(x = param_normalized, y = output_normalized)) +
        geom_line(color = "steelblue", linewidth = 1) +
        geom_point(color = "darkblue", size = 1.5) +
        geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
        geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
        facet_wrap(~parameter, scales = "free", ncol = 3) +
        labs(
          title = paste("AgNP-BSA Parameter Response Curves -", target_label, "(", coating_label, ")"),
          subtitle = paste("Impact of parameter changes on size evolution with", coating, "coating"),
          x = "Parameter Change (%)",
          y = y_label
        ) +
        theme_minimal() +
        theme(
          plot.title = element_text(size = 14, face = "bold"),
          strip.text = element_text(size = 10, face = "bold")
        )
      
      return(p)
    },
    
    # Generate comprehensive sensitivity report
    generate_report = function(target_variable = "size") {
      cat("\n=== AgNP-BSA SENSITIVITY ANALYSIS REPORT ===\n")
      cat("Target Variable:", toupper(target_variable), "\n")
      cat("Analysis Date:", Sys.Date(), "\n\n")
      
      indices <- self$calculate_sensitivity_indices()
      
      cat("TOP 5 MOST INFLUENTIAL PARAMETERS:\n")
      cat("===================================\n")
      for (i in 1:min(5, nrow(indices))) {
        param <- indices[i, ]
        cat(sprintf("%d. %s\n", i, param$parameter))
        cat(sprintf("   Total Sensitivity Range: %.1f%%\n", param$total_range))
        cat(sprintf("   Low Effect (-50%%): %.1f%%\n", param$low_effect))
        cat(sprintf("   High Effect (+50%%): %.1f%%\n", param$high_effect))
        cat("\n")
      }
      
      cat("PARAMETER CATEGORIES:\n")
      cat("====================\n")
      kinetic_params <- c("k_on_base", "k_off_base", "c_init")
      physical_params <- c("L_base", "layer_adjust_factor", "k_off_adjust")
      coating_params <- c("coating_factor_citrate", "coating_factor_PVP")
      size_params <- c("size_scaling_small", "size_scaling_medium", "size_scaling_large", 
                       "size_effect_exp_small", "size_effect_exp_medium", "size_effect_exp_large")
      model_params <- c("cooperative_factor", "surface_exp", "crowding_factor", "time_factor", "time_offset")
      
      for (category in list(
        list(name = "Kinetic Parameters", params = kinetic_params),
        list(name = "Physical Parameters", params = physical_params),
        list(name = "Coating Parameters", params = coating_params),
        list(name = "Size-Dependent Parameters", params = size_params),
        list(name = "Model Parameters", params = model_params)
      )) {
        cat(sprintf("%s:\n", category$name))
        cat_params <- indices[indices$parameter %in% category$params, ]
        if (nrow(cat_params) > 0) {
          for (j in 1:nrow(cat_params)) {
            cat(sprintf("  - %s (Range: %.1f%%, Rank: %d)\n", 
                        cat_params$parameter[j], 
                        cat_params$total_range[j],
                        cat_params$importance_rank[j]))
          }
        }
        cat("\n")
      }
      
      return(indices)
    }
  )
)

# =============================================================================
# 5. COMPREHENSIVE EXECUTION FUNCTION
# =============================================================================

run_complete_nanoparticle_sensitivity_analysis <- function() {
  cat("=============================================================================================================\n")
  cat("COMPREHENSIVE NANOPARTICLE PROTEIN CORONA SIZE and ZETA POTENTIAL EVOLUTION MODEL SENSITIVITY ANALYSIS\n")
  cat("=============================================================================================================\n\n")
  
  # Store all results
  all_results <- list()
  
  # 1. AgNP-HSA Analysis
  cat("1. RUNNING AgNP-HSA SENSITIVITY ANALYSIS\n")
  cat("==========================================\n")
  
  tryCatch({
    # AgNP-HSA Size Analysis
    cat("Running AgNP-HSA Size Analysis...\n")
    agnp_hsa_size_analyzer <- AgNP_HSA_SensitivityAnalyzer$new()
    agnp_hsa_size_results <- agnp_hsa_size_analyzer$run_oat_sensitivity(
      target_variable = "size", 
      time_point = 1.0, 
      initial_size = 23.0
    )
    agnp_hsa_size_indices <- agnp_hsa_size_analyzer$generate_report("size")
    agnp_hsa_size_tornado <- agnp_hsa_size_analyzer$plot_tornado("Size Evolution")
    agnp_hsa_size_response <- agnp_hsa_size_analyzer$plot_response_curves(top_n = 6, target_variable = "size")
    
    cat("AgNP-HSA Size Analysis completed successfully!\n\n")
  }, error = function(e) {
    cat("Error in AgNP-HSA Size Analysis:", e$message, "\n\n")
    agnp_hsa_size_analyzer <- NULL
    agnp_hsa_size_indices <- NULL
    agnp_hsa_size_tornado <- NULL
    agnp_hsa_size_response <- NULL
  })
  
  tryCatch({
    # AgNP-HSA Zeta Analysis
    cat("Running AgNP-HSA Zeta Analysis...\n")
    agnp_hsa_zeta_analyzer <- AgNP_HSA_SensitivityAnalyzer$new()
    agnp_hsa_zeta_results <- agnp_hsa_zeta_analyzer$run_oat_sensitivity(
      target_variable = "zeta", 
      time_point = 1.0, 
      initial_size = 23.0
    )
    agnp_hsa_zeta_indices <- agnp_hsa_zeta_analyzer$generate_report("zeta potential")
    agnp_hsa_zeta_tornado <- agnp_hsa_zeta_analyzer$plot_tornado("Zeta Potential Evolution")
    agnp_hsa_zeta_response <- agnp_hsa_zeta_analyzer$plot_response_curves(top_n = 6, target_variable = "zeta")
    
    cat("AgNP-HSA Zeta Analysis completed successfully!\n\n")
  }, error = function(e) {
    cat("Error in AgNP-HSA Zeta Analysis:", e$message, "\n\n")
    agnp_hsa_zeta_analyzer <- NULL
    agnp_hsa_zeta_indices <- NULL
    agnp_hsa_zeta_tornado <- NULL
    agnp_hsa_zeta_response <- NULL
  })
  
  all_results$agnp_hsa <- list(
    size_analyzer = agnp_hsa_size_analyzer,
    zeta_analyzer = agnp_hsa_zeta_analyzer,
    size_indices = agnp_hsa_size_indices,
    zeta_indices = agnp_hsa_zeta_indices,
    size_tornado = agnp_hsa_size_tornado,
    zeta_tornado = agnp_hsa_zeta_tornado,
    size_response = agnp_hsa_size_response,
    zeta_response = agnp_hsa_zeta_response
  )
  
  # 2. AuNP-HSA Analysis
  cat("\n2. RUNNING AuNP-HSA SENSITIVITY ANALYSIS\n")
  cat("==========================================\n")
  
  tryCatch({
    # AuNP-HSA Size Analysis
    cat("Running AuNP-HSA Size Analysis...\n")
    aunp_hsa_size_analyzer <- AuNP_HSA_SensitivityAnalyzer$new()
    aunp_hsa_size_results <- aunp_hsa_size_analyzer$run_oat_sensitivity(
      target_variable = "size", 
      time_point = 1.0, 
      initial_size = 10.0
    )
    aunp_hsa_size_indices <- aunp_hsa_size_analyzer$generate_report("size")
    aunp_hsa_size_tornado <- aunp_hsa_size_analyzer$plot_tornado("Size Evolution")
    aunp_hsa_size_response <- aunp_hsa_size_analyzer$plot_response_curves(top_n = 6, target_variable = "size")
    
    cat("AuNP-HSA Size Analysis completed successfully!\n\n")
  }, error = function(e) {
    cat("Error in AuNP-HSA Size Analysis:", e$message, "\n\n")
    aunp_hsa_size_analyzer <- NULL
    aunp_hsa_size_indices <- NULL
    aunp_hsa_size_tornado <- NULL
    aunp_hsa_size_response <- NULL
  })
  
  tryCatch({
    # AuNP-HSA Zeta Analysis
    cat("Running AuNP-HSA Zeta Analysis...\n")
    aunp_hsa_zeta_analyzer <- AuNP_HSA_SensitivityAnalyzer$new()
    aunp_hsa_zeta_results <- aunp_hsa_zeta_analyzer$run_oat_sensitivity(
      target_variable = "zeta", 
      time_point = 24.0, 
      initial_size = 10.0,
      initial_zeta = -48.93
    )
    aunp_hsa_zeta_indices <- aunp_hsa_zeta_analyzer$generate_report("zeta potential")
    aunp_hsa_zeta_tornado <- aunp_hsa_zeta_analyzer$plot_tornado("Zeta Potential Evolution")
    aunp_hsa_zeta_response <- aunp_hsa_zeta_analyzer$plot_response_curves(top_n = 6, target_variable = "zeta")
    
    cat("AuNP-HSA Zeta Analysis completed successfully!\n\n")
  }, error = function(e) {
    cat("Error in AuNP-HSA Zeta Analysis:", e$message, "\n\n")
    aunp_hsa_zeta_analyzer <- NULL
    aunp_hsa_zeta_indices <- NULL
    aunp_hsa_zeta_tornado <- NULL
    aunp_hsa_zeta_response <- NULL
  })
  
  all_results$aunp_hsa <- list(
    size_analyzer = aunp_hsa_size_analyzer,
    zeta_analyzer = aunp_hsa_zeta_analyzer,
    size_indices = aunp_hsa_size_indices,
    zeta_indices = aunp_hsa_zeta_indices,
    size_tornado = aunp_hsa_size_tornado,
    zeta_tornado = aunp_hsa_zeta_tornado,
    size_response = aunp_hsa_size_response,
    zeta_response = aunp_hsa_zeta_response
  )
  
  # 3. AuNP-BSA Analysis
  cat("\n3. RUNNING AuNP-BSA SENSITIVITY ANALYSIS\n")
  cat("==========================================\n")
  
  tryCatch({
    cat("Running AuNP-BSA Size Analysis...\n")
    aunp_bsa_analyzer <- AuNP_BSA_SensitivityAnalyzer$new()
    aunp_bsa_results <- aunp_bsa_analyzer$run_oat_sensitivity(
      target_variable = "size", 
      time_point = 1.5, 
      initial_size = 9.4
    )
    aunp_bsa_indices <- aunp_bsa_analyzer$generate_report("size")
    aunp_bsa_tornado <- aunp_bsa_analyzer$plot_tornado("Size Evolution")
    aunp_bsa_response <- aunp_bsa_analyzer$plot_response_curves(top_n = 6, target_variable = "size")
    
    cat("AuNP-BSA Analysis completed successfully!\n\n")
  }, error = function(e) {
    cat("Error in AuNP-BSA Analysis:", e$message, "\n\n")
    aunp_bsa_analyzer <- NULL
    aunp_bsa_indices <- NULL
    aunp_bsa_tornado <- NULL
    aunp_bsa_response <- NULL
  })
  
  all_results$aunp_bsa <- list(
    analyzer = if(exists("aunp_bsa_analyzer")) aunp_bsa_analyzer else NULL,
    indices = if(exists("aunp_bsa_indices")) aunp_bsa_indices else NULL,
    tornado = if(exists("aunp_bsa_tornado")) aunp_bsa_tornado else NULL,
    response = if(exists("aunp_bsa_response")) aunp_bsa_response else NULL
  )
  
  # 4. AgNP-BSA Analysis (Citrate)
  cat("\n4. RUNNING AgNP-BSA SENSITIVITY ANALYSIS (CITRATE)\n")
  cat("====================================================\n")
  
  tryCatch({
    cat("Running AgNP-BSA Citrate Analysis...\n")
    agnp_bsa_citrate_analyzer <- AgNP_BSA_SensitivityAnalyzer$new()
    agnp_bsa_citrate_results <- agnp_bsa_citrate_analyzer$run_oat_sensitivity(
      target_variable = "size", 
      time_point = 1.0, 
      initial_size = 38.2,
      coating = "citrate"
    )
    agnp_bsa_citrate_indices <- agnp_bsa_citrate_analyzer$generate_report("size (citrate)")
    agnp_bsa_citrate_tornado <- agnp_bsa_citrate_analyzer$plot_tornado("Size Evolution (Citrate)")
    agnp_bsa_citrate_response <- agnp_bsa_citrate_analyzer$plot_response_curves(top_n = 6, target_variable = "size", coating = "citrate")
    
    cat("AgNP-BSA Citrate Analysis completed successfully!\n\n")
  }, error = function(e) {
    cat("Error in AgNP-BSA Citrate Analysis:", e$message, "\n\n")
    agnp_bsa_citrate_analyzer <- NULL
    agnp_bsa_citrate_indices <- NULL
    agnp_bsa_citrate_tornado <- NULL
    agnp_bsa_citrate_response <- NULL
  })
  
  # 5. AgNP-BSA Analysis (PVP)
  cat("\n5. RUNNING AgNP-BSA SENSITIVITY ANALYSIS (PVP)\n")
  cat("================================================\n")
  
  tryCatch({
    cat("Running AgNP-BSA PVP Analysis...\n")
    agnp_bsa_pvp_analyzer <- AgNP_BSA_SensitivityAnalyzer$new()
    agnp_bsa_pvp_results <- agnp_bsa_pvp_analyzer$run_oat_sensitivity(
      target_variable = "size", 
      time_point = 1.0, 
      initial_size = 40.5,
      coating = "PVP"
    )
    agnp_bsa_pvp_indices <- agnp_bsa_pvp_analyzer$generate_report("size (PVP)")
    agnp_bsa_pvp_tornado <- agnp_bsa_pvp_analyzer$plot_tornado("Size Evolution (PVP)")
    agnp_bsa_pvp_response <- agnp_bsa_pvp_analyzer$plot_response_curves(top_n = 6, target_variable = "size", coating = "PVP")
    
    cat("AgNP-BSA PVP Analysis completed successfully!\n\n")
  }, error = function(e) {
    cat("Error in AgNP-BSA PVP Analysis:", e$message, "\n\n")
    agnp_bsa_pvp_analyzer <- NULL
    agnp_bsa_pvp_indices <- NULL
    agnp_bsa_pvp_tornado <- NULL
    agnp_bsa_pvp_response <- NULL
  })
  
  all_results$agnp_bsa <- list(
    citrate_analyzer = agnp_bsa_citrate_analyzer,
    pvp_analyzer = agnp_bsa_pvp_analyzer,
    citrate_indices = agnp_bsa_citrate_indices,
    pvp_indices = agnp_bsa_pvp_indices,
    citrate_tornado = agnp_bsa_citrate_tornado,
    pvp_tornado = agnp_bsa_pvp_tornado,
    citrate_response = agnp_bsa_citrate_response,
    pvp_response = agnp_bsa_pvp_response
  )
  
  # Display plots - only if they exist
  cat("\n=========================================================================\n")
  cat("DISPLAYING TORNADO PLOTS FOR ALL SYSTEMS\n")
  cat("=========================================================================\n")
  
  if (!is.null(all_results$agnp_hsa$size_tornado)) {
    cat("AgNP-HSA Size Tornado Plot:\n")
    print(all_results$agnp_hsa$size_tornado)
  }
  
  if (!is.null(all_results$agnp_hsa$zeta_tornado)) {
    cat("AgNP-HSA Zeta Tornado Plot:\n")
    print(all_results$agnp_hsa$zeta_tornado)
  }
  
  if (!is.null(all_results$aunp_hsa$size_tornado)) {
    cat("AuNP-HSA Size Tornado Plot:\n")
    print(all_results$aunp_hsa$size_tornado)
  }
  
  if (!is.null(all_results$aunp_hsa$zeta_tornado)) {
    cat("AuNP-HSA Zeta Tornado Plot:\n")
    print(all_results$aunp_hsa$zeta_tornado)
  }
  
  if (!is.null(all_results$aunp_bsa$tornado)) {
    cat("AuNP-BSA Tornado Plot:\n")
    print(all_results$aunp_bsa$tornado)
  }
  
  if (!is.null(all_results$agnp_bsa$citrate_tornado)) {
    cat("AgNP-BSA Citrate Tornado Plot:\n")
    print(all_results$agnp_bsa$citrate_tornado)
  }
  
  if (!is.null(all_results$agnp_bsa$pvp_tornado)) {
    cat("AgNP-BSA PVP Tornado Plot:\n")
    print(all_results$agnp_bsa$pvp_tornado)
  }
  
  cat("\n=========================================================================\n")
  cat("DISPLAYING RESPONSE CURVES FOR ALL SYSTEMS\n")
  cat("=========================================================================\n")
  
  if (!is.null(all_results$agnp_hsa$size_response)) {
    cat("AgNP-HSA Size Response Curves:\n")
    print(all_results$agnp_hsa$size_response)
  }
  
  if (!is.null(all_results$agnp_hsa$zeta_response)) {
    cat("AgNP-HSA Zeta Response Curves:\n")
    print(all_results$agnp_hsa$zeta_response)
  }
  
  if (!is.null(all_results$aunp_hsa$size_response)) {
    cat("AuNP-HSA Size Response Curves:\n")
    print(all_results$aunp_hsa$size_response)
  }
  
  if (!is.null(all_results$aunp_hsa$zeta_response)) {
    cat("AuNP-HSA Zeta Response Curves:\n")
    print(all_results$aunp_hsa$zeta_response)
  }
  
  if (!is.null(all_results$aunp_bsa$response)) {
    cat("AuNP-BSA Response Curves:\n")
    print(all_results$aunp_bsa$response)
  }
  
  if (!is.null(all_results$agnp_bsa$citrate_response)) {
    cat("AgNP-BSA Citrate Response Curves:\n")
    print(all_results$agnp_bsa$citrate_response)
  }
  
  if (!is.null(all_results$agnp_bsa$pvp_response)) {
    cat("AgNP-BSA PVP Response Curves:\n")
    print(all_results$agnp_bsa$pvp_response)
  }
  
  cat("\n=========================================================================\n")
  cat("SENSITIVITY ANALYSIS COMPLETE - ALL SYSTEMS ANALYZED\n")
  cat("=========================================================================\n")
  
  return(all_results)
}

# =============================================================================
# 6. EXECUTE COMPLETE ANALYSIS
# =============================================================================

  # Execute the complete comprehensive sensitivity analysis
  complete_sensitivity_results <- run_complete_nanoparticle_sensitivity_analysis()
