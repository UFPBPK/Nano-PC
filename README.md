# Nano-PC
# Nanoparticle Size and Zeta Potential Evolution Models
**📅 Date**: May 2025
## 📋 Project Description

This repository contains comprehensive computational models for studying the size and zeta potential evolution of silver (AgNP) and gold (AuNP) nanoparticles in the presence of different protein coronas. The models analyze how nanoparticle properties change over time when interacting with Bovine Serum Albumin (BSA) and Human Serum Albumin (HSA) using advanced differential equation systems and object-oriented programming.

## 📁 Repository Contents

### 🥈 Silver Nanoparticle (AgNP) Models

#### 🔧 Parameter Optimization
- **`AgNP-HSA ParameterOptimization_May22.R`**
  - **Advanced multi-stage optimization** for crowding factor parameters (a, b)
  - Three-stage approach: broad parameter sweep → focused search → gradient descent
  - Uses Zhao et al. (2021) AgNP20 experimental data for calibration
  - Includes visualization tools (error heatmaps, convergence plots)
  - Comprehensive validation with RMSE and error percentage calculations

#### 📏 Size Evolution Models
- **`AgNP-BSA Size Evolution Model_May22.R`**
  - **R6 object-oriented model** for BSA corona formation on silver nanoparticles
  - Uses differential equations with coating-specific effects (citrate vs PVP)
  - Includes calibration and validation datasets with multiple size ranges
  - Experimental data sources: Shannahan et al. (2015a), Shannahan et al. (2015b), Zhang et al. (2023)

- **`AgNP-HSA Size Evolution Model_May22.R`**
  - **R6 object-oriented model** for HSA corona formation on silver nanoparticles
  - Advanced differential equations with size-dependent kinetics
  - Calibration: Zhao et al. (2021) AgNP20 data
  - Validation: Zhao et al. (2021) AgNP40/AgNP80, Shannahan et al. (2015a) AgNP19 data
  - Comprehensive validation reporting with error metrics

#### ⚡ Zeta Potential Evolution Models
- **`AgNP-HSA Zeta Potential Evolution Model_May22.R`**
  - **R6 object-oriented model** for zeta potential evolution of AgNP-HSA systems
  - Complex charge dynamics with size, charge, and incubation time factors
  - Calibration: Zhao et al. (2021) SNP20 data
  - Validation: Zhao et al. (2021) SNP40/SNP80, Shannahan et al. (2015a)
  - Includes sophisticated performance metrics and charge-type analysis

### 🥇 Gold Nanoparticle (AuNP) Models

#### 📏 Size Evolution Models
- **`AuNP-BSA Size Evolution Model_May22.R`**
  - **R6 object-oriented model** for BSA corona formation on gold nanoparticles
  - Calibration: Hanigan-Diebel et al. (2024) AuNP-MPTMA/MEEE/MHA data
  - Validation: Kumar et al. (2024), Mishra & Das (2022) multiple size series
  - Advanced size-dependent kinetics and protein layer calculations

- **`AuNP-HSA Size Evolution Model_May22.R`**
  - **R6 object-oriented model** for HSA corona formation on gold nanoparticles
  - Calibration: Guglielmelli et al. (2023) 10 nm AuNP data
  - Validation: Extensive dataset from Goy-Lopez et al. (2012) (5-100 nm series), Capomaccio et al. (2015) 18.5 nm data
  - Size-dependent protein layer thickness calculations

#### ⚡ Zeta Potential Evolution Models
- **`AuNP-HSA Zeta Potential Evolution Model_May22.R`**
  - **R6 object-oriented model** for zeta potential evolution with multiple coating types
  - Calibration: Dai et al. (2023) CTAB-coated AuNPs
  - Validation: Dell'Aglio et al. (2021), Goy-Lopez et al. (2012) citrate-coated series
  - Coating-specific effects (CTAB vs citrate), charge-dependent dynamics
  - Advanced validation reporting with coating-type and charge-type stratification

### 📊 Comprehensive Sensitivity Analysis Framework

- **`Complete_Sensitivity_Analysis_Framework_Jun25.R`**
  - **Comprehensive parameter sensitivity analysis** for all NP-PC systems
  - **R6 object-oriented framework** with four specialized analyzer classes:
    - `AgNP_HSA_SensitivityAnalyzer`: Silver-HSA size and zeta potential sensitivity
    - `AuNP_HSA_SensitivityAnalyzer`: Gold-HSA size and zeta potential sensitivity  
    - `AuNP_BSA_SensitivityAnalyzer`: Gold-BSA size evolution sensitivity
    - `AgNP_BSA_SensitivityAnalyzer`: Silver-BSA size evolution sensitivity (citrate and PVP coatings)
  - **One-at-a-time (OAT) sensitivity analysis** with ±50% parameter variation
  - **Advanced visualization capabilities**: tornado plots, and parameter response curves
  - **Automated execution function** that runs complete analysis across all systems
  - **Comprehensive reporting** with parameter ranking and categorical analysis



## ⚙️ Prerequisites

### 📦 Required R Packages
```r
# Install required packages if not already installed
install.packages(c(
  "deSolve",      # Differential equations (REQUIRED for ALL scripts)
  "dplyr",        # Data manipulation (REQUIRED for ALL scripts)
  "R6",           # Object-oriented programming (REQUIRED for ALL scripts EXCEPT Parameter Optimization)
  "ggplot2",      # Data visualization (Parameter Optimization only)
  "tidyr",        # Data reshaping (Parameter Optimization only)
  "gridExtra",    # Multiple plots (Parameter Optimization only)
  "grid"          # Graphics operations (Parameter Optimization only)
))

# Note: 'stats' package is included with base R installation
```

### 💻 System Requirements
- R version 4.0.0 or higher
- RStudio (recommended for interactive use)
- Sufficient memory for differential equation solving (especially for long time simulations)

## 🚀 Usage

### ▶️ Running Individual Models
Each script is **completely self-contained** and **runs automatically** when sourced:

```r
# ALL scripts execute immediately when sourced - no manual steps needed
source("AgNP-BSA Size Evolution Model_May22.R")        # Runs automatically
source("AgNP-HSA Size Evolution Model_May22.R")        # Runs automatically  
source("AuNP-BSA Size Evolution Model_May22.R")        # Runs automatically
source("AuNP-HSA Size Evolution Model_May22.R")        # Runs automatically
source("AgNP-HSA Zeta Potential Evolution Model_May22.R")  # Runs automatically
source("AuNP-HSA Zeta Potential Evolution Model_May22.R")  # Runs automatically
source("AgNP-HSA ParameterOptimization_May22.R")       # Runs automatically
```

**What happens automatically:**
1. Libraries are loaded
2. Models are defined and calibrated
3. Simulations are executed
4. Results are printed to console
5. Data is stored in `simulation_output` variable

### 🔍 Accessing Results After Sourcing
```r
# After sourcing any script, results are automatically available:
results <- simulation_output$results
validation <- simulation_output$validation
summary_table <- simulation_output$summary_table  # (if available)

# You can also create custom simulations:
simulator <- AgNP_HSASizeSimulator$new()
custom_results <- simulator$simulate_all_data()
```

## 📈 Sensitivity Analysis

### 🎯 Running Complete Sensitivity Analysis
The sensitivity analysis framework provides comprehensive parameter sensitivity assessment across all nanoparticle-protein systems:

```r
# Execute complete sensitivity analysis for all systems
source("Complete_Sensitivity_Analysis_Framework_May28.R")

# Results are automatically stored in complete_sensitivity_results
tornado_plots <- complete_sensitivity_results$agnp_hsa$size_tornado
response_curves <- complete_sensitivity_results$aunp_hsa$zeta_response
sensitivity_indices <- complete_sensitivity_results$aunp_bsa$indices

# Analyze specific systems independently
agnp_analyzer <- AgNP_HSA_SensitivityAnalyzer$new()
size_sensitivity <- agnp_analyzer$run_oat_sensitivity("size", time_point = 1.0, initial_size = 23.0)
tornado_plot <- agnp_analyzer$plot_tornado("Size Evolution")
response_curves <- agnp_analyzer$plot_response_curves(top_n = 6)
sensitivity_report <- agnp_analyzer$generate_report("size")
```

### ✨ Sensitivity Analysis Features
- **📊 Parameter Variation Range**: ±50% around calibrated base values
- **🔍 Analysis Method**: One-at-a-time (OAT) sensitivity analysis
- **🎯 Target Variables**: Size evolution and zeta potential evolution
- **📈 Visualization Options**: Tornado plots, parameter response curves, sensitivity rankings
- **📐 Statistical Metrics**: Total sensitivity range, individual parameter effects, importance rankings

## 📊 Data Requirements

The models use experimental data from published literature:

### 🥈 Silver Nanoparticle Data Sources:
- **Zhao et al. (2021)**: Calibration & internal validation dataset for AgNP-HSA size evolution model, and calibration & internal validation dataset for AgNP-HSA zeta potential evolution model
- **Shannahan et al. (2015a)**: Validation dataset for AgNP-HSA & AgNP-BSA size evolution model, and external validation dataset for AgNP-HSA zeta potential evolution model
- **Shannahan et al. (2015b)**: Calibration & internal validation dataset for AgNP-BSA size evolution model
- **Zhang et al. (2023)**: External validation dataset for AgNP-BSA size evolution model
- **Hanigan-Diebel et al. (2024)**: 
- ***Shannahan et al. (2015a): Shannahan JH, Podila R, Aldossari AA, Emerson H, Powell BA, Ke PC, et al. Formation of a Protein Corona on Silver Nanoparticles Mediates Cellular Toxicity via Scavenger Receptors. Toxicological Sciences. 2015;143(1):136-46. ***
- ***Shannahan et al. (2015b): Shannahan JH, Podila R, Brown JM. A hyperspectral and toxicological analysis of protein corona impact on silver nanoparticle properties, intracellular modifications, and macrophage activation. Int J Nanomed. 2015;10:6509-21.***

### 🥇 Gold Nanoparticle Data Sources:
- **Guglielmelli et al. (2023)**: Calibration dataset for AuNP-HSA size evolution model
- **Capomaccio et al. (2015)**: External validation dataset for AuNP-HSA size evolution model
- **Goy-Lopez et al. (2012)**: External validation dataset for AuNP-HSA size and zeta potential evolution model
- **Hanigan-Diebel et al. (2024)**: Calibration dataset for AuNP-BSA size evolution model
- **Kumar et al. (2024)**: External validation dataset for AuNP-BSA size evolution model
- **Mishra & Das (2022)**: External validation dataset for AuNP-BSA size evolution model
- **Dai et al. (2023)**: Calibration dataset for AuNP-HSA zeta potential evolution model
- **Dell'Aglio et al. (2021)**: External validation dataset for AuNP-HSA zeta potential evolution model




## 📈 Model Outputs

Each script generates:

### 📏 Size Evolution Models:
- **Time-series plots**: Nanoparticle diameter vs. time
- **Validation tables**: Comparison of simulated vs. experimental final sizes
- **Error metrics**: Percentage errors and RMSE calculations
- **Surface coverage dynamics**: Protein binding kinetics over time

### ⚡ Zeta Potential Models:
- **Charge evolution plots**: Zeta potential vs. time curves  
- **Comprehensive validation reports**: Multi-dataset performance analysis
- **Stratified error analysis**: Performance by coating type, charge polarity, size range
- **Surface coverage correlations**: Relationship between binding and charge changes

### Parameter Optimization:
- **Error heatmaps**: Parameter space exploration visualization
- **Convergence plots**: Optimization algorithm performance
- **Final parameter values**: Optimized crowding factors with uncertainty estimates

## 🔑 Key Model Features

### 🧮 Advanced Mathematical Framework:
- **Differential equation systems**: Complex ODE models with multiple modifying factors
- **Size-dependent kinetics**: Protein binding rates vary with nanoparticle diameter
- **Cooperative binding effects**: Protein-protein interactions on nanoparticle surfaces
- **Crowding effects**: Non-linear protein packing dynamics
- **Time-dependent processes**: Evolving association/dissociation rates

### Coating-Specific Modeling:
- **Surface chemistry effects**: Different behavior for citrate, PVP, CTAB coatings
- **Charge-dependent interactions**: Positive vs. negative initial zeta potentials
- **Material-specific parameters**: Distinct models for silver vs. gold nanoparticles

### ✅ Comprehensive Validation:
- **Multi-dataset validation**: Testing against multiple independent experimental studies
- **Error stratification**: Performance analysis by size, charge, coating, time
- **Statistical metrics**: RMSE, MAE, percentage errors across all datasets

### Sensitivity Analysis Outputs:
- **Tornado plots**: Parameter importance ranking with directional effects visualization
- **Response curves**: Detailed parameter-output relationships across variation ranges
- **Sensitivity indices**: Quantitative parameter importance metrics and rankings
- **Comprehensive reports**: Statistical analysis with parameter categorization and top influential factors
- **Interactive visualizations**: Publication-ready plots with customizable aesthetics
- **Cross-system comparisons**: Comparative sensitivity analysis across different NP-PC systems

## 📂 File Organization

```
├── Silver Nanoparticle (AgNP) Models/
│   ├── AgNP-HSA ParameterOptimization_May22.R
│   ├── AgNP-BSA Size Evolution Model_May22.R
│   ├── AgNP-HSA Size Evolution Model_May22.R
│   └── AgNP-HSA Zeta Potential Evolution Model_May22.R
├── Gold Nanoparticle (AuNP) Models/
│   ├── AuNP-BSA Size Evolution Model_May22.R
│   ├── AuNP-HSA Size Evolution Model_May22.R
│   └── AuNP-HSA Zeta Potential Evolution Model_May22.R
└── Sensitivity Analysis Framework/
    └── Complete_Sensitivity_Analysis_Framework_May28.R
```

## 🧮 Mathematical Model Details

### ⚙️ Core Differential Equation:
```
dθ/dt = k_on × C_protein × S × F_coop × F_crowd - k_off × θ^n
```

Where:
- θ represents the protein surface coverage (ranging from 0 to 1) (unitless)
- dθ/dt is the rate of change in surface coverage over time (unit: s⁻¹ (per second))
- k_on is the association rate constant (unit: M⁻¹·s⁻¹ (inverse molar per second))
- C_protein is the initial protein concentration (unit: mol/L or M (molar))
- S is the surface availability term (unitless)
- F_coop is the cooperative binding term (unitless)
- F_crowd is the crowding effect term (unitless)
- k_off is the dissociation rate constant (unit: s⁻¹ (per second))
- n is an exponent that varies depending on the specific NP-PC system (unitless)


### 📏 Size Calculation:
```
d(t) = d_initial + 2L × θ(t)
```

Where L is the protein layer thickness (varies by protein type and nanoparticle size)

### ⚡ Zeta Potential Calculation:
```
ζ(t) = ζ_initial + Δζ_max × F_charge × F_size × F_coating × θ(t)
```

Where θ(t) represents the coverage-dependent charge relationship

## 📊 Sensitivity Analysis Methodology

### 🎯 One-at-a-Time (OAT) Analysis Framework
The sensitivity analysis employs a systematic OAT methodology to assess parameter influence on model outputs:
Parameter Sensitivity = (Output_high - Output_baseline) / Output_baseline × 100%
Total Sensitivity Range = |Effect_high - Effect_low|

### 📋 Parameter Categories Analyzed:
- **Kinetic Parameters**: Association/dissociation rate constants, protein concentrations
- **Physical Parameters**: Layer thickness, surface coating factors, size-dependent adjustments
- **Model-Specific Parameters**: Cooperative binding factors, crowding effects, surface availability terms
- **System-Specific Parameters**: Charge-dependent factors, coating-specific modifications, time-dependent terms

### 📊 Tornado Plot Interpretation:
- **Bar Length**: Indicates parameter importance (longer bars = higher sensitivity)
- **Bar Direction**: Shows effect polarity (left = negative effect, right = positive effect)
- **Parameter Ranking**: Ordered by total sensitivity range for clear importance hierarchy

## 📝 Citation

If you use these models in your research, please cite:
```
Chen, X., & Lin, Z. (2025). Mechanistic Modeling of Size and Surface Charge Changes of Gold and Silver Nanoparticles During In Vitro Nanoparticle Protein Corona Formation. (under review)
```

## 📞 Contact

**👥 Authors**: Xinyue Chen, Zhoumeng Lin  
**📧 Email**: linzhoumeng@ufl.edu  
**🏛️ Institution**: University of Florida  
- Department of Environmental and Global Health, College of Public Health and Health Professions
- Center for Environmental and Human Toxicology  
- Center for Pharmacometrics and Systems Pharmacology  



## 🙏 Acknowledgments

- Funding sources: The study was supported by the National Institute of Biomedical Imaging and Bioengineering of National Institutes of Health (R01EB031022 and R03EB026045).
- Experimental data providers: Multiple research groups as cited in individual model files

---

**Note**: All models are pre-calibrated with experimental data and **execute automatically when sourced**. Simply source any script file and the complete simulation will run immediately, displaying results in the console and storing data in the `simulation_output` variable. No manual function calls are required.
