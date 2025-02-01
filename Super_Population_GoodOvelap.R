# Required R packages for data generation and analysis
library(survey)         # Survey data analysis
library(doParallel)     # Parallel processing
library(foreach)        # Looping in parallel
library(xtable)         # Table formatting for LaTeX
library(data.table)     # Data manipulation
library(sampling)       # Sampling techniques
library(MatchIt)        # Propensity score matching
library(Hmisc)          # Data utility functions
library(nloptr)         # Nonlinear optimization
library(numDeriv)       # Numerical derivatives
library(MASS)           # Statistical methods
library(patchwork)      # Plot layout organization
library(boot)           # Bootstrap resampling
library(magrittr)       # Pipe operator (%>%)

# Set the working directory
setwd("/Volumes/T7 Shield/R_SW_Final_Data_0118_2025/TestReformSimulation/")

# Define global parameters for population and sample generation
N           = 1000000            # Total population size
n           = N * 0.005          # Sample size (0.5% of the population)
n_covar     = 7                  # Number of covariates (6 regular + 1 associated with missingness if needed)
n_strata    = 10                 # Number of strata in the population
n_cluster   = 20                 # Number of clusters within each stratum
n_scenarios = 1                  # Number of scenarios to consider
tau_strata  = 0.35               # Standard deviation of covariate means across strata
min_prob    = 0.000              # Minimum probability threshold to avoid near-zero values
tau_cluster = c(0.15)            # Standard deviation of covariate means across clusters
s.m         = 1                  # Scaling factor for the simulation
# n.reps      = 5000             # Number of simulation replicates

# Set the seed for reproducibility
set.seed(20220930)


############################################################################################################################
# Generate Super Population Data
############################################################################################################################

# Initialize list to store datasets across scenarios
data_sets = list()


for(scenario in 1:n_scenarios){
  # Means in Matrix
  mu_strata  = matrix(rnorm(n_strata * n_covar, mean = 0, sd = tau_strata), nrow = n_strata, ncol = n_covar)  # Mean structure for strata
  mu_cluster = matrix(rnorm(n_cluster * n_covar, mean = 0, sd = tau_cluster[scenario]), nrow = n_cluster, ncol = n_covar)  # Mean structure for clusters
  
  # Identifying the STRATA EFFECT by stratum and covariates
  # Define row names (Strata)
  strata_row = rep(NA, n_strata)
  for (strata in 1:n_strata) {
    strata_row[strata] = paste("Strata", strata, sep = "_")
  }
  rownames(mu_strata) = strata_row
  
  # Define column names (Covariates)
  strata_col = rep(NA, n_covar)
  for (covariate in 1:n_covar) {
    strata_col[covariate] = paste("Covariate", covariate, sep = "_")
  }
  colnames(mu_strata) = strata_col
  
  # Identifying the CLUSTER EFFECT by cluster and covariates
  # Define row names (Clusters)
  cluster_row = rep(NA, n_cluster)
  for (cluster in 1:n_cluster) {
    cluster_row[cluster] = paste("Cluster", cluster, sep = "_")
  }
  rownames(mu_cluster) = cluster_row
  
  # Define column names (Covariates)
  cluster_col = rep(NA, n_covar)
  for (covariate in 1:n_covar) {
    cluster_col[covariate] = paste("Covariate", covariate, sep = "_")
  }
  colnames(mu_cluster) = cluster_col
  
  ############################################################################################################################
  # CREATING THE COVARIATES
  ############################################################################################################################
  
  cluster_l  = list()  # Temporary storage for cluster-level data
  strata_l   = list()  # Temporary storage for stratum-level data
  variable_l = list()  # Temporary storage for each covariate
  
  # Generate covariates for all subjects across clusters and strata
  for (variable in 1:n_covar) {  # Iterate over covariates
    for (strata in 1:n_strata) {  # Iterate over strata
      for (cluster in 1:n_cluster) {  # Iterate over clusters
        
        # Generate data for subjects in a specific cluster within a stratum
        data = cbind(
          rep((strata - 1) * n_cluster + cluster, N / (n_strata * n_cluster)),  # Cluster index
          rep(strata, N / (n_strata * n_cluster)),  # Stratum index
          rnorm(N / (n_strata * n_cluster), mean = mu_cluster[cluster, variable] + mu_strata[strata, variable])  # Covariate values
        )
        
        colnames(data) = c("Cluster", "Strata", paste("X", variable, sep = ""))  # Define column names
        cluster_l[[cluster]] = data
      }
      strata_l[[strata]] = do.call("rbind", cluster_l)  # Combine cluster data within each stratum
    }
    variable_l[[variable]] = data.frame(do.call("rbind", strata_l))  # Combine stratum data for each covariate
  }
  
  # Combine all covariate data into a single dataset
  data = variable_l[[1]]  # Initialize with the first covariate
  for (variable in 2:n_covar) {
    data = cbind(data, variable_l[[variable]][, 3])  # Add covariate values (3rd column corresponds to covariate values Xi)
  }
  
  # Define column names for the covariate dataset
  var_col = rep(NA, n_covar)
  for (covariate in 1:n_covar) {
    var_col[covariate] = paste("x", covariate, sep = "")  # Define covariate names
  }
  
  colnames(data) = c("Cluster", "Strata", var_col)  # Assign final column names

  
  ############################################################################################################################
  # PROPENSITY SCORE MODEL
  ############################################################################################################################
  
  # Logistic (expit) function for probability calculation
  expit = function(x) { exp(x) / (1 + exp(x)) }
  
  # Define coefficients for the logistic regression model
  a0 = log(35 / 80)  # Intercept
  a1 = log(1.10)     # Coefficient for x1
  a2 = log(1.25)     # Coefficient for x2
  a3 = log(1.50)     # Coefficient for x3
  a4 = log(1.75)     # Coefficient for x4
  a5 = log(2.00)     # Coefficient for x5
  a6 = log(2.50)     # Coefficient for x6
  a7 = log(1.10)     # Interaction term coefficient for x1 * x2
  
  delta_z = 0.6  # Scaling factor for treatment assignment effects
  
  # Calculate the propensity score formula using the logistic model
  formula = a0 + delta_z * (
    a1 * data$x1 + a2 * data$x2 + a3 * data$x3 + a4 * data$x4 + 
      a5 * data$x5 + a6 * data$x6 + a7 * data$x1 * data$x2
  )
  
  # Calculate propensity scores
  prob_treat = expit(formula)  # Treatment probabilities for each individual
  
  # Assign propensity scores to the dataset
  data$prob_treat = prob_treat
  
  # Simulate treatment assignment based on propensity scores
  data$z = rbinom(N, 1, prob_treat)  # Treatment indicator: z ~ Bernoulli(prob_treat)
  
  ############################################################################################################################
  # VISUALIZE PROPENSITY SCORE DISTRIBUTION
  ############################################################################################################################
  
  # Load required package for visualization
  library(ggplot2)
  
  # Convert treatment indicator into a factor for labeling in the plot
  data$TreatmentGroup <- factor(data$z, 
                                levels = c(0, 1), 
                                labels = c("Control (z=0)", "Treatment (z=1)"))
  
  # Create a histogram of propensity scores by treatment group
  p_good_overlap <- ggplot(data, aes(x = prob_treat, fill = TreatmentGroup)) +
    geom_histogram(alpha = 0.6, position = "identity", bins = 50, color = "black") +
    scale_fill_manual(values = c("Control (z=0)" = "pink", "Treatment (z=1)" = "purple")) +
    labs(x = "Propensity Score", y = "Frequency", fill = "Treatment Group") +
    theme_minimal(base_size = 16) +
    theme(
      # Customize axis titles and text for readability
      axis.title = element_text(face = "bold", size = 18),
      axis.text = element_text(size = 16),
      # Customize legend for better interpretability
      legend.title = element_text(face = "bold", size = 16),
      legend.text = element_text(size = 14),
      # Optionally, adjust plot title formatting if used
      plot.title = element_text(face = "bold", size = 20)
    )
  
  # Display the plot
  print(p_good_overlap)
  
  ############################################################################################################################
  # OUTCOME MODEL CONTINUOUS (Y)
  ############################################################################################################################
  
  # Define coefficients for the outcome model
  b0  = 0.00   # Intercept
  b1  = 2.50   # Coefficient for x1
  b2  = -2.00  # Coefficient for x2
  b3  = 1.75   # Coefficient for x3
  b4  = -1.25  # Coefficient for x4
  b5  = 1.50   # Coefficient for x5
  b6  = 1.10   # Coefficient for x6
  b70 = 2.50   # Interaction term coefficient for x1 * x2 in y0
  b71 = 1.50   # Interaction term coefficient for x1 * x2 in y1
  
  # Define treatment effect scaling factors
  delta0 = 0.3  # Scaling factor for control group (y0)
  delta1 = 1    # Constant treatment effect
  delta2 = 0.2  # Interaction scaling factor for treatment effect
  
  # Generate counterfactual outcomes for the control group
  data$y0 = b0 + delta0 * (
    b1 * data$x1 + b2 * data$x2 + b3 * data$x3 + 
      b4 * data$x4 + b5 * data$x5 + b6 * data$x6 + 
      b70 * data$x1 * data$x2
  ) + rnorm(N)  # Random noise with default sd = 1
  
  # Generate counterfactual outcomes for the treatment group
  data$y1 = data$y0 + delta1 + delta2 * (
    b1 * data$x1 + b2 * data$x2 + b3 * data$x3 + 
      b4 * data$x4 + b5 * data$x5 + b6 * data$x6 + 
      b71 * data$x1 * data$x2
  ) + rnorm(N)  # Random noise with default sd = 1
  
  # Calculate observed outcome based on treatment assignment
  data$y = data$z * data$y1 + (1 - data$z) * data$y0  # Observed outcome for each subject
  
  
  data_sets[[scenario]] = data

}





SuperPopGen_GoodOverlap_AllinOne_0118_2025<- data_sets

save(SuperPopGen_GoodOverlap_AllinOne_0118_2025, file = "SuperPopGen_GoodOverlap_AllinOne_0118_2025.RData")


