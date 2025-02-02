# Load the pre-generated superpopulation data
load("SuperPopGen_GoodOverlap_AllinOne_0118_2025.RData")

n.reps = 5000  # Number of simulation replicates

# Set up parallel processing with 10 cores
cl <- parallel::makeCluster(10)  # Create a cluster with 10 cores
registerDoParallel(cl)           # Register the cluster for parallel execution

# Initialize a list to store results for each scenario
results <- vector("list", n_scenarios * s.m)

# Loop over all scenarios
for (scenario in 1:(n_scenarios * s.m)) {
  
  cat("Scenario", scenario, "\n")  # Log the current scenario
  
  # Extract the dataset for the current scenario
  data = data_sets[[scenario]]
  data$ones = 1                 # Indicator variable for sampling
  data$id = 1:nrow(data)        # Assign unique IDs to each subject
  ########################################################################################################################
  # Calculate Population-Level Parameters
  ########################################################################################################################
  
  # Extract propensity scores
  es_sp = data$prob_treat
  
  # Calculate superpopulation-level treatment effects
  ate_sp = mean(data$y1) - mean(data$y0)  # ATE
  att_sp = mean(data$y1[data$z == 1]) - mean(data$y0[data$z == 1])  # ATT
  ato_sp = sum(es_sp * (1 - es_sp) * data$y1) / sum(es_sp * (1 - es_sp)) - sum(es_sp * (1 - es_sp) * data$y0) / sum(es_sp * (1 - es_sp))  # ATO
  
  # Store superpopulation treatment effects
  results[[scenario]]$ate_sp <- ate_sp
  results[[scenario]]$att_sp <- att_sp
  results[[scenario]]$ato_sp <- ato_sp
  
  # ########################################################################################################################
  # # Define Multi-Stage Sampling Parameters
  # ########################################################################################################################
  # 
  # size1 = rep(100000, 10)  # Population size per stratum
  # size2 = rep(5, 10)       # Number of clusters per stratum
  # size3 = c(rep(850 / 5, 5), rep(750 / 5, 5), rep(700 / 5, 5), rep(650 / 5, 5), 
  #           rep(600 / 5, 5), rep(400 / 5, 5), rep(350 / 5, 5), rep(300 / 5, 5), 
  #           rep(250 / 5, 5), rep(150 / 5, 5))  # Cluster-level sample sizes
  
  ########################################################################################################################
  # Parallel Sampling and Estimation
  ########################################################################################################################
  
  # Perform sampling and estimation in parallel
  ests.list <- foreach(i.rep = 1:n.reps, .packages = c("survey", "MatchIt", "sampling", "Hmisc")) %dopar% {
    print(paste(scenario, "-", i.rep))
    
    ####################################################################################################################
    # Create the Sample 
    ####################################################################################################################
    
    # # Multistage sampling (stratified by strata, cluster, and individual levels)
    # s = mstage(data, stage = list("stratified", "cluster", ""), 
    #            varnames = list("Strata", "Cluster", "ones"), 
    #            size = list(size1, size2, size3), 
    #            method = list("", "srswor", "srswor"))
    # 
    # sample = getdata(data, s)[[3]]  # Extract the final sample
    # save(sample, i.rep, file = paste(scenario, "-", i.rep, ".RData"))  # Save the sample
    
    # Coefficients for sample selection model
    c0 = log(0.005 / (1 - 0.005)) 
    c1 = log(1.05)  
    c2 = log(1.10)  
    c3 = log(1.15)  
    c4 = log(1.10) 
    c5 = log(1.05) 
    c6 = log(1.10) 
    delta_z_in_s = log(0.9) 
    
    formula_S = c0 + delta_z_in_s * data$z + c1 * data$x1 + c2 * data$x2 +
      c3 * data$x3 + c4 * data$x4 + c5 * data$x5 + c6 * data$x6
    
    prob_S = expit(formula_S)
    data$Prob = prob_S
    data$S = rbinom(N, 1, prob_S)
    sample = data[data$S == 1, ]

    save(sample,i.rep, file = paste(scenario,"-",i.rep,".RData"))
    
    
    ####################################################################################################################
    # Calculate Sample-Level Treatment Effects
    ####################################################################################################################
    
    att_fp = mean(sample$y1[sample$z == 1]) - mean(sample$y0[sample$z == 1])  # ATT in sample
    ate_fp = mean(sample$y1) - mean(sample$y0)  # ATE in sample
    es_fp = sample$prob_treat  # Propensity scores in sample
    ato_fp = mean(es_fp * (1 - es_fp) * sample$y1) / mean(es_fp * (1 - es_fp)) - 
      mean(es_fp * (1 - es_fp) * sample$y0) / mean(es_fp * (1 - es_fp))  # ATO in sample
    
    # Add survey weights
    sample$s.wt = 1 / sample$Prob  # Survey weights
    
    ####################################################################################################################
    # Perform Estimation Under Different Model Specifications
    ####################################################################################################################
    
    # Define covariate specifications
    covariatesOM_Correct = c("x1", "x2", "x3", "x4", "x5", "x6", "x1:x2")  # Correct outcome model
    covariatesPS_Correct = c("x1", "x2", "x3", "x4", "x5", "x6", "x1:x2")  # Correct propensity score model
    covariatesOM_MisInter = c("x1", "x2", "x3", "x4", "x5", "x6")          # Misspecified outcome model
    covariatesPS_MisInter = c("x1", "x2", "x3", "x4", "x5", "x6")          # Misspecified propensity score model
    covariatesSMD = c("x1", "x2", "x3", "x4", "x5", "x6")
    
    # Generate results for all combinations of model correctness
    
    # Cor|Cor
    results_scenario_BothCorrect=c(
      EstSOD(sample, "s.wt", "z", "y", "U", "W", covariatesPS_Correct, covariatesOM_Correct, covariatesSMD, min_prob, "Cluster", "Strata"),
      EstSOD(sample, "s.wt", "z", "y", "W", "W", covariatesPS_Correct, covariatesOM_Correct, covariatesSMD, min_prob, "Cluster", "Strata"),
      EstSOD(sample, "s.wt", "z", "y", "C", "W", covariatesPS_Correct, covariatesOM_Correct, covariatesSMD, min_prob, "Cluster", "Strata"),
      EstSOD(sample, "s.wt", "z", "y", "CW", "W", covariatesPS_Correct, covariatesOM_Correct, covariatesSMD, min_prob, "Cluster", "Strata")
      
    )
    
    
    # Mis|Cor
    results_scenario_MisPSInter=c(
      EstSOD(sample, "s.wt", "z", "y", "U", "W", covariatesPS_MisInter, covariatesOM_Correct, covariatesSMD, min_prob, "Cluster", "Strata"),
      EstSOD(sample, "s.wt", "z", "y", "W", "W", covariatesPS_MisInter, covariatesOM_Correct, covariatesSMD, min_prob, "Cluster", "Strata"),
      EstSOD(sample, "s.wt", "z", "y", "C", "W", covariatesPS_MisInter, covariatesOM_Correct, covariatesSMD, min_prob, "Cluster", "Strata"),
      EstSOD(sample, "s.wt", "z", "y", "CW", "W", covariatesPS_MisInter, covariatesOM_Correct, covariatesSMD, min_prob, "Cluster", "Strata")
      
    )
    
    
    # Cor|Mis
    results_scenario_MisOMInter=c(
      EstSOD(sample, "s.wt", "z", "y", "U", "W", covariatesPS_Correct, covariatesOM_MisInter , covariatesSMD, min_prob, "Cluster", "Strata"),
      EstSOD(sample, "s.wt", "z", "y", "W", "W", covariatesPS_Correct, covariatesOM_MisInter, covariatesSMD, min_prob, "Cluster", "Strata"),
      EstSOD(sample, "s.wt", "z", "y", "C", "W", covariatesPS_Correct, covariatesOM_MisInter, covariatesSMD, min_prob, "Cluster", "Strata"),
      EstSOD(sample, "s.wt", "z", "y", "CW", "W", covariatesPS_Correct, covariatesOM_MisInter, covariatesSMD, min_prob, "Cluster", "Strata")
    )
    
    
    # Mis|Mis
    results_scenario_MisBoth=c(
      EstSOD(sample, "s.wt", "z", "y", "U", "W", covariatesPS_MisInter, covariatesOM_MisInter, covariatesSMD, min_prob, "Cluster", "Strata"),
      EstSOD(sample, "s.wt", "z", "y", "W", "W", covariatesPS_MisInter, covariatesOM_MisInter, covariatesSMD, min_prob, "Cluster", "Strata"),
      EstSOD(sample, "s.wt", "z", "y", "C", "W", covariatesPS_MisInter, covariatesOM_MisInter, covariatesSMD, min_prob, "Cluster", "Strata"),
      EstSOD(sample, "s.wt", "z", "y", "CW", "W", covariatesPS_MisInter, covariatesOM_MisInter, covariatesSMD, min_prob, "Cluster", "Strata")
      
    )
    
    # Return results for the current replicate
    return(list(
      BothCorrect = results_scenario_BothCorrect,
      MisPSInter = results_scenario_MisPSInter,
      MisOMInter = results_scenario_MisOMInter,
      MisBoth = results_scenario_MisBoth
    ))
  }
  
  # Aggregate results across replicates for each scenario
  for (estimate in names(ests.list[[1]])) {
    ests.estimate <- lapply(ests.list, `[[`, estimate)
    results[[scenario]][[paste0("ests_", estimate)]] <- do.call(rbind, ests.estimate)
  }
}

# Save the complete results
results_5000_Success_GoodOverlap_AllinOne_0118_2025 <- results
save(results_5000_Success_GoodOverlap_AllinOne_0118_2025, file = "results_5000_Success_GoodOverlap_AllinOne_0118_2025.RData")

# Stop the parallel processing cluster
stopCluster(cl)
