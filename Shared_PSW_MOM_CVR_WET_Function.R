
##############################################################################
# Propensity Score Weighting (PSW) Estimator
##############################################################################

PSW_Estimator <- function(data, 
                          target = c("ATE","ATT","ATO")) {
  
  target <- match.arg(target)
  
  # Extract outcome and treatment from data
  y <- data$y
  z <- data$z
  
  # Compute weighted means for mu1 and mu0 depending on target
  if (target == "ATE") {
    mu1 <- sum(z * y * data$ate.wt) / sum(z * data$ate.wt)
    mu0 <- sum((1 - z) * y * data$ate.wt) / sum((1 - z) * data$ate.wt)
    
  } else if (target == "ATT") {
    mu1 <- sum(z * y * data$att.wt) / sum(z * data$att.wt)
    mu0 <- sum((1 - z) * y * data$att.wt) / sum((1 - z) * data$att.wt)
    
  } else if (target == "ATO") {
    mu1 <- sum(z * y * data$ato.wt) / sum(z * data$ato.wt)
    mu0 <- sum((1 - z) * y * data$ato.wt) / sum((1 - z) * data$ato.wt)
    
  } else {
    warning("PSW_Estimator only supports 'ATE', 'ATT', or 'ATO'. Returning NA.")
    return(list(mu1=NA, mu0=NA, est_psw=NA))
  }
  
  est_psw <- mu1 - mu0
  
  list(
    mu1     = mu1,
    mu0     = mu0,
    est_psw = est_psw
  )
}



##############################################################################
# Sandwich Variance Estimator for PSW
##############################################################################

PSW_SandwichVariance <- function(
    data, 
    target       = c("ATE","ATT","ATO"), 
    mu1_hat,             # The estimated mu1 (from PSW_Estimator)
    mu0_hat,             # The estimated mu0 (from PSW_Estimator)
    est_psw,             # The PSW estimate = mu1_hat - mu0_hat
    true_value,          # e.g. ate_sp, att_sp, or ato_sp (for coverage)
    # Additional arguments needed for the sandwich variance:
    beta_ps,             # the esp PS model coefs (weighted or unweighted)
    beta_ps_fp,          # the efp PS model coefs
    covM,                # model.matrix(...) for esp
    covM_fp,             # model.matrix(...) for efp
    weights_ps,          # e.g. data$s.wt or 1
    min_prob = 0.00,            # min_prob to keep esp from 0 or 1
    n                    # total sample size in data
) {

  target <- match.arg(target)
  
  y <- data$y
  z <- data$z
  
  mu1 <- mu1_hat
  mu0 <- mu0_hat
  
  W    <- covM
  p    <- ncol(W)
  W_fp <- covM_fp
  p_fp <- ncol(W_fp)
  
  phi <- function(theta) {
    local_mu1  <- theta[1]
    local_mu0  <- theta[2]
    local_beta <- theta[3:(p+2)]
    
    esp <- plogis(drop(W %*% local_beta))
    esp <- pmax(pmin(esp, 1 - min_prob), min_prob)
    
    efp <- plogis(drop(W_fp %*% beta_ps_fp))
    efp <- pmax(pmin(efp, 1 - min_prob), min_prob)
    
    # Define the score components for the chosen estimand
    if (target == "ATE") {
      f1 <- (y - local_mu1) * z       * (1/esp)       * data$s.wt
      f2 <- (y - local_mu0) * (1 - z) * (1/(1 - esp)) * data$s.wt
      
    } else if (target == "ATT") {
      f1 <- (y - local_mu1) * z       * 1             * data$s.wt
      f2 <- (y - local_mu0) * (1 - z) * (esp/(1-esp)) * data$s.wt
      
    } else if (target == "ATO") {
      f1 <- (y - local_mu1) * z       * (1 - esp) * data$s.wt
      f2 <- (y - local_mu0) * (1 - z) * esp       * data$s.wt
      
    } else {
      warning("PSW_SandwichVariance only supports 'ATE', 'ATT', or 'ATO'. Returning NA.")
      # Return a zero matrix
      out_dim <- p + 2
      return(matrix(0, nrow=out_dim, ncol=nrow(data)))
    }
    
    # f3: esp logistic PS score
    f3 <- weights_ps * W * (z - esp)
    
    rbind(f1, f2, t(f3))
  }
  
  mphi <- function(theta) {
    rowMeans( phi(theta) )
  }
  
  # Vector of parameters: (mu1, mu0, beta_ps)
  theta_hat <- c(mu1_hat, mu0_hat, beta_ps)
  
  # Derivative wrt all parameters
  Atheta <- numDeriv::jacobian(mphi, theta_hat)
  invAtheta <- tryCatch(
    solve(Atheta),
    error = function(e) {
      message("Matrix is singular, using pseudo-inverse via MASS::ginv")
      MASS::ginv(Atheta)
    }
  )
  
  # Meat operator => crossprod of phi
  phis_hat <- phi(theta_hat)            
  B_mat    <- tcrossprod(phis_hat) / n  
  
  # Final sandwich variance
  Var_mat <- invAtheta %*% B_mat %*% t(invAtheta) / n
  
  a <- matrix(0, nrow=2, ncol=length(theta_hat))
  a[1,1] <- 1
  a[2,2] <- 1
  covmu <- a %*% Var_mat %*% t(a)
  
  contrast <- matrix(c(1, -1), nrow=1)
  var_psw  <- drop(contrast %*% covmu %*% t(contrast))
  sd_psw   <- sqrt(var_psw)
  
  coverage_psw <- as.numeric(
    (est_psw - 1.96*sd_psw < true_value) &&
      (est_psw + 1.96*sd_psw > true_value)
  )
  
  list(
    var_psw       = var_psw,
    sd_psw        = sd_psw,
    coverage_psw  = coverage_psw
  )
}


##############################################################################
# Moment (MOM) Estimator
##############################################################################

MOM_Estimator <- function(
    data,
    target = c("ATE", "ATT", "ATO"),  
    cluster_formula = ~1,            # Default if no clusters
    strata_formula  = NULL,          # Default if no strata
    covariatesOM    = NULL
) {
  target <- match.arg(target)
  
  # (a) Build a survey design with weights = ~1 
  #     because we fit outcome models unweighted in the MOM approach
  design_OM_mom <- survey::svydesign(
    ids     = cluster_formula,
    strata  = strata_formula,
    weights = ~1,
    data    = data
  )
  
  # (b) Subset the design to z=1 and z=0
  design_OM_mom_1 <- subset(design_OM_mom, z == 1)
  design_OM_mom_0 <- subset(design_OM_mom, z == 0)
  
  # (c) Build formula for outcome models
  #     e.g. "y ~ x1 + x2 + x3"
  formula_OM_1 <- stats::as.formula(
    paste("y", paste(covariatesOM, collapse=" + "), sep=" ~ ")
  )
  formula_OM_0 <- stats::as.formula(
    paste("y", paste(covariatesOM, collapse=" + "), sep=" ~ ")
  )
  
  # (d) Fit these outcome models
  fit_1_mom <- survey::svyglm(formula_OM_1, design = design_OM_mom_1)
  fit_0_mom <- survey::svyglm(formula_OM_0, design = design_OM_mom_0)
  
  # (e) Predict potential outcomes for ALL rows of the data
  pred_1_mom <- stats::predict(fit_1_mom, newdata = data, type = "response")
  pred_0_mom <- stats::predict(fit_0_mom, newdata = data, type = "response")
  
  # (f) Compute the final MOM estimate: depends on target ∈ {ATE, ATT, ATO}
  y <- data$y
  z <- data$z
  
  if (target == "ATE") {
    mu1_mom <-
      sum(data$h_ate * pred_1_mom) / sum(data$h_ate) +
      sum(z * data$ate.wt * (y - pred_1_mom)) / sum(z * data$ate.wt)
    
    mu0_mom <-
      sum(data$h_ate * pred_0_mom) / sum(data$h_ate) +
      sum((1 - z) * data$ate.wt * (y - pred_0_mom)) / sum((1 - z) * data$ate.wt)
    
  } else if (target == "ATT") {
    mu1_mom <-
      sum(data$h_att * pred_1_mom) / sum(data$h_att) +
      sum(z * data$att.wt * (y - pred_1_mom)) / sum(z * data$att.wt)
    
    mu0_mom <-
      sum(data$h_att * pred_0_mom) / sum(data$h_att) +
      sum((1 - z) * data$att.wt * (y - pred_0_mom)) / sum((1 - z) * data$att.wt)
    
  } else if (target == "ATO") {
    mu1_mom <-
      sum(data$h_ato * pred_1_mom) / sum(data$h_ato) +
      sum(z * data$ato.wt * (y - pred_1_mom)) / sum(z * data$ato.wt)
    
    mu0_mom <-
      sum(data$h_ato * pred_0_mom) / sum(data$h_ato) +
      sum((1 - z) * data$ato.wt * (y - pred_0_mom)) / sum((1 - z) * data$ato.wt)
    
  } else {
    warning("MOM_Estimator: only supports ATE, ATT, or ATO. Returning NA.")
    return(list(mu1_mom=NA, mu0_mom=NA, est_mom=NA,
                fit_1_mom=NULL, fit_0_mom=NULL))
  }
  
  est_mom <- mu1_mom - mu0_mom
  
  # Return everything needed for subsequent usage
  list(
    mu1_mom   = mu1_mom,      
    mu0_mom   = mu0_mom,      
    est_mom   = est_mom,      
    
    # The fitted outcome models
    fit_1_mom = fit_1_mom,
    fit_0_mom = fit_0_mom,
    
    # The predicted potential outcomes
    pred_1_mom = pred_1_mom,
    pred_0_mom = pred_0_mom
  )
}


##############################################################################
# Clever Covariate (CVR) Estimator
##############################################################################

CVR_Estimator <- function(
    data,
    target = c("ATE","ATT","ATO"),
    cluster_formula = ~1,
    strata_formula  = NULL,
    covariatesOM    = NULL
) {
  target <- match.arg(target)
  
  # 1) Build a survey design with weights=~1 for outcome regression
  design_OM_CVR <- survey::svydesign(
    ids     = cluster_formula,
    strata  = strata_formula,
    weights = ~1,
    data    = data
  )
  
  # 2) Subset by z == 1 and z == 0
  design_OM_CVR_1 <- subset(design_OM_CVR, z == 1)
  design_OM_CVR_0 <- subset(design_OM_CVR, z == 0)
  
  # 3) Depending on the target, add "ate.wt", "att.wt", or "ato.wt" as a predictor
  if (target == "ATE") {
    formula_1_cvr <- as.formula(
      paste("y ~", paste(c(covariatesOM, "ate.wt"), collapse=" + "))
    )
    formula_0_cvr <- formula_1_cvr
    
  } else if (target == "ATT") {
    formula_1_cvr <- as.formula(
      paste("y ~", paste(c(covariatesOM, "att.wt"), collapse=" + "))
    )
    formula_0_cvr <- formula_1_cvr
    
  } else if (target == "ATO") {
    formula_1_cvr <- as.formula(
      paste("y ~", paste(c(covariatesOM, "ato.wt"), collapse=" + "))
    )
    formula_0_cvr <- formula_1_cvr
    
  } else {
    warning("CVR_Estimator only supports ATE, ATT, ATO. Returning NA.")
    return(list(mu1=NA, mu0=NA, est_cvr=NA, fit_1_cvr=NULL, fit_0_cvr=NULL))
  }
  
  # 4) Fit the two outcome models
  fit_1_cvr <- survey::svyglm(formula_1_cvr, design=design_OM_CVR_1)
  fit_0_cvr <- survey::svyglm(formula_0_cvr, design=design_OM_CVR_0)
  
  # 5) Predict potential outcomes for everyone in data
  pred_1_cvr <- predict(fit_1_cvr, newdata=data, type="response")
  pred_0_cvr <- predict(fit_0_cvr, newdata=data, type="response")
  
  # 6) Compute mu1_cvr and mu0_cvr
  if (target == "ATE") {
    mu1_cvr <- sum(data$h_ate * pred_1_cvr) / sum(data$h_ate)
    mu0_cvr <- sum(data$h_ate * pred_0_cvr) / sum(data$h_ate)
  } else if (target == "ATT") {
    mu1_cvr <- sum(data$h_att * pred_1_cvr) / sum(data$h_att)
    mu0_cvr <- sum(data$h_att * pred_0_cvr) / sum(data$h_att)
  } else {  # ATO
    mu1_cvr <- sum(data$h_ato * pred_1_cvr) / sum(data$h_ato)
    mu0_cvr <- sum(data$h_ato * pred_0_cvr) / sum(data$h_ato)
  }
  
  est_cvr <- mu1_cvr - mu0_cvr
  
  return(list(
    mu1        = mu1_cvr,
    mu0        = mu0_cvr,
    est_cvr    = est_cvr,
    fit_1_cvr  = fit_1_cvr,
    fit_0_cvr  = fit_0_cvr,
    pred_1_cvr = pred_1_cvr,
    pred_0_cvr = pred_0_cvr
  ))
}


##############################################################################
# Weighted Regression (WET) Estimator
##############################################################################

WET_Estimator <- function(
    data,
    target          = c("ATE", "ATT", "ATO"),
    cluster_formula = ~1,          
    strata_formula  = NULL,
    covariatesOM    = NULL
) {
  target <- match.arg(target)
  
  formula_within_group_WET <- as.formula(
    paste("y ~", paste(covariatesOM, collapse=" + "))
  )
  
  # 2) Choose which weight to use in the survey design
  #    - ATE => data$ate.wt
  #    - ATT => data$att.wt
  #    - ATO => data$ato.wt
  if (target == "ATE") {
    design_OM_WET <- survey::svydesign(
      ids     = cluster_formula,
      strata  = strata_formula,
      weights = ~data$ate.wt,
      data    = data
    )
    
    design_OM_WET_1 <- subset(design_OM_WET, z == 1)
    design_OM_WET_0 <- subset(design_OM_WET, z == 0)
    
    fit_1_wet <- survey::svyglm(formula_within_group_WET, design = design_OM_WET_1)
    fit_0_wet <- survey::svyglm(formula_within_group_WET, design = design_OM_WET_0)
    
    pred_1_wet <- predict(fit_1_wet, newdata=data, type="response")
    pred_0_wet <- predict(fit_0_wet, newdata=data, type="response")
    
    # Combine via tilting function h_ate
    mu1_wet <- sum(data$h_ate * pred_1_wet) / sum(data$h_ate)
    mu0_wet <- sum(data$h_ate * pred_0_wet) / sum(data$h_ate)
    
  } else if (target == "ATT") {
    design_OM_WET <- survey::svydesign(
      ids     = cluster_formula,
      strata  = strata_formula,
      weights = ~data$att.wt,
      data    = data
    )
    
    design_OM_WET_1 <- subset(design_OM_WET, z == 1)
    design_OM_WET_0 <- subset(design_OM_WET, z == 0)
    
    fit_1_wet <- survey::svyglm(formula_within_group_WET, design = design_OM_WET_1)
    fit_0_wet <- survey::svyglm(formula_within_group_WET, design = design_OM_WET_0)
    
    pred_1_wet <- predict(fit_1_wet, newdata=data, type="response")
    pred_0_wet <- predict(fit_0_wet, newdata=data, type="response")
    
    mu1_wet <- sum(data$h_att * pred_1_wet) / sum(data$h_att)
    mu0_wet <- sum(data$h_att * pred_0_wet) / sum(data$h_att)
    
  } else { # target == "ATO"
    design_OM_WET <- survey::svydesign(
      ids     = cluster_formula,
      strata  = strata_formula,
      weights = ~data$ato.wt,
      data    = data
    )
    
    design_OM_WET_1 <- subset(design_OM_WET, z == 1)
    design_OM_WET_0 <- subset(design_OM_WET, z == 0)
    
    fit_1_wet <- survey::svyglm(formula_within_group_WET, design = design_OM_WET_1)
    fit_0_wet <- survey::svyglm(formula_within_group_WET, design = design_OM_WET_0)
    
    pred_1_wet <- predict(fit_1_wet, newdata=data, type="response")
    pred_0_wet <- predict(fit_0_wet, newdata=data, type="response")
    
    mu1_wet <- sum(data$h_ato * pred_1_wet) / sum(data$h_ato)
    mu0_wet <- sum(data$h_ato * pred_0_wet) / sum(data$h_ato)
  }
  
  est_wet <- mu1_wet - mu0_wet
  
  list(
    mu1         = mu1_wet,
    mu0         = mu0_wet,
    est_wet     = est_wet,
    fit_1_wet   = fit_1_wet,   
    fit_0_wet   = fit_0_wet,
    pred_1_wet  = pred_1_wet,
    pred_0_wet  = pred_0_wet
  )
}


##############################################################################
# Sandwich Variance Estimator for MOM
##############################################################################

MOM_SandwichVariance <- function(
    data,
    target        = c("ATE", "ATT", "ATO"),
    
    # The PSW-based mu1, mu0 estimates => e.g. mu1_pate, mu0_pate
    mu1_hat,            
    mu0_hat,            
    
    # The final MOM difference
    est_mom,            
    
    # For coverage calculation
    true_value,         # e.g. ate_sp, att_sp, or ato_sp
    
    # The fitted outcome models
    fit_1_mom,          # e.g. the svyglm fit for z=1
    fit_0_mom,          # e.g. the svyglm fit for z=0
    
    # PS arguments
    beta_ps,            
    beta_ps_fp,         
    covM,               
    covM_fp,            
    weights_ps,         
    min_prob   = 0.00,  
    n          = nrow(data)
) {

  target <- match.arg(target)
  
  y <- data$y
  z <- data$z

  # Extract outcome-model coefs (gamma1, gamma0)
  gamma1.h <- as.numeric(coef(fit_1_mom))
  gamma0.h <- as.numeric(coef(fit_0_mom))
  
  # Build the same design matrix used in fit_1_mom and fit_0_mom
  XY <- model.matrix(stats::formula(fit_1_mom), data = data)
  
  # Dimensions
  W    <- covM
  p    <- ncol(W)
  W_fp <- covM_fp
  p_fp <- ncol(W_fp)
  q    <- ncol(XY)    
  
  # Re-extract the predicted potential outcomes from data
  m1.h <- data$pred_1_mom
  m0.h <- data$pred_0_mom
  
  # Compute the "aug1h.h" and "aug1z.h" as well as "aug0h.h" and "aug0z.h"
  if (target == "ATE") {
    aug1h.h <- sum(data$h_ate * m1.h) / sum(data$h_ate)
    aug1z.h <- sum(z * m1.h * data$ate.wt) / sum(z * data$ate.wt)
    
    aug0h.h <- sum(data$h_ate * m0.h) / sum(data$h_ate)
    aug0z.h <- sum((1 - z) * m0.h * data$ate.wt) / sum((1 - z) * data$ate.wt)
    
  } else if (target == "ATT") {
    aug1h.h <- sum(data$h_att * m1.h) / sum(data$h_att)
    aug1z.h <- sum(z * m1.h * data$att.wt) / sum(z * data$att.wt)
    
    aug0h.h <- sum(data$h_att * m0.h) / sum(data$h_att)
    aug0z.h <- sum((1 - z) * m0.h * data$att.wt) / sum((1 - z) * data$att.wt)
    
  } else if (target == "ATO") {
    aug1h.h <- sum(data$h_ato * m1.h) / sum(data$h_ato)
    aug1z.h <- sum(z * m1.h * data$ato.wt) / sum(z * data$ato.wt)
    
    aug0h.h <- sum(data$h_ato * m0.h) / sum(data$h_ato)
    aug0z.h <- sum((1 - z) * m0.h * data$ato.wt) / sum((1 - z) * data$ato.wt)
    
  } else {
    warning("MOM_SandwichVariance only supports ATE, ATT, ATO. Returning NA.")
    return(list(var_mom=NA, sd_mom=NA, coverage_mom=NA))
  }
  
  # Define v1, v2, v3
  v1 <- aug1h.h - aug0h.h
  v2 <- mu1_hat - aug1z.h
  v3 <- mu0_hat - aug0z.h
  
  # Define phi(...) for the augmented approach
  phi <- function(theta) {
    loc_v1      <- theta[1]
    loc_v2      <- theta[2]
    loc_v3      <- theta[3]
    loc_beta    <- theta[4:(3 + p)]
    loc_beta_fp <- theta[(4 + p):(3 + p + p_fp)]
    loc_gamma1  <- theta[(4 + p + p_fp):(3 + p + p_fp + q)]
    loc_gamma0  <- theta[(4 + p + p_fp + q):(3 + p + p_fp + 2*q)]
    
    esp <- plogis(drop(W    %*% loc_beta))
    esp <- pmax(pmin(esp, 1 - min_prob), min_prob)
    
    efp <- plogis(drop(W_fp %*% loc_beta_fp))
    efp <- pmax(pmin(efp, 1 - min_prob), min_prob)
    
    r_z <- ifelse(z == 1, esp/efp, (1 - esp)/(1 - efp))
    
    m1_local <- drop(XY %*% loc_gamma1)
    m0_local <- drop(XY %*% loc_gamma0)
    
    # f1, f2, f3 differ by target
    if (target == "ATE") {
      f1 <- data$s.wt*(1/r_z)*(m1_local - m0_local - loc_v1)
      f2 <- z * (1/esp)*data$s.wt * (y - m1_local - loc_v2)
      f3 <- (1 - z)*(1/(1-esp))*data$s.wt * (y - m0_local - loc_v3)
    } else if (target == "ATT") {
      f1 <- esp*data$s.wt*(1/r_z)*(m1_local - m0_local - loc_v1)
      f2 <- z * 1 * data$s.wt * (y - m1_local - loc_v2)
      f3 <- (1 - z)*(esp/(1-esp))*data$s.wt*(y - m0_local - loc_v3)
    } else { # "ATO"
      f1 <- esp*(1-esp)*data$s.wt*(1/r_z)*(m1_local - m0_local - loc_v1)
      f2 <- z*(1-esp)*data$s.wt*(y - m1_local - loc_v2)
      f3 <- (1 - z)*esp*data$s.wt*(y - m0_local - loc_v3)
    }
    
    # f7, f8 => outcome-model score wrt gamma1, gamma0
    f7 <- XY * ((y - m1_local) * z)
    f8 <- XY * ((y - m0_local) * (1 - z))
    
    # f9, f10 => logistic PS components
    f9  <- weights_ps * W    * (z - esp)
    f10 <-               W_fp * (z - efp)
    
    rbind(
      f1,
      f2,
      f3,
      t(f7),
      t(f8),
      t(f9),
      t(f10)
    )
  }
  
  # Helper: rowMeans of phi
  mphi <- function(theta) rowMeans(phi(theta))
  
  #  param vector = (v1, v2, v3, beta_ps, beta_ps_fp, gamma1, gamma0)
  theta_hat <- c(
    v1, v2, v3,
    beta_ps,
    beta_ps_fp,
    gamma1.h,
    gamma0.h
  )

  # Derivative (Jacobian) and Meat (crossprod of phi)
  Atheta <- numDeriv::jacobian(mphi, theta_hat)
  invAtheta <- tryCatch(
    solve(Atheta),
    error = function(e) {
      message("MOM_SandwichVariance: Atheta singular => using MASS::ginv")
      MASS::ginv(Atheta)
    }
  )
  
  phis_hat <- phi(theta_hat)           
  B_mat    <- tcrossprod(phis_hat) / n  
  
  Var_mat <- invAtheta %*% B_mat %*% t(invAtheta) / n
  
  dims_total <- length(theta_hat)
  a <- matrix(0, nrow=3, ncol=dims_total)
  a[1,1] <- 1  # picks out v1
  a[2,2] <- 1  # picks out v2
  a[3,3] <- 1  # picks out v3
  
  cov_v123 <- a %*% Var_mat %*% t(a)  # 3x3
  
  # (v1 + v2 - v3)
  contrast <- matrix(c(1,1,-1), nrow=1) 
  var_mom  <- drop(contrast %*% cov_v123 %*% t(contrast))
  sd_mom   <- sqrt(var_mom)
  
  coverage_mom <- as.numeric(
    (est_mom - 1.96*sd_mom < true_value) &&
      (est_mom + 1.96*sd_mom > true_value)
  )
  
  list(
    var_mom       = var_mom,
    sd_mom        = sd_mom,
    coverage_mom  = coverage_mom
  )
}


##############################################################################
# Sandwich Variance Estimator for CVR
##############################################################################

CVR_SandwichVariance <- function(
    data,
    target = c("ATE","ATT","ATO"),
    
    # The PSW-based mu1, mu0 estimates => e.g. mu1_pate, mu0_pate
    mu1_hat,
    mu0_hat,
    
    # The final CVR difference
    est_cvr,
    
    # Optional: True value for coverage
    true_value,       # e.g. ate_sp, att_sp, or ato_sp
    
    # Fitted outcome models for "clever" approach
    fit_1_cvr,
    fit_0_cvr,
    
    # PS arguments
    beta_ps,
    beta_ps_fp,
    covM,
    covM_fp,
    weights_ps,
    
    min_prob = 0.00,
    n        = nrow(data)
) {

  target <- match.arg(target)

  y <- data$y
  z <- data$z
  
  # 2) Extract the outcome-model coefs from fit_1_cvr, fit_0_cvr
  gamma1.h <- as.numeric(coef(fit_1_cvr))
  gamma0.h <- as.numeric(coef(fit_0_cvr))
  
  # Build the design matrix XY that matches fit_1_cvr
  XY <- model.matrix(stats::formula(fit_1_cvr), data=data)

  W    <- covM
  p    <- ncol(W)
  W_fp <- covM_fp
  p_fp <- ncol(W_fp)
  q    <- ncol(XY)
  
  # Re-extract the predicted potential outcomes
  m1.h <- data$pred_1_cvr  
  m0.h <- data$pred_0_cvr
  
  # 4) Compute the "augmentation" pieces => v1, v2, v3
  if (target == "ATE") {
    clever1h.h <- sum(data$h_ate * m1.h) / sum(data$h_ate)
    clever1z.h <- sum(z * m1.h * data$ate.wt) / sum(z * data$ate.wt)
    clever0h.h <- sum(data$h_ate * m0.h) / sum(data$h_ate)
    clever0z.h <- sum((1 - z) * m0.h * data$ate.wt) / sum((1 - z) * data$ate.wt)
    
  } else if (target == "ATT") {
    clever1h.h <- sum(data$h_att * m1.h) / sum(data$h_att)
    clever1z.h <- sum(z * m1.h * data$att.wt) / sum(z * data$att.wt)
    clever0h.h <- sum(data$h_att * m0.h) / sum(data$h_att)
    clever0z.h <- sum((1 - z) * m0.h * data$att.wt) / sum((1 - z) * data$att.wt)
    
  } else if (target == "ATO") {
    clever1h.h <- sum(data$h_ato * m1.h) / sum(data$h_ato)
    clever1z.h <- sum(z * m1.h * data$ato.wt) / sum(z * data$ato.wt)
    clever0h.h <- sum(data$h_ato * m0.h) / sum(data$h_ato)
    clever0z.h <- sum((1 - z) * m0.h * data$ato.wt) / sum((1 - z) * data$ato.wt)
    
  } else {
    warning("CVR_SandwichVariance only supports ATE, ATT, ATO. Returning NA.")
    return(list(var_cvr=NA, sd_cvr=NA, coverage_cvr=NA))
  }
  
  v1 <- (clever1h.h - clever0h.h)
  v2 <- (mu1_hat - clever1z.h)
  v3 <- (mu0_hat - clever0z.h)
  
  # 5) Build phi(...) function
  phi <- function(theta) {
    loc_v1      <- theta[1]
    loc_v2      <- theta[2]
    loc_v3      <- theta[3]
    loc_beta    <- theta[4:(3 + p)]
    loc_beta_fp <- theta[(4 + p):(3 + p + p_fp)]
    loc_gamma1  <- theta[(4 + p + p_fp)           : (3 + p + p_fp +    q)]
    loc_gamma0  <- theta[(4 + p + p_fp +    q)    : (3 + p + p_fp + 2*q)]
    
    esp <- plogis(drop(W    %*% loc_beta))
    esp <- pmax(pmin(esp, 1 - min_prob), min_prob)
    
    efp <- plogis(drop(W_fp %*% loc_beta_fp))
    efp <- pmax(pmin(efp, 1 - min_prob), min_prob)
    
    r_z <- ifelse(z==1, esp/efp, (1-esp)/(1-efp))
    
    m1_local <- drop(XY %*% loc_gamma1)
    m0_local <- drop(XY %*% loc_gamma0)
    
    # f1, f2, f3 => differ by target
    if (target == "ATE") {
      f1 <- data$s.wt*(1/r_z)*(m1_local - m0_local - loc_v1)
      f2 <- z*(1/esp)*data$s.wt * (y - m1_local - loc_v2)
      f3 <- (1 - z)*(1/(1-esp))*data$s.wt * (y - m0_local - loc_v3)
    } else if (target == "ATT") {
      f1 <- esp*data$s.wt*(1/r_z)*(m1_local - m0_local - loc_v1)
      f2 <- z*1*data$s.wt * (y - m1_local - loc_v2)
      f3 <- (1 - z)*(esp/(1-esp))*data$s.wt * (y - m0_local - loc_v3)
    } else {
      # ATO
      f1 <- esp*(1-esp)*data$s.wt*(1/r_z)*(m1_local - m0_local - loc_v1)
      f2 <- z*(1-esp)*data$s.wt * (y - m1_local - loc_v2)
      f3 <- (1 - z)*esp*data$s.wt * (y - m0_local - loc_v3)
    }
    
    # f7, f8 => outcome-model score
    f7 <- XY * ((y - m1_local)*z)
    f8 <- XY * ((y - m0_local)*(1-z))
    
    # f9 => esp PS score
    f9  <- weights_ps * W    * (z - esp)
    # f10 => efp PS score
    f10 <-              W_fp * (z - efp)
    
    rbind(
      f1,
      f2,
      f3,
      t(f7),
      t(f8),
      t(f9),
      t(f10)
    )
  }
  
  mphi <- function(theta) rowMeans(phi(theta))
  
  theta_hat <- c(
    v1, v2, v3,
    beta_ps,
    beta_ps_fp,
    gamma1.h,
    gamma0.h
  )
  
  Atheta <- numDeriv::jacobian(mphi, theta_hat)
  invAtheta <- tryCatch(
    solve(Atheta),
    error = function(e) {
      message("CVR_SandwichVariance: Atheta singular => using MASS::ginv")
      MASS::ginv(Atheta)
    }
  )
  
  phis_hat <- phi(theta_hat)              
  B_mat    <- tcrossprod(phis_hat) / n    
  
  Var_mat <- invAtheta %*% B_mat %*% t(invAtheta) / n
  
  # 8) Extract var of (v1, v2, v3) => contrast(1,1,-1)
  dims_total <- length(theta_hat)
  a <- matrix(0, nrow=3, ncol=dims_total)
  a[1,1] <- 1
  a[2,2] <- 1
  a[3,3] <- 1
  
  cov_v123 <- a %*% Var_mat %*% t(a)
  
  contrast <- matrix(c(1,1,-1), nrow=1)  
  var_cvr  <- drop(contrast %*% cov_v123 %*% t(contrast))
  sd_cvr   <- sqrt(var_cvr)
  
  coverage_cvr <- as.numeric(
    (est_cvr - 1.96*sd_cvr < true_value) &&
      (est_cvr + 1.96*sd_cvr > true_value)
  )
  
  list(
    var_cvr       = var_cvr,
    sd_cvr        = sd_cvr,
    coverage_cvr  = coverage_cvr
  )
}



##############################################################################
# Sandwich Variance Estimator for WET
##############################################################################

WET_SandwichVariance <- function(
    data,
    target         = c("ATE", "ATT", "ATO"),
    
    # The PSW-based mu1, mu0 estimates => e.g. mu1_pate, mu0_pate
    mu1_hat,              
    mu0_hat,
    
    # The final WET difference
    est_wet,         
    
    # Optional: True value for coverage
    true_value,          # e.g. ate_sp, att_sp, or ato_sp
    
    # Fitted outcome models for Weighted Regression
    fit_1_wet,
    fit_0_wet,
    
    # PS arguments
    beta_ps,
    beta_ps_fp,
    covM,
    covM_fp,
    weights_ps,
    min_prob = 0.00,
    n        = nrow(data)
){
  target <- match.arg(target)
  
  y <- data$y
  z <- data$z
  
  # 1) Extract outcome-model coefs (gamma1, gamma0)
  gamma1.h <- as.numeric(coef(fit_1_wet))
  gamma0.h <- as.numeric(coef(fit_0_wet))
  
  # Build the design matrix used by fit_1_wet
  XY <- model.matrix(stats::formula(fit_1_wet), data=data)
  
  # Dimensions
  W    <- covM
  p    <- ncol(W)
  W_fp <- covM_fp
  p_fp <- ncol(W_fp)
  q    <- ncol(XY)
  
  # 2) Re-extract predicted potential outcomes from WET_Estimator
  m1.h <- data$pred_1_wet
  m0.h <- data$pred_0_wet
  
  # 3) Compute the augmentation pieces => (v1, v2, v3)
  if (target == "ATE") {
    weighted1h.h <- sum(data$h_ate * m1.h) / sum(data$h_ate)
    weighted1z.h <- sum(z * m1.h * data$ate.wt) / sum(z * data$ate.wt)
    
    weighted0h.h <- sum(data$h_ate * m0.h) / sum(data$h_ate)
    weighted0z.h <- sum((1 - z)* m0.h * data$ate.wt) / sum((1 - z)* data$ate.wt)
    
  } else if (target == "ATT") {
    weighted1h.h <- sum(data$h_att * m1.h) / sum(data$h_att)
    weighted1z.h <- sum(z * m1.h * data$att.wt) / sum(z * data$att.wt)
    
    weighted0h.h <- sum(data$h_att * m0.h) / sum(data$h_att)
    weighted0z.h <- sum((1 - z)* m0.h * data$att.wt) / sum((1 - z)* data$att.wt)
    
  } else { # target == "ATO"
    weighted1h.h <- sum(data$h_ato * m1.h) / sum(data$h_ato)
    weighted1z.h <- sum(z * m1.h * data$ato.wt) / sum(z * data$ato.wt)
    
    weighted0h.h <- sum(data$h_ato * m0.h) / sum(data$h_ato)
    weighted0z.h <- sum((1 - z)* m0.h * data$ato.wt) / sum((1 - z)* data$ato.wt)
  }
  
  v1 <- weighted1h.h - weighted0h.h
  v2 <- mu1_hat - weighted1z.h
  v3 <- mu0_hat - weighted0z.h
  
  # 4) Define the phi(...) function for Weighted Regression approach
  phi <- function(theta){
    loc_v1      <- theta[1]
    loc_v2      <- theta[2]
    loc_v3      <- theta[3]
    loc_beta    <- theta[4:(3+p)]
    loc_beta_fp <- theta[(4+p):(3+p+p_fp)]
    loc_gamma1  <- theta[(4+p+p_fp)          : (3+p+p_fp +    q)]
    loc_gamma0  <- theta[(4+p+p_fp +    q)   : (3+p+p_fp + 2*q)]
    
    esp <- plogis(drop(W    %*% loc_beta))
    esp <- pmax(pmin(esp, 1 - min_prob), min_prob)
    
    efp <- plogis(drop(W_fp %*% loc_beta_fp))
    efp <- pmax(pmin(efp, 1 - min_prob), min_prob)
    
    r_z <- ifelse(z==1, esp/efp, (1-esp)/(1-efp))
    
    m1_local <- drop(XY %*% loc_gamma1)
    m0_local <- drop(XY %*% loc_gamma0)

    if (target == "ATE") {
      f1 <- data$s.wt*(1/r_z)*(m1_local - m0_local - loc_v1)
      f2 <- z*(1/esp)*data$s.wt * (y - m1_local - loc_v2)
      f3 <- (1 - z)*(1/(1-esp))*data$s.wt * (y - m0_local - loc_v3)
      
    } else if (target == "ATT") {
      f1 <- esp*data$s.wt*(1/r_z)*(m1_local - m0_local - loc_v1)
      f2 <- z*(1)*data$s.wt * (y - m1_local - loc_v2)
      f3 <- (1 - z)*(esp/(1-esp))*data$s.wt * (y - m0_local - loc_v3)
      
    } else {
      f1 <- esp*(1-esp)*data$s.wt*(1/r_z)*(m1_local - m0_local - loc_v1)
      f2 <- z*(1-esp)*data$s.wt * (y - m1_local - loc_v2)
      f3 <- (1 - z)*esp*data$s.wt * (y - m0_local - loc_v3)
    }
    
    # f7, f8 => from outcome model wrt gamma1, gamma0
    if (target == "ATE") {
      f7 <- XY * ((y - m1_local)*z)        * data$ate.wt
      f8 <- XY * ((y - m0_local)*(1 - z))  * data$ate.wt
    } else if (target == "ATT") {
      f7 <- XY * ((y - m1_local)*z)        * data$att.wt
      f8 <- XY * ((y - m0_local)*(1 - z))  * data$att.wt
    } else { # ATO
      f7 <- XY * ((y - m1_local)*z)        * data$ato.wt
      f8 <- XY * ((y - m0_local)*(1 - z))  * data$ato.wt
    }

    f9  <- weights_ps * W    * (z - esp)
    f10 <- W_fp * (z - efp)
    
    rbind(
      f1,
      f2,
      f3,
      t(f7),
      t(f8),
      t(f9),
      t(f10)
    )
  }
  
  mphi <- function(theta) rowMeans(phi(theta))

  theta_hat <- c(
    v1, v2, v3,
    beta_ps,
    beta_ps_fp,
    gamma1.h,
    gamma0.h
  )
  
  Atheta <- numDeriv::jacobian(mphi, theta_hat)
  invAtheta <- tryCatch(
    solve(Atheta),
    error=function(e) {
      message("WET_SandwichVariance: Atheta singular => using MASS::ginv")
      MASS::ginv(Atheta)
    }
  )
  
  phis_hat <- phi(theta_hat)            
  B_mat    <- tcrossprod(phis_hat) / n  
  
  Var_mat <- invAtheta %*% B_mat %*% t(invAtheta) / n
  
  dims_total <- length(theta_hat)
  a <- matrix(0, nrow=3, ncol=dims_total)
  a[1,1] <- 1 # picks v1
  a[2,2] <- 1 # picks v2
  a[3,3] <- 1 # picks v3
  
  cov_v123 <- a %*% Var_mat %*% t(a)
  
  # The estimand is v1 + v2 - v3 = (1,1,-1)
  contrast <- matrix(c(1,1,-1), nrow=1)
  var_wet  <- drop(contrast %*% cov_v123 %*% t(contrast))
  sd_wet   <- sqrt(var_wet)
  
  coverage_wet <- as.numeric(
    (est_wet - 1.96*sd_wet < true_value) &&
      (est_wet + 1.96*sd_wet > true_value)
  )
  
  list(
    var_wet       = var_wet,
    sd_wet        = sd_wet,
    coverage_wet  = coverage_wet
  )
}





##############################################################################
# Pooled Standard Deviation Calculator for SMD
##############################################################################

calc_pooled_sd <- function(data, variable, wt_name) {
  # Weighted means in each group z=0, z=1
  # Use "tapply" with the user-specified weight column
  w_mean <- tapply(data[[variable]] * data[[wt_name]],
                   data[["z"]], sum) /
    tapply(data[[wt_name]], data[["z"]], sum)
  
  data[["z"]] <- factor(data[["z"]])
  
  # Match each observation to its group's weighted mean
  matched_mean <- w_mean[ as.numeric(data[["z"]]) ]
  
  # Weighted variance for each group
  squared_diff <- (data[[variable]] - matched_mean)^2
  group_var <- tapply(squared_diff * data[[wt_name]],
                      data[["z"]], sum) /
    tapply(data[[wt_name]], data[["z"]], sum)
  
  # Weighted group sizes
  group_sizes <- tapply(data[[wt_name]], data[["z"]], sum)
  
  # Pooled SD
  pooled_sd <- sqrt( sum((group_sizes - 1)*group_var) /
                       sum(group_sizes - 1) )
  return(pooled_sd)
}


##############################################################################
# Standardized Mean Difference (SMD) Calculator
##############################################################################

SMD_Calculator <- function(data, varlist, target = c("ATE","ATT","ATO"), name) {
  target <- match.arg(target)
  
  if (target == "ATE") {
    wt_name <- "ate.wt"
    prefix  <- "SMD_ATE_SW"
  } else if (target == "ATT") {
    wt_name <- "att.wt"
    prefix  <- "SMD_ATT_SW"
  } else {
    wt_name <- "ato.wt"
    prefix  <- "SMD_ATO_SW"
  }
  
  # Store the SMD for each variable in `varlist`
  smd_values <- vector("numeric", length(varlist))
  names_smd  <- character(length(varlist))
  
  # Weighted sums in each group => Weighted means => difference => divide by pooled SD
  wts <- data[[wt_name]]  # The actual numeric weight vector
  
  for (i in seq_along(varlist)) {
    var_i <- varlist[i]
    
    w_mean_diff <- diff(
      tapply(data[[var_i]] * wts, data[["z"]], sum) /
        tapply(wts, data[["z"]], sum)
    )
    
    # Pooled SD
    psd <- calc_pooled_sd(data, var_i, wt_name)
    
    smd_val <- abs(w_mean_diff) / psd
    smd_values[i] <- smd_val
    
    # e.g. "SMD_ATE_SW.x1.<name>" 
    names_smd[i]  <- paste(prefix, var_i, name, sep=".")
  }
  
  names(smd_values) <- names_smd
  return(smd_values)
}





