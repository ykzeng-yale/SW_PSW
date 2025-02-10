# data=sample
# sample.wt="s.wt"
# treatment="z"
# outcome="y"
# method_PS="W"
# SW_Design = "Independent"
# 
# covariatesPS = c("x1","x2","x3","x4","x5","x6","x1:x2")
# covariatesOM = c("x1","x2","x3","x4","x5","x6","x1:x2")
# covariatesSMD = c("x1","x2","x3","x4","x5","x6")
# 
# covariatesPS=covariatesPS
# covariatesOM=covariatesOM
# covariatesSMD=covariatesSMD
# 
# cluster_name = "Cluster"
# strata_name  = "Strata"
# 
# min_prob = min_prob
# 
# 


##############################################################################
# Updated EstSOD (Survey Observational Data) Using Shared PSW/MOM/CVR/WET Estimator & Sandwich Functions
##############################################################################

EstSOD <- function(
    data, 
    sample.wt    = "s.wt",    # name of sampling weight column
    treatment    = "z",
    outcome      = "y",
    method_PS,                # e.g. "U","W","C","CW"
    SW_Design, # e.g. "Retrospective", "Independent", "Prospective"
    covariatesPS,
    covariatesOM,
    covariatesSMD,
    min_prob,
    cluster_name = "Cluster",
    strata_name  = "Strata"
) {
  ###################################################################
  # 1) Handle cluster/strata columns and define formulas
  ###################################################################
  
  if (! (cluster_name %in% colnames(data)) ) {
    data[[cluster_name]] <- 1
    cluster_formula <- ~1
  } else {
    cluster_formula <- stats::as.formula(paste0("~", cluster_name))
  }
  
  if (! (strata_name %in% colnames(data)) ) {
    data[[strata_name]] <- 1
    strata_formula <- NULL
  } else {
    strata_formula <- stats::as.formula(paste0("~", strata_name))
  }
  
  # For naming
  name <- paste0(method_PS, ".PS|")
  
  ###################################################################
  # 2) Fit the logistic PS model 
  #    (mimicking your original approach but using method_PS)
  ###################################################################
  
  formula_ps_main <- as.formula(
    paste(treatment, paste(covariatesPS, collapse=" + "), sep=" ~ ")
  )
  
  if (method_PS == "U") {
    # Unweighted logistic
    design_ps <- survey::svydesign(ids=~1, weights=~1, data=data)
    weights_ps <- 1
    mod_ps <- survey::svyglm(formula_ps_main, design=design_ps, family=quasibinomial())
    
  } else if (method_PS == "W") {
    # Weighted logistic using sample.wt
    design_ps <- survey::svydesign(
      ids     = cluster_formula,
      strata  = strata_formula,
      weights = ~s.wt, 
      data    = data
    )
    weights_ps <- data$s.wt
    mod_ps <- survey::svyglm(formula_ps_main, design=design_ps, family=quasibinomial())
    
  } else if (method_PS == "C") {
    # Unweighted logistic but "s.wt" included as a covariate
    formula_ps_c <- as.formula(
      paste(treatment, paste(c(covariatesPS, "s.wt"), collapse=" + "), sep=" ~ ")
    )
    design_ps <- survey::svydesign(ids=~1, weights=~1, data=data)
    weights_ps <- 1
    mod_ps <- survey::svyglm(formula_ps_c, design=design_ps, family=quasibinomial())
    
  } else if (method_PS == "CW") {
    # Weighted logistic + "s.wt" as a covariate
    formula_ps_c <- as.formula(
      paste(treatment, paste(c(covariatesPS, "s.wt"), collapse=" + "), sep=" ~ ")
    )
    design_ps <- survey::svydesign(
      ids     = cluster_formula,
      strata  = strata_formula,
      weights = ~s.wt,
      data    = data
    )
    weights_ps <- data$s.wt
    mod_ps <- survey::svyglm(formula_ps_c, design=design_ps, family=quasibinomial())
    
  } else {
    stop("method_PS must be one of 'U','W','C','CW'.")
  }
  
  # Extract esp logistic PS coefs
  coef_ps <- as.numeric(coef(mod_ps))
  data$ps.w <- predict(mod_ps, newdata=data, type="response")
  
  # Extract efp logistic PS coefs
  formula_ps_fp <- as.formula(
    paste(treatment, paste(covariatesPS, collapse=" + "), sep=" ~ ")
  )
  design_ps_fp <- survey::svydesign(ids=~1, weights=~1, data=data)
  mod_ps_fp    <- survey::svyglm(formula_ps_fp, design=design_ps_fp, family=quasibinomial())
  coef_ps_fp   <- as.numeric(coef(mod_ps_fp))
  
  data$ps.w_fp <- predict(mod_ps_fp, newdata=data, type="response")

  min_prob <- min_prob
  data$ps.w    <- pmax(pmin(data$ps.w,    1 - min_prob), min_prob)
  data$ps.w_fp <- pmax(pmin(data$ps.w_fp, 1 - min_prob), min_prob)
  
  covM    <- model.matrix(mod_ps)
  covM_fp <- model.matrix(mod_ps_fp)
  
  esp    <- data$ps.w
  efp    <- data$ps.w_fp
  
  data$z <- data[[treatment]]
  r_z <- ifelse(data$z==1, esp/efp, (1-esp)/(1-efp))
  
  data$esp <- esp
  data$efp <- efp
  data$r_z <- r_z
  
  
  if (SW_Design == "Retrospective") { # Survey weight: P(S=1 | Z, X)
    
    data$ate.wt <- ifelse(data$z==1, 1/esp, 1/(1-esp)) * data$s.wt
    data$att.wt <- ifelse(data$z==1, 1, esp/(1-esp))   * data$s.wt
    data$ato.wt <- ifelse(data$z==1, 1-esp, esp)       * data$s.wt
    
    data$h_ate <- data$s.wt * (1/r_z)
    data$h_att <- esp*data$s.wt * (1/r_z)
    data$h_ato <- esp*(1-esp)*data$s.wt * (1/r_z)
    
  } else if (SW_Design == "Independent") { # Survey weight: P(S=1 | X), Z and S are independent, esp = efp
    
    data$ate.wt <- ifelse(data$z==1, 1/efp, 1/(1-efp)) * data$s.wt
    data$att.wt <- ifelse(data$z==1, 1, efp/(1-efp))   * data$s.wt
    data$ato.wt <- ifelse(data$z==1, 1-efp, efp)       * data$s.wt
    
    data$h_ate <- data$s.wt 
    data$h_att <- efp*data$s.wt
    data$h_ato <- efp*(1-efp)*data$s.wt
    
  } else if (SW_Design == "Prospective") { # Survey weight: P(S=1 | X), Z still depends on S.
    
    data$ate.wt <- ifelse(data$z==1, 1/efp, 1/(1-efp)) * data$s.wt
    data$att.wt <- ifelse(data$z==1, esp/efp, esp/(1-efp))   * data$s.wt
    data$ato.wt <- ifelse(data$z==1, esp*(1-esp)/efp, esp*(1-esp)/(1-efp))   * data$s.wt
    
    data$h_ate <- data$s.wt 
    data$h_att <- esp*data$s.wt
    data$h_ato <- esp*(1-esp)*data$s.wt
    
  } else {
    stop("SW_Design must be 'Retrospective', 'Independent' or 'Prospective'.")
  }
  
  
  ###################################################################
  # 4) Now call the shared estimator/variance functions
  ###################################################################
  
  ###########################
  # 4A) PSW for ATE
  ###########################
  psw_ate <- PSW_Estimator(data, target="ATE")
  mu1_pate <- psw_ate$mu1
  mu0_pate <- psw_ate$mu0
  pate     <- psw_ate$est_psw
  
  # Sandwich
  sv_ate <- PSW_SandwichVariance(
    data     = data,
    target     = "ATE",
    SW_Design  = SW_Design,
    mu1_hat    = mu1_pate,
    mu0_hat    = mu0_pate,
    est_psw    = pate,
    true_value = ate_sp,  
    beta_ps    = coef_ps,
    beta_ps_fp = coef_ps_fp,
    covM       = covM,
    covM_fp    = covM_fp,
    weights_ps = weights_ps,
    min_prob   = min_prob,
    n          = nrow(data)
  )
  ate_san_var      <- sv_ate$var_psw
  ate_san_sd       <- sv_ate$sd_psw
  ate_san_coverage <- sv_ate$coverage_psw
  
  ###########################
  # 4A) PSW for ATT
  ###########################
  psw_att <- PSW_Estimator(data, target="ATT")
  mu1_patt <- psw_att$mu1
  mu0_patt <- psw_att$mu0
  patt     <- psw_att$est_psw
  
  sv_att <- PSW_SandwichVariance(
    data     = data,
    target     = "ATT",
    SW_Design  = SW_Design,
    mu1_hat    = mu1_patt,
    mu0_hat    = mu0_patt,
    est_psw    = patt,
    true_value = att_sp,
    beta_ps    = coef_ps,
    beta_ps_fp = coef_ps_fp,
    covM       = covM,
    covM_fp    = covM_fp,
    weights_ps = weights_ps,
    min_prob   = min_prob,
    n          = nrow(data)
  )
  att_san_var      <- sv_att$var_psw
  att_san_sd       <- sv_att$sd_psw
  att_san_coverage <- sv_att$coverage_psw
  
  ###########################
  # 4A) PSW for ATO
  ###########################
  psw_ato <- PSW_Estimator(data, target="ATO")
  mu1_pato <- psw_ato$mu1
  mu0_pato <- psw_ato$mu0
  pato     <- psw_ato$est_psw
  
  sv_ato <- PSW_SandwichVariance(
    data     = data,
    target     = "ATO",
    SW_Design  = SW_Design,
    mu1_hat    = mu1_pato,
    mu0_hat    = mu0_pato,
    est_psw    = pato,
    true_value = ato_sp,
    beta_ps    = coef_ps,
    beta_ps_fp = coef_ps_fp,
    covM       = covM,
    covM_fp    = covM_fp,
    weights_ps = weights_ps,
    min_prob   = min_prob,
    n          = nrow(data)
  )
  ato_san_var      <- sv_ato$var_psw
  ato_san_sd       <- sv_ato$sd_psw
  ato_san_coverage <- sv_ato$coverage_psw
  
  ###################################################################
  # 4B) MOM for ATE, ATT, ATO
  ###################################################################
  # ATE (MOM)
  res_ate_mom <- MOM_Estimator(
    data          = data,
    target          = "ATE",
    cluster_formula = cluster_formula,
    strata_formula  = strata_formula,
    covariatesOM    = covariatesOM
  )

  data$pred_1_mom <- res_ate_mom$pred_1_mom
  data$pred_0_mom <- res_ate_mom$pred_0_mom
  
  mu1_ate_mom <- res_ate_mom$mu1_mom
  mu0_ate_mom <- res_ate_mom$mu0_mom
  pate_mom    <- res_ate_mom$est_mom
  
  # Sandwich for MOM (ATE)
  sv_ate_mom <- MOM_SandwichVariance(
    data     = data,
    target     = "ATE",
    SW_Design  = SW_Design,
    mu1_hat    = mu1_pate, 
    mu0_hat    = mu0_pate,
    est_mom    = pate_mom,
    true_value = ate_sp,
    fit_1_mom  = res_ate_mom$fit_1_mom,
    fit_0_mom  = res_ate_mom$fit_0_mom,
    beta_ps    = coef_ps,
    beta_ps_fp = coef_ps_fp,
    covM       = covM,
    covM_fp    = covM_fp,
    weights_ps = weights_ps,
    min_prob   = min_prob,
    n          = nrow(data)
  )
  ate_san_var_mom      <- sv_ate_mom$var_mom
  ate_san_sd_mom       <- sv_ate_mom$sd_mom
  ate_san_coverage_mom <- sv_ate_mom$coverage_mom

  # Do similarly for ATT (MOM)
  res_att_mom <- MOM_Estimator(data, target="ATT",
                               cluster_formula=cluster_formula,
                               strata_formula =strata_formula,
                               covariatesOM   =covariatesOM
  )
  data$pred_1_mom <- res_att_mom$pred_1_mom
  data$pred_0_mom <- res_att_mom$pred_0_mom
  mu1_att_mom <- res_att_mom$mu1_mom
  mu0_att_mom <- res_att_mom$mu0_mom
  patt_mom    <- res_att_mom$est_mom
  
  sv_att_mom <- MOM_SandwichVariance(
    data     = data,
    target     = "ATT",
    SW_Design  = SW_Design,
    mu1_hat    = mu1_patt,
    mu0_hat    = mu0_patt,
    est_mom    = patt_mom,
    true_value = att_sp,
    fit_1_mom  = res_att_mom$fit_1_mom,
    fit_0_mom  = res_att_mom$fit_0_mom,
    beta_ps    = coef_ps,
    beta_ps_fp = coef_ps_fp,
    covM       = covM,
    covM_fp    = covM_fp,
    weights_ps = weights_ps,
    min_prob   = min_prob,
    n          = nrow(data)
  )
  att_san_var_mom      <- sv_att_mom$var_mom
  att_san_sd_mom       <- sv_att_mom$sd_mom
  att_san_coverage_mom <- sv_att_mom$coverage_mom
  
  # And ATO (MOM)
  res_ato_mom <- MOM_Estimator(data, target="ATO",
                               cluster_formula=cluster_formula,
                               strata_formula =strata_formula,
                               covariatesOM   =covariatesOM
  )
  data$pred_1_mom <- res_ato_mom$pred_1_mom
  data$pred_0_mom <- res_ato_mom$pred_0_mom
  mu1_ato_mom <- res_ato_mom$mu1_mom
  mu0_ato_mom <- res_ato_mom$mu0_mom
  pato_mom    <- res_ato_mom$est_mom
  
  sv_ato_mom <- MOM_SandwichVariance(
    data     = data,
    target     = "ATO",
    SW_Design  = SW_Design,
    mu1_hat    = mu1_pato,  
    mu0_hat    = mu0_pato,
    est_mom    = pato_mom,
    true_value = ato_sp,
    fit_1_mom  = res_ato_mom$fit_1_mom,
    fit_0_mom  = res_ato_mom$fit_0_mom,
    beta_ps    = coef_ps,
    beta_ps_fp = coef_ps_fp,
    covM       = covM,
    covM_fp    = covM_fp,
    weights_ps = weights_ps,
    min_prob   = min_prob,
    n          = nrow(data)
  )
  ato_san_var_mom      <- sv_ato_mom$var_mom
  ato_san_sd_mom       <- sv_ato_mom$sd_mom
  ato_san_coverage_mom <- sv_ato_mom$coverage_mom
  
  ###################################################################
  # 4C) CVR for ATE, ATT, ATO
  ###################################################################
  
  # ATE (CVR)
  res_ate_cvr <- CVR_Estimator(
    data          = data,
    target          = "ATE",
    cluster_formula = cluster_formula,
    strata_formula  = strata_formula,
    covariatesOM    = covariatesOM
  )
  data$pred_1_cvr <- res_ate_cvr$pred_1_cvr
  data$pred_0_cvr <- res_ate_cvr$pred_0_cvr
  mu1_pate_cvr <- res_ate_cvr$mu1
  mu0_pate_cvr <- res_ate_cvr$mu0
  pate_cvr  <- res_ate_cvr$est_cvr
  
  sv_ate_cvr <- CVR_SandwichVariance(
    data     = data,
    target     = "ATE",
    SW_Design  = SW_Design,
    mu1_hat    = mu1_pate, 
    mu0_hat    = mu0_pate,
    est_cvr    = pate_cvr,
    true_value = ate_sp,
    fit_1_cvr  = res_ate_cvr$fit_1_cvr,
    fit_0_cvr  = res_ate_cvr$fit_0_cvr,
    beta_ps    = coef_ps,
    beta_ps_fp = coef_ps_fp,
    covM       = covM,
    covM_fp    = covM_fp,
    weights_ps = weights_ps,
    min_prob   = min_prob,
    n          = nrow(data)
  )
  ate_san_var_cvr      <- sv_ate_cvr$var_cvr
  ate_san_sd_cvr       <- sv_ate_cvr$sd_cvr
  ate_san_coverage_cvr <- sv_ate_cvr$coverage_cvr
  
  # ATT (CVR)
  res_att_cvr <- CVR_Estimator(
    data          = data,
    target          = "ATT",
    cluster_formula = cluster_formula,
    strata_formula  = strata_formula,
    covariatesOM    = covariatesOM
  )
  data$pred_1_cvr <- res_att_cvr$pred_1_cvr
  data$pred_0_cvr <- res_att_cvr$pred_0_cvr
  mu1_patt_cvr <- res_att_cvr$mu1
  mu0_patt_cvr <- res_att_cvr$mu0
  patt_cvr  <- res_att_cvr$est_cvr
  
  sv_att_cvr <- CVR_SandwichVariance(
    data     = data,
    target     = "ATT",
    SW_Design  = SW_Design,
    mu1_hat    = mu1_patt,  
    mu0_hat    = mu0_patt,
    est_cvr    = patt_cvr,
    true_value = att_sp,
    fit_1_cvr  = res_att_cvr$fit_1_cvr,
    fit_0_cvr  = res_att_cvr$fit_0_cvr,
    beta_ps    = coef_ps,
    beta_ps_fp = coef_ps_fp,
    covM       = covM,
    covM_fp    = covM_fp,
    weights_ps = weights_ps,
    min_prob   = min_prob,
    n          = nrow(data)
  )
  att_san_var_cvr      <- sv_att_cvr$var_cvr
  att_san_sd_cvr       <- sv_att_cvr$sd_cvr
  att_san_coverage_cvr <- sv_att_cvr$coverage_cvr
  
  # ATO (CVR)
  res_ato_cvr <- CVR_Estimator(
    data          = data,
    target          = "ATO",
    cluster_formula = cluster_formula,
    strata_formula  = strata_formula,
    covariatesOM    = covariatesOM
  )
  data$pred_1_cvr <- res_ato_cvr$pred_1_cvr
  data$pred_0_cvr <- res_ato_cvr$pred_0_cvr
  mu1_pato_cvr <- res_ato_cvr$mu1
  mu0_pato_cvr <- res_ato_cvr$mu0
  pato_cvr  <- res_ato_cvr$est_cvr
  
  sv_ato_cvr <- CVR_SandwichVariance(
    data     = data,
    target     = "ATO",
    SW_Design  = SW_Design,
    mu1_hat    = mu1_pato,  
    mu0_hat    = mu0_pato,
    est_cvr    = pato_cvr,
    true_value = ato_sp,
    fit_1_cvr  = res_ato_cvr$fit_1_cvr,
    fit_0_cvr  = res_ato_cvr$fit_0_cvr,
    beta_ps    = coef_ps,
    beta_ps_fp = coef_ps_fp,
    covM       = covM,
    covM_fp    = covM_fp,
    weights_ps = weights_ps,
    min_prob   = min_prob,
    n          = nrow(data)
  )
  ato_san_var_cvr      <- sv_ato_cvr$var_cvr
  ato_san_sd_cvr       <- sv_ato_cvr$sd_cvr
  ato_san_coverage_cvr <- sv_ato_cvr$coverage_cvr
  
  
  ###################################################################
  # 4D) WET for ATE, ATT, ATO
  ###################################################################
  
  # ATE (WET)
  res_ate_wet <- WET_Estimator(
    data          = data,
    target          = "ATE",
    covariatesOM    = covariatesOM,
    cluster_formula = cluster_formula,
    strata_formula  = strata_formula
  )
  data$pred_1_wet <- res_ate_wet$pred_1_wet
  data$pred_0_wet <- res_ate_wet$pred_0_wet
  mu1_pate_wet <- res_ate_wet$mu1
  mu0_pate_wet <- res_ate_wet$mu0
  pate_wet<- res_ate_wet$est_wet
  
  sv_ate_wet <- WET_SandwichVariance(
    data     = data,
    target     = "ATE",
    SW_Design  = SW_Design,
    mu1_hat    = mu1_pate,  
    mu0_hat    = mu0_pate,
    est_wet    = pate_wet,
    true_value = ate_sp,
    fit_1_wet  = res_ate_wet$fit_1_wet,
    fit_0_wet  = res_ate_wet$fit_0_wet,
    beta_ps    = coef_ps,
    beta_ps_fp = coef_ps_fp,
    covM       = covM,
    covM_fp    = covM_fp,
    weights_ps = weights_ps,
    min_prob   = min_prob,
    n          = nrow(data)
  )
  ate_san_var_wet      <- sv_ate_wet$var_wet
  ate_san_sd_wet       <- sv_ate_wet$sd_wet
  ate_san_coverage_wet <- sv_ate_wet$coverage_wet
  
  # ATT (WET)
  res_att_wet <- WET_Estimator(
    data          = data,
    target          = "ATT",
    covariatesOM    = covariatesOM,
    cluster_formula = cluster_formula,
    strata_formula  = strata_formula
  )
  data$pred_1_wet <- res_att_wet$pred_1_wet
  data$pred_0_wet <- res_att_wet$pred_0_wet
  mu1_patt_wet <- res_att_wet$mu1
  mu0_patt_wet <- res_att_wet$mu0
  patt_wet<- res_att_wet$est_wet
  
  sv_att_wet <- WET_SandwichVariance(
    data     = data,
    target     = "ATT",
    SW_Design  = SW_Design,
    mu1_hat    = mu1_patt,
    mu0_hat    = mu0_patt,
    est_wet    = patt_wet,
    true_value = att_sp,
    fit_1_wet  = res_att_wet$fit_1_wet,
    fit_0_wet  = res_att_wet$fit_0_wet,
    beta_ps    = coef_ps,
    beta_ps_fp = coef_ps_fp,
    covM       = covM,
    covM_fp    = covM_fp,
    weights_ps = weights_ps,
    min_prob   = min_prob,
    n          = nrow(data)
  )
  att_san_var_wet      <- sv_att_wet$var_wet
  att_san_sd_wet       <- sv_att_wet$sd_wet
  att_san_coverage_wet <- sv_att_wet$coverage_wet
  
  # ATO (WET)
  res_ato_wet <- WET_Estimator(
    data          = data,
    target          = "ATO",
    covariatesOM    = covariatesOM,
    cluster_formula = cluster_formula,
    strata_formula  = strata_formula
  )
  data$pred_1_wet <- res_ato_wet$pred_1_wet
  data$pred_0_wet <- res_ato_wet$pred_0_wet
  mu1_pato_wet <- res_ato_wet$mu1
  mu0_pato_wet <- res_ato_wet$mu0
  pato_wet<- res_ato_wet$est_wet
  
  sv_ato_wet <- WET_SandwichVariance(
    data     = data,
    target     = "ATO",
    SW_Design  = SW_Design,
    mu1_hat    = mu1_pato,
    mu0_hat    = mu0_pato,
    est_wet    = pato_wet,
    true_value = ato_sp,
    fit_1_wet  = res_ato_wet$fit_1_wet,
    fit_0_wet  = res_ato_wet$fit_0_wet,
    beta_ps    = coef_ps,
    beta_ps_fp = coef_ps_fp,
    covM       = covM,
    covM_fp    = covM_fp,
    weights_ps = weights_ps,
    min_prob   = min_prob,
    n          = nrow(data)
  )
  ato_san_var_wet      <- sv_ato_wet$var_wet
  ato_san_sd_wet       <- sv_ato_wet$sd_wet
  ato_san_coverage_wet <- sv_ato_wet$coverage_wet
  
  ###################################################################
  # 5) Compute biases 
  ###################################################################

  bias.ate <- pate - ate_sp
  bias.att <- patt - att_sp
  bias.ato <- pato - ato_sp
  names(bias.ate) <- "bias.ate"
  names(bias.att) <- "bias.att"
  names(bias.ato) <- "bias.ato"
  
  bias.ate_mom <- pate_mom - ate_sp
  bias.att_mom <- patt_mom - att_sp
  bias.ato_mom <- pato_mom - ato_sp
  
  names(bias.ate_mom) <- "bias.ate_mom"
  names(bias.att_mom) <- "bias.att_mom"
  names(bias.ato_mom) <- "bias.ato_mom"
  
  bias.ate_cvr <- pate_cvr - ate_sp
  bias.att_cvr <- patt_cvr - att_sp
  bias.ato_cvr <- pato_cvr - ato_sp
  
  names(bias.ate_cvr) <- "bias.ate_cvr"
  names(bias.att_cvr) <- "bias.att_cvr"
  names(bias.ato_cvr) <- "bias.ato_cvr"
  
  bias.ate_wet <- pate_wet - ate_sp
  bias.att_wet <- patt_wet - att_sp
  bias.ato_wet <- pato_wet - ato_sp
  
  names(bias.ate_wet) <- "bias.ate_wet"
  names(bias.att_wet) <- "bias.att_wet"
  names(bias.ato_wet) <- "bias.ato_wet"
  
  ###################################################################
  # 6) Build final named result vector
  ###################################################################
  
  name_result <- c(
    paste("PATE", name, sep="|"), 
    paste("PATT", name, sep="|"), 
    paste("PATO", name, sep="|"),
    
    paste("Bias|PATE", name, sep="|"), 
    paste("Bias|PATT", name, sep="|"), 
    paste("Bias|PATO", name, sep="|"),
    
    paste("SanVar|PATE", name, sep="|"), 
    paste("SanVar|PATT", name, sep="|" ),
    paste("SanVar|PATO", name, sep="|" ),
    
    paste("SanCov|PATE", name, sep="|"), 
    paste("SanCov|PATT", name, sep="|" ),
    paste("SanCov|PATO", name, sep="|" ),
    
    paste("PATE_MOM", name, sep="|"), 
    paste("PATT_MOM", name, sep="|"), 
    paste("PATO_MOM", name, sep="|"),
    
    paste("Bias|PATE_MOM", name, sep="|"), 
    paste("Bias|PATT_MOM", name, sep="|"), 
    paste("Bias|PATO_MOM", name, sep="|"),
    
    paste("SanVar|PATE_MOM", name, sep="|"), 
    paste("SanVar|PATT_MOM", name, sep="|"),
    paste("SanVar|PATO_MOM", name, sep="|"),
    
    paste("SanCov|PATE_MOM", name, sep="|"),
    paste("SanCov|PATT_MOM", name, sep="|"),
    paste("SanCov|PATO_MOM", name, sep="|"),
    
    paste("PATE_CVR", name, sep="|"), 
    paste("PATT_CVR", name, sep="|"),
    paste("PATO_CVR", name, sep="|"),
    
    paste("Bias|PATE_CVR", name, sep="|"),
    paste("Bias|PATT_CVR", name, sep="|"),
    paste("Bias|PATO_CVR", name, sep="|"),
    
    paste("SanVar|PATE_CVR", name, sep="|"),
    paste("SanVar|PATT_CVR", name, sep="|"),
    paste("SanVar|PATO_CVR", name, sep="|"),
    
    paste("SanCov|PATE_CVR", name, sep="|"),
    paste("SanCov|PATT_CVR", name, sep="|"),
    paste("SanCov|PATO_CVR", name, sep="|"),
    
    paste("PATE_WET", name, sep="|"),
    paste("PATT_WET", name, sep="|"),
    paste("PATO_WET", name, sep="|"),
    
    paste("Bias|PATE_WET", name, sep="|"),
    paste("Bias|PATT_WET", name, sep="|"),
    paste("Bias|PATO_WET", name, sep="|"),
    
    paste("SanVar|PATE_WET", name, sep="|"),
    paste("SanVar|PATT_WET", name, sep="|"),
    paste("SanVar|PATO_WET", name, sep="|"),
    
    paste("SanCov|PATE_WET", name, sep="|"),
    paste("SanCov|PATT_WET", name, sep="|"),
    paste("SanCov|PATO_WET", name, sep="|")
  )
  
  # Build the final results, in the same order:
  part_of_result <- c(
    # PSW
    pate, patt, pato,
    bias.ate, bias.att, bias.ato,
    ate_san_var, att_san_var, ato_san_var,
    ate_san_coverage, att_san_coverage, ato_san_coverage,
    
    # MOM
    pate_mom, patt_mom, pato_mom,
    bias.ate_mom, bias.att_mom, bias.ato_mom,
    ate_san_var_mom, att_san_var_mom, ato_san_var_mom,
    ate_san_coverage_mom, att_san_coverage_mom, ato_san_coverage_mom,
    
    # CVR
    pate_cvr, patt_cvr, pato_cvr,
    bias.ate_cvr, bias.att_cvr, bias.ato_cvr,
    ate_san_var_cvr, att_san_var_cvr, ato_san_var_cvr,
    ate_san_coverage_cvr, att_san_coverage_cvr, ato_san_coverage_cvr,
    
    # WET
    pate_wet, patt_wet, pato_wet,
    bias.ate_wet, bias.att_wet, bias.ato_wet,
    ate_san_var_wet, att_san_var_wet, ato_san_var_wet,
    ate_san_coverage_wet, att_san_coverage_wet, ato_san_coverage_wet
  )
  
  names(part_of_result) <- name_result
  

  # Add SMD 
  smd_varlist <- covariatesSMD
  
  smd_ate_sw <- SMD_Calculator(
    data   = data,
    varlist  = smd_varlist,
    target   = "ATE",
    name     = name
  )
  smd_att_sw <- SMD_Calculator(
    data   = data,
    varlist  = smd_varlist,
    target   = "ATT",
    name     = name
  )
  smd_ato_sw <- SMD_Calculator(
    data   = data,
    varlist  = smd_varlist,
    target   = "ATO",
    name     = name
  )
  
  ###################################################################
  # 7) Return
  ###################################################################
  
  results_sample <- c(
    part_of_result,
    smd_ate_sw,
    smd_att_sw,
    smd_ato_sw
  )
  
  return(results_sample)
}
