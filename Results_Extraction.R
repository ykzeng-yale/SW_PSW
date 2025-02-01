##############################################################################
# Results Extract
##############################################################################


load("/Volumes/T7 Shield/R_SW_Final_Data_0118_2025/GoodOverlap_Multi_Stage_Update_Variance/results_5000_Success_GoodOverlap_AllinOne_0118_2025")
results <- results_5000_Success_GoodOverlap_AllinOne_0118_2025

# Cor|Cor

results_sum <- results
results_sum[[1]]$ests<-results_sum[[1]]$ests_BothCorrect

# Mis|Cor

results_sum <- results
results_sum[[1]]$ests<-results_sum[[1]]$ests_MisPSInter

# Cor|Mis

results_sum <- results
results_sum[[1]]$ests<-results_sum[[1]]$ests_MisOMInter

# Mis|Mis

results_sum <- results
results_sum[[1]]$ests<-results_sum[[1]]$ests_MisBoth


###########################################################################
# Relative Bias
###########################################################################


ATE.ests <- c("PATE|U.PS|W.OM", "PATE|W.PS|W.OM", "PATE|C.PS|W.OM", "PATE|CW.PS|W.OM")
ATT.ests <- c("PATT|U.PS|W.OM", "PATT|W.PS|W.OM", "PATT|C.PS|W.OM", "PATT|CW.PS|W.OM")
ATO.ests <- c("PATO|U.PS|W.OM", "PATO|W.PS|W.OM", "PATO|C.PS|W.OM", "PATO|CW.PS|W.OM")

ATE_BIAS <- sapply(results_sum, function(x) {colMeans(x$ests[,ATE.ests] - x$ate_sp)/x$ate_sp}*100)
ATT_BIAS <- sapply(results_sum, function(x) {colMeans(x$ests[,ATT.ests] - x$att_sp)/x$att_sp}*100)
ATO_BIAS <- sapply(results_sum, function(x) {colMeans(x$ests[,ATO.ests] - x$ato_sp)/x$ato_sp}*100)

xtable(rbind(ATE_BIAS,ATT_BIAS,ATO_BIAS))
round(rbind(ATE_BIAS,ATT_BIAS,ATO_BIAS),3)


ATE_MOM.ests <- c("PATE_MOM|U.PS|W.OM", "PATE_MOM|W.PS|W.OM", "PATE_MOM|C.PS|W.OM", "PATE_MOM|CW.PS|W.OM")
ATT_MOM.ests <- c("PATT_MOM|U.PS|W.OM", "PATT_MOM|W.PS|W.OM", "PATT_MOM|C.PS|W.OM", "PATT_MOM|CW.PS|W.OM")
ATO_MOM.ests <- c("PATO_MOM|U.PS|W.OM", "PATO_MOM|W.PS|W.OM", "PATO_MOM|C.PS|W.OM", "PATO_MOM|CW.PS|W.OM")

ATE_MOM_BIAS <- sapply(results_sum, function(x) {colMeans(x$ests[,ATE_MOM.ests] - x$ate_sp)/x$ate_sp}*100)
ATT_MOM_BIAS <- sapply(results_sum, function(x) {colMeans(x$ests[,ATT_MOM.ests] - x$att_sp)/x$att_sp}*100)
ATO_MOM_BIAS <- sapply(results_sum, function(x) {colMeans(x$ests[,ATO_MOM.ests] - x$ato_sp)/x$ato_sp}*100)

xtable(rbind(ATE_MOM_BIAS,ATT_MOM_BIAS,ATO_MOM_BIAS))
round(rbind(ATE_MOM_BIAS,ATT_MOM_BIAS,ATO_MOM_BIAS),3)


ATE_CVR.ests <- c("PATE_CVR|U.PS|W.OM", "PATE_CVR|W.PS|W.OM", "PATE_CVR|C.PS|W.OM", "PATE_CVR|CW.PS|W.OM")
ATT_CVR.ests <- c("PATT_CVR|U.PS|W.OM", "PATT_CVR|W.PS|W.OM", "PATT_CVR|C.PS|W.OM", "PATT_CVR|CW.PS|W.OM")
ATO_CVR.ests <- c("PATO_CVR|U.PS|W.OM", "PATO_CVR|W.PS|W.OM", "PATO_CVR|C.PS|W.OM", "PATO_CVR|CW.PS|W.OM")

ATE_CVR_BIAS <- sapply(results_sum, function(x) {colMeans(x$ests[,ATE_CVR.ests] - x$ate_sp)/x$ate_sp}*100)
ATT_CVR_BIAS <- sapply(results_sum, function(x) {colMeans(x$ests[,ATT_CVR.ests] - x$att_sp)/x$att_sp}*100)
ATO_CVR_BIAS <- sapply(results_sum, function(x) {colMeans(x$ests[,ATO_CVR.ests] - x$ato_sp)/x$ato_sp}*100)


xtable(rbind(ATE_CVR_BIAS,ATT_CVR_BIAS,ATO_CVR_BIAS))
round(rbind(ATE_CVR_BIAS,ATT_CVR_BIAS,ATO_CVR_BIAS),3)


ATE_WET.ests <- c("PATE_WET|U.PS|W.OM", "PATE_WET|W.PS|W.OM", "PATE_WET|C.PS|W.OM", "PATE_WET|CW.PS|W.OM")
ATT_WET.ests <- c("PATT_WET|U.PS|W.OM", "PATT_WET|W.PS|W.OM", "PATT_WET|C.PS|W.OM", "PATT_WET|CW.PS|W.OM")
ATO_WET.ests <- c("PATO_WET|U.PS|W.OM", "PATO_WET|W.PS|W.OM", "PATO_WET|C.PS|W.OM", "PATO_WET|CW.PS|W.OM")

ATE_WET_BIAS <- sapply(results_sum, function(x) {colMeans(x$ests[,ATE_WET.ests] - x$ate_sp)/x$ate_sp}*100)
ATT_WET_BIAS <- sapply(results_sum, function(x) {colMeans(x$ests[,ATT_WET.ests] - x$att_sp)/x$att_sp}*100)
ATO_WET_BIAS <- sapply(results_sum, function(x) {colMeans(x$ests[,ATO_WET.ests] - x$ato_sp)/x$ato_sp}*100)

xtable(rbind(ATE_WET_BIAS,ATT_WET_BIAS,ATO_WET_BIAS))
round(rbind(ATE_WET_BIAS,ATT_WET_BIAS,ATO_WET_BIAS),3)


round(rbind(ATE_BIAS,ATE_MOM_BIAS,ATE_CVR_BIAS,ATE_WET_BIAS),3)
round(rbind(ATT_BIAS,ATT_MOM_BIAS,ATT_CVR_BIAS,ATT_WET_BIAS),3)
round(rbind(ATO_BIAS,ATO_MOM_BIAS,ATO_CVR_BIAS,ATO_WET_BIAS),3)




###########################################################################
# Monte Carlo Variance and Relative Efficiency
###########################################################################

#######################################
# Compute Monte Carlo Variances
#######################################


ATE.ests <- c("PATE|U.PS|W.OM", "PATE|W.PS|W.OM", "PATE|C.PS|W.OM", "PATE|CW.PS|W.OM")
ATT.ests <- c("PATT|U.PS|W.OM", "PATT|W.PS|W.OM", "PATT|C.PS|W.OM", "PATT|CW.PS|W.OM")
ATO.ests <- c("PATO|U.PS|W.OM", "PATO|W.PS|W.OM", "PATO|C.PS|W.OM", "PATO|CW.PS|W.OM")

ATE_MOM.ests <- c("PATE_MOM|U.PS|W.OM", "PATE_MOM|W.PS|W.OM", "PATE_MOM|C.PS|W.OM", "PATE_MOM|CW.PS|W.OM")
ATT_MOM.ests <- c("PATT_MOM|U.PS|W.OM", "PATT_MOM|W.PS|W.OM", "PATT_MOM|C.PS|W.OM", "PATT_MOM|CW.PS|W.OM")
ATO_MOM.ests <- c("PATO_MOM|U.PS|W.OM", "PATO_MOM|W.PS|W.OM", "PATO_MOM|C.PS|W.OM", "PATO_MOM|CW.PS|W.OM")

ATE_CVR.ests <- c("PATE_CVR|U.PS|W.OM", "PATE_CVR|W.PS|W.OM", "PATE_CVR|C.PS|W.OM", "PATE_CVR|CW.PS|W.OM")
ATT_CVR.ests <- c("PATT_CVR|U.PS|W.OM", "PATT_CVR|W.PS|W.OM", "PATT_CVR|C.PS|W.OM", "PATT_CVR|CW.PS|W.OM")
ATO_CVR.ests <- c("PATO_CVR|U.PS|W.OM", "PATO_CVR|W.PS|W.OM", "PATO_CVR|C.PS|W.OM", "PATO_CVR|CW.PS|W.OM")

ATE_WET.ests <- c("PATE_WET|U.PS|W.OM", "PATE_WET|W.PS|W.OM", "PATE_WET|C.PS|W.OM", "PATE_WET|CW.PS|W.OM")
ATT_WET.ests <- c("PATT_WET|U.PS|W.OM", "PATT_WET|W.PS|W.OM", "PATT_WET|C.PS|W.OM", "PATT_WET|CW.PS|W.OM")
ATO_WET.ests <- c("PATO_WET|U.PS|W.OM", "PATO_WET|W.PS|W.OM", "PATO_WET|C.PS|W.OM", "PATO_WET|CW.PS|W.OM")


ATE_MC_VAR <- sapply(results_sum, function(x) { apply(x$ests[, ATE.ests], 2, var) })
ATT_MC_VAR <- sapply(results_sum, function(x) { apply(x$ests[, ATT.ests], 2, var) })
ATO_MC_VAR <- sapply(results_sum, function(x) { apply(x$ests[, ATO.ests], 2, var) })

ATE_MOM_MC_VAR <- sapply(results_sum, function(x) { apply(x$ests[, ATE_MOM.ests], 2, var) })
ATT_MOM_MC_VAR <- sapply(results_sum, function(x) { apply(x$ests[, ATT_MOM.ests], 2, var) })
ATO_MOM_MC_VAR <- sapply(results_sum, function(x) { apply(x$ests[, ATO_MOM.ests], 2, var) })

ATE_CVR_MC_VAR <- sapply(results_sum, function(x) { apply(x$ests[, ATE_CVR.ests], 2, var) })
ATT_CVR_MC_VAR <- sapply(results_sum, function(x) { apply(x$ests[, ATT_CVR.ests], 2, var) })
ATO_CVR_MC_VAR <- sapply(results_sum, function(x) { apply(x$ests[, ATO_CVR.ests], 2, var) })

ATE_WET_MC_VAR <- sapply(results_sum, function(x) { apply(x$ests[, ATE_WET.ests], 2, var) })
ATT_WET_MC_VAR <- sapply(results_sum, function(x) { apply(x$ests[, ATT_WET.ests], 2, var) })
ATO_WET_MC_VAR <- sapply(results_sum, function(x) { apply(x$ests[, ATO_WET.ests], 2, var) })

round(ATE_MC_VAR, 12)
round(ATT_MC_VAR, 12)
round(ATO_MC_VAR, 12)

round(ATE_MOM_MC_VAR, 12)
round(ATT_MOM_MC_VAR, 12)
round(ATO_MOM_MC_VAR, 12)

round(ATE_CVR_MC_VAR, 12)
round(ATT_CVR_MC_VAR, 12)
round(ATO_CVR_MC_VAR, 12)

round(ATE_WET_MC_VAR, 12)
round(ATT_WET_MC_VAR, 12)
round(ATO_WET_MC_VAR, 12)

#######################################
# Relative Efficiency Computation
#######################################

PATE_W.PS_Correct_MC <- 0.005139163 # Baseline MC Variance for PATE|W.PS|W.OM estimator (from PATE|W.PS|W.OM ATE_MC_VAR)

#####################################
# Relative Efficiency based on MC variance
#####################################

Rel_Eff_ATE_MC <- PATE_W.PS_Correct_MC / ATE_MC_VAR
Rel_Eff_ATT_MC <- PATE_W.PS_Correct_MC / ATT_MC_VAR
Rel_Eff_ATO_MC <- PATE_W.PS_Correct_MC / ATO_MC_VAR

Rel_Eff_ATE_MOM_MC <- PATE_W.PS_Correct_MC / ATE_MOM_MC_VAR
Rel_Eff_ATT_MOM_MC <- PATE_W.PS_Correct_MC / ATT_MOM_MC_VAR
Rel_Eff_ATO_MOM_MC <- PATE_W.PS_Correct_MC / ATO_MOM_MC_VAR

Rel_Eff_ATE_CVR_MC <- PATE_W.PS_Correct_MC / ATE_CVR_MC_VAR
Rel_Eff_ATT_CVR_MC <- PATE_W.PS_Correct_MC / ATT_CVR_MC_VAR
Rel_Eff_ATO_CVR_MC <- PATE_W.PS_Correct_MC / ATO_CVR_MC_VAR

Rel_Eff_ATE_WET_MC <- PATE_W.PS_Correct_MC / ATE_WET_MC_VAR
Rel_Eff_ATT_WET_MC <- PATE_W.PS_Correct_MC / ATT_WET_MC_VAR
Rel_Eff_ATO_WET_MC <- PATE_W.PS_Correct_MC / ATO_WET_MC_VAR

xtable(rbind(Rel_Eff_ATE_MC, Rel_Eff_ATT_MC, Rel_Eff_ATO_MC))
round(rbind(Rel_Eff_ATE_MC, Rel_Eff_ATT_MC, Rel_Eff_ATO_MC), 3)

xtable(rbind(Rel_Eff_ATE_MOM_MC, Rel_Eff_ATT_MOM_MC, Rel_Eff_ATO_MOM_MC))
round(rbind(Rel_Eff_ATE_MOM_MC, Rel_Eff_ATT_MOM_MC, Rel_Eff_ATO_MOM_MC), 3)

xtable(rbind(Rel_Eff_ATE_CVR_MC, Rel_Eff_ATT_CVR_MC, Rel_Eff_ATO_CVR_MC))
round(rbind(Rel_Eff_ATE_CVR_MC, Rel_Eff_ATT_CVR_MC, Rel_Eff_ATO_CVR_MC), 3)

xtable(rbind(Rel_Eff_ATE_WET_MC, Rel_Eff_ATT_WET_MC, Rel_Eff_ATO_WET_MC))
round(rbind(Rel_Eff_ATE_WET_MC, Rel_Eff_ATT_WET_MC, Rel_Eff_ATO_WET_MC), 3)


round(rbind(Rel_Eff_ATE_MC, Rel_Eff_ATE_MOM_MC, Rel_Eff_ATE_CVR_MC, Rel_Eff_ATE_WET_MC), 3)
round(rbind(Rel_Eff_ATT_MC, Rel_Eff_ATT_MOM_MC, Rel_Eff_ATT_CVR_MC, Rel_Eff_ATT_WET_MC), 3)
round(rbind(Rel_Eff_ATO_MC, Rel_Eff_ATO_MOM_MC, Rel_Eff_ATO_CVR_MC, Rel_Eff_ATO_WET_MC), 3)





###########################################################################
# Coverage
###########################################################################

San_ATE.cov <- c("SanCov|PATE|U.PS|W.OM", "SanCov|PATE|W.PS|W.OM", "SanCov|PATE|C.PS|W.OM", "SanCov|PATE|CW.PS|W.OM")
San_ATT.cov <- c("SanCov|PATT|U.PS|W.OM", "SanCov|PATT|W.PS|W.OM", "SanCov|PATT|C.PS|W.OM", "SanCov|PATT|CW.PS|W.OM")
San_ATO.cov <- c("SanCov|PATO|U.PS|W.OM", "SanCov|PATO|W.PS|W.OM", "SanCov|PATO|C.PS|W.OM", "SanCov|PATO|CW.PS|W.OM")

San_ATE_COV <- sapply(results_sum, function(x) {colMeans(x$ests[,San_ATE.cov])})
San_ATT_COV <- sapply(results_sum, function(x) {colMeans(x$ests[,San_ATT.cov])})
San_ATO_COV <- sapply(results_sum, function(x) {colMeans(x$ests[,San_ATO.cov])})

xtable(rbind(San_ATE_COV,San_ATT_COV,San_ATO_COV))
round(rbind(San_ATE_COV,San_ATT_COV,San_ATO_COV),3)


San_ATE_MOM.cov <- c("SanCov|PATE_MOM|U.PS|W.OM", "SanCov|PATE_MOM|W.PS|W.OM", "SanCov|PATE_MOM|C.PS|W.OM", "SanCov|PATE_MOM|CW.PS|W.OM")
San_ATT_MOM.cov <- c("SanCov|PATT_MOM|U.PS|W.OM", "SanCov|PATT_MOM|W.PS|W.OM", "SanCov|PATT_MOM|C.PS|W.OM", "SanCov|PATT_MOM|CW.PS|W.OM")
San_ATO_MOM.cov <- c("SanCov|PATO_MOM|U.PS|W.OM", "SanCov|PATO_MOM|W.PS|W.OM", "SanCov|PATO_MOM|C.PS|W.OM", "SanCov|PATO_MOM|CW.PS|W.OM")

San_ATE_COV_MOM <- sapply(results_sum, function(x) {colMeans(x$ests[,San_ATE_MOM.cov])})
San_ATT_COV_MOM <- sapply(results_sum, function(x) {colMeans(x$ests[,San_ATT_MOM.cov])})
San_ATO_COV_MOM <- sapply(results_sum, function(x) {colMeans(x$ests[,San_ATO_MOM.cov])})

xtable(rbind(San_ATE_COV_MOM,San_ATT_COV_MOM,San_ATO_COV_MOM))
round(rbind(San_ATE_COV_MOM,San_ATT_COV_MOM,San_ATO_COV_MOM),3)


San_ATE_CVR.cov <- c("SanCov|PATE_CVR|U.PS|W.OM", "SanCov|PATE_CVR|W.PS|W.OM", "SanCov|PATE_CVR|C.PS|W.OM", "SanCov|PATE_CVR|CW.PS|W.OM")
San_ATT_CVR.cov <- c("SanCov|PATT_CVR|U.PS|W.OM", "SanCov|PATT_CVR|W.PS|W.OM", "SanCov|PATT_CVR|C.PS|W.OM", "SanCov|PATT_CVR|CW.PS|W.OM")
San_ATO_CVR.cov <- c("SanCov|PATO_CVR|U.PS|W.OM", "SanCov|PATO_CVR|W.PS|W.OM", "SanCov|PATO_CVR|C.PS|W.OM", "SanCov|PATO_CVR|CW.PS|W.OM")

San_ATE_COV_CVR <- sapply(results_sum, function(x) {colMeans(x$ests[,San_ATE_CVR.cov])})
San_ATT_COV_CVR <- sapply(results_sum, function(x) {colMeans(x$ests[,San_ATT_CVR.cov])})
San_ATO_COV_CVR <- sapply(results_sum, function(x) {colMeans(x$ests[,San_ATO_CVR.cov])})

xtable(rbind(San_ATE_COV_CVR,San_ATT_COV_CVR,San_ATO_COV_CVR))
round(rbind(San_ATE_COV_CVR,San_ATT_COV_CVR,San_ATO_COV_CVR),3)


San_ATE_WET.cov <- c("SanCov|PATE_WET|U.PS|W.OM", "SanCov|PATE_WET|W.PS|W.OM", "SanCov|PATE_WET|C.PS|W.OM", "SanCov|PATE_WET|CW.PS|W.OM")
San_ATT_WET.cov <- c("SanCov|PATT_WET|U.PS|W.OM", "SanCov|PATT_WET|W.PS|W.OM", "SanCov|PATT_WET|C.PS|W.OM", "SanCov|PATT_WET|CW.PS|W.OM")
San_ATO_WET.cov <- c("SanCov|PATO_WET|U.PS|W.OM", "SanCov|PATO_WET|W.PS|W.OM", "SanCov|PATO_WET|C.PS|W.OM", "SanCov|PATO_WET|CW.PS|W.OM")

San_ATE_COV_WET <- sapply(results_sum, function(x) {colMeans(x$ests[,San_ATE_WET.cov])})
San_ATT_COV_WET <- sapply(results_sum, function(x) {colMeans(x$ests[,San_ATT_WET.cov])})
San_ATO_COV_WET <- sapply(results_sum, function(x) {colMeans(x$ests[,San_ATO_WET.cov])})

xtable(rbind(San_ATE_COV_WET,San_ATT_COV_WET,San_ATO_COV_WET))
round(rbind(San_ATE_COV_WET,San_ATT_COV_WET,San_ATO_COV_WET),3)

round(rbind(San_ATE_COV,San_ATE_COV_MOM,San_ATE_COV_CVR,San_ATE_COV_WET),3)
round(rbind(San_ATT_COV,San_ATT_COV_MOM,San_ATT_COV_CVR,San_ATT_COV_WET),3)
round(rbind(San_ATO_COV,San_ATO_COV_MOM,San_ATO_COV_CVR,San_ATO_COV_WET),3)








