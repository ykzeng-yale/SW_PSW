# SW_PSW

**Authors:** Yukang Zeng, Fan Li, and Guangyu Tong

## Overview
This repository contains simulation code and examples accompanying the paper:

> **Moving Toward Best Practice When Using Propensity Score Weighting in Survey Observational Studies**  
> Yukang Zeng, Fan Li, and Guangyu Tong, _arXiv:2501.16156 [stat.ME] (2025)_  
> [https://arxiv.org/abs/2501.16156](https://arxiv.org/abs/2501.16156)

The code demonstrates how to generate super-populations under both good and poor overlap conditions; implement multi-stage and treatment‐dependent sampling designs; and estimate treatment effects (ATE, ATT, ATO) using four methods:
- **PSW** (Propensity Score Weighting)
- **MOM** (Moment-based Augmented Estimator)
- **CVR** (Clever Covariate)
- **WET** (Weighted Regression)

Robust sandwich variance estimators are provided for three sampling scenarios:
- **Retrospective:** Sampling probabilities depend on treatment (or group) and covariates.
- **Independent:** Sampling is independent of treatment assignment.
- **Prospective:** Sampling occurs before treatment assignment (though treatment may be influenced by sampling).

Covariate balance is assessed via standardized mean differences (SMD).

## Files
- **Super_Population_GoodOverlap.R**  
  Generates a super-population with good overlap, including covariates, a logistic PS model, simulated treatment assignment, counterfactual outcomes, and visualizes the PS distribution.
  
- **Super_Population_PoorOverlap.R**  
  Similar to the above but generates a super-population with poor overlap conditions.
  
- **Shared_PSW_MOM_CVR_WET_Function.R**  
  Contains unified functions for the four estimators and their sandwich variance estimators:
  - `PSW_Estimator` & `PSW_SandwichVariance`
  - `MOM_Estimator` & `MOM_SandwichVariance`
  - `CVR_Estimator` & `CVR_SandwichVariance`
  - `WET_Estimator` & `WET_SandwichVariance`  
  Also includes `SMD_Calculator` for computing standardized mean differences.
  
- **EstSOD_Function.R**  
  A wrapper function `EstSOD()` that fits logistic PS models (with modes “U”, “W”, “C”, “CW”), constructs survey-based weights for ATE, ATT, and ATO under the three design scenarios (Retrospective, Independent, Prospective), calls the four estimators, and computes their sandwich variances, biases, coverage probabilities, and SMDs.
  
- **Sampling_MultiStage_and_Estimation.R**  
  Demonstrates multi-stage stratified and clustered sampling from a super-population and applies the estimation functions.
  
- **Sampling_DependOnTreatment_and_Estimation.R**  
  Implements sampling when the sampling process depends on treatment assignment.
  
- **Results_Extraction.R**  
  Provides scripts to load and summarize simulation results (e.g., bias, variance, coverage, relative efficiency).

## Simulation Scenarios
Simulations are conducted under varying overlap conditions and sampling designs. Three survey design scenarios are supported:
- **Retrospective:** Sampling probabilities depend on treatment (or group) and covariates.
- **Independent:** Sampling is independent of treatment.
- **Prospective:** Sampling occurs prior to treatment assignment, though treatment may be influenced by the sampling process.

Estimation methods are applied for targets ATE, ATT, and ATO under different model specifications.

## Citation
If you use or adapt these scripts, please cite:

*Moving Toward Best Practice When Using Propensity Score Weighting in Survey Observational Studies*  
Yukang Zeng, Fan Li, and Guangyu Tong  
_arXiv:2501.16156 [stat.ME], 2025_  
[https://arxiv.org/abs/2501.16156](https://arxiv.org/abs/2501.16156)

## Contact
For questions, bug reports, or suggestions, please open an issue on GitHub or contact the authors.
