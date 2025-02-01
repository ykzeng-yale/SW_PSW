```rmd
---
title: "SW_PSW README"
output: github_document
---

# SW_PSW

**SW_PSW** is a collection of R scripts and functions to illustrate simulation-based analyses for **Propensity Score Weighting in Survey Observational Studies**, corresponding to the paper:

> **Moving toward best practice when using propensity score weighting in survey observational studies**  
> by **Yukang Zeng, Fan Li, Guangyu Tong**.  
> _arXiv:2501.16156 \[stat.ME\] (2025)_  
> [arXiv link](https://arxiv.org/abs/2501.16156)

---

## Overview

In this repository, we demonstrate how to:

1. **Generate super-populations** that exhibit “good overlap” between treatment and control groups.
2. **Implement multi-stage sampling** to obtain survey data from the super-population.
3. **Estimate causal effects** under different target populations (ATE, ATT, ATO) using:
   - **PSW** (Propensity Score Weighting),
   - **MOM** (Moment-based Augmented Weighted Estimator),
   - **CVR** (Clever Covariate Adjustment),
   - **WET** (Weighted Regression).
4. **Compute robust variance estimates** (sandwich variance) for each estimator.
5. **Assess performance** via bias, variance, coverage in simulation studies.
6. **Extract, summarize, and visualize** simulation results (e.g., relative bias, coverage, standardized mean differences).

---

## Repository Files

```
.
├─ EstSOD_Function.R
├─ README.md
├─ Results_Extraction.R
├─ Sampling_MultiStage_and_Estimation.R
├─ Shared_PSW_MOM_CVR_WET_Function.R
├─ Super_Population_GoodOvelap.R
└─ ...
```

1. **Super_Population_GoodOvelap.R**  
   - Generates a large super-population with defined covariates, logistic model for treatment assignment, and counterfactual outcomes.  
   - Saves `SuperPopGen_GoodOverlap_AllinOne_0118_2025.RData`.

2. **Shared_PSW_MOM_CVR_WET_Function.R**  
   - Contains unified functions for the four estimators (PSW, MOM, CVR, WET) and their sandwich variances:
     - `PSW_Estimator`, `PSW_SandwichVariance`
     - `MOM_Estimator`, `MOM_SandwichVariance`
     - `CVR_Estimator`, `CVR_SandwichVariance`
     - `WET_Estimator`, `WET_SandwichVariance`

3. **EstSOD_Function.R**  
   - High-level function **`EstSOD()`** that:
     - Fits a logistic propensity model (`U`, `W`, `C`, `CW`).
     - Constructs weights (`ate.wt`, `att.wt`, `ato.wt`) and tilting functions.
     - Computes **PSW**, **MOM**, **CVR**, **WET** estimates for ATE/ATT/ATO.
     - Computes sandwich variances, biases, standardized mean differences.
     - Returns one named vector with all results.

4. **Sampling_MultiStage_and_Estimation.R**  
   - Illustrates how to perform **multi-stage sampling** from the super-population:
     - Defines sample sizes (strata, clusters, individuals).
     - Uses the `sampling` package (`mstage()`) to draw complex survey samples.
     - Repeats sampling + calling `EstSOD()` in parallel (using `foreach` + `doParallel`).
     - Aggregates the results, saves to `.RData`.

5. **Results_Extraction.R**  
   - Shows how to load `.RData` results, compute and display summary metrics (bias, coverage, relative efficiency, etc.).

---

## Usage

1. **Set up environment**  
   - Ensure R (≥ 4.0) is installed.
   - Install needed packages:

```r
install.packages(c(
  "survey", "MatchIt", "sampling", "doParallel", "foreach", 
  "data.table", "numDeriv", "MASS", "Hmisc", "ggplot2"
))
```

2. **Generate super-population**  
```r
source("Super_Population_GoodOvelap.R")
# This saves 'SuperPopGen_GoodOverlap_AllinOne_0118_2025.RData'
```

3. **Run sampling & estimation**  
```r
source("Sampling_MultiStage_and_Estimation.R")
# Defines scenario(s), runs replicate sampling, calls 'EstSOD()', saves .RData
```

4. **Analyze or extract results**  
```r
source("Results_Extraction.R")
# Summarizes simulation outcomes: bias, coverage, relative efficiency, etc.
```

---

## Corresponding Paper

This repository supports analyses in the paper:

> **Moving toward best practice when using propensity score weighting in survey observational studies**  
> by **Yukang Zeng, Fan Li, Guangyu Tong**.  
> _arXiv:2501.16156 \[stat.ME\] (2025)_  
> [arXiv link](https://arxiv.org/abs/2501.16156)

We provide:

- Balancing-weights framework for complex surveys.  
- M-estimator derivations & asymptotic properties.  
- Extensive simulation studies.  
- Empirical examples with real survey data.

---

## Contributing / Contact

- **Issues / Pull Requests** are welcome for suggestions, improvements, or questions.  
- For additional details, please see the paper or contact the authors as listed there.

---

## License

MIT License or similar—see [LICENSE](https://opensource.org/licenses/MIT) (if provided).  
Feel free to reuse and adapt for research and educational purposes.
```
