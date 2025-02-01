## Title: SW_PSW

## Author
**Yukang Zeng, Fan Li and Guangyu Tong**

## Description
This folder contains simulation code and examples for the paper:

> **Moving toward best practice when using propensity score weighting in survey observational studies**  
> by **Yukang Zeng, Fan Li, Guangyu Tong**.  
> _arXiv:2501.16156 [stat.ME] (2025)_  
> [arXiv link](https://arxiv.org/abs/2501.16156)

The code demonstrates how to generate a super-population with good overlap, perform multi-stage sampling, and estimate treatment effects under different targets (ATE, ATT, ATO) using four methods:

- **PSW (Propensity Score Weighting)**
- **MOM (Moment-based Augmented Weighted Estimator)**
- **CVR (Clever Covariate)**
- **WET (Weighted Regression)**

We also provide sandwich variance formulas and code to compute standardized mean differences (SMD).

---

## Files

1. **Super_Population_GoodOvelap.R**  
   - Generates a large super-population with user-defined covariates, logistic PS model, counterfactual outcomes, and an example histogram.  
   - Outputs `SuperPopGen_GoodOverlap_AllinOne_0118_2025.RData`.

2. **Shared_PSW_MOM_CVR_WET_Function.R**  
   - Defines **unified** functions for each estimator and their robust (sandwich) variance:
     - `PSW_Estimator`, `PSW_SandwichVariance`
     - `MOM_Estimator`, `MOM_SandwichVariance`
     - `CVR_Estimator`, `CVR_SandwichVariance`
     - `WET_Estimator`, `WET_SandwichVariance`
   - Also includes `SMD_Calculator` for standardized mean differences.

3. **EstSOD_Function.R**  
   - A high-level wrapper function named `EstSOD()` that:
     1. Fits logistic PS (four modes: "U","W","C","CW").
     2. Constructs survey-based weights for ATE, ATT, ATO.
     3. Calls **PSW**, **MOM**, **CVR**, **WET** estimators.
     4. Computes sandwich variances, biases, coverage, and SMD.
     5. Returns a single, named result vector of all quantities.

4. **Sampling_MultiStage_and_Estimation.R**  
   - Example multi-stage sampling:
     1. Loads the super-population.
     2. Performs multi-stage stratified+cluster sampling using the `sampling` package.
     3. Repeats sampling + `EstSOD()` calls in parallel.
     4. Stores results in `.RData` format.

5. **Results_Extraction.R**  
   - Demonstrates how to read the saved `.RData`, extract or summarize results:
     - Monte Carlo bias, variance, coverage, relative efficiency, etc.
     - Example code for printing or tabulating with `xtable`.

---

## Usage

1. **Generate Super-Population**  
   ```{r gen_pop}
   source("Super_Population_GoodOvelap.R")
   # Outputs SuperPopGen_GoodOverlap_AllinOne_0118_2025.RData
   ```

2. **Call Multi-Stage Sampling & Estimation**
   ```{r sampling}
   source("Sampling_MultiStage_and_Estimation.R")
   # Repeatedly samples from super-pop, calls `EstSOD()`, saves results
   ```

3. **Extract & Summarize**
   ```{r extract}
   source("Results_Extraction.R")
   # Summarize bias, coverage, etc. from the stored results
   ```

## Example Simulation Scenarios

* **Scenario**: We control the number of strata (10) and clusters per stratum (20), then define sample sizes at each stage.
* **Replicates**: We typically use `n.reps=5000` for stable results.
* **Covariates**: 6 baseline covariates, plus an optional 7th for missingness scenarios.
* **Treatment**: logistic model with moderate overlap.
* **Outcome**: continuous outcome with a known super-population ATE, ATT, ATO.

## Citation

If you use or adapt these scripts for your research, please cite:

**Moving toward best practice when using propensity score weighting in survey observational studies**
Yukang Zeng, Fan Li, Guangyu Tong
arXiv:2501.16156 [stat.ME], 2025.
https://arxiv.org/abs/2501.16156

## Contact

For questions, bug reports, or requests, please open a GitHub Issue or contact the authors listed in the manuscript.
