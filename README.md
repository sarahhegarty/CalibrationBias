# CalibrationBias
Investigating group fairness through estimating differential calibration bias with varying coefficient (VC) models

# Contents - Overview
1. vc-funcs-2025-07-21.R: R script containing functions for implementing the varying coefficient model approach
2. simulation-code: folder containing R scripts for simulation data generation and implementing Monte Carlo experiments in a range of settings
3. simulation-summary-results: folder containing the Rmd files for summarizing the Monte Carlo experiment results and the resulting html reports
4. synthetic-data-example: folder containing the Rmd files needed to reproduce the synthetic data example as well as resulting html reports and figures. Note that the synthetic dataset used is available at https://gitlab.com/labsysmed/dissecting-bias.

# Details: vc-funcs
- The primary function which implements the VC model is vcfit.bvs().
- Function for calculating the residual squares criterion (RSC), rsc(), and mean squared error (MSE) of the primary curve, zheng.mse(), and the bias curve, zheng.mse.cubics(), are provided.
- Supporting functions such as colMax() which returns the maximum value in a column, kernel() which calculates the kernel are also included.
  
  
