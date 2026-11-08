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
  
# Details: Simulations
- dataGen.R contains the function for generating data in the simulations. It relies on cal.base() and cal.curve() for the specification of the calibration curve for each group. The specification of those curves used in the current work can be found in:
    - cal.base.R: the overall calibration curve
    - delta.shift.R: Delta taking the form of a constant shift relative to the overall curve, corresponds to the null when delta = 0
    - delta.bspline0.R: corresponds to the increasing-trend curve
    - delta.bspline1.R: corresponds to the zero-centered curve
    - delta.bspline2.R: corresponds to the down-shifted curve
- 13.04.03.04-runOnce.R contains the code for running one interaction of the simulation, including data generation, selecting the bandwidths, fitting the varying coefficient model and estimating the test statistics.
- 13.04.03.04-runMC-[setting]-[delta]-n[N].R: implement the Monte Carlo experiments for the simulation setting defined by the functional form of the calibration curve (setting = bs0/increasing-treand, bs1/zero-centered, bs2/down-shifted), effect size at tau=0.5 (delta=0.1, 0.2), and group sample size (N=500,1000).
