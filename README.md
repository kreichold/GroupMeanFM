# GroupMeanFM
Group-Mean Fully Modified OLS Estimation and Inference in Panel Cointegrating Regressions

## Introduction
This repository contains MATLAB code to estimate panel cointegrating polynomial regressions using a group-mean fully modified OLS estimator as proposed in Wagner and Reichold (2023). The code also allows to perform standard and cross-section dependence robust inference based upon the group-mean fully modified OLS estimator. No knowledge about the presence or absence of deterministic drifts in the stochastic regressors is required.

## Usage
Download the files and move them into your current working directory, `pwd`.

## Main Function

+ **GroupMeanFMOLS.m**
This is the only function that needs to be executed by the practitioner. It returns group-mean fully modified OLS estimation results. 

## Auxiliary Functions

+ **And_HAC91.m**
This function implements the data-dependent bandwidth selection rule of Andrews (1991).

+ **demean_detrend.m**
This function allows to demean and detrend time series.

+ **lr_varmod.m**
This function estimates long-run variance matrices.

+ **lr_weights.m**
This function computes kernel weights used in **lr_varmod.m**.

## Reference
Andrews, D.W.K. (1991). Heteroskedasticity and Autocorrelation Consistent Covariance Matrix Estimation. *Econometrica* **59**, 817--858.

Wagner M., Reichold, K. (2023). [Panel Cointegrating Polynomial Regressions: Group-Mean Fully Modified OLS Estimation and Inference](https://doi.org/10.48550/arXiv.2204.01373](https://doi.org/10.1080/07474938.2023.2178141). *Econometric Reviews*. Forthcoming.
