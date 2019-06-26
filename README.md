# Structural GARCH Code

This MATLAB package is used to estimate the Structural GJR-GARCH (SGJR) model from the Engle and Siriwardane paper [Structural GARCH: The Volatility-Leverage Connection](https://academic.oup.com/rfs/article-abstract/31/2/449/4139801?redirectedFrom=fulltext).

## Running the Code
Clone this repository into your MATLAB working directory:

```bash
cd /home/user/projects/project/using/SGJR/
git clone [repository URI]  # creates an 'sgjr' directory in your working directory
```

Then, in your MATLAB code, add the SGJR directory to your path and reference it using the `sgjr.` prefix.

```MATLAB
addpath(fullfile(pwd, 'sgjr'));
load(fullfile(pwd, 'sgjr', testData));
AIG_cell = struct2cell(AIG);
result = sgjr(AIG_cell{:}, zerocurve);
```

The call to the `sgjr` function will return an SGJR result object containing:

* parameters - A vector of parameters Structural-GARCH parameters of the form (omega, alpha, gamma, beta, phi)
* vcv - 5X5 QMLE Robust Variance-Covariance Matrix
* stdErrorsMLE - MLE standard errors
* stdErrorsQMLE - QMLE standard errors
* equityVariance - daily equity variance
* assetVariance - (T - burn_index) x 1 series of asset variances            # TODO: What is burn_index?
* leverageMultiplier - (T - burn_index) x 1 series of leverage multiplier   # TODO: What is burn_index?
* loglikelihood - log-likelihood at optimum
* BIC - Bayes Information Criterion
* numRuns - number of different starting parameters
* smoothedBvDebt - (T - burn_index) x 1 series of smoothed bvdebt          # TODO: What is burn_index?
* mvEquity - (T - burn_index) x 1 series of market value of equity         # TODO: What is burn_index?
* asset2Debt - (T - burn_index) x 1 series of B-S asset to debt series     # TODO: What is burn_index?
* deltaImplied - (T - burn_index) x 1 series of B-S delta                  # TODO: What is burn_index?
* volatilityForecast - (T - burn_index) x 1 series of asset forecast vol   # TODO: What is burn_index?
* tau - TODO - define
* gjr - Parameters from estimation of a GJR-GARCH regression on the returns provided

#### Options

You can also pass in a `forecastType` option to specify either a 'Constant' or 'Dynamic' forecast ('Constant' is the default).

```MATLAB
result = sgjr(AIG_cell{:}, zerocurve, forecastType, 'Dynamic');
```

* Dataset required

## Theory
