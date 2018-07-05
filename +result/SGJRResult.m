classdef SGJRResult
% A result object that holds the estimation results for an SGJR estimation.  An SJRResult object
% contains the following properties:
%
%    parameters              - A vector of parameters Structural-GARCH parameters of the
%                              form (omega, alpha, gamma, beta,phi)
%    vcv                     - 5X5 QMLE Robust Variance-Covariance Matrix
%    stdErrorsMLE            - MLE standard errors
%    stdErrorsQMLE           - QMLE standard errors
%    equityVariance          - daily equity variance
%    assetVariance           - (T-burn_index)x1 series of asset variances
%    leverageMultiplier      - (T-burn_index)x1 series of leverage multiplier
%    loglikelihood           - log-likelihood at optimum
%    BIC                     - Bayes Information Criterion
%    numRuns                 - number of different starting parameters
%    smoothedBvDebt          - (T-burn_index)x1 series of smoothed bvdebt
%    mvEquity                - (T-burn_index)x1 series of market value of equity
%    asset2Debt              - (T-burn_index)x1 series of B-S asset to debt series
%    deltaImplied            - (T-burn_index)x1 series of B-S delta
%    volatilityForecast      - (T-burn_index)x1 series of asset forecast vol
%    tau                     - TODO - define
%    gjr                     - Parameters from estimation of a GJR-GARCH regression on the returns provided

  properties
    parameters
    vcv
    stdErrorsMLE
    stdErrorsQMLE
    equityVariance
    assetVariance
    leverageMultiplier
    logLikelihood
    smoothedBvDebt
    mvEquity
    asset2Debt
    deltaImplied
    volatilityForecast
    tau
    gjr
  end

  properties (Dependent)
    BIC
    tStatistics
    tStatisticsRobust
  end

  methods
    function ic = get.BIC(obj)
      T = length(obj.equityVariance);
      ic = -2 * obj.logLikelihood / T + (length(obj.parameters)) * log(T) / T;
    end

    function tstat = get.tStatistics(obj)
      tstat = obj.parameters ./ obj.stdErrorsMLE;
    end

    function tstat = get.tStatisticsRobust(obj)
      tstat = obj.parameters ./ obj.stdErrorsQMLE;
    end
  end
end