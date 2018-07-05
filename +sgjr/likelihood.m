function [ ...
  llf, dLL_dtheta, lls, assetVariance, equityVolatility, smoothedBvDebt, mvEquityOut, ...
  garchLeverageMultiplier, asset2Debt, deltaImplied, annualizedForecastVol, tauOut] = ...
  likelihood(parameters, ret, equity, debt, riskFree, tau, forecastType)
% PURPOSE:
%        Generate recursive asset volatility, equity volatility, and delta
%        for Structural GARCH
%
%
% INPUTS:
%    parameters      - A 5x1 vector of GJR parameters [omega; alpha; gamma;
%                       beta; phi].
%    inputData       - A struct with the following attributes
%        return          - A T by 1 vector of zero mean returns
%        equity          - A T by 1 vector of market value of equity
%        debt            - A T by 1 vector of book value of debt
%        riskFree        - A T by 1 vector of annualized risk-free rates
%        forecastType    - Either 'DF' for dynamic forecast or 'CF' for
%                          constant forecast.  Details in ES (2014)
%

% TODO - determine which outputs we actually need
% OUTPUTS:
%    parameters              - A vector of parameters Structural-GARCH parameters of the
%                              form (omega, alpha, gamma, beta, phi)
%    lls                     - time series of scores from ML
%    assetVariance           - (T-burnIndex)x1 series of asset variances
%    equityVolatility        - (T-burnIndex)x1 series of annualized equity vol
%    smoothedBvDebt          - (T-burnIndex)x1 series of smoothed bvDebt
%    mvEquityOut             - (T-burnIndex)x1 series of market value of equity
%    garchLeverageMultiplier - (T-burnIndex)x1 series of leverage multiplier
%    asset2Debt              - (T-burnIndex)x1 series of B-S asset to debt series
%    deltaImplied            - (T-burnIndex)x1 series of B-S delta
%    annualizedForecastVol   - (T-burnIndex)x1 series of asset forecast vol
%    tauOut                  - (T-burnIndex)x1 series of debt maturities

% NOTES:
%   -   The outputs for smoothed debt, market value of equity, and garchLeverageMultiplier
%       are truncated and lagged (so in T-1 info set) by the burn index, which
%       is a burnout phase for the model
%
% Author: Emil Siriwardane
% esiriwar@stern.nyu.edu
% Date: 3/15/2014

  validateInput(parameters, ret, equity, debt, riskFree, tau, forecastType)
  % validateInput(parameters, inputData);

  if ~isreal(parameters)
    assetVariance = [];
    llf = Inf;
    lls = [];
    equityVolatility = [];

    return;
  end

  % TODO - what's the significance of the number 21?  Will it always be 21?
  % Since we smooth the debt, we discard the first #burnIndex-1 observations
  % i.e. we take the data from t = burnIndex:end
  burnIndex = 21;

  %%% Smooth the debt %%%
  phi = 0.01;
  % smoothedDebt = util.exponentiallySmooth(inputData.debt, phi);
  smoothedDebt = util.exponentiallySmooth(debt, phi);
  equity2debt = equity ./ smoothedDebt;
  % equity2debt = inputData.equity ./ smoothedDebt;

  %%% Truncate debt and equity for output - lagged so in T-1 info set %%%
  smoothedBvDebt = smoothedDebt(burnIndex - 1:end - 1);
  mvEquityOut = equity(burnIndex - 1:end - 1);
  tauOut = tau(burnIndex - 1:end - 1);
  % mvEquityOut = inputData.equity(burnIndex - 1:end - 1);
  % tauOut = inputData.tau(burnIndex - 1:end - 1);

  %%% Generate equity and asset vol series %%%
  [assetVariance, garchLeverageMultiplier, asset2Debt, deltaImplied,annualizedForecastVol, ...
   dLM_dtheta, dha_dtheta] = sgjr.core(parameters, ret, equity, debt, riskFree, tau, forecastType, equity2debt);
  % dLM_dtheta, dha_dtheta] = sgjr.core(parameters, inputData, equity2debt);

  % This will be T-1 units long since t-1 garchLeverageMultiplier is used for time t volatility.
  equityVariance = (garchLeverageMultiplier(1:end - 1).^2) .* assetVariance(2:end);

  %Corresponds to time 2,...,.T.  Set equityVariance(1) = unconditional var of returns

  %%% Throw away burnIndex in observations %%%
  equityVariance = equityVariance(burnIndex - 1:end);

  %%% Likelihood calculation %%%
  [lls, llf] = computeLikelihood(equityVariance, ret(burnIndex:end));

  dLL_dtheta = computeDerivatives( ...
    assetVariance, equityVariance, ret, ...
    garchLeverageMultiplier, dLM_dtheta, dha_dtheta, burnIndex ...
  );

  %Other outputs
  assetVariance = assetVariance(burnIndex:end);
  asset2Debt = asset2Debt(burnIndex-1:end-1);
  deltaImplied = deltaImplied(burnIndex-1:end-1);
  annualizedForecastVol = annualizedForecastVol(burnIndex-1:end-1);
  garchLeverageMultiplier = garchLeverageMultiplier(burnIndex-1:end-1);

  if isnan(llf) || sum(assetVariance < 0) > 0
    llf = Inf;
  end

  equityVolatility = sqrt(252 * equityVariance);
end

function validateInput(parameters, ret, equity, debt, riskFree, tau, forecastType)
% function validateInput(parameters, inputData)
  numParam = size(parameters, 1);

  if numParam ~= 5
    error('The number of input parameters should be 5');
  end

  [numObs, numReturnCols] = size(ret);
  % [numObs, numReturnCols] = size(inputData.return);
  if numReturnCols ~= 1
    error(['The return series should have col dim 1, but has col dim ' num2str(tmp)])
  end

  numEquityObs = size(equity, 1);
  % numEquityObs = size(inputData.equity, 1);
  if numEquityObs ~= numObs
    error('The market value of equity series and the return series should be the same dim');
  end

  numDebtObs = size(debt);
  if numDebtObs ~= numObs
    error('The book value of debt series and the return series should be the same dim');
  end

  try
    [ret equity debt riskFree]; %#ok
  catch
    error([ ...
      'The return, market value of equity, book value of debt, ' ...
      'and risk-free rate series must all be same size' ...
    ]);
  end
end

function [lls, llf] = computeLikelihood(equityVariance, returnData)
  lls = -0.5 * log(2 * pi) - 0.5 * log(equityVariance) - 0.5 * (returnData.^2) ./ equityVariance;
  lls = -lls;

  llf = sum(lls);
end

function dLL_dtheta = computeDerivatives( ...
  assetVariance, equityVariance, returnData, ...
  garchLeverageMultiplier, dLM_dtheta, dha_dtheta, burnIndex ...
)
  dhe_dtheta = ones(length(assetVariance) - 1, 5) * NaN;

  % TODO - vectorize calculation?
  dhe_dtheta(:, 1) = ...
      2 * garchLeverageMultiplier(1:end - 1) .* dLM_dtheta(1:end - 1, 1) .* assetVariance(2:end) + ...
      (garchLeverageMultiplier(1:end-1).^2) .* dha_dtheta(2:end, 1);

  dhe_dtheta(:, 2) = ...
      2 * garchLeverageMultiplier(1:end - 1) .* dLM_dtheta(1:end - 1, 2) .* assetVariance(2:end) + ...
      (garchLeverageMultiplier(1:end - 1).^2) .* dha_dtheta(2:end, 2);

  dhe_dtheta(:, 3) = ...
      2 * garchLeverageMultiplier(1:end - 1) .* dLM_dtheta(1:end - 1, 3) .* assetVariance(2:end) + ...
      (garchLeverageMultiplier(1:end - 1).^2) .* dha_dtheta(2:end, 3);

  dhe_dtheta(:, 4) = ...
      2 * garchLeverageMultiplier(1:end - 1) .* dLM_dtheta(1:end - 1, 4) .* assetVariance(2:end) + ...
      (garchLeverageMultiplier(1:end - 1).^2) .* dha_dtheta(2:end, 4);

  dhe_dtheta(:, 5) = ...
      2 * garchLeverageMultiplier(1:end - 1) .* dLM_dtheta(1:end - 1, 5) .* assetVariance(2:end) + ...
      (garchLeverageMultiplier(1:end - 1).^2) .* dha_dtheta(2:end, 5);

  LL_sum_helper = 1 ./ equityVariance - returnData(burnIndex:end).^2 ./ equityVariance.^2;

  % TODO - vectorize calculation?
  dLL_dtheta = ones(5,1);
  dLL_dtheta(1) = 0.5 * sum(LL_sum_helper .* dhe_dtheta(burnIndex-1:end, 1));
  dLL_dtheta(2) = 0.5 * sum(LL_sum_helper .* dhe_dtheta(burnIndex-1:end, 2));
  dLL_dtheta(3) = 0.5 * sum(LL_sum_helper .* dhe_dtheta(burnIndex-1:end, 3));
  dLL_dtheta(4) = 0.5 * sum(LL_sum_helper .* dhe_dtheta(burnIndex-1:end, 4));
  dLL_dtheta(5) = 0.5 * sum(LL_sum_helper .* dhe_dtheta(burnIndex-1:end, 5));
end