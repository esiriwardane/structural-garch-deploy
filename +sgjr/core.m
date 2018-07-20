function [ ...
  assetVariance, garchLeverageMultiplier, asset2Debt, deltaImplied, annualizedForecastVol, dLM_dtheta, ...
  dhaDTheta] = core(param, inputData, equity2Debt)
%%%%%%%% HELPER FUNCTION FOR STRUCTURAL GARCH %%%%%%%%%

  T = length(inputData.return);

  % Initialize parameters and construct matrices
  omega = param(1);
  archp = param(2);
  asymp = param(3);
  garchp = param(4);
  phi = param(5);

  [asset2Debt, deltaImplied, garchLeverageMultiplier, ...
   annualizedForecastVol, ra, lrv] = deal(NaN * ones(T, 1));

  [dLM_dtheta, dhaDTheta, dsigmaA_dtheta] = deal(NaN * ones(T, 5));

  %%% Compute unconditional daily variance and unconditional variance over
  %%% the life of the option %%%
  longRunVar = omega / (1 - archp - 0.5 * asymp - garchp);
  lrVarOption = longRunVar * 252 * inputData.tau;
  assetVariance = longRunVar * ones(T, 1);

  if longRunVar < 0 || isreal(longRunVar)==0
    return
  end

  %%% Run the recursion
  for I = 1:T - 1
    theta = archp + 0.5 * asymp + garchp;
    k = inputData.tau(I) * 252;

    if I == 1
      ra(I) = inputData.return(I);   %LM(0) = 1
    else
      ra(I) = inputData.return(I) / garchLeverageMultiplier(I - 1);
    end

    assetVariance(I + 1) = omega + archp * (ra(I)^2) + ...
        asymp * (ra(I)^2) * (ra(I) < 0) + garchp * assetVariance(I);   %Tomorrow's asset variance forecast

    if isDynamicForecast(inputData.forecastType)
      sigma_tk2 = (omega / (1 - theta)) * k - (omega / (1 - theta)) * ((1 - theta.^k) / (1 - theta)) + ...
          assetVariance(I + 1) * ((1 - theta.^k) / (1 - theta));
      if sigma_tk2 < 0
        sigma_tk2 = lrVarOption(I);
      end
    else % Constant Forecast forecast (Default) %%%
      sigma_tk2 = lrVarOption(I);
    end

    %%% Compute annualized asset volatility forecast %%%
    annualizedForecastVol(I) = sqrt(sigma_tk2 / inputData.tau(I));

    %%% Invert Black-Scholes to get asset to debt ratio %%%
    asset2Debt(I) = sgjr.assetFixedpoint( ...
      1, inputData.riskFree(I), inputData.tau(I), sigma_tk2, ...
      equity2Debt(I), (equity2Debt(I) + 1)/equity2Debt(I) ...
    );
    d1Contemp = (log(asset2Debt(I)) + ...
                 (inputData.riskFree(I) * inputData.tau(I) + sigma_tk2 / 2)) / sqrt(sigma_tk2);

    deltaImplied(I) = 0.5 * erfc(-d1Contemp ./ sqrt(2));

    garchLeverageMultiplier(I) = (deltaImplied(I)^phi) * ((asset2Debt(I) / equity2Debt(I))^phi);

    %% Compute analytical derivatives
    %Helpers for dLM
    vega = sgjr.fastBsVega( ...
      asset2Debt(I), 1, inputData.riskFree(I), annualizedForecastVol(I), inputData.tau(I) ...
    );
    assetSigmaA = -vega / deltaImplied(I);

    d2 = d1Contemp - annualizedForecastVol(I) * sqrt(inputData.tau(I));
    bsVanna = -sgjr.fastNormPdf(d1Contemp, 0, 1) * d2 / annualizedForecastVol(I);
    bsGamma = sgjr.fastNormPdf(d1Contemp, 0 , 1) / ...
              (asset2Debt(I) * annualizedForecastVol(I) * sqrt(inputData.tau(I)));
    deltaSigmaA = bsGamma * assetSigmaA + bsVanna;

    preMultiply = phi * garchLeverageMultiplier(I) * ...
        (deltaSigmaA / deltaImplied(I) + assetSigmaA / asset2Debt(I));
    deltaAE = deltaImplied(I) * asset2Debt(I) / equity2Debt(I);

    if I == 1
      dhaDTheta(I + 1, 1) = 1 + garchp / (1 - theta);
      dhaDTheta(I + 1, 2) = ra(I)^2 + garchp * omega / (1 - theta)^2;
      dhaDTheta(I + 1, 3) = (ra(I) < 0) * ra(I)^2 + 0.5 * garchp * omega / (1 - theta)^2;
      dhaDTheta(I + 1, 4) = garchp * omega / (1 - theta)^2 + longRunVar;
      dhaDTheta(I + 1, 5) = 0;
    else
      %Helper constants
      ralmA = archp * ra(I)^2 / garchLeverageMultiplier(I - 1);
      ralmG = asymp * (ra(I) < 0) * ra(I)^2 / garchLeverageMultiplier(I - 1);

      dhaDTheta(I + 1, 1) = 1 - 2 * (ralmA + ralmG) * dLM_dtheta(I - 1, 1) + ...
          garchp * dhaDTheta(I, 1);
      dhaDTheta(I + 1, 2) = -2 * (ralmA + ralmG) * dLM_dtheta(I - 1, 2) + ...
          ra(I)^2 + garchp * dhaDTheta(I, 2);
      dhaDTheta(I + 1, 3) = -2*(ralmA + ralmG)*dLM_dtheta(I - 1, 3) + ...
          (ra(I) < 0) * ra(I)^2 + garchp * dhaDTheta(I, 3);
      dhaDTheta(I + 1, 4) = -2 * (ralmA + ralmG) * dLM_dtheta(I - 1, 4) + ...
          garchp * dhaDTheta(I, 4) + assetVariance(I);
      dhaDTheta(I + 1, 5) = -2 * (ralmA + ralmG) * dLM_dtheta(I - 1, 5) + ...
          garchp * dhaDTheta(I, 5);
    end

    if isDynamicForecast(inputData.forecastType)
      %Compute long-run variance
      lrv(I) = annualizedForecastVol(I)^2 * inputData.tau(I);

      %Define some helper constants
      mtheta = 1 - theta;
      mthetak = 1 - theta^k;
      mratiok = mthetak / mtheta;
      helper = (-k * theta^(k - 1)) / (mtheta) + mthetak / (mtheta^2);

      commonTerm = (k * omega) / mtheta^2 - (omega / mtheta^2) * mratiok - ...
          (omega / mtheta) * helper + helper * assetVariance(I + 1);

      %Compute derivatives
      dlrvDOmega = (k / mtheta) - mratiok / mtheta + (mthetak / mtheta) * dhaDTheta(I + 1, 1);
      dlrvDAlpha =  commonTerm + mratiok * dhaDTheta(I + 1, 2);
      dlrvDGamma =  0.5 * commonTerm + mratiok * dhaDTheta(I + 1, 3);
      dlrvDBeta =  commonTerm + mratiok * dhaDTheta(I + 1, 4);
      dlrvDPhi = mratiok * dhaDTheta(I + 1, 5);

      lrv_mult = 0.5 / sqrt(inputData.tau(I) * lrv(I));
      dsigmaA_dtheta(I, 1) = lrv_mult * dlrvDOmega;
      dsigmaA_dtheta(I, 2) = lrv_mult * dlrvDAlpha;
      dsigmaA_dtheta(I, 3) = lrv_mult * dlrvDGamma;
      dsigmaA_dtheta(I, 4) = lrv_mult * dlrvDBeta;
      dsigmaA_dtheta(I, 5) = lrv_mult * dlrvDPhi;
    else
      %Compute long-run asset vol derivatives
      aFac = sqrt(252) * 0.5;

      dsigmaA_dtheta(I, 1) = aFac * (omega * (1-theta))^(-1/2);
      dsigmaA_dtheta(I, 2) = aFac * sqrt(omega) / (1-theta)^(3/2);
      dsigmaA_dtheta(I, 3) = aFac * 0.5 * sqrt(omega) / (1-theta)^(3/2);
      dsigmaA_dtheta(I, 4) = aFac * sqrt(omega) / (1-theta)^(3/2);
      dsigmaA_dtheta(I, 5) = 0;
    end

    dLM_dtheta(I, 1) = preMultiply * dsigmaA_dtheta(I, 1);
    dLM_dtheta(I, 2) = preMultiply * dsigmaA_dtheta(I, 2);
    dLM_dtheta(I, 3) = preMultiply * dsigmaA_dtheta(I, 3);
    dLM_dtheta(I, 4) = preMultiply * dsigmaA_dtheta(I, 4);
    dLM_dtheta(I, 5) = log(deltaAE) * garchLeverageMultiplier(I) +  dsigmaA_dtheta(I, 5) * preMultiply;
  end
end

function TF = isDynamicForecast(forecastType)
  TF = strcmp(forecastType, 'DF');
end