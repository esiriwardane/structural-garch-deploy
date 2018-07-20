function result = estimate(returnData, mvEquity, bvDebt, riskFreeRate, tau, forecastType)
% PURPOSE:
%        Estimates the Structural-GJR Model of Engle and Siriwardane
%        (2014), which we refer to as ES (2014) henceforth
%
%
% INPUTS:
%    returnData      - A T by 1 vector of zero mean returns
%    mvEquity        - A T by 1 vector of market value of equity
%    bvDebt          - A T by 1 vector of book value of debt
%    riskFreeRate    - A T by 1 vector of annualized risk-free rates
%    tau             - A scalar denoting the time to maturity of the debt
%                      Can also be a T by 1 vector of maturities
%    forecastType    - Either 'DF' for dynamic forecast or 'CF' for
%                      constant forecast.  Details in ES (2014)
%
% OUTPUTS:
%    result          - an SGJRResult object
%
% SEE ALSO: sgjr.likelihood, result.SGJRResult
%
% Author: Emil Siriwardane
% esiriwar@stern.nyu.edu
% Date: 3/15/2014

  validateInput(returnData, mvEquity, bvDebt, tau);

  if isscalar(tau)
    tau = repmat(tau, T, 1);
  end

  inputData = struct( ...
    'return', returnData, 'equity', mvEquity, 'debt', bvDebt, ...
    'riskFree', riskFreeRate, 'tau', tau, 'forecastType', forecastType ...
  );

  optimizationOptions = createOptimizationOptions();
  startingParameters = getGarchStartingParameters(inputData.return);

  % TODO -- populate the GJR coefficients in the result object
  problem = createOptimization(inputData, optimizationOptions, startingParameters);
  startingPoints = createStartingPoints();
  parameters = estimateSgjr(problem, startingPoints);

  result = collectResults(parameters, inputData);
end

function validateInput(returnData, mvEquity, bvDebt, tau)
  [numObs, returnColumns] = size(returnData);

  if returnColumns ~= 1
    error(['The return series should have col dim 1, but has col dim ' num2str(tmp)]);
  end

  numEquityObs = size(mvEquity, 1);

  if numEquityObs ~= numObs
    error('The market value of equity series and the return series should be the same dim');
  end

  numDebtObs = size(bvDebt, 1);
  if numDebtObs ~= numObs
    error('The book value of debt series and the return series should be the same dim');
  end

  if ~isscalar(tau) && length(tau) ~= numObs
    error('Tau and the return series should be the same dim');
  end
end

function options = createOptimizationOptions()
  options.optim = optimset('fmincon');

  options.optim = optimset(options.optim , 'TolFun'      , 1e-8);
  options.optim = optimset(options.optim , 'TolX'        , 1e-8);
  options.optim = optimset(options.optim , 'TolCon'      , 1e-12);
  options.optim = optimset(options.optim , 'Display'     , 'off');
  options.optim = optimset(options.optim , 'Algorithm'   ,'sqp');
  options.optim = optimset(options.optim , 'MaxFunEvals' , 2000*(1+1+1+1));
  options.optim = optimset(options.optim , 'GradObj'     , 'off');
  options.optim = optimset(options.optim , 'GradObj'     , 'on');
  options.optim = optimset(options.optim , 'HessUpdate'  , 'bfgs');
  options.optim = optimset(options.optim , 'MaxIter'     , 10000);

  %Optimization inequality constraints
  options.sumA = [ -1  0  0.0  0  0 ;            % omega > 0
                    0 -1  0.0  0  0; ...              % alpha > -0.01 (slightly relaxed)
                    0  0  0.0 -1  0; ...              % beta > 0
                    0 -1 -1.0  0  0; ...             % alpha + gamma > 0
                    0  1  0.5  1  0;  ...             % alpha + 0.5*gamma + beta < 1
                    zeros(1,4) -1; ...            % phi>0
                    zeros(1,4) 1];                % phi < 10

  options.sumB = [0 -0.01 0 0 0.995 0 10]';
end

%% TODO -- refactor this so that it works with new AND old MATLAB
function starting = getGarchStartingParameters(returnData)
  spec = garchset( 'Distribution' , 'Gaussian'  , 'P', 1, 'Q', 1,...
                   'VarianceModel', 'GJR', 'Display','off');

  [spec, ~, ~, ~, ~] = garchfit(spec, returnData * 100);


  starting_GARCH = [spec.K / 100^2; spec.ARCH; spec.Leverage; spec.GARCH];
  theta = spec.ARCH + 0.5 * spec.Leverage + spec.GARCH;

  starting = [starting_GARCH; 1];
  if theta == 1
    starting(4) = 0.9;
  end
end

function problem = createOptimization(inputData, options, startingParameters)
  problem = createOptimProblem( ...
    'fmincon', 'x0', startingParameters, ...
    'objective', @(x) sgjr.likelihood(x, inputData), ...
    'options', options ...
  );

  problem.Aineq = options.sumA;
  problem.bineq = options.sumB;
end

function startingPoints = createStartingPoints()
  omega_s = [0.01; 0.05];
  alpha_s = [0.01; 0.03];
  gamma_s = [0.06;0.03];
  beta_s = [0.8; 0.9; 0.6];
  c_1s = 1;

  [omega_s, alpha_s, gamma_s, beta_s, c_1s ] = ndgrid(omega_s, alpha_s, gamma_s, beta_s, c_1s);
  pmatrix = [omega_s(:) alpha_s(:) gamma_s(:) beta_s(:) c_1s(:) ];

  startingPoints = CustomStartPointSet(pmatrix);
end

function parameters = estimateSgjr(problem, startingPoints)
  ms = MultiStart;
  ms.StartPointsToRun='bounds-ineqs';
  ms.UseParallel='always';
  ms.Display='off';

  s = matlabpool('size');
  if s==0
    matlabpool open
    parameters = run(ms, problem, startingPoints);
    matlabpool close
  else
    parameters = run(ms, problem, startingPoints);
  end
end

function results = collectResults(parameters, inputData)
  results = result.SGJRResult;

  results.parameters = parameters;
  [ ...
    LL, ~, ~, results.assetVariance, sigma_eAnn, results.smoothedBvDebt, ...
    results.mvEquity, results.leverageMultiplier, results.asset2Debt, ...
    results.deltaImplied, results.volatilityForecast, results.tau ...
  ] = sgjr.likelihood(parameters, inputData);

  results.logLikelihood = LL * -1;
  results.equityVariance = (sigma_eAnn.^2) / 252;

  % Calculate QMLE covariance matrix and robust t-stats
  [results.vcv, VC_mle] = sgjr.robustvcv(@sgjr.likelihood, parameters, inputData);

  results.stdErrorsQMLE = diag(sqrt(results.vcv));

  % Calculate MLE covariance matrix and MLE standard errors
  results.stdErrorsMLE = diag(sqrt(VC_mle));
end