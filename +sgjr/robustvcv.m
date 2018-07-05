function [VCV, VC, A, B, scores, grossScores] = robustvcv(fun, theta, varargin)
% Compute Robust Variance Covariance matrix numerically, using centered
% differentiation.  Avoids taking second derivatives using results from
% Bollerslev-Wooldridge (1991)
%
%
% USAGE:
%     [VCV,A,B,SCORES,HESS,GROSSSCORES]=robustvcv(FUN,THETA,NW,VARARGIN)
%
% INPUTS:
%     FUN           - Function name ('fun') or function handle (@fun) which will
%                       return the sum of the log-likelihood (scalar) as the 1st output, the individual
%                       log likelihoods (T by 1 vector) as the third
%                       output, and the ANNUALIZED VOLATILITY series as the
%                       fifth output
%     THETA         - Parameter estimates at the optimum, usually from fmin*
%     VARARGIN      - Other inputs to the log-likelihood function, such as data
%
% OUTPUTS:
%     VCV           - Estimated robust covariance matrix (see Bollerslev and Wooldridge 1991)
%     VC            - Estimated MLE covariance matrix
%     A             - A portion of robust covariance;
%     B             - B portion of robust covariance, outer product of
%                     gradients
%     SCORES        - T x num_parameters matrix of scores
%     GROSSSCORES  - Numerical scores (1 by num_parameters) of the objective function, usually for diagnostics
%
% COMMENTS:
%     This function simplifies calculating sandwich covariance estimators for (Q)MLE estimation
%     Additionally, instead of using numerical Hessians, it uses only first
%     derivatives.  It also avoids using the noisy OPG estimator.  For full
%     exposition of the theory, refer to the original Bollerslev and
%     Wooldridge paper

%Theory:
%  For a QMLE Estimator, the asymptotic covariance matrix is given by
%  V = A^(-1)*B*A^(-1) where
%
%  A = Hessian of the log likelihood function.  In practice, this can be
%  alternatively computed using derivatives of the volatility series
%  B = outer product of the scores, or t * cov(scores)
%
%  For an MLE Estimator, the asymptotic covariance matrix is given by
%  V=A^(-1). This is typically more stable than using the numerical Hessian or the OPG

% Copyright: Emil Siriwardane
% esiriwar@stern.nyu.edu
% Date: 10/17/2011

  % convert theta to a column vector, if it isn't
  theta = theta(:);

  k = length(theta);
  h = max(abs(theta * eps^(1/3)), 1e-13);
  h = diag(h);

  %Likelihood function, likelihood series, and annualized volatility series
  %evaluated at the optimum
  [~, ~, like, ~, hopt] = feval(fun, theta, varargin{:});

  t = length(like);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Calculate Scores and Derivatives of Variance Series
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  LLFp = zeros(k, 1);
  LLFm = zeros(k, 1);
  likep = zeros(t, k);
  likem = zeros(t, k);

  hp = zeros(t, k);
  hm = zeros(t, k);

  for i = 1:k
    thetaph = theta + h(:, i);
    [LLFp(i), ~, likep(:, i), ~, hp(:, i)] = feval(fun, thetaph, varargin{:});
    thetamh = theta - h(:, i);
    [LLFm(i), ~, likem(:, i), ~, hm(:, i)]=feval(fun, thetamh, varargin{:});
  end

  %%% De-annualize and convert to a variance series
  hp = hp.^2 / 252;
  hm = hm.^2 / 252;

  %%% Compute numerical derivatives
  scores = zeros(t, k);
  grossScores = zeros(k, 1);
  hchanges = zeros(t, k);

  h = diag(h);
  for i = 1:k
    scores(:, i) = (likep(:, i) - likem(:, i)) ./ (2 * h(i));
    grossScores(i) = (LLFp(i) - LLFm(i)) ./ (2 * h(i));
    hchanges(:, i) = (hp(:, i) - hm(:, i)) ./ (2 * h(i));
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Compute QMLE and MLE covariance matrix
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %%% Calculate A-matrix
  hdivide = sqrt(2) * (hopt.^2 / 252);
  a = hchanges ./ repmat(hdivide, 1, k);
  A = a' * a;

  %%% Calculate B-matrix
  B = scores' * scores;       %Use the OPG estimator for the middle of the sandwhich

  %%% QMLE Covariance matrix
  VCV = A \ B / A;

  %%%% MLE Covariance matrix
  VC = inv(A);
end