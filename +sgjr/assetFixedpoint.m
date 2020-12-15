function [assetFp, fval] = assetFixedpoint( ...
  bvDebtMA, riskFreeRate, time, holdingPeriodVariance, mvEquityAug, x ...
)
%%%%% Helper function to solve fixed point problem in Structural GARCH model of Engle and Siriwardane (2014) %%%%%
% VOLATLITY NEEDS TO BE OVER LIFE OF OPTION ALREADY

  %%% Adjusted inputs %%%%
  volatility = sqrt(holdingPeriodVariance);
  subtract = mvEquityAug;
  strike = bvDebtMA;

  evaluateFunctionAtPrice = ...
    @(price) evaluateFunction(price, strike, riskFreeRate, time, volatility, subtract);

  %%%% Initialization %%%%
  tolerance = 1e-6;

  %%%%%% Evalaute starting guess %%%%%%
  fx = evaluateFunctionAtPrice(x);
  if fx==0
    assetFp=x;
    fval = fx;
    return;
  end

  [assetFp, a, fa, fb] = findChangeOfSign(x, fx, evaluateFunctionAtPrice);
  if isnan(assetFp)
    assetFp = NaN;
    fval = NaN;
    return;
  end

  fc = fb;
  [c, e, d] = deal(Inf);

  % Main loop, exit from middle of the loop
  while fb ~= 0 && a ~= assetFp
    % Insure that assetFp is the best result so far, a is the previous
    % value of assetFp, and c is on the opposite side of the zero from assetFp.
    if sign(fb) == sign(fc)
      c = a;
      fc = fa;
      d = assetFp - a;
      e = d;
    end

    if abs(fc) < abs(fb)
      a = assetFp;
      assetFp = c;
      c = a;
      fa = fb;
      fb = fc;
      fc = fa;
    end

    % Convergence test and possible exit
    toler = 2.0 * tolerance * max(abs(assetFp), 1.0);
    m = 0.5 * (c - assetFp);
    if (convergedWithinTolerance(m, fb, toler))
      break
    end

    % Choose bisection or interpolation
    if (abs(e) < toler) || (abs(fa) <= abs(fb))
      [d, e] = bisectValue(m);
    else
      [d, e] = interpolateValue(m, assetFp, a, c, d, e, fa, fb, fc, toler);
    end

    % Next point
    a = assetFp;
    fa = fb;
    if abs(d) > toler
      assetFp = assetFp + d;
    elseif assetFp > c
      assetFp = assetFp - toler;
    else
      assetFp = assetFp + toler;
    end

    %%%%% Function evaluation %%%%%
    fb = evaluateFunctionAtPrice(assetFp);
  end % Main loop

  fval = fb; % assetFp is the best value
end

% TODO -- name this better -- evaluateBlackScholes?
function fx = evaluateFunction(price, strike, riskFreeRate, time, volatility, subtract)
  d1 = (log(price/strike) + (riskFreeRate*time+volatility^2/2))/sqrt(volatility^2);
  d2 = d1 - sqrt(volatility^2);

  normcdf_d1 = 0.5 * erfc(-d1 ./ sqrt(2));
  normcdf_d2 = 0.5 * erfc(-d2 ./ sqrt(2));

  fx =  price * normcdf_d1 - normcdf_d2 * strike * exp(-riskFreeRate * time) - subtract;
end

function dx = determineDx(x)
  if x ~= 0
    dx = x/50;
  else
    dx = 1/50;
  end
end

% TODO - find a better name for this function?
function [assetFp, a, fa, fb] = findChangeOfSign(x, fx, evaluateFunctionAtPrice)
  dx = determineDx(x);

  twosqrt = sqrt(2);
  a = x;
  fa = fx;
  fb = fx;
  assetFp = x;

  while sign(fa) == sign(fb)
    dx = twosqrt * dx;
    a = max([x - dx, 0]);

    fa = evaluateFunctionAtPrice(a);
    if ~isfinite(fa) || ~isreal(fa) || ~isfinite(a)
      assetFp = NaN;
      return
    end

    assetFp = x + dx;

    fb = evaluateFunctionAtPrice(assetFp);
    if ~isfinite(fb) || ~isreal(fb) || ~isfinite(assetFp)
      assetFp = NaN;
      return
    end
  end
end

function TF = convergedWithinTolerance(m, fb, toler)
  TF = (abs(m) <= toler) || (fb == 0.0);
end

function [d, e] = bisectValue(m)
  [d, e] = deal(m);
end

function [p, q] = linearInterpolation(s, m)
  p = 2.0 * m * s;
  q = 1.0 - s;
end

function [p, q] = inverseQuadraticInterpolation(s, m, assetFp, a, fa, fb, fc)
  q = fa / fc;
  r = fb / fc;

  p = s * (2.0 * m * q * (q - r) - (assetFp - a) * (r - 1.0));
  q = (q - 1.0) * (r - 1.0) * (s - 1.0);
end

function TF = interpolatedPointIsAcceptable(m, p, q, e, toler)
  TF = (2.0 * p < 3.0 * m * q - abs(toler * q)) && (p < abs(0.5 * e * q));
end

function [d, e] = interpolateValue(m, assetFp, a, c, d, e, fa, fb, fc, toler)
  s = fb / fa;

  if (a == c)
    [p, q] = linearInterpolation(s, m);
  else
    [p, q] = inverseQuadraticInterpolation(s, m, assetFp, a, fa, fb, fc);
  end

  if p > 0
    q = -q;
  else
    p = -p;
  end

  if (interpolatedPointIsAcceptable(m, p, q, e, toler))
    e = d;
    d = p / q;
  else
    [d, e] = deal(m);
  end
end
