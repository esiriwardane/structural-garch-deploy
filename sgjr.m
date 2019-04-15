function result = sgjr( ...
  dates, returnData, mvEquity, insuranceReserves, totalDeposits, shortTermDebt, longTermDebt, ...
  otherLiabilities, zeroCurve, varargin ...
)
  % SGJR Structural GARCH code
  %
  % INPUTS:
  %    dates             - A T by 1 vector of dates (either a cell array, or a vector of datenums)
  %    returnData        - A T by 1 vector of zero mean returns
  %    mvEquity          - A T by 1 vector of market value of equity
  %    insuranceReserves - A T by 1 vector of insurance reserves (if the firm is an insurance
  %                        company, otherwise a vector of zeros)
  %    totalDeposits     - A T by 1 vector of insurance reserves (if the firm is a bank,
  %                        otherwise a vector of zeros
  %    shortTermDebt     - A T by 1 vector of the short term and current debt of the firm
  %    longTermDebt      - A T by 1 vector of the long term debt of the firm
  %    otherLiabilities  - A T by 1 vector of the other relevant liabilities of the firm
  %    zerocurve         - Zero curve data, where the first column is a column of dates
  %                        (YYYYMMDD), the second column is a column of maturities (in days), and
  %                        the third column is a column of rates
  % TODO -- check with Emil on the defition of the zerocurve argument and if it's the best format
  % to work with
  %
  % OPTIONS:
  %    forecastType      - 'Constant' (default) || 'Dynamic'
  %
  % OUTPUTS:
  %    result            - An SGJRResult object
  %
  % SEE ALSO:  result/SGJRResult

  % TODO - ask Emil the format of zero curve data
  options = processOptions(varargin{:});

  % TOOD ensure all data are in column vector format
  if iscell(dates)
    dates = datenum(dates);
  end

  debtMaturities = [30 1 2 8 3];

  bvDebtComponents = [ ...
    insuranceReserves, totalDeposits, shortTermDebt, longTermDebt, otherLiabilities ...
  ];

  [dates, returnData, mvEquity, bvDebtComponents] = ...
    removeDaysWithZeroBv(dates, returnData, mvEquity, bvDebtComponents);

  bvDebt = sum(bvDebtComponents, 2);

  % TODO create a better name for this function
  smoothedTau = createSmoothedTau(bvDebtComponents, debtMaturities);
  riskFreeRate = util.getRiskFreeRates(dates, zeroCurve, smoothedTau * 252);

  result = sgjr.estimate(returnData, mvEquity, bvDebt, riskFreeRate, smoothedTau, options.forecastType);
end

function options = processOptions(varargin)
  options = struct(varargin{:});

  if ~isfield(options, 'forecastType')
    options.forecastType = 'Constant';
  elseif ~validateForecastType(options.forecastType)
    error( ...
      'sgjr:invalidForecastType', ...
      'Invalid forecast type - %s - valid forecast types are either ''Constant'' or ''Dynamic''', ...
      options.forecastType ...
    );
  end
end

function [dates, returnData, mvEquity, bvDebt] = ...
      removeDaysWithZeroBv(dates, returnData, mvEquity, bvDebt)
  % Remove days in which book value of debt == 0 from the dataset
  zeroBvIndices = sum(bvDebt, 2) == 0;

  dates(zeroBvIndices) = [];
  returnData(zeroBvIndices) = [];
  mvEquity(zeroBvIndices) = [];
  bvDebt(zeroBvIndices, :) = [];
end

function smoothedTau = createSmoothedTau(bvDebtComponents, debtMaturities)
  bvDebt = sum(bvDebtComponents, 2);
  dWeights = bvDebtComponents ./ repmat(bvDebt, 1, 5);
  T = length(dWeights);

  rawTau = sum(repmat(debtMaturities, T, 1) .* dWeights, 2);
  smoothedTau = util.exponentiallySmooth(rawTau, 0.01);
end
