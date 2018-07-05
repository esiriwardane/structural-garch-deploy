function riskFreeRates = getRiskFreeRates(targetDates, zeroCurve, maturity)
  if isscalar(maturity)
    maturity = repmat(maturity, length(targetDates), 1);
  end

  %Convert dates to serial dates
  % TODO -- ask Emil if we want to allow users to MEX files?
  % coder.extrinsic ('datenum','datestr','str2num','unique','num2str')
  riskFreeRates = NaN * ones(length(targetDates), 1);

  uniqueDates = datenum(num2str(unique(zeroCurve(:, 1))), 'yyyymmdd');
  defaultClosestDate = datenum(num2str(zeroCurve(1, 1)), 'yyyymmdd');

  for t = 1:length(targetDates)
    closestDate = findClosestDate(uniqueDates, targetDates(t), defaultClosestDate);

    [nearestBelow, nearestAbove] = findNearestCurvePoints(zeroCurve, closestDate, maturity(t));

    if nearestBelow.maturity == maturity(t)
        riskFreeRates(t, 1) = nearestBelow.rate;
    elseif isempty(nearestAbove.maturity)
        nearestAbove.index = find((zeroCurve(:, 1) == closestDate), 1, 'last');
        riskFreeRates(t, 1) = zeroCurve(nearestAbove.index, 3) / 100;
    else
        riskFreeRates(t, 1) = interp1( ...
          [nearestBelow.maturity; nearestAbove.maturity], ...
          [nearestBelow.rate; nearestAbove.rate], ...
          maturity(t) ...
        );
    end

  end
end

function closestDate= findClosestDate(dates, targetDate, default)
  closestDate = dates(find(dates <= targetDate, 1, 'last'));

  if isempty(closestDate)
    closestDate = default;
  end

  closestDate = str2double(datestr(closestDate, 'yyyymmdd'));
end

function [nearestBelow, nearestAbove] = findNearestCurvePoints(zeroCurve, closestDate, maturity)
  nearestBelow = ...
      findNearest(zeroCurve, closestDate, maturity, @(curve, maturity) curve <= maturity, 'last');

  nearestAbove = ...
      findNearest(zeroCurve, closestDate, maturity, @(curve, maturity) curve >= maturity, 'first');
end

function nearest = findNearest(zeroCurve, date, maturity, compareMaturities, firstOrLast)
  nearest.index = find( ...
    (zeroCurve(:, 1) == date) & (compareMaturities(zeroCurve(:, 2), maturity)), 1, firstOrLast ...
  );

  nearest.maturity = zeroCurve(nearest.index, 2);
  nearest.rate = zeroCurve(nearest.index, 3) / 100;
end
