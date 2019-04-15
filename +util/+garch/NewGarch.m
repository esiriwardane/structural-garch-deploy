classdef NewGarch < util.garch.GarchInterface
  properties
    fitResults
  end

  methods
    function fit(obj, data)
      obj.fitResults = gjr(1, 1).estimate(data, 'Display', 'off');
    end

    function [params, theta] = getStartingParams(obj, data)
      params = [ ...
        obj.fitResults.Constant / 100^2; ...
        cell2mat(obj.fitResults.ARCH)'; ...
        cell2mat(obj.fitResults.Leverage)'; ...
        cell2mat(obj.fitResults.GARCH)' ...
      ];

      theta = sum(cell2mat(obj.fitResults.ARCH)) + ...
              0.5 * sum(cell2mat(obj.fitResults.Leverage)) + ...
              sum(cell2mat(obj.fitResults.GARCH));
    end
  end
end