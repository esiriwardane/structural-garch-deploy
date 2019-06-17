classdef NewGarch < util.garch.GarchInterface
  properties
    fitResults
  end

  methods (Static)
    function result = getParam(param)
      if (isempty(param))
        result = [0];
      else
        result = cell2mat(param)';
      end
    end
  end

  methods
    function fit(obj, data)
      obj.fitResults = gjr(1, 1).estimate(data, 'Display', 'off');
    end

    function [params, theta] = getStartingParams(obj, data)
      params = [ ...
        obj.fitResults.Constant / 100^2; ...
        obj.getParam(obj.fitResults.ARCH); ...
        obj.getParam(obj.fitResults.Leverage); ...
        obj.getParam(obj.fitResults.GARCH) ...
      ];

      theta = sum(cell2mat(obj.fitResults.ARCH)) + ...
              0.5 * sum(cell2mat(obj.fitResults.Leverage)) + ...
              sum(cell2mat(obj.fitResults.GARCH));
    end
  end
end