classdef OldSchoolGarch < util.garch.GarchInterface
  properties
    fitResults
  end

  methods
    function fit(obj, data)
      spec = garchset( ...
        'Distribution' , 'Gaussian'  , 'P', 1, 'Q', 1, 'VarianceModel', 'GJR', 'Display', 'off' ...
      );

      [obj.fitResults, ~, ~, ~, ~] = garchfit(spec, data);
    end

    function [params, theta] = getStartingParams(obj, data)
      params = [ ...
        obj.fitResults.K / 100^2; ...
        obj.fitResults.ARCH; ...
        obj.fitResults.Leverage; ...
        obj.fitResults.GARCH ...
      ];

      theta = obj.fitResults.ARCH + 0.5 * obj.fitResults.Leverage + obj.fitResults.GARCH;
    end
  end
end