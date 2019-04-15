classdef GarchInterface < handle
  properties (Abstract)
    fitResults
  end

  methods (Abstract)
    fit(obj, data)
    % Fit the GARCH model

    [params, theta] = getStartingParams(obj, data)
    % get the starting parameters for the SGJR model
  end

  methods
    function [params, theta] = getStartingParameters(obj, data)
      if (nargin < 2 && isempty(obj.fitResults))
        error( ...
          'SGJR:StartingParamters:NeedsEstimate', ...
          'You have to estimate the GARCH model before getting the starting parameters' ...
        );
      elseif isempty(obj.fitResults)
        obj.fit(data);
      end

      [params, theta] = obj.getStartingParams(data);
    end
  end
end