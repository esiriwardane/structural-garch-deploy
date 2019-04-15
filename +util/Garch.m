classdef Garch
  properties
    garchObj
    originalWarnings
  end

  methods
    function obj = Garch
      obj.originalWarnings = warning ('off','econ:garchfit:FunctionToBeRemoved');
      warning ('off','econ:garchset:FunctionToBeRemoved');
      warning ('off','econ:garchinfer:FunctionToBeRemoved');
      warning ('off','econ:garchpred:FunctionToBeRemoved');

      if obj.hasOldGarch
        obj.garchObj = util.garch.OldSchoolGarch();
      else
        obj.garchObj = util.garch.NewGarch();
      end
    end

    function fit(obj, data)
      obj.garchObj.fit(data);
    end

    function TF = hasOldGarch(obj)
      try
        garchspec();
        TF= true;
      catch ME
        if strcmp(ME.message, 'Undefined function or variable ''garchspec''.')
          TF = false;
        else
          rethrow(ME);
        end
      end
    end

    function [params, theta] = getStartingParameters(obj, data)
      [params, theta] = obj.garchObj.getStartingParameters(data);
    end

    function delete(obj)
      warning(obj.originalWarnings);
    end
  end
end