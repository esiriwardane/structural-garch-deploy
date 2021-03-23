classdef Pool
  properties
    poolObject
  end

  methods
    function obj = Pool
      if obj.hasMatlabPool
        obj.poolObject = util.pool.MatlabPool();
      else
        obj.poolObject = util.pool.ParPool();
      end
    end

    function TF = hasMatlabPool(obj)
      try
        matlabpool('size');
        TF = true;
      catch ME
        if strcmp(ME.message, 'Undefined function or variable ''matlabpool''.') || ...
            strcmp(ME.message, 'Undefined function ''matlabpool'' for input arguments of type ''char''.')
          TF = false;
        else
          rethrow(ME);
        end
      end
    end

    function close(obj)
      obj.poolObject.close();
    end

    function TF = isOpen(obj)
      TF = obj.poolObject.isOpen();
    end

    function open(obj)
      obj.poolObject.open();
    end
  end
end