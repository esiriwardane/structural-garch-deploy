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
        matlabpool;
        matlabpool close;
        TF = true;
      catch ME
        if strcmp(ME.message, 'Undefined function or variable ''matlabpool''.')
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