classdef ParPool < util.pool.PoolInterface
  properties
    poolObj
  end

  methods
    function close(obj)
      try
        delete(gcp('nocreate'));
      catch ME
      end

      delete(obj.poolObj);
    end

    function TF = isOpen(obj)
      TF = ~isempty(gcp('nocreate'));
    end

    function open(obj)
      obj.poolObj = parpool;
    end
  end
end