classdef MatlabPool < util.pool.PoolInterface
  methods
    function close(obj)
      matlabpool close;
    end

    function TF = isOpen(obj)
      TF = ~isempty(gcp);
    end

    function open(obj)
      matlabpool open;
    end
  end
end