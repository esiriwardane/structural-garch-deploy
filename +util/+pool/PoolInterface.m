classdef PoolInterface < handle
  methods (Abstract)
    close(obj)
    % Close a pool

    TF = isOpen(obj)
    % Determine if a pool is open

    open(obj)
    % Open a pool
  end
end