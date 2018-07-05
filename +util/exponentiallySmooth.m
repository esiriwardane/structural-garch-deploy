% TODO - Ask Emil if we want to allow uwers to MEX files
function s = exponentiallySmooth(x, phi)
  s = NaN(length(x), 1);
  s(1) = x(1);

  for T = 2:length(x)
    s(T) = phi * x(T) + (1 - phi) * s(T - 1);
  end
end
