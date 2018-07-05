function [v] = fastBsVega(S, X, r, sig, T)
  d1 = (log(S ./ X) + (r + 0.5 * sig.^2) * T) / (sig * sqrt(T));
  v = S * sqrt(T) * sgjr.fastNormPdf(d1, 0, 1);
end