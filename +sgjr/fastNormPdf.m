function [v] = fastNormPdf(x, mu, sd)
  var = sd.^2;
  pi = 3.1415926;
  denom = (2 * pi * var).^.5;
  num = exp(-(x - mu).^2 ./ (2 * var));
  v = num ./ denom;
end