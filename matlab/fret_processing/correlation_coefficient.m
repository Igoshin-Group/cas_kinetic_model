function c_xy = correlation_coefficient(x, y)
% correlation_coefficient: Calculate the Pearson correlation between
%                          two signals, x and y.
%   c_xy = correlation_coefficient(x, y)
% Inputs:
%   x, y = signal vectors
% Outputs:
%   c_xy = Pearson correlation

N = length(x);

% Check if the signals are the same length
if N ~= length(y)
    error('Signals must be the same length');
end

% Calculate the Pearson correlation
mu_x = mean(x);
mu_y = mean(y);
std_x = std(x);
std_y = std(y);

c_xy = (1 / (N - 1)) * sum(((x - mu_x) / std_x) .* ((y - mu_y) / std_y));
end
