function r = log_random(lb,ub,n)
% log_random: generate a vector of log-random numbers
%   r = log_random(lb,ub,n)
% inputs:
%   lb = lower bound
%   ub = upper bound
%   n = number of samples
% outputs:
%   r = log-random numbers (Nx1 vector)

% create random vector
r = 10 .^ (lb + (ub - lb) .* rand(n,1));
end