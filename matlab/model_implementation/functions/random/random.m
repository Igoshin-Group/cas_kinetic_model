function r = random(lb,ub,n)
% random: generate a vector of random numbers
%   r = random(lb,ub,n)
% inputs:
%   lb = lower bound
%   ub = upper bound
%   n = number of samples
% outputs:
%   r = random numbers (Nx1 vector)

% create random vector
r = (lb + (ub - lb) .* rand(n,1));
end