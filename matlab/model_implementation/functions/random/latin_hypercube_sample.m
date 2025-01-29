function s = latin_hypercube_sample(num_samples,num_parameters,lb,ub,use_log)
% latin_hypercube_sample: generate an LHS sample
%   s = latin_hypercube_sample(num_samples,num_parameters,lb,ub)
% inputs:
%   num_samples = number of sample points
%   num_parameters = number of parameters
%   lb = parameter lower bounds (num_parameter x 1 vector, or scalar)
%   ub = parameter upper bounds (num_parameter x 1 vector, or scalar)
%   use_log = use logarithmic parameter bounds, [10^lb, 10^ub] (default: true)
% outputs:
%   s = LHS sample (num_samples x num_parameters matrix)

% set default parameters
if ~exist("use_log","var"); use_log = true; end

% process the parameter bounds
if length(lb) == 1  % generate a num_parameters-length vector if scalar
    lb = ones(num_parameters,1) .* lb;
else % convert to a column vector otherwise
    lb = lb(:);
end

if length(ub) == 1 % generate a num_parameters-length vector if scalar
    ub = ones(num_parameters,1) .* ub;
else % convert to a column vector otherwise
    ub = ub(:);
end

% generate the LHS sample
s = lhsdesign(num_samples,num_parameters);

% scale to [lb, ub] for each parameter
for ii = 1:num_parameters
    if use_log % logarithmic sampling
        s(:,ii) = 10 .^ (lb(ii) + ((ub(ii) - lb(ii)) .* s(:,ii)));
    else % linear sampling
        s(:,ii) = (lb(ii) + ((ub(ii) - lb(ii)) .* s(:,ii)));
    end
end
end