function psets = on_target_parameter_generator(x,n)
% parameter_generator: create *n* model parameter sets
%   psets = parameter_generator(x,n)
% inputs:
%   x = model parameters (vector)
%   n = number of parameter sets
% outputs:
%   psets = model parameter sets (struct array)

% model description:
%   update description for model variant
%   create a new branch for each model to keep things organized

% extract shared parameters
num_shared = 7;
if num_shared > 0
    shared_parameters = x(1:num_shared); x(1:num_shared) = [];    
else
    shared_parameters = [];
end

% reshape the parameter vector
x = reshape(x,[],n);

% create parameter structures
for ii = 1:n
    psets(ii) = parameter_generator_helper(shared_parameters, x(:,ii));
end
end

function p = parameter_generator_helper(shared_parameters, free_parameters)
% parameter_generator_helper: creates a model parameter set
%   p = parameter_generator_helper(x)
% inputs:
%   shared_parameters = model parameters shared across conditions (vector)
%   free_parameters = model parameters that vary across conditions (vector)
% outputs:
%   p = model parameter set (struct)

% create the parameter structure
p = struct();

% set shared parameters (which are fixed across conditions)
p.k5 = shared_parameters(1);
p.k55 = shared_parameters(2);
p.k6 = shared_parameters(3);
p.k66 = shared_parameters(4);

p.k22 = shared_parameters(5);

p.k1 = shared_parameters(6);
p.k11 = shared_parameters(7);

p.k9 = 1e-10; % set to small value (see Dagdas et al. 2017)
p.k10 = 1e-10; % set to small value (see Dagdas et al. 2017)

% set free parameters (which vary across conditions)
p.k2 = free_parameters(1);
p.k7 = free_parameters(2);
p.k8 = free_parameters(3);

p.k3 = free_parameters(4);
p.k33 = free_parameters(5);
p.k4 = free_parameters(6);
p.k44 = free_parameters(7);
p.kclv = free_parameters(8);

% set parameters that are defined *a priori*
% note: check that the residence time is reasonable
% p.k1 = 1;   % PAM-binding rate (1/sec)
% p.k11 = log(2)/0.75;    % PAM-unbinding rate (1/sec)


% set discrimination factors (for on-target substrate)
p.f1 = 1;   % PAM-binding discrimination factor (dimensionless)
p.f11 = 1;  % PAM-unbinding discrimination factor (dimensionless)
p.f2 = 1;   %
p.f22 = 1;
p.f3 = 1;
p.f33 = 1;
p.f4 = 1;
p.f44 = 1;
p.f5 = 1;   
p.f55 = 1;
p.f6 = 1;
p.f66 = 1;
p.f7 = 1;
p.f77 = 1;    % set by detailed balance
p.f8 = 1;
p.f88 = 1;    % set by detailed balance
p.f9 = 1;       % PAM-binding disc. factor, checkpoint conf. (dimensionless)
p.f99 = 1;    % set by detailed balance
p.f10 = 1;      % PAM-binding disc. factor, closed conformation (dimensionless)
p.f1010 = 1;  % set by detailed balance
p.fclv = 1;

% calculate model parameters set by detailed balance
p.k77 = p.k7 * (p.k33*p.k22*p.k5) / (p.k3*p.k2*p.k11);
p.k88 = p.k8 * (p.k44*p.k33*p.k22*p.k5*p.k6) / (p.k4*p.k3*p.k2*p.k55*p.k66);
p.k99 = (p.k11 * p.k55 * p.k9) / (p.k5 * p.k1);
p.k1010 = (p.k10 * p.k66 * p.k55 * p.k11) / (p.k6 * p.k5 * p.k1);

p.f77 = p.f7 * ((p.f33 * p.f22 * p.f5) / (p.f2 * p.f3 * p.f55));
p.f88 = p.f8 * ((p.f44 * p.f33 * p.f22 * p.f5 * p.f6)/(p.f4 * p.f3 * p.f2 * p.f55 * p.f66));
p.f99 = p.f9 * ((p.f11 * p.f55) / (p.f5 * p.f1));
p.f1010 = p.f10 * ((p.f66 * p.f55 * p.f11) / (p.f6 * p.f5 * p.f1));

% p.k9 = p.k99 * ((p.k1 * p.k5) / (p.k55 * p.k11));
% p.k10 = p.k1010 * ((p.k1 * p.k5 * p.k6) / (p.k66 * p.k55 * p.k11));
% p.f9 = p.f99 * ((p.f1 * p.f5) / (p.f55 * p.f11));
% p.f10 = p.f1010 * ((p.f1 * p.f5 * p.f6) / (p.f66 * p.f55 * p.f11));
end