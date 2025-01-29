function psets = off_target_parameter_generator(x,on_target_model)
% parameter_generator: create *n* model parameter sets
%   psets = parameter_generator(x,n)
% inputs:
%   x = model parameters (vector)
%   n = number of parameter sets or conditions (int)
%   on_target_model (struct)
% outputs:
%   psets = model parameter sets (struct array)

% model description:
%   update description for model variant
%   create a new branch for each model to keep things organized

% create parameter structure
psets = parameter_generator_helper(x, on_target_model);
end

function p = parameter_generator_helper(free_parameters, on_target_model)
% parameter_generator_helper: creates a model parameter set
%   p = parameter_generator_helper(x)
% inputs:
%   shared_parameters = model parameters shared across conditions (vector)
%   free_parameters = model parameters that vary across conditions (vector)
% outputs:
%   p = model parameter set (struct)

% create the parameter structure (based on on-target model)
p = on_target_model;

% set discrimination factors for *a priori* defined parameters
p.f1 = 1;

p.f9 = 1;
p.f10 = 1;

% set discrimination factors for shared parameters
p.f5 = 1;   
p.f55 = 1;
p.f6 = 1;
p.f66 = 1;

p.f2 = 1;

% set free discrimination factors (for off-target substrate)

p.f22 = free_parameters(1) / p.k22;

p.f33 = free_parameters(2) / p.k33;

p.f44 = free_parameters(3) / p.k44;

p.f3 = free_parameters(4) / p.k3;

p.f4 = free_parameters(5) / p.k4;

p.fclv = free_parameters(6) / p.kclv;

p.f7 = free_parameters(7) / p.k7;

p.f8 = free_parameters(8) / p.k8;

p.f11 = free_parameters(9) / p.k11;

% calculate model parameters set by detailed balance
p.k77 = p.k7 * (p.k33*p.k22*p.k5) / (p.k3*p.k2*p.k11);
p.k88 = p.k8 * (p.k44*p.k33*p.k22*p.k5*p.k6) / (p.k4*p.k3*p.k2*p.k55*p.k66);
p.k99 = (p.k11 * p.k55 * p.k9) / (p.k5 * p.k1);
p.k1010 = (p.k10 * p.k66 * p.k55 * p.k11) / (p.k6 * p.k5 * p.k1);

p.f77 = p.f7 * ((p.f33 * p.f22 * p.f5) / (p.f2 * p.f3 * p.f55));
p.f88 = p.f8 * ((p.f44 * p.f33 * p.f22 * p.f5 * p.f6)/(p.f4 * p.f3 * p.f2 * p.f55 * p.f66));
p.f99 = p.f9 * ((p.f11 * p.f55) / (p.f5 * p.f1));
p.f1010 = p.f10 * ((p.f66 * p.f55 * p.f11) / (p.f6 * p.f5 * p.f1));

%p.k9 = p.k99 * ((p.k1 * p.k5) / (p.k55 * p.k11));
%p.k10 = p.k1010 * ((p.k1 * p.k5 * p.k6) / (p.k66 * p.k55 * p.k11));
%p.f9 = p.f99 * ((p.f1 * p.f5) / (p.f55 * p.f11));
%p.f10 = p.f1010 * ((p.f1 * p.f5 * p.f6) / (p.f66 * p.f55 * p.f11));
end