function [psets_on, psets_off] = create_variant_parameter_set(x, model)
% create parameter set for mixed on-target and off-target data

% num_variants = model.num_variants;
num_substrates = model.num_substrates;
num_shared = model.num_shared_parameters;
num_free = model.num_free_parameters;
num_parameters = num_shared + (num_free * num_substrates);

if length(x) ~= num_parameters
    error("There is a problem somewhere: num_parameters should equal length(x)");
end

x_on = x(1:num_shared+num_free);
x_off = x(num_shared+num_free+1:end);

on_model = model; on_model.num_substrates = 1;
off_model = model; off_model.num_substrates = off_model.num_substrates - 1;

pset_on = create_on_target_parameters(x_on,on_model);
psets_off = create_off_target_parameters(x_off,off_model,pset_on);

psets = [pset_on, psets_off];
end