function [y, mfpt, psets_off] = off_target_model(x,base_model,psets_on,sub_model, balance_function)
% on_target_model: compute Cas9 model outputs for on-target substrates
%   y = on_target_model(x,n,sub_model)
% inputs:
%   x = model parameters (vector)
%   off_target_base_model = off-target model definition (struct)
%   on_target_pset = on-target model parameters (struct)
%   sub_model = sub-model function (function handle)
% outputs:
%   y = model output (vector)

% create parameter set structures
psets_off = create_off_target_parameter_sets(x, base_model, psets_on);

% compute model outputs
y = []; mfpt = [];
for ii = 1:base_model.num_variants
    pset = psets_off(ii);
    y = [y; sub_model(pset)];
    mfpt = [mfpt; compute_off_target_mfpt(pset, balance_function)];
end

% flatten the output vector
y = reshape(y.',1,[]);
end

