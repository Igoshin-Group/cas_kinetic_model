function y = cas_target_model(x, base_model, on_sub_model, off_sub_model)
% cas_target_model: compute Cas9 model outputs for combined on-target and
% off-target substrates
%   y = cas_target_model(x, base_model, on_sub_model, off_sub_model)
% inputs:
%   x = model parameters (vector)
%   base_model =
%   on_sub_model =
%   off_sub_model =
% outputs:
%   y = model output (vector)

num_substrates = base_model.num_substrates; % fitting multiple substrates at a time

% create parameter set structures
psets = create_variant_parameter_set(x, base_model);
pset_on = psets(1);
pset_off = psets(2:end);

% compute on-target model
y = on_sub_model(pset_on);

% compute off-target model
for ii = 1:num_substrates-1
    y = [y; off_sub_model(pset_off(ii))];
end

% flatten the output vector
y = reshape(y.',1,[]);
end