function [y, mfpt, psets_on] = on_target_model(x, model, sub_model, balance_function)
% on_target_model: compute Cas9 model outputs for on-target substrates
%   y = on_target_model(x,n,sub_model)
% inputs:
%   x = model parameters (vector)
%   n = number of parameter sets or conditions
%   sub_model
% outputs:
%   y = model output (vector)

num_variants = model.num_variants;
num_substrates = model.num_substrates;

mfpt = zeros(num_variants,1);

% create parameter set structures
% create parameter set structures
[psets_on, ~] = create_model_parameter_sets(x, model);

y = []; % container for model outputs

% compute on-target model
for ii = 1:num_variants
    y = [y; sub_model(psets_on(ii))];
    mfpt(ii) = compute_off_target_mfpt(psets_on(ii), balance_function);
end

% predictions are loaded into 'y' as [VAR1_ON, VAR2_ON,
%                                     VAR1_OFF1, VAR1_OFF2,
%                                     VAR2_OFF1, VAR2_OFF2]

% % flatten the output vector
% y = reshape(y.',2,[]);
end
