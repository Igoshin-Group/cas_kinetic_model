function [y, mfpt, psets_on, psets_off] = simultaneous_model(x, model, on_sub_model, off_sub_model, balance_function)
% cas_target_model: compute Cas9 model outputs for combined on-target and
% off-target substrates
%   y = cas_target_model(x, base_model, on_sub_model, off_sub_model)
% inputs:
%   x = model parameters (vector)
%   model = model structural information (struct)
%   on_sub_model = function handle for computing model observables
%   off_sub_model = function handle for computing model observables
% outputs:
%   y = model output (vector)

num_variants = model.num_variants;
num_substrates = model.num_substrates;

mfpts_on = zeros(num_variants,1);
mfpts_off = zeros(num_variants,num_substrates-1);

% create parameter set structures
[psets_on, psets_off] = create_model_parameter_sets(x, model);

y = []; % container for model outputs

% compute on-target model
for ii = 1:num_variants
    y = [y; on_sub_model(psets_on(ii))];
    mfpts_on(ii) = compute_off_target_mfpt(psets_on(ii), balance_function);
end

% compute off-target model
for ii = 1:num_variants
    for jj = 1:num_substrates-1
        y = [y; off_sub_model(psets_off(ii,jj))];
        mfpts_off(ii,jj) = compute_off_target_mfpt(psets_off(ii,jj), balance_function);
    end
end

% predictions are loaded into 'y' as [VAR1_ON, VAR2_ON,
%                                     VAR1_OFF1, VAR1_OFF2,
%                                     VAR2_OFF1, VAR2_OFF2]

% re-format mfpt data
mfpt= [mfpts_on];
for ii = 1:num_variants
    for jj = 1:num_substrates-1
        mfpt = [mfpt; mfpts_off(ii,jj)];
    end
end

% global index
% if mod(index,50) == 0 % display only every 50 iterations
%     fprintf("MFPT (calc): %d %d %d %d %d %d\n",mfpt);
% end
% index = index + 1; % increment global counter
end