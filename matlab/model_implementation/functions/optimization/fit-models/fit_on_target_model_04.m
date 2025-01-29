function [x, fx] = fit_on_target_model_04(seed, alpha)
% function_fit_on_target_model: fit the on-target model

% load and pre-process data for fitting
dataset = load_all_data();

% extract on-target substrate data
on_target_wt = dataset.wt(1);
on_target_hf1 = dataset.hf1(1);
on_target_enh = dataset.enh(1);

% add dissociation constants (from Chen et al. 2017, see Figure 1b)
on_target_wt.dissociation_constant = 14.23;
on_target_hf1.dissociation_constant = 17.80;
on_target_enh.dissociation_constant = 12.70;

% create on-target data structure
on_target_data(1) = on_target_wt;
on_target_data(2) = on_target_hf1;
on_target_data(3) = on_target_enh;

num_variants = length(on_target_data);
num_substrates = num_variants;

base_model = load_model('../models/on-target/on-target-model-04.json');

% optimization settings
fitting_options.dissociation_constant = true;
fitting_options.normalized_binding = false;
fitting_options.mean_first_passage_time = true;
fitting_options.state_occupancy = true;
fitting_options.dwell_mean = true;
fitting_options.dwell_variance = true;
fitting_options.dwell_skewness = false;
fitting_options.transition_bias = false;
fitting_options.relative_cleavage_rate = false;
fitting_options.non_productive_state = false;

weights = generate_weight_vector(fitting_options);
weights = repmat(weights,num_substrates,1); % duplicate for each substrate

sub_model_func = @(pset)sub_model(pset,fitting_options);
data_vector_func = @(dataset)create_data_vector(dataset,fitting_options);

% create model function handle
model_func = @(x)on_target_model(x,base_model,sub_model_func);

% optimization settings
[lb, ub] = create_optimization_bounds(base_model);

fit_function = @(p,data)on_target_model_fit(p,data,weights,model_func,alpha);

% create output folder and delete old file
output_folder = "../output/fits/on-target/model-04/"+"alpha_"+num2str(alpha)+"/";
output_file = output_folder+"model_parameters_"+num2str(seed)+".csv"; % delete old file for testing purposes
[status, msg] = mkdir(output_folder);

% get on-target substrate data
fitting_data = on_target_data;
mfpt_obs = zeros(num_variants,1);
for ii = 1:num_variants
    mfpt_obs(ii) = on_target_data(ii).cleavage_mfpt;
end


% create y_obs vectors
y_obs = data_vector_func(fitting_data);
y_norm = y_obs.^2 ; y_norm(y_norm == 0) = 1;

data.y_obs = y_obs;
data.y_norm = y_norm;
data.mfpt_obs = mfpt_obs;

% generate the model fits (pass in [] for MFPTs here)
[x, fx] = generate_fits(fit_function,data,lb,ub,seed);

% save the fit
save_fit_to_file(x, fx, output_file);
end