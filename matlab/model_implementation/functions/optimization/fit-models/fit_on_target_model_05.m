function [x, fx] = fit_on_target_model_05(seed)
% function_fit_on_target_model: fit the on-target model
% uses state occupancies from FRET histograms rather from the fitted HMM traces

% initialize the path
init_path();

% load and pre-process data for fitting
mfpt_threshold = inf;
use_fret_histograms = false;
dataset = load_all_data_updated(mfpt_threshold,use_fret_histograms);

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

base_model = load_model('../models/on_target_model_05.json');
base_model.num_substrates = num_substrates;

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
model_func = @(x,n)on_target_model(x,n,sub_model_func,@(x)create_on_target_parameters(x,base_model));

% optimization settings
[lb, ub] = create_optimization_bounds(base_model);

fit_function = @(p,n,y_obs,y_norm)model_fit(p,n,y_obs,y_norm,weights,model_func);
num_runs = 1; % number of fitting attempts within the generate_fits function
num_fitting_runs = 1; % number of PSO runs to perform

% create output folder and delete old file
output_folder = '../output/on_target_model_fits/model_05/';
output_file = join([output_folder, 'model_parameters_', num2str(seed), '.csv']); % delete old file for testing purposes
[status, msg] = mkdir(output_folder);

% get on-target substrate data
fitting_data = on_target_data;

% create y_obs vectors
y_obs = data_vector_func(fitting_data);
y_norm = y_obs.^2 ; y_norm(y_norm == 0) = 1;

% generate the model fits (pass in [] for MFPTs here)
[x, fx] = generate_fits(fit_function,y_obs,y_norm,num_substrates,lb,ub,seed);

% save the fit
save_fit_to_file(x, fx, output_file);
end