function [x, fx] = fit_simultaneous_model_02d(seed, substrate, alpha)

% load and pre-process data for fitting
dataset = load_distal_data(inf);

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
% on_target_data(3) = on_target_enh;

% create off-target data structure
off_idx = [substrate]; % indices of off-target substrates for fitting (substrate 48 has 4bp distal mismatch)
off_target_data(1,length(off_idx)) = dataset.wt(off_idx);
off_target_data(2,length(off_idx)) = dataset.hf1(off_idx);
% off_target_data(3,length(off_idx)) = dataset.enh(off_idx);


% load the model
base_model = load_model('../models/simultaneous/cas-simultaneous-model-02d.json');
num_substrates = base_model.num_substrates;
num_variants = base_model.num_variants;

% generate weights for fitting data points
[on_sub_model, on_data_vector, weight_vector] = create_on_target_sub_model();
[off_sub_model, off_data_vector] = create_off_target_sub_model();

weight_vector = repmat(weight_vector,num_substrates*num_variants,1); % duplicate for each substrate


% create model function handle
model_func = @(x)simultaneous_model(x,base_model,on_sub_model,off_sub_model);

% optimization settings
[lb, ub] = create_optimization_bounds(base_model);


% create output folder and delete old file
output_folder = "../output/fits/simultaneous/model-02d-alt-mfpt/sub_"+num2str(substrate)+"/"+"alpha_"+num2str(alpha)+"/";
output_file = output_folder+"model_parameters_"+num2str(seed)+".csv"; % delete old file for testing purposes
[status, msg] = mkdir(output_folder);

% create data structure
y_obs_on_target = on_data_vector(on_target_data);
on_target_mfpts = zeros(num_variants,1);
for ii = 1:num_variants
    on_target_mfpts(ii) = on_target_data(ii).cleavage_mfpt;
end

y_obs_off_target = [];
off_target_mfpts = zeros(num_variants,num_substrates-1);

for ii = 1:num_variants
    for jj = 1:num_substrates-1
        y_obs_off_target = [y_obs_off_target, off_data_vector(off_target_data(ii,jj))];
        off_target_mfpts(ii,jj) = off_target_data(ii,jj).cleavage_mfpt;
    end
end

% re-format data for fitting
y_obs = [y_obs_on_target, y_obs_off_target];
mfpt_obs = [on_target_mfpts];
for ii = 1:num_variants
    for jj = 1:num_substrates-1
        mfpt_obs = [mfpt_obs; off_target_mfpts(ii,jj)];
    end
end

data.y_obs = y_obs;
data.y_norm = y_obs.^2 ; data.y_norm(data.y_norm == 0) = 1;
data.mfpt_obs = mfpt_obs;

options.lb = lb;
options.ub = ub;
options.seed = seed;
options.algorithm = 'particleswarm';
options.model = base_model;

% generate the model fits (pass in [] for MFPTs here)
fit_function = @(p,data)simultaneous_model_fit(p,data,weight_vector,model_func,alpha);
[x, fx] = generate_fits(fit_function,data,options);

% save the fit
save_fit_to_file(x, fx, output_file);
end

%% Local Functions
function [sub_model_func, data_vector, weight_vector] = create_on_target_sub_model()
% optimization settings
fitting_options.dissociation_constant = true;
fitting_options.normalized_binding = false;
fitting_options.mean_first_passage_time = false;
fitting_options.state_occupancy = true;
fitting_options.dwell_mean = true;
fitting_options.dwell_variance = true;
fitting_options.dwell_skewness = false;
fitting_options.transition_bias = false;
fitting_options.relative_cleavage_rate = false;
fitting_options.non_productive_state = false;
% create sub_model function and data_vector function
sub_model_func = @(pset)sub_model(pset,fitting_options);
data_vector = @(dataset)create_data_vector(dataset,fitting_options);
weight_vector = generate_weight_vector(fitting_options);
end

function [sub_model_func, data_vector] = create_off_target_sub_model()
% optimization settings
fitting_options.dissociation_constant = false;
fitting_options.normalized_binding = true;
fitting_options.mean_first_passage_time = false;
fitting_options.state_occupancy = true;
fitting_options.dwell_mean = true;
fitting_options.dwell_variance = true;
fitting_options.dwell_skewness = false;
fitting_options.transition_bias = false;
fitting_options.relative_cleavage_rate = false;
fitting_options.non_productive_state = false;
% create sub_model function and data_vector function
sub_model_func = @(pset)sub_model(pset,fitting_options);
data_vector = @(dataset)create_data_vector(dataset,fitting_options);
end