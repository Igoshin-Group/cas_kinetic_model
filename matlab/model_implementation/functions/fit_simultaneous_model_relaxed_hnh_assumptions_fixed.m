function [x, fx] = fit_simultaneous_model_relaxed_hnh_assumptions_fixed(seed, substrate, lb, ub, x0, folder)

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

fit_function = @(x) objective_function_relaxed_hnh_assumptions_fixed(x, on_target_data, off_target_data);

% create output folder and delete old file
output_folder = folder+"sub_"+num2str(substrate)+"/";
output_file = output_folder+"model_parameters_"+num2str(seed)+".csv"; % delete old file for testing purposes
[status, msg] = mkdir(output_folder);


% % optimizer options
% options.lb = lb;
% options.ub = ub;
% options.seed = seed;
% options.algorithm = 'simulannealbnd';
% options.x0 = x0;

% optimizer options
options.lb = lb;
options.ub = ub;
options.seed = seed;
options.algorithm = 'particleswarm';
options.x0 = x0;
options.initial_swarm_matrix = x0;
options.max_time = (5 * 60 * 60) - (4 * 60); % 5 hour running time with 4 min margin (to save)

% fit the model
[x, fx] = generate_fits(fit_function,options);

% save the fit
save_fit_to_file(x, fx, output_file);
end