clc;
clear;
% close all;

init_path();

%% Load data
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
off_idx = [48]; % indices of off-target substrates for fitting (substrate 48 has 4bp distal mismatch)
off_target_data(1,length(off_idx)) = dataset.wt(off_idx);
off_target_data(2,length(off_idx)) = dataset.hf1(off_idx);
% off_target_data(3,length(off_idx)) = dataset.enh(off_idx);

% compute off-target dissociation constants
num_variants = 2;
for ii = 1:num_variants
    for jj = 1:length(off_idx)
        off_target_data(ii,jj).dissociation_constant = off_target_data(ii,jj).normalized_binding * on_target_data(ii).dissociation_constant;
    end
end

num_substrates = 2;
%% Fitting settings
% generate weights for fitting data points
[on_sub_model, on_data_vector, weight_vector] = create_on_target_sub_model();
[off_sub_model, off_data_vector] = create_off_target_sub_model();

weight_vector = repmat(weight_vector,num_substrates*num_variants,1); % duplicate for each substrate

% NOTE: verify that the fitting options are the same as used to fit

%% Load fits
substrates = off_idx;
ranks = [1];
[models, energies, scores] = load_model_ensemble(substrates, ranks);

model = models{1};
psets_on = model(:,1);
psets_off = model(:,2);

pset_on_wt = psets_on(1);
pset_on_hf1 = psets_on(2);

save("./2d_model_parameters.mat", "psets_on");


%% Generate on-target fit figure
close all;
figure("Name","On-Target");
hold on;
colors = ["#527722", "#A2142F", "#0072BD"];
for ii=1:num_variants
    % get predicted and observed data points
    y_pred = on_sub_model(psets_on(ii));
    y_obs = on_data_vector(on_target_data(ii));
    % plot
    plot(y_obs, y_pred, "o", "MarkerFaceColor", colors(ii), "MarkerEdgeColor","k", "LineWidth", 1, "MarkerSize", 6);
end

xlabel("Observed Data", "FontSize", 14);
ylabel("Fitted Data", "FontSize", 14);
title("On-Target Fit Quality", "FontSize", 18);
legend("WT","HF1");
legend("Location", "NorthWest");

set(gca,"XScale", "log", "YScale", "log");


% add reference line for perfect correlation
x_min = 0; x_max = 1e2;
ref = fplot(@(x) 1.25.*x, [x_min, x_max]);
ref.Color = 'k';
ref.LineStyle = '--';
ref.DisplayName = "Reference Line";
ref.HandleVisibility = "off";

ref = fplot(@(x) 0.75.*x, [x_min, x_max]);
ref.Color = 'k';
ref.LineStyle = '--';
ref.DisplayName = "Reference Line";
ref.HandleVisibility = "off";

ref = fplot(@(x) x, [x_min, x_max]);
ref.Color = 'k';
ref.LineStyle = '--';
ref.LineWidth = 2;
ref.DisplayName = "Reference Line";
ref.HandleVisibility = "off";
hold off;

xlim([x_min, x_max]);
ylim([x_min, x_max]);
pbaspect([1,1,1]);

%save_figure(gcf, "../figures/", "qc/", "model-on-target-quality");
create_quality_control_figures_publication(psets_on,on_target_data,"On-Target");

%% Generate off-target fit figure
figure("Name","Off-Target");
hold on;
colors = ["#527722", "#A2142F", "#0072BD"];
for ii=1:num_variants
    % get predicted and observed data points
    y_pred = off_sub_model(psets_off(ii));
    y_obs = off_data_vector(off_target_data(ii));
    % plot
    plot(y_obs, y_pred, "o", "MarkerFaceColor", colors(ii), "MarkerEdgeColor","k", "LineWidth", 1, "MarkerSize", 6);
end

xlabel("Observed Data", "FontSize", 14);
ylabel("Fitted Data", "FontSize", 14);
title("Off-Target Fit Quality", "FontSize", 18);
legend("WT","HF1");
legend("Location", "NorthWest");

set(gca,"XScale", "log", "YScale", "log");


% add reference line for perfect correlation
x_min = 0; x_max = 1e3;
ref = fplot(@(x) 1.25.*x, [x_min, x_max]);
ref.Color = 'k';
ref.LineStyle = '--';
ref.DisplayName = "Reference Line";
ref.HandleVisibility = "off";

ref = fplot(@(x) 0.75.*x, [x_min, x_max]);
ref.Color = 'k';
ref.LineStyle = '--';
ref.DisplayName = "Reference Line";
ref.HandleVisibility = "off";

ref = fplot(@(x) x, [x_min, x_max]);
ref.Color = 'k';
ref.LineStyle = '--';
ref.LineWidth = 2;
ref.DisplayName = "Reference Line";
ref.HandleVisibility = "off";
hold off;

xlim([x_min, x_max]);
ylim([x_min, x_max]);
pbaspect([1,1,1]);

%save_figure(gcf, "../figures/", "qc/", "model-off-target-quality");

% create_quality_control_figures(psets_off,off_target_data,"Off-Target");

%% Clean-up
% switch back to the working directory
cd(working_directory);

%% Local Functions
function [sub_model_func, data_vector, weight_vector] = create_on_target_sub_model()
% optimization settings
fitting_options.dissociation_constant = true;
fitting_options.normalized_binding = false;
fitting_options.mean_first_passage_time = true;
fitting_options.state_occupancy = true;
fitting_options.dwell_mean = true;
fitting_options.dwell_variance = false;
fitting_options.dwell_skewness = false;
fitting_options.dwell_cv = true;
fitting_options.transition_bias = false;
fitting_options.relative_cleavage_rate = false;
fitting_options.non_productive_state = false;
fitting_options.transition_path_cv = false;
% create sub_model function and data_vector function
sub_model_func = @(pset)sub_model(pset,fitting_options);
data_vector = @(dataset)create_data_vector(dataset,fitting_options);
weight_vector = generate_weight_vector(fitting_options);
end

function [sub_model_func, data_vector] = create_off_target_sub_model()
% optimization settings
fitting_options.dissociation_constant = false;
fitting_options.normalized_binding = true;
fitting_options.mean_first_passage_time = true;
fitting_options.state_occupancy = true;
fitting_options.dwell_mean = true;
fitting_options.dwell_variance = false;
fitting_options.dwell_skewness = false;
fitting_options.dwell_cv = true;
fitting_options.transition_bias = false;
fitting_options.relative_cleavage_rate = false;
fitting_options.non_productive_state = false;
fitting_options.transition_path_cv = false;
% create sub_model function and data_vector function
sub_model_func = @(pset)sub_model(pset,fitting_options);
data_vector = @(dataset)create_data_vector(dataset,fitting_options);
end