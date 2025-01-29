clc;
clear;
close all;

working_directory = pwd(); % get the script directory
switch_to_project_directory(); % switch to the code directory

%% Load the fits
% Load the distal substrate indices
load("../data/matlab/distal_substrate_indices.mat");
substrates = off_idx;
num_substrates = length(substrates);

C = 1; % [dCas9] (nM)

ranks = [1, 2, 3, 4, 5]; % We are using the top 3 on-target models

% Load the model ensemble
[models, energies, scores] = load_model_ensemble(substrates, ranks);
num_models = length(models);
num_variants = size(models{end},1);

score_threshold = 2;

rank = 1;

best_model = models{rank};
best_energies = energies{rank};

psets_on = best_model(:,1);
psets_off = best_model(:,[2:end]);

% pset_on_wt = psets_on(1);
% pset_on_hf1 = psets_on(2);

%save("../model_parameters.mat", "pset_on_wt", "pset_on_hf1");


%% Load dataset
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
for ii = 1:length(substrates)
    off_target_data(1,ii) = dataset.wt(substrates(ii));
    off_target_data(2,ii) = dataset.hf1(substrates(ii));
    % off_target_data(3,ii) = dataset.enh(off_idx(ii));
end

% compute off-target dissociation constants
for ii = 1:num_variants
    for jj = 1:length(off_idx)
        off_target_data(ii,jj).dissociation_constant = (1/off_target_data(ii,jj).normalized_binding) * on_target_data(ii).dissociation_constant;
    end
end

dataset = [on_target_data', off_target_data]; % should be equivalent to the fitted dataset


%% Compute sensitivities

parameters = ["k1a","k11a","k2", "k22", "k3", "k33", "k4", "k44", "k5", "k55", "k6", "k66", "k7", "k8", "kclv"];

V_outputs = zeros(length(parameters), 2);
for ii = 1:length(parameters)
    cur_parameter = parameters(ii);
    for jj = 1:2
        V_outputs(ii,jj) = logarithmic_gains(psets_on(jj), cur_parameter, @compute_on_target_speed);
    end
end

figure();
x_labels = parameters;
X = categorical(x_labels);
X = reordercats(X, x_labels);
bar(X,V_outputs);

% heatmap(["WT","HF1"],parameters,outputs);

xlabel("Parameter");
ylabel("Sensitivity (log-gain)");
legend("WT","HF1");
title("On-Target Speed Sensitivity Analysis");
pbaspect([2,1,1]);
ylim([-1,1]);
saveas(gcf,"../best_model_sensitivity_on_target_speed.eps","epsc");

%%
S_outputs = zeros(length(parameters), 2);
for ii = 1:length(parameters)
    cur_parameter = parameters(ii);
    for jj = 1:2
        S_outputs(ii,jj) = specificity_logarithmic_gains(psets_off(jj,:), cur_parameter);
    end
end

figure();
x_labels = parameters;
X = categorical(x_labels);
X = reordercats(X, x_labels);
bar(X,S_outputs);

% heatmap(["WT","HF1"],parameters,outputs);

xlabel("Parameter");
ylabel("Sensitivity (log-gain)");
legend("WT","HF1");
title("Off-Target Specificity Sensitivity Analysis");
%% Local Functions
function average_error = compute_average_error(psets)
    errors = nan(length(psets),1);
    for ii = 1:length(psets)
        errors(ii) = compute_cleavage_error(psets(ii));
    end
    average_error = mean(errors);
end

function S = compute_average_specificity(psets)
    S = 1 - compute_average_error(psets);
end

function V = compute_on_target_speed(pset)
    V = 1 ./ compute_off_target_mfpt(pset);
end

function lg = specificity_logarithmic_gains(models, parameter)
% determine the initial points
x0 = models(1).(parameter);
property = @compute_average_specificity;
y0 = property(models);
N = length(y0);
lg = zeros(N,1);

for ii = 1:N
    % compute the derivative at (x0, y0)
    h = 0.01 * x0; % set the step-size
    
    for jj = 1:length(models)
        % create the parameter sets for each point (central difference, 8th
        % order accuracy)
        x1 = models; x1(jj).(parameter) = x1(jj).(parameter) - 4*h;
        x2 = models; x2(jj).(parameter) = x2(jj).(parameter) - 3*h;
        x3 = models; x3(jj).(parameter) = x3(jj).(parameter) - 2*h;
        x4 = models; x4(jj).(parameter) = x4(jj).(parameter) - 1*h;
        x5 = models; x5(jj).(parameter) = x5(jj).(parameter) + 1*h;
        x6 = models; x6(jj).(parameter) = x6(jj).(parameter) + 2*h;
        x7 = models; x7(jj).(parameter) = x7(jj).(parameter) + 3*h;
        x8 = models; x8(jj).(parameter) = x8(jj).(parameter) + 4*h;        
    end

    % compute the output for each point
    y1 = property(x1);
    y2 = property(x2);
    y3 = property(x3);
    y4 = property(x4);
    y5 = property(x5);
    y6 = property(x6);
    y7 = property(x7);
    y8 = property(x8);
    
    % compute the approximate derivative using finite differences
    D1 = (1/280) * y1(ii);
    D2 = -(4/105) * y2(ii);
    D3 = (1/5) * y3(ii);
    D4 = -(4/5) * y4(ii);
    D5 = (4/5) * y5(ii);
    D6 = -(1/5) * y6(ii);
    D7 = (4/105) * y7(ii);
    D8 = -(1/280) * y8(ii);
    D = (D1 + D2 + D3 + D4 + D5 + D6 + D7 + D8) / (h);

    % compute log-gain (dLogP/dLogx = x/P * dP/dx)
    lg(ii) = (x0 ./ y0(ii)) * D;   
end
end

function lg = logarithmic_gains(model, parameter, property)
% determine the initial points
x0 = model.(parameter);
y0 = property(model);
N = length(y0);
lg = zeros(N,1);
for ii = 1:N
    % compute the derivative at (x0, y0)
    h = 0.01 * x0; % set the step-size
    
    % create the parameter sets for each point (central difference, 8th
    % order accuracy)
    x1 = model; x1.(parameter) = x1.(parameter) - 4*h;
    x2 = model; x2.(parameter) = x2.(parameter) - 3*h;
    x3 = model; x3.(parameter) = x3.(parameter) - 2*h;
    x4 = model; x4.(parameter) = x4.(parameter) - 1*h;
    x5 = model; x5.(parameter) = x5.(parameter) + 1*h;
    x6 = model; x6.(parameter) = x6.(parameter) + 2*h;
    x7 = model; x7.(parameter) = x7.(parameter) + 3*h;
    x8 = model; x8.(parameter) = x8.(parameter) + 4*h;
    
    % compute the output for each point
    y1 = property(x1);
    y2 = property(x2);
    y3 = property(x3);
    y4 = property(x4);
    y5 = property(x5);
    y6 = property(x6);
    y7 = property(x7);
    y8 = property(x8);
    
    % compute the approximate derivative using finite differences
    D1 = (1/280) * y1(ii);
    D2 = -(4/105) * y2(ii);
    D3 = (1/5) * y3(ii);
    D4 = -(4/5) * y4(ii);
    D5 = (4/5) * y5(ii);
    D6 = -(1/5) * y6(ii);
    D7 = (4/105) * y7(ii);
    D8 = -(1/280) * y8(ii);
    D = (D1 + D2 + D3 + D4 + D5 + D6 + D7 + D8) / (h);

    % compute log-gain (dLogP/dLogx = x/P * dP/dx)
    lg(ii) = (x0 ./ y0(ii)) * D;   
end
end