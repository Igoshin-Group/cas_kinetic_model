clc;
clear;
close all;

%% Load the fits
% Load the distal substrate indices (off_idx)
load("./data/matlab/distal_substrate_indices.mat");
substrates = off_idx;
num_substrates = length(off_idx) + 1; % add 1 for the on-target

C = 1; % [dCas9] (nM)

ranks = [1,2,3,4,5]; % We are using the top 3 on-target models

% Load the model ensemble
[models, energies, scores] = load_model_ensemble(substrates, ranks);
num_models = length(models);
num_variants = size(models{end},1);

score_threshold = 2;


mu = zeros(5,4,3);
sigma = zeros(5,4,3);

for rank_index = 1:length(ranks)
rank = ranks(rank_index);
model = models{rank_index};
psets_on = model(:,1); model(:,1) = [];
psets_off = model;

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
% 5, 6, 7
% off_idx = [2, 5, 6, 7, 36, 37, 38, 39, 48]; % smFRET substrates % 39 must be added back
for ii = 1:length(off_idx)
    off_target_data(1,ii) = dataset.wt(off_idx(ii));
    off_target_data(2,ii) = dataset.hf1(off_idx(ii));
    % off_target_data(3,ii) = dataset.enh(off_idx(ii));
end

% compute off-target dissociation constants
for ii = 1:num_variants
    for jj = 1:length(off_idx)
        off_target_data(ii,jj).dissociation_constant = (1/off_target_data(ii,jj).normalized_binding) * on_target_data(ii).dissociation_constant;
    end
end

dataset = [on_target_data', off_target_data]; % should be equivalent to the fitted dataset

pset_on_wt = psets_on(1);
pset_on_hf1 = psets_on(2);

psets_wt = [pset_on_wt, psets_off(1,:)];
psets_hf1 = [pset_on_hf1, psets_off(2,:)];

[pset_on_clv_swap,~] = swap_cleavage_rate(pset_on_wt, pset_on_hf1);
[pset_on_conf_swap,~] = swap_conformational_transitions(pset_on_wt, pset_on_hf1);

%% Panel A: WT/WT
mismatch_counts = [0, off_target_data(1,:).distal_mismatches];
distal_mismatches = [off_target_data(2,:).distal_mismatches];
figure('Renderer', 'painters', 'Position', [10 10 900 600])
subplot(2,2,1);
z_pred = zeros(num_substrates,1); z_pred(1) = 1;
for ii = 1:num_substrates
    z_pred(ii) = compute_numerical_cleavage_error(psets_wt(ii));
end

err1 = z_pred;
boxplot(z_pred,mismatch_counts);
title("WT Conf. / WT R-Loop","FontSize",18);
ylabel("Error, \eta (unitless)","FontSize",14);
ylim([1e-2,1.5]);
set(gca,"YScale","log");
%% Panel B: WT Conf /HF1 R-Loop
subplot(2,2,2);
z_pred = zeros(num_substrates,1); z_pred(1) = 1;
for ii = 1:num_substrates
    % WT-Cas9 with HF1 R-Loop Transitions
    [pset_swap,tmp] = swap_r_loop_transitions(psets_wt(ii),psets_hf1(ii));
    z_pred(ii) = compute_numerical_cleavage_error(pset_swap);
end

err2 = z_pred;

boxplot(z_pred,mismatch_counts);
title("WT Conf. / HF1 R-Loop","FontSize",18);
ylim([1e-2,1.5]);
set(gca,"YScale","log");
%% Panel C: HF1 Conf /WT R-Loop Transitions
subplot(2,2,3);
z_pred = zeros(num_substrates,1); z_pred(1) = 1;
for ii = 1:num_substrates
    % WT-Cas9 with HF1 HNH dynamics
    [pset_swap,tmp] = swap_conformational_transitions(psets_wt(ii),psets_hf1(ii));
    z_pred(ii) = compute_numerical_cleavage_error(pset_swap);
end

boxplot(z_pred,mismatch_counts);
title("HF1 Conf. / WT R-Loop","FontSize",18);
ylabel("Error, \eta (unitless)","FontSize",14);
xlabel("Num. Distal Mismatches (#)","FontSize",14);
ylim([1e-2,1.5]);
set(gca,"YScale","log");

err3 = z_pred;

%% Panel D: HF1/HF1

subplot(2,2,4);
z_pred = zeros(num_substrates,1); z_pred(1) = 1;
for ii = 1:num_substrates
    z_pred(ii) = compute_numerical_cleavage_error(psets_hf1(ii));
end

boxplot(z_pred,mismatch_counts);
title("HF1 Conf. / HF1 R-Loop","FontSize",18);
xlabel("Num. Distal Mismatches (#)","FontSize",14);

ylim([-inf,1.5]);
set(gca,"YScale","log");

% save_figure(gcf, "../figures/", "pred/", "parameter-swap-boxplot-top-model-"+num2str(rank));
err4 = z_pred;

%% Generate bar plot
errors = [err1(1:end), err2(1:end), err3(1:end), err4(1:end)];
mismatch_counts = [0, off_target_data(1,:).distal_mismatches];

temp_errors = cell(5,4);
error_sum = zeros(5,4);
bin_counts = zeros(5,1);
for ii = 1:length(mismatch_counts)
    num_mismatches = mismatch_counts(ii)+1; % add 1 to mismatch count for proper indexing
    for jj = 1:4
        error_sum(num_mismatches,jj) = error_sum(num_mismatches,jj) + errors(ii,jj);
        temp_errors{num_mismatches,jj}(end+1) = errors(ii,jj);
    end
    bin_counts(num_mismatches) = bin_counts(num_mismatches) + 1;
end

mean_error = error_sum ./ bin_counts;


% compute statistics for ensemble
for count = 1:5
    for cond = 1:4
        mu(count,cond,rank_index) = mean(temp_errors{count, cond});
        sigma(count,cond,rank_index) = std(temp_errors{count, cond});
    end
end

% set errors for on-target, by definition 1 for both WT-Cas9 and HF1-Cas9
% mean_error(1,1) = 1;
% mean_error(1,2) = 1;
% mean_error(1,3) = 1;
% mean_error(1,4) = 1;
% mean_error(1,:) = [];
figure();
bar([0:4],mean_error);
set(gca,"YScale","log");
xlabel("Distal Mismatches (#)");
ylabel("Error");
title("Effects of Conformational Transitions on Error");
legend("WT-Cas9 (native)", "WT HNH / HF1 R-Loop", "HF1 HNH / WT R-Loop", "HF1-Cas9 (native)");
ylim([-inf,1e1]);

fprintf("Min. error for (%i) is %f\n", rank, min(mean_error, [], "all"));

%save_figure(gcf, "../figures/", "pred/", "parameter-swap-bar-plot-"+num2str(rank));

end

%% Compute average statistics for ensemble
mean_errors = zeros(4,4);
sigma_errors = zeros(4,4);
for ii = 2:5
    for jj = 1:4
        values = squeeze(mu(ii,jj,:));
        values(values<0) = nan;
        %values = round(values, 3); % round to 3 decimal
        mean_errors(ii-1,jj) = mean(values,"omitnan");
        sigma_errors(ii-1,jj) = std(values,"omitnan") ./ sqrt(3);
    end
end

%% Plot bargraph with error bars
X = [2, 2, 2, 2; 3, 3, 3, 3; 4, 4, 4, 4; 5, 5, 5, 5] - 1;
Y = mean_errors;
CI = sigma_errors;

figure();
barplot(X,Y,CI);

xlabel("Distal Mismatches (#)");
ylabel("Error");
title("Effects of Conformational Transitions on Error");
legend("WT-Cas9 (native)", "WT HNH / HF1 R-Loop", "HF1 HNH / WT R-Loop", "HF1-Cas9 (native)");
set(gca,"YScale","log");
pbaspect([0.5,1,1]);
save_figure(gcf, "../../figures/", "pred/", "wt-hf1-parameter-swap");