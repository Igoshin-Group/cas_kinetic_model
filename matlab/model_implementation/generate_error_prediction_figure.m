clc;
clear;
close all;


%
% This is the correct version to use to generate error predictions.
% 2024-09-19
%

init_path();

%% Load the fits
% Load the distal substrate indices
load("./data/matlab/distal_substrate_indices.mat");
substrates = off_idx;
num_substrates = length(substrates);

C = 1; % [dCas9] (nM)

ranks = [1,2,3,4,5]; % We are using the top 3 on-target models

% Load the model ensemble
[models, energies, scores] = load_model_ensemble(substrates, ranks);
num_models = length(models);
num_variants = size(models{end},1);

score_threshold = inf;


best_model = models{1};
on_target_wt = best_model(1,1);


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

%% Get smFRET substrate mask
mask = zeros(num_substrates+1,1);

for ii = 1:num_substrates+1
    mismatch_pattern = dataset(1,ii).mismatch_pattern;
    mismatch_pattern = mismatch_pattern{1};
    distal = mismatch_pattern(1:4);
    % look for the correct mismatch patterns
    if contains(distal, '0000') || contains(distal, '1000') || contains(distal, '1100') || contains(distal, '1110') || contains(distal, '1111')
        mask(ii) = 1;
        fprintf("%s\n", distal);
    else
        continue;
    end
end

%% Compute cleavage error predictions
pred_means = cell(5,2);
substrate_counts = zeros(5,1);

for rank = 1:num_models
    psets = models{rank};
    cur_scores = scores{rank};

clv_speed = @(p) 1 ./ compute_off_target_mfpt(p);
clv_error = @(p) compute_numerical_cleavage_error(p); %@(p) compute_forward_error(p);
mismatch_counts = [dataset(1,:).total_mismatches];
num_variants = 2;
num_conditions = 5;%length(unique(mismatch_counts));
[mismatch_mode,mismatch_frequency] = mode(mismatch_counts);
num_samples = max(mismatch_frequency);
data = nan(num_conditions,num_samples,num_variants);
observation_array = ones(num_conditions,1);

wt_errors = zeros(length(mismatch_counts),2);
hf1_errors = zeros(length(mismatch_counts),2);
for ii = 1:num_substrates+1
    
    pset_wt = psets(1,ii);
    pset_hf1 = psets(2,ii);
    
%     % test kclv similar to on-target
%     pset_wt.fclv = 1;
%     pset_hf1.fclv = 1;

    % get the current mismatch count
    num_mismatches = mismatch_counts(ii);
    
    % check the mask (optional)
%     fprintf("Checking smFRET mask (only use substrates matching the smFRET dataset sequence pattern).\n");
%     if mask(ii) == 0
%         continue;
%     end

    % check if the score is above the threshold
    cur_score = cur_scores(ii);
    if cur_score > score_threshold
        wt_error = nan;
        hf1_error = nan;
        
        wt_errors(ii,:) = nan;
        hf1_errors(ii,:) = nan;
            % get the observation index
        sample_index = observation_array(num_mismatches+1);
        % increment the observation index
        observation_array(num_mismatches+1) = sample_index + 1;
        substrate_counts(num_mismatches+1) = substrate_counts(num_mismatches+1) + 1;
        % store the cleavage errors
        data(num_mismatches+1, sample_index, 1) = wt_error;
        data(num_mismatches+1, sample_index, 2) = hf1_error;
    else 
        % compute the cleavage errors
        wt_error = clv_error(pset_wt);
        hf1_error = clv_error(pset_hf1);

        wt_errors(ii,1) = clv_error(pset_wt);
        wt_errors(ii,2) = compute_cleavage_error(pset_wt);

        hf1_errors(ii,1) = clv_error(pset_hf1);
        hf1_errors(ii,2) = compute_cleavage_error(pset_hf1);

        % get the observation index
        sample_index = observation_array(num_mismatches+1);
        % increment the observation index
        observation_array(num_mismatches+1) = sample_index + 1;
        substrate_counts(num_mismatches+1) = substrate_counts(num_mismatches+1) + 1;
        % store the cleavage errors
        data(num_mismatches+1, sample_index, 1) = wt_error;
        data(num_mismatches+1, sample_index, 2) = hf1_error;
    end
end


pred_mean = mean(data,2,'omitnan'); % average over the samples (axis 2)
pred_mean = reshape(pred_mean,[num_conditions,num_variants]);

for variant = 1:2
    for count = 1:5
        pred_means{count,variant}(end+1) = pred_mean(count,variant);
    end
end

end

pred_mean = zeros(5,2);
pred_std = zeros(5,2);
for variant = 1:2
    for count = 1:5
        pred_mean(count,variant) = (mean(pred_means{count,variant}, "omitnan"));
        pred_std(count,variant) = (std(pred_means{count,variant}, "omitnan"));
    end
end




%% Load the Kim et al. 2020 dataset
% From Kim et al. 2020 (obtained from supplemental tables)
kim_dataset = readtable("./data/kim_2020_processed_distal_data.csv");

kim_dataset_size = height(kim_dataset);
variants = 2;

%% Compute ranges for each mismatch count
cleavage_errors = cell(5,2);
counts = zeros(5,1);
for ii = 1:kim_dataset_size
    row = kim_dataset(ii,:);
    distal_mismatches = row.distal_mismatches;
    index = distal_mismatches + 1;
    cleavage_errors{index,1}(end+1) = row.wt;
    cleavage_errors{index,2}(end+1) = row.hf1;
    counts(index) = counts(index) + 1;
end

mu = zeros(5,2);
sigma = zeros(5,2);
for ii = 1:5
    mu(ii,1) = (mean(cleavage_errors{ii,1}, "omitnan"));
    sigma(ii,1) = (std(cleavage_errors{ii,1}, "omitnan"))./ sqrt(length(cleavage_errors{ii,1}));
   
    mu(ii,2) = (mean(cleavage_errors{ii,2}, "omitnan"));
    sigma(ii,2) = (std(cleavage_errors{ii,2}, "omitnan")) ./ sqrt(length(cleavage_errors{ii,2}));

end
%% Plot test range
width = 0.25;
spacing = 0.15;
figure();
hold on;

face_colors = ['b', 'c'];

for ii = 2:4
    x0 = ii-1;
    y0 = mu(ii);
    
    r = draw_region(x0 - spacing, mu(ii,1), width, sigma(ii,1), [], []);
    r.FaceColor = face_colors(1);      
    r.EdgeColor = 'none';
    
    r = draw_region(x0 + spacing, mu(ii,2), width, sigma(ii,2), [], []);
    r.FaceColor = face_colors(2);      
    r.EdgeColor = 'none';
    
    errorbar(x0 - spacing, pred_mean(ii,1), pred_std(ii,1) ./ sqrt(3), "Marker", "o", "Color", "g", "LineStyle","none");
    errorbar(x0 + spacing, pred_mean(ii,2), pred_std(ii,2) ./ sqrt(3), "Marker", "o", "Color", "r", "LineStyle","none");
end


%% Use Jones data instead
cleavage_errors = cell(5,2);
counts = zeros(5,1);

wt_on_target_mfpt = dataset(1,1).cleavage_mfpt;
hf1_on_target_mfpt = dataset(2,1).cleavage_mfpt;

for ii = 1:length(dataset)
    row = dataset(:,ii);
    wt = row(1);
    hf1 = row(2);
    distal_mismatches = wt.distal_mismatches;
    index = distal_mismatches + 1;
    cleavage_errors{index,1}(end+1) = wt_on_target_mfpt ./ wt.cleavage_mfpt;
    cleavage_errors{index,2}(end+1) = hf1_on_target_mfpt ./ hf1.cleavage_mfpt;
    counts(index) = counts(index) + 1;    
end

%% Finish figure
mu = zeros(5,2);
sigma = zeros(5,2);
for ii = 1:5
    mu(ii,1) = (mean(cleavage_errors{ii,1}, "omitnan"));
    sigma(ii,1) = (std(cleavage_errors{ii,1}, "omitnan")) ./ sqrt(length(cleavage_errors{ii,1}));
    
    mu(ii,2) = (mean(cleavage_errors{ii,2}, "omitnan"));
    sigma(ii,2) = (std(cleavage_errors{ii,2}, "omitnan")) ./ sqrt(length(cleavage_errors{ii,2}));
end

face_colors = ['m', 'y'];
for ii = 2:4
    x0 = ii-1;
    y0 = mu(ii);
    
    r = draw_region(x0 - spacing, mu(ii,1), width, sigma(ii,1), [], []);
    r.FaceColor = face_colors(1);      
    r.EdgeColor = 'none';
    
    r = draw_region(x0 + spacing, mu(ii,2), width, sigma(ii,2), [], []);
    r.FaceColor = face_colors(2);      
    r.EdgeColor = 'none';
    
    %errorbar(x0 - spacing, pred_mean(ii,1), pred_std(ii,1) ./ sqrt(3), "Marker", "*", "Color", "g", "LineStyle","none");
    %errorbar(x0 + spacing, pred_mean(ii,2), pred_std(ii,2) ./ sqrt(3), "Marker", "*", "Color", "r", "LineStyle","none");
end



hold off;

xlabel("Distal Mismatches (#)", "FontSize", 14);
ylabel("Error (unitless)", "FontSize", 14);
%title("Kim et al. 2020 - Cleavage Errors");
legend("WT (Jones 2022)","HF1 (Jones 2022)", "WT (Model)", "HF1 (Model)");
legend("Location","SouthWest");

save_figure(gcf, "../../figures/", "pred/", "error-prediction-linear");

set(gca,"YScale","log");
pbaspect([1,1,1]);

save_figure(gcf, "../../figures/", "pred/", "error-prediction-logarithmic");

%% Local Functions
function r = draw_region(x0,y0,width,height,hatch,angle)
half_width = width/2; 
half_height = height; % should be +/- 1 SD, not +/- 1/2 SD (in log-space, 1/2 SD would be SD-log10(2))
x = [x0 - half_width, x0 - half_width, x0 + half_width, x0 + half_width];
y = [y0 - half_height, y0 + half_height, y0 + half_height, y0 - half_height];

% idx = [x<0];
% x(idx) = 0.01;
% idx = [y<0];
% y(idx) = 0.01;

r = fill(x,y,"g");

if ~isempty(hatch)
    external.hatchfill(r,hatch,angle,5);
end
%plot(x0, y0, "k*");
end

function r = draw_region_ci(x0,y0,width,height)
half_width = width/2; 
x = [x0 - half_width, x0 - half_width, x0 + half_width, x0 + half_width];
y = [height(1), height(2), height(2), height(1)];

% idx = [x<0];
% x(idx) = 0.01;
% idx = [y<0];
% y(idx) = 0.01;

r = fill(x,y,"g");
external.hatchfill(r);

plot(x0, y0, "k*");
end