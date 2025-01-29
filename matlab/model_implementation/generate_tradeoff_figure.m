clc;
clear;
close all;

working_directory = pwd(); % get the script directory
switch_to_project_directory(); % switch to the code directory

%% Load the fits
% Load the distal substrate indices
load("../data/matlab/distal_substrate_indices.mat");
substrates = off_idx;

C = 1; % [dCas9] (nM)

ranks = [1,2,3,4,5]; % We are using the top 3 on-target models

% Load the model ensemble
[models, energies, scores] = load_model_ensemble(substrates, ranks);
num_models = length(models);
num_variants = size(models{end},1);

score_threshold = 2;

%% Compute native specificity and speed for each substrate
native_specificity = zeros(num_variants, num_models);
native_speed = zeros(num_variants, num_models);

for ii = 1:num_models
    native_specificity(:,ii) = compute_average_specificity(models{ii});
    native_speed(:,ii) = compute_speed(models{ii});
end

mean_specificity = mean(native_specificity, 2);
mean_speed = mean(native_speed, 2);


%% Compute speed-specificity for perturbed parameters
parameters = ["E2R_dag","E3R_dag","E4R_dag","E5R_dag","E6R_dag","E7R_dag","E8R_dag","E9R_dag","EclvR_dag"];
num_parameters = length(parameters);

%range = linspace(0.8, 1.2, 200); % range 75% to 125% of fitted energy
range = linspace(-2, 2, 200); % range +/- 2 KT

num_points = length(range);

specificity = cell(num_models,1);
speed = cell(num_models,1);

error_function = @compute_numerical_cleavage_error;
mfpt_function = @compute_numerical_mfpt;

for ii = 1:num_models
    psets = models{ii};
    cur_scores = scores{ii};
    psets_on = psets(:,1); psets(:,1) = [];
    psets_off = psets;
    
    num_substrates = size(psets_off,2);
    
    model_energies = energies{ii};
    
    specificities = nan(num_variants, num_substrates, num_parameters, num_points);
    speeds = nan(num_variants, num_parameters, num_points);
    
    mean_specificity = zeros(num_variants, num_parameters, num_points);
    mean_speed = zeros(num_variants, num_parameters, num_points);
    
    for jj = 1:num_parameters
        parameter = parameters(jj);
        for kk = 1:num_variants
            for ll = 1:num_substrates
               
%                 % skip poorly fitted substrates for now
%                 cur_score = cur_scores(ll);
%                 if cur_score > score_threshold; continue; end
                
                % extract parameter sets and energies for the current
                % substrate
                model_on = psets_on(kk);
                model_off = psets_off(kk,ll);
                cur_energies = model_energies(kk,ll);
                base_energy = cur_energies.(parameter); % the fitted energy value
               
                for qq = 1:num_points
                    value = range(qq);
                    cur_energies.(parameter) = base_energy + value; % add the perturbation to the base energy
                    
                    % re-compute model parameters based on the new energy
                    % landscape
                    model_on = adjust_model_energies(model_on, cur_energies, C);
                    model_off = adjust_model_energies(model_off, cur_energies, C);
                    
                    % compute model outputs
                    specificities(kk, ll, jj, qq) = 1 - error_function(model_off);
                    speeds(kk, jj, qq) = 1 ./ mfpt_function(model_on);
               end
            end  
            mean_specificity(kk, jj, :) = mean(specificities(kk, :, jj, :), "omitnan");
            mean_speed(kk, jj, :) = speeds(kk, jj, :);            
        end
        
    end
    speed{ii} = mean_speed;
    specificity{ii} = mean_specificity;
end

% compute mean and SEM for the ensemble
mu_speed = zeros(num_variants, num_parameters, num_points);
mu_specificity = zeros(num_variants, num_parameters, num_points);
CI_speed = zeros(num_variants, num_parameters, num_points);
CI_specificity = zeros(num_variants, num_parameters, num_points);

for ii = 1:num_variants
    for jj = 1:num_parameters
        for kk = 1:num_points
            model_speeds = zeros(num_models, num_variants, num_parameters, num_points);
            model_specificities = zeros(num_models, num_variants, num_parameters, num_points);
            for ll = 1:num_models
                model_speeds(ll, :, :, :) = speed{ll};
                model_specificities(ll, :, :, :) = specificity{ll};
            end
            mu_speed(ii,jj,kk) = mean(model_speeds(:,ii,jj,kk));
            mu_specificity(ii,jj,kk) = mean(model_specificities(:,ii,jj,kk));
            CI_speed(ii,jj,kk) = std(model_speeds(:,ii,jj,kk)) ./ sqrt(num_models);
            CI_specificity(ii,jj,kk) = std(model_specificities(:,ii,jj,kk)) ./ sqrt(num_models);
        end
    end
end

%% Generate figures
for ii = 1:num_parameters
    figure();
    hold on;
    
    parameter = parameters(ii);
    
    title_str = "Parameter: " + parameter;
    title(title_str, "FontSize", 14);
    
    colors = ["g", "r"];
    markers = ["o", "square", "^"];
    
    for jj = 1:num_variants

        X = squeeze(mu_speed(jj, ii, :));
        Y = squeeze(mu_specificity(jj, ii, :));
        
        y_ci = squeeze(CI_specificity(jj, ii, :));
        
        xconf = [X; X(end:-1:1)];
        yconf = [Y+y_ci; Y(end:-1:1)-y_ci];
        
        ci_bounds = fill(xconf,yconf,'red');
        if jj == 1
            ci_bounds.FaceColor = [0.8 1 0.8];
            ci_bounds.EdgeColor = 'none';
        else
            ci_bounds.FaceColor = [1 0.8 0.8];
            ci_bounds.EdgeColor = 'none';
        end

        color = colors(jj);
        marker = "o";

        %scatter(X, Y, 10, range, "filled", marker);
        plot(X,Y,"k","LineWidth",1.5);
        
        plot(mean(native_speed(jj,:)), mean(native_specificity(jj,:)), "MarkerFaceColor", color, "Marker", "o", "LineWidth", 2);
    end
    %colormap("copper");
    %cb = colorbar();
    %cb.Label.String = "Free-Energy Perturbation (KT)";
    hold off;
    
    ylabel("Specificity (unitless)", "FontSize", 12);
    xlabel("Speed (sec)", "FontSize", 12);
    
    figure_file = "ensemble-tradeoff-combined-fixed"+parameter;
    
    ylim([0,1.5]); % set lower bound to 0 for ease of discussion; in theory we can reach
                   % regimes in which the off-target is more favorable
                   % (leading to Specificity < 0)
    
    pbaspect([1,1,1]);
    %save_figure(gcf, "../figures/", "tradeoffs-ensemble-ci-fixed-arrow-lower-bound/", figure_file);
    %close;
        
end

%% Clean-up
% switch back to the working directory
cd(working_directory);

%% Local Functions
function err = compute_average_error(models)
num_variants = size(models,1);
num_substrates = size(models,2);
% need to make sure that we do not count the on-target substrate!
raw_err = zeros(num_variants, num_substrates-1);
for ii = 1:num_variants
    for jj = 2:num_substrates
        raw_err(ii,jj-1) = compute_numerical_cleavage_error(models(ii,jj));
    end
end
err = mean(raw_err,2);
end

function S = compute_average_specificity(models)
err = compute_average_error(models);
S = 1 - err;
end

function V = compute_speed(models)
V = zeros(2,1);
V(1) = 1 ./ compute_numerical_mfpt(models(1,1));
V(2) = 1 ./ compute_numerical_mfpt(models(2,1));
end