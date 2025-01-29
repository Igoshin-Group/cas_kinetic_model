function [models, energies, scores] = load_model_ensemble(substrates, ranks)
models = cell(size(ranks));
energies = cell(size(ranks));
scores = cell(size(ranks));
for ii = 1:length(ranks)
    [models{ii}, energies{ii}, scores{ii}] = load_parameter_sets(substrates, ranks(ii));
end
end

function [models, energies, scores] = load_parameter_sets(substrates, rank)

% load on-target model

folder = "../../output/fitting/on_target_model_fit/sub_48/";
[fits, scores] = load_parameter_table(folder);
[on_model_fit, on_score] = get_parameter_set_by_rank(fits, scores, rank);

C = 1; % [Cas9] (nM)
[model, on_target_energies] = generate_model_from_energies(on_model_fit, C);

on_target_model = model(:,1);

num_substrates = 2;
num_variants = 2;

% load off-target model
psets_off = [];
energies_off = [];
folder_base = "../../output/fitting/off_target_model_fit/Ensemble-"+...
    num2str(rank)+"/";
num_substrates = length(substrates); % add one for the on-target
best_scores = zeros(num_substrates+1,1);
best_scores(1) = on_score;
for ii = 1:num_substrates
    folder = folder_base + "sub_" + num2str(substrates(ii)) + "/";
    
    [fits, scores] = load_parameter_table(folder);

    % find the minimum score
    [min_score, idx] = min(scores);

    % get the best parameter set
    best_fit = fits(idx,:);
    best_scores(ii+1) = min_score; % shift by 1 for the on-target
    
    C = 1; % [Cas9] (nM)
    [model, energies] = generate_off_target_model_from_energies(best_fit, on_target_energies, C);
    energies = energies(:);
    
    psets_off_cur = model(:,2);
    
%     for jj = 1:2
%         psets_off_cur(jj) = rate_floor(psets_off_cur(jj), 1e-2); % set a floor
%         psets_off_cur(jj) = rate_ceiling(psets_off_cur(jj), 1e2); % set a ceiling
%     end
    if isempty(psets_off)
        psets_off = psets_off_cur;
    else
        psets_off(:,end+1) = psets_off_cur;
    end
    
    if isempty(energies_off)
        energies_off = energies;
    else
        energies_off(:,end+1) = energies;
    end
end

models = [on_target_model, psets_off];
energies = energies_off;
scores = best_scores;
end

function adjusted = rate_floor(model, threshold)
    fields = fieldnames(model);
    for ii = 1:length(fields)
        field = fields{ii};
        rate = model.(field);
        if rate < threshold; model.(field) = threshold; end
    end
    adjusted = model;
end

function adjusted = rate_ceiling(model, threshold)
    fields = fieldnames(model);
    for ii = 1:length(fields)
        field = fields{ii};
        rate = model.(field);
        if rate > threshold; model.(field) = threshold; end
    end
    adjusted = model;
end