function run_off_target_fit_cluster(seed, substrate_idx, on_target_model_idx)
    % load the indices for the distal substrates
    load("./data/matlab/distal_substrate_indices.mat");
    substrate = off_idx(substrate_idx);
    
    % load the on-target model
    folder = "../../output/fitting/on_target_model_fit/sub_48/";
    [fits, scores] = load_parameter_table(folder);
    model = get_parameter_set_by_rank(fits, scores, on_target_model_idx);
    
    C = 1; % [Cas9] (nM)
    [~, energies] = generate_model_from_energies(model, C);
    x0 = create_off_target_initial_condition(energies);

    % generate optimization bounds
    lb_val = -20;
    ub_val = 20;
    lb = lb_val * ones(32,1);
    ub = ub_val * ones(32,1);
    
    folder = "../outputs/fitting/off_target_model_fit/"+ ...
             "Ensemble-"+num2str(on_target_model_idx)+"/"; % output folder
    
    % run off-target model fit
    [x, fx] = fit_simultaneous_model_off_target(seed, substrate, lb, ub, x0, folder, model);
    fprintf("Score: %f\n", fx);
    for ii = 1:length(x)
        fprintf("Parameter %i: %f\n", ii, x(ii));
    end
end