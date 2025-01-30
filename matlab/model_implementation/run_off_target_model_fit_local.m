%
%   Cas9 Project
%   File: run_simultaneous_model_fits_local.m
%   Description: Runs the fitting code for simultaneous
%                on-/off-target model fits on the local PC.
%

init_path();

load("./data/matlab/distal_substrate_indices.mat");

rank = 1; % use the best on-target model

num_substrates = length(substrates);

for ii = 1:num_substrates
    fprintf("Substrate ID: %i\n", substrates(ii));

    folder = "../../output/fitting/off-target-model-fit/test/";

    num_fits = 1;
    for seed = 1:num_fits
        run_off_target_fit_cluster(seed, ii, rank);
    end
end