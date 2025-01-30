%
%   Cas9 Project
%   File: run_simultaneous_model_fits_local.m
%   Description: Runs the fitting code for simultaneous
%                on-/off-target model fits on the local PC.
%

init_path();

substrate = 48;

fprintf("Substrate ID: %i\n", substrate);

% load energies from prior best fit model
load("./data/matlab/model-02g-wt-48.mat");
load("./data/matlab/model-02g-hf1-48.mat");

x0 = [Ewt; Ehf1];    

% generate optimization bounds
lb_val = -20;
ub_val = 20;
lb = lb_val * ones(68,1);
ub = ub_val * ones(68,1);

folder = "../../output/fitting/on-target-model-fit/test/";

num_fits = 1;
for seed = 1:num_fits
    % run on-target model fit
    [x, fx] = fit_simultaneous_model_relaxed_hnh_assumptions_fixed(seed, substrate, lb, ub, x0, folder);
    fprintf("Score: %f\n", fx);
    for ii = 1:length(x)
        fprintf("Parameter %i: %f\n", ii, x(ii));
    end
end