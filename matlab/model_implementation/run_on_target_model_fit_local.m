%
%   Cas9 Project
%   File: run_simultaneous_model_fits_local.m
%   Description: Runs the fitting code for simultaneous
%                on-/off-target model fits on the local PC.
%

working_directory = pwd(); % get the script directory
switch_to_project_directory(); % switch to the code directory

clc;

% load energies from prior best fit model
load("model-02g-wt-48.mat");
load("model-02g-hf1-48.mat");

x0 = [Ewt; Ehf1];

for ii = 1:length(smFRET_substrates)
    substrate = smFRET_substrates(ii);
    for seed = 1:1
        [x, fx] = fit_simultaneous_model(seed, substrate, x0);
    end
end

% switch back to the working directory
cd(working_directory);