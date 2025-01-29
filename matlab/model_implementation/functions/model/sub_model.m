function y = sub_model(pset, options)
% sub_model: compute Cas9 model statistics
%   y = sub_model(pset,use_bind,use_mfpt,use_state_prob,use_dwell_mean,use_split,use_error)
% inputs:
%   pset = parameter set (struct)
%   use_bind = compute normalized binding (ndABA) for each substrate
%   use_kd = compute dissociation constant (Kd) for each substrate
%   use_mfpt = compute cleavage MFPT for each substrate
%   use_state_prob = compute steady-state HNH occupancies for each substrate
%   use_dwell_mean = compute HNH state dwell means for each substrate
%   use_dwell_variance = compute HNH state dwell variances for each substrate
%   use dwell_skewness = compute HNH state dwell skewness for each substrate
%   use_split = compute HNH domain transition splitting probability for each substrate
%   use_error = compute the cleavage error for each substrate
% outputs:
%   y = model output (vector)
% note: observables should always be in the same order

% pre-assign observables
normalized_binding = [];
dissociation_constant = [];
cleavage_mfpt = [];
fret_state_occupancies = [];
fret_dwell_means = [];
fret_dwell_variance = [];
fret_dwell_skewness = [];
fret_dwell_cv = [];
fret_transition_split = [];
cleavage_error = [];
non_productive_state = [];
transition_cv = [];

% compute observables given the parameter set
if options.normalized_binding
    normalized_binding = compute_normalized_binding_affinity(pset); % nan
end

if options.dissociation_constant
    dissociation_constant = analytical_dissociation_constant(pset);
end

if options.mean_first_passage_time
    cleavage_mfpt = compute_off_target_mfpt(pset); % should be valid for on-target as well
end

if options.state_occupancy
    fret_state_occupancies = compute_state_occupancies(pset);
end

if options.dwell_mean
    fret_dwell_means = compute_dwell_means(pset);
end

if options.dwell_variance
    fret_dwell_variance = compute_dwell_variance(pset);
end

if options.dwell_skewness
    fret_dwell_skewness = compute_dwell_skewness(pset);
end

if options.dwell_cv
    fret_dwell_cv = compute_dwell_cv(pset);
end

if options.transition_bias
    fret_transition_split = compute_transition_splitting_probability(pset);
end

if options.relative_cleavage_rate
    cleavage_error = compute_cleavage_error(pset);
end

if options.non_productive_state
    non_productive_state = compute_non_productive_fraction(pset);
end

if options.transition_path_cv % transition path time moment
    transition_cv = transition_path_cv(pset);
end

% package output vector; empty observables will not cause errors
y = [normalized_binding; ...
     dissociation_constant; ...
     cleavage_mfpt; ...
     fret_state_occupancies; ...
     fret_dwell_means; ...
     fret_dwell_variance; ...
     fret_dwell_skewness; ...
     fret_dwell_cv; ...
     fret_transition_split; ...
     cleavage_error; ...
     non_productive_state; ...
     transition_cv];
end
