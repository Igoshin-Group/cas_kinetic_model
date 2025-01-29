function y_obs = create_data_vector(dataset,options)
% create_data_vector: create the y_obs vector for fitting
%   y_obs = create_data_vector(dataset,use_bind,use_mfpt,use_state_prob,use_dwell_mean,use_split,use_error)
% input:
%   dataset = Cas9 high-throughput and FRET data (struct)
%   use_bind = include normalized binding (ndABA) for each substrate
%   use_mfpt = include cleavage MFPT for each substrate
%   use_state_prob = include steady-state HNH occupancies for each substrate
%   use_dwell_mean = include HNH state dwell means for each substrate
%   use_dwell_variance = include HNH state dwell variance for each substrate
%   use_dwell_skewness = include HNH state dwell skewness for each substrate
%   use_split = include HNH domain transition splitting probability for each substrate
%   use_error = include the cleavage error for each substrate
% output:
%   y_obs = fitting data (vector)

% initialization
y_obs = [];
num_targets = length(dataset);

% create vector
for ii=1:num_targets
    % extract the data for the current substrate
    cur_data = dataset(ii);
    
    % pre-assign data variables
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
    
    % collect data for the current substrate (ii)
    if options.normalized_binding
        normalized_binding = cur_data.normalized_binding;
    end
    
    if options.dissociation_constant
        dissociation_constant = cur_data.dissociation_constant;
    end
    
    if options.mean_first_passage_time
        cleavage_mfpt = cur_data.cleavage_mfpt;
    end
    
    if options.state_occupancy
        fret_state_occupancies = cur_data.fret.state_occupancies;
    end
    
    if options.dwell_mean
        fret_dwell_means = cur_data.fret.dwell_means;
    end

    if options.dwell_variance
        fret_dwell_variance = cur_data.fret.dwell_variance;
    end

    if options.dwell_skewness
        fret_dwell_skewness = cur_data.fret.dwell_skewness;
    end
    
    if options.dwell_cv
        fret_dwell_cv = cur_data.fret.dwell_cv;
    end
    
    if options.transition_bias
        fret_transition_split = cur_data.fret.transition_split;
    end
    
    if options.relative_cleavage_rate
        cleavage_error = cur_data.cleavage_error;
    end

    if options.non_productive_state
        non_productive_state = cur_data.non_productive_state;
    end
    
    if options.transition_path_cv
        transition_cv = cur_data.fret.transition_path_cv;
    end
    
    % package data for current substrate (ii)
    y = [normalized_binding; ...
         dissociation_constant; ...
         cleavage_mfpt; ...
         fret_state_occupancies; ...
         fret_dwell_means; ...
         fret_dwell_variance; ...
         fret_dwell_skewness; ...
         fret_dwell_cv; ...
         fret_transition_split;
         cleavage_error; ...
         non_productive_state; ...
         transition_cv];

    % add packaged data to the overall output vector
    y_obs = [y_obs; y];
end

% flatten vector
y_obs = reshape(y_obs.',1,[]);
end
