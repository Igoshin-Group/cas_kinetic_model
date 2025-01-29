function fret_data = process_fret_data(file_paths,num_states,dt)
% process_fret_data: process empirical FRET data and create data structure
%   fret_data = process_fret_data(file_paths,num_states,dt)
% inputs:
%   file_paths = file paths to FRET data (previously processed by ebFRET)
%   num_states = number of FRET states (default: 3, see Chen et al. 2017)
%   dt = acquisition time-step (default: 0.1 sec, see Chen et al. 2017)
% outputs:
%   fret_data = FRET data structure (struct array, idx = num. mismatches)
% note: file paths should be in order of ascending PAM-distal mismatches
% note: for WT-Cas9, we assume 1,2-bp mismatches are equivalent to 0bp

% initialization
num_conditions = length(file_paths);
fret_data = []; % empty placeholder for the output struct array

% process each data file
for ii=1:num_conditions
    % extract the FRET trace data
    file_path = file_paths{ii};
    traces = extract_traces(file_path);

    % process traces
    dwell_means = get_empirical_dwell_means(traces,num_states,dt);
    state_occupancies = get_empirical_state_occupancies(traces,num_states);
    transition_split = get_empirical_transition_split(traces);

    % create data structure
    temp_struct = struct();
    temp_struct.dwell_means = dwell_means;
    temp_struct.state_occupancies = state_occupancies;
    temp_struct.transition_split = transition_split;
    temp_struct.num_mismatches = ii-1; % num. PAM-distal mismatches
    temp_struct.traces = traces;
    
    % store data in fret_data structure
    if isempty(fret_data)
        fret_data = temp_struct;
    else
        fret_data(ii) = temp_struct;
    end
end
end