function cv = get_transition_path_cv(traces,num_states,dt)

% set num_states if not provided
if ~exist("num_states","var")
    num_states = 3; % corresponds with Chen et al. 2017 paper
end

% set time_step if not provided
if ~exist("dt","var")
    dt = 0.1; % acquisition time-step (sec; see Chen et al. 2017)
end

num_traces = length(traces);

transition_paths = [];
for ii = 1:num_traces
    
    % get the state path for the current trace
    state_path = get_fret_state_path(traces(ii));
    
    % put the state path and time into a struct for processing
    trajectory.state = state_path;
    trajectory.time = (0:length(state_path)-1)*dt;
    
    % calculate the transition paths
    transitions = determine_transition_path_distribution(trajectory);
    transition_paths = [transition_paths, transitions];
end

mu = mean(transition_paths);
sigma = std(transition_paths);
cv = sigma / mu;
end