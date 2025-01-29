function state_path = get_fret_state_path(trace, state_values)
% get_fret_state_path: retrieve the viterbi path for a fret trace
%   state_path = get_fret_state_path(trace)
% inputs:
%   trace = FRET trace
%   state_values = FRET state values
% outputs:
%   state_path = FRET state path (vector)

% set default state_values based on Dagdas et al. 2017, Sci. Adv.
if ~exist("state_values","var")
    state_values = [0.19, 0.34, 0.97];
end

% initialization
trace_length = length(trace.state_mean);
num_states = length(state_values);
thresholds = zeros(num_states-1,1);

% check that state_values is ascending
if ~issorted(state_values)
    error("The values in state_values must be ascending");
end

% compute thresholds between states
for ii=1:num_states-1
    thresholds(ii) = (state_values(ii) + state_values(ii+1)) / 2;
end

% calculate the state path given the thresholds
state_path = zeros(trace_length,1);
for ii=1:trace_length
    current_value = trace.state_mean(ii);
    tf = current_value > thresholds; % tf = [0, 0, 0, ...] if state 1
    state_path(ii) = sum(tf,1) + 1;
end
end