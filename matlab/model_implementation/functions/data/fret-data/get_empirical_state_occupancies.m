function P = get_empirical_state_occupancies(traces,num_states)
% get_empirical_state_occupancies: compute FRET state probabilities
%   P = get_empirical_state_occupancies(traces,num_states)
% inputs:
%   traces = FRET trace data (struct array)
%   num_states = number of FRET states (default: 3; see Chen et al. 2017)
% outputs:
%   P = steady-state occupancies for each FRET state (num_statesx1 vector)

% set num_states if not passed in
if ~exist("num_states","var")
    num_states = 3; % consistent with Chen et al. 2017
end

% initialization
P = zeros(num_states,1);
num_traces = length(traces);

for ii=1:num_traces
    state_path = get_fret_state_path(traces(ii));
    fret_state = state_path(end); % assume final state == steady-state
    P(fret_state) = P(fret_state) + 1; % increment counter
end

% normalize to probability vector
Z = sum(P);
P = P ./ Z; % final state occupancy probability vector

% replace P=0 with small epsilon and re-normalize
P(P==0) = 1E-4;
Z = sum(P);
P = P ./ Z;
end