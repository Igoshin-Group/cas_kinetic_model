function split = get_empirical_transition_split(traces)
% get_empirical_transition_split: get state 2 transition splitting prob.
%   split = get_empirical_transition_split(traces)
% inputs:
%   traces = FRET trace data (struct array)
% outputs:
%   split = transition splitting probability (state 2 to 1)
% note: here, we assume three (3) FRET states exist

% initialize
num_traces = length(traces);
N21 = 0;
N23 = 0;

% determine splitting probability from data
for ii=1:num_traces
    state_path = get_fret_state_path(traces(ii));
    for jj=1:length(state_path)-2
        % check if a transition has occured
        if state_path(jj) ~= state_path(jj+1)
            % transition has occured, check if it is from state 2
            if state_path(jj) == 2
                next_state = state_path(jj+1);
                if next_state == 1
                    N21 = N21 + 1;
                elseif next_state == 3                    N23 = N23 + 1;
                end
            end
        end
    end
end
% compute splitting probability from transition frequencies
Ntot = N21 + N23;
split = N21 / Ntot; % splitting probability (state 2 to 1)
end