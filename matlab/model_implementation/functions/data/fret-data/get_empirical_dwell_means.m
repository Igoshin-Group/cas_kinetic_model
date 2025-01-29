function dwell_means = get_empirical_dwell_means(traces,num_states,dt)
% get_empirical_dwell_means: compute FRET state dwell times from data
%   dwell_times = get_empirical_dwell_times(traces,num_states,dt)
% inputs:
%   traces = FRET trace data (struct array)
%   num_states = number of FRET states possible (default=3)
%   dt = acquisition time-step (sec, default=0.1)
% outputs:
%   dwell_means = state dwell time means (vector)

% set num_states if not provided
if ~exist("num_states","var")
    num_states = 3; % corresponds with Chen et al. 2017 paper
end

% set time_step if not provided
if ~exist("dt","var")
    dt = 0.1; % acquisition time-step (sec; see Chen et al. 2017)
end

num_traces = length(traces);
dwell_means = cell(num_states,1); % container to store dwell times
for ii=1:num_traces
    state_path = get_fret_state_path(traces(ii));
    current_dwells = cell(num_states,1);
    
    % determine dwell times
    counter = 0;
    for jj=1:length(state_path)-2
        counter = counter + 1;
        % check if a transition has occured
        if state_path(jj) ~= state_path(jj+1)
            % transition has occured
            current_dwells{state_path(jj)}(end+1) = counter;
            counter = 0; % reset the counter
        end
    end

    % drop the first and last transition (we do not know start and end)
    for jj=1:num_states
        if length(current_dwells{jj}) > 2
            current_dwells{jj}(1) = []; % drop first
            current_dwells{jj}(end) = []; % drop last
        end
        % save current_dwells to dwell_times structure
        dwell_means{jj} = [dwell_means{jj}, current_dwells{jj}]; % shift by one unit of time, dt, to account for minimum measurement frequency
    end
end

% convert dwell times to units of time (not acquisition time-steps)
for ii=1:num_states
    dwell_means{ii} = dwell_means{ii}*dt;
    
    % compute mean dwell times
    dwell_means{ii} = mean(dwell_means{ii});
    

end
% convert to vector
dwell_means = cell2mat(dwell_means);
end
