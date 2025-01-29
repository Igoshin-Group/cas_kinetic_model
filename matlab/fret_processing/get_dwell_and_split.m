function [dwell, split] = get_dwell_and_split(path)
% calculate the dwell times and splitting probabilities given a dataset
% read the data
data = readmatrix(path);
traces = extract_fitted_traces(data);

state_paths = struct();
for ii = 1:length(traces)
    state_paths(ii).path = determine_state_path(traces(ii));
end

% variables to calculate splitting probabilities
N21 = 0;
N23 = 0;

% compute dwell time histogram
state1_dwells = [];
state2_dwells = [];
state3_dwells = [];
for ii = 1:length(state_paths)
    path = state_paths(ii).path;
    
    state1_dwell = [];
    state2_dwell = [];
    state3_dwell = [];
    
    counter = 0;
    for jj = 1:length(path)-2
        counter = counter + 1;
        if path(jj) ~= path(jj+1)
            % transition has occured
            switch path(jj)
                case 1
                    state1_dwell(end+1) = counter;
                case 2
                    state2_dwell(end+1) = counter;
                    if path(jj+1) == 1
                        N21 = N21 + 1;
                    elseif path(jj+1) == 3
                        N23 = N23 + 1;
                    end
                case 3
                    state3_dwell(end+1) = counter;
            end
            counter = 0;
        end
    end
    
    % drop first and last
    if length(state1_dwell) > 2
        state1_dwell(1) = [];
        state1_dwell(end) = [];
        state1_dwells = [state1_dwells state1_dwell];        
    end
    
    if length(state2_dwell) > 2
        state2_dwell(1) = [];
        state2_dwell(end) = [];
        state2_dwells = [state2_dwells state2_dwell];        
        end
    
    if length(state3_dwell) > 2
        state3_dwell(1) = [];
        state3_dwell(end) = [];
        state3_dwells = [state3_dwells state3_dwell];        
    end
end

% convert to real time (not time-steps)
dt = 0.1; % timestep (s)
state1_dwells = state1_dwells*dt;
state2_dwells = state2_dwells*dt;
state3_dwells = state3_dwells*dt;

dwell = {state1_dwells, state2_dwells, state3_dwells};

% calculate splitting probability
NTot = N21 + N23;
pi21 = N21 / NTot;
pi23 = N23 / NTot;
split= [pi21 pi23];
end

function traces = extract_fitted_traces(data)
num_traces = max(data(:,1));
traces = struct();
for ii = 1:num_traces
    idx = find(data(:,1) == ii);
    cur_trace = data(idx,:);
    traces(ii).donor = cur_trace(:,2);
    traces(ii).acceptor = cur_trace(:,3);
    traces(ii).fret = cur_trace(:,4);
    traces(ii).state = cur_trace(:,5);
    traces(ii).state_mean = cur_trace(:,6);
end
end

function path = determine_state_path(trace)
% state FRET values (from Dagdas et al. 2017)
state1 = 0.19;
state2 = 0.34;
state3 = 0.97;

% compute threshold between states
threshold1 = (state1 + state2) / 2;
threshold2 = (state2 + state3) / 2;

% calculate the state path
path = [];
for ii = 1:length(trace.state_mean)
    fret_val = trace.state_mean(ii);
    if fret_val <= threshold1
        path(ii) = 1;
    elseif fret_val > threshold1 && fret_val <= threshold2
        path(ii) = 2;
    else
        path(ii) = 3;
    end
end
end