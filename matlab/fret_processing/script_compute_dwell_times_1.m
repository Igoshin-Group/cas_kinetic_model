% Script to compute the dwell times for each state
% Date: 2021-07-07
% Updated: 2021-07-07

clc;
clear;
close all;

calculate_state_dwell_times("../Output/Processed Data/ebFRET-Fits/wt-0bp_traces.dat", "../Output/Dwell-Times/wt-0bp/");
calculate_state_dwell_times("../Output/Processed Data/ebFRET-Fits/wt-3bp_traces.dat", "../Output/Dwell-Times/wt-3bp/");
calculate_state_dwell_times("../Output/Processed Data/ebFRET-Fits/wt-4bp_traces.dat", "../Output/Dwell-Times/wt-4bp/");

calculate_state_dwell_times("../Output/Processed Data/ebFRET-Fits/hf1-0bp_traces.dat", "../Output/Dwell-Times/hf1-0bp/");
calculate_state_dwell_times("../Output/Processed Data/ebFRET-Fits/hf1-1bp_traces.dat", "../Output/Dwell-Times/hf1-1bp/");
calculate_state_dwell_times("../Output/Processed Data/ebFRET-Fits/hf1-2bp_traces.dat", "../Output/Dwell-Times/hf1-2bp/");
calculate_state_dwell_times("../Output/Processed Data/ebFRET-Fits/hf1-3bp_traces.dat", "../Output/Dwell-Times/hf1-3bp/");
calculate_state_dwell_times("../Output/Processed Data/ebFRET-Fits/hf1-4bp_traces.dat", "../Output/Dwell-Times/hf1-4bp/");

calculate_state_dwell_times("../Output/Processed Data/ebFRET-Fits/eSpy-0bp_traces.dat", "../Output/Dwell-Times/eSpy-0bp/");
calculate_state_dwell_times("../Output/Processed Data/ebFRET-Fits/eSpy-1bp_traces.dat", "../Output/Dwell-Times/eSpy-1bp/");
calculate_state_dwell_times("../Output/Processed Data/ebFRET-Fits/eSpy-2bp_traces.dat", "../Output/Dwell-Times/eSpy-2bp/");
calculate_state_dwell_times("../Output/Processed Data/ebFRET-Fits/eSpy-3bp_traces.dat", "../Output/Dwell-Times/eSpy-3bp/");
calculate_state_dwell_times("../Output/Processed Data/ebFRET-Fits/eSpy-4bp_traces.dat", "../Output/Dwell-Times/eSpy-4bp/");


function calculate_state_dwell_times(path, output_folder)
% compute the dwell times for each state

[dwell, split] = get_dwell_and_split(path);

state1_dwells = dwell{1};
state2_dwells = dwell{2};
state3_dwells = dwell{3};



% create outputs
mkdir(output_folder);
cwd = pwd();

cd(output_folder);

figure();
histogram(state1_dwells);
xlabel("Dwell Time (s)");
ylabel("Frequency (#)");
title("Open State (State 1) Dwell Time Distribution");
saveas(gca,"./State1_Dwell.png");

figure();
histogram(state2_dwells);
xlabel("Dwell Time (s)");
ylabel("Frequency (#)");
title("Open State (State 2) Dwell Time Distribution");
saveas(gca,"./State2_Dwell.png");

figure();
histogram(state3_dwells);
xlabel("Dwell Time (s)");
ylabel("Frequency (#)");
title("Open State (State 3) Dwell Time Distribution");
saveas(gca,"./State3_Dwell.png");

save("./dwell_times.mat","dwell", "split");

cd(cwd); % return to the working directory
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