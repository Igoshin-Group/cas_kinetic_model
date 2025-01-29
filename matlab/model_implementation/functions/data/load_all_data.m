function dataset = load_all_data(threshold, state_occupancies_from_fret_histograms)
% load_all_data: load the project data into a dataset struct
%   dataset = load_all_data()
% inputs:
%   use_pam_proximal = if false, filter out PAM-proximal mismatched substrates
% outputs:
%   dataset = project data (struct)

if ~exist('threshold', 'var')
    threshold = 60000; % 60000 (sec) MFPT threshold (approximate time of Jones experiment)
end

if ~exist('state_occupancies_from_fret_histograms', 'var')
    state_occupancies_from_fret_histograms = false;
end

% high-throughput file paths
wt_file_path = "./data/processed/high-throughput/wt_dataset.xlsx";
hf1_file_path = "./data/processed/high-throughput/hf1_dataset.xlsx";
enh_file_path = "./data/processed/high-throughput/enh_dataset.xlsx";

% load corrected state occupancies (calculated from FRET traces without photobleach detection)
if state_occupancies_from_fret_histograms
state_occupancy_data = load("./data/processed/fret-data/state_occupancy_data.mat");
end

% FRET data file paths
% NOTE: we are missing WT-Cas9 1bp and 2bp data, so we assume similar
% dynamics to 0bp PAM-distal mismatches; the FRET distributions for these
% mismatch variants appear to be similar (see Chen et al. 2017)\
wt_fret_file_paths = ["./data/processed/fret-data/wt-cas9/wt-0bp_traces.dat",...
                      "./data/processed/fret-data/wt-cas9/wt-0bp_traces.dat",...
                      "./data/processed/fret-data/wt-cas9/wt-0bp_traces.dat",...
                      "./data/processed/fret-data/wt-cas9/wt-3bp_traces.dat",...
                      "./data/processed/fret-data/wt-cas9/wt-4bp_traces.dat"];
hf1_fret_file_paths = ["./data/processed/fret-data/hf1-cas9/hf1-0bp_traces.dat",...
                      "./data/processed/fret-data/hf1-cas9/hf1-1bp_traces.dat",...
                      "./data/processed/fret-data/hf1-cas9/hf1-2bp_traces.dat",...
                      "./data/processed/fret-data/hf1-cas9/hf1-3bp_traces.dat",...
                      "./data/processed/fret-data/hf1-cas9/hf1-4bp_traces.dat"];
enh_fret_file_paths = ["./data/processed/fret-data/enh-cas9/enh-0bp_traces.dat",...
                      "./data/processed/fret-data/enh-cas9/enh-1bp_traces.dat",...
                      "./data/processed/fret-data/enh-cas9/enh-2bp_traces.dat",...
                      "./data/processed/fret-data/enh-cas9/enh-3bp_traces.dat",...
                      "./data/processed/fret-data/enh-cas9/enh-4bp_traces.dat"];

% load high-throughput datasets
mfpt_threshold = threshold; % set MFPT threshold based on input

% WT-Cas9
[wt_dataset,tf] = extract_high_throughput_data(wt_file_path,mfpt_threshold);
invalid_idx = tf;

% HF1-Cas9
[hf1_dataset,tf] = extract_high_throughput_data(hf1_file_path,mfpt_threshold);
invalid_idx = invalid_idx | tf;

% eSpy-Cas9 (Enh-Cas9)
[enh_dataset,tf] = extract_high_throughput_data(enh_file_path,mfpt_threshold);
invalid_idx = invalid_idx | tf;

%remove invalid datapoints (MFPT > threshold; see loading function)
wt_dataset(invalid_idx) = [];
hf1_dataset(invalid_idx) = [];
enh_dataset(invalid_idx) = [];

% load FRET datasets
num_states = 3; % number of discrete FRET states (see Chen et al. 2017)
dt = 0.1; % acquisition time-step (see Chen et al. 2017)

wt_fret_data = process_fret_data(wt_fret_file_paths,num_states,dt);
hf1_fret_data = process_fret_data(hf1_fret_file_paths,num_states,dt);
enh_fret_data = process_fret_data(enh_fret_file_paths,num_states,dt);

% combine datasets
wt_dataset = combine_fret_and_high_throughput_data(wt_dataset,wt_fret_data);
hf1_dataset = combine_fret_and_high_throughput_data(hf1_dataset,hf1_fret_data);
enh_dataset = combine_fret_and_high_throughput_data(enh_dataset,enh_fret_data);

% update state occupancies
if state_occupancies_from_fret_histograms
state_occupancies = state_occupancy_data.data;

for ii = 1:length(wt_dataset)
    mismatch_idx = wt_dataset(ii).distal_mismatches + 1;
    wt_dataset(ii).fret.state_occupancies = state_occupancies.wt(mismatch_idx).state_occupancies;
    hf1_dataset(ii).fret.state_occupancies = state_occupancies.hf1(mismatch_idx).state_occupancies;
    enh_dataset(ii).fret.state_occupancies = state_occupancies.enh(mismatch_idx).state_occupancies;
end
end

% package the data
dataset = struct();
dataset.wt = wt_dataset;
dataset.hf1 = hf1_dataset;
dataset.enh = enh_dataset;

% filter out substrates where the error is > 1
idx = find([dataset.wt.cleavage_error] > 1);
dataset.wt(idx) = [];
dataset.hf1(idx) = [];
dataset.enh(idx) = [];
end
