function [dataset,tf] = extract_high_throughput_data_updated(file_path, threshold)
% extract_high_throughput_data: extract a high-throughput dataset
%   [dataset,tf] = extract_high_throughput_data(file_path, threshold)
% inputs:
%   file_path = file path to data (.xlsx)
%   threshold = MFPT cut-off threshold (sec)
% outputs:
%   dataset = high-throughput data (struct array, for one Cas9 variant)
%   tf = indices of datapoints with MFPT above the threshold

% set threshold if not supplied
if ~exist("threshold","var")
    warning("Using default MFPT threshold");
    threshold = 60000;
end

% load data
data_table = readtable(file_path);
num_substrates = height(data_table); % returns number of rows

% create dataset structure
dataset = [];

for ii=1:num_substrates
    data = struct();

    data.cleavage_error = data_table(ii,:).cleavage_error;
    data.cleavage_error_uncertainty = data_table(ii,:).cleavage_error_sd;
    data.normalized_binding = data_table(ii,:).ndABA;
    data.normalized_binding_uncertainty = data_table(ii,:).ndABA_sd;

    data.cleavage_mfpt = data_table(ii,:).cleavage_mfpt;
    data.cleavage_rate = data_table(ii,:).kclv;
    data.cleavage_rate_uncertainty = data_table(ii,:).kclv_sd;

    data.distal_mismatches = data_table(ii,:).distal_mismatches;
    data.total_mismatches = data_table(ii,:).total_mismatches;
    data.proximal_mismatches = data.total_mismatches - data.distal_mismatches;
    data.non_productive_state = [];
    data.dissociation_constant = [];

    data.valid_pam = data_table(ii,:).valid_pam;
    
    data.mismatch_pattern = data_table(ii,:).mismatch_pattern;
    data.substrate_id = ii;

    % store data in dataset
    if isempty(dataset)
        dataset = data; % create the first struct if necessary
    else
        dataset(ii) = data;
    end
end
% get indices of invalid datapoints
tf = [dataset.cleavage_mfpt] > threshold;
%tf = or(tf,[dataset.cleavage_rate] < 0.0001*dataset(1).cleavage_rate); % cleavage_rate < 0.01% of on-target
%tf = or(tf,[dataset.normalized_binding] < 0); % negative normalized binding
end

