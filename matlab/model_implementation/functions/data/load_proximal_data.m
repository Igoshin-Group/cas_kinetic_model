function dataset = load_proximal_data(threshold)
% load_proximal_data: load the project data into a dataset struct
%   dataset = load_proximal_data()
% inputs:
%   none
% outputs:
%   dataset = project data (struct)

if ~exist('threshold', 'var')
    threshold = 60000; % 60000 (sec) MFPT threshold (approximate time of Jones experiment)
end

% load all the data
dataset = load_all_data(threshold);

% remove substrates possessing PAM-proximal mismatches
idx = find([dataset.wt.distal_mismatches] ~= 0);
dataset.wt(idx) = [];
dataset.hf1(idx) = [];
dataset.enh(idx) = [];

% filter out substrates where the error is > 1
idx = find([dataset.wt.cleavage_error] > 1);
dataset.wt(idx) = [];
dataset.hf1(idx) = [];
dataset.enh(idx) = [];
end
