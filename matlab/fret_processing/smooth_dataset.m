function ds = smooth_dataset(dataset, order, window_length)
% smooth_dataset: denoise an smFRET trace using a SG filter
% ds = smooth_dataset(dataset, order, window_length)
% inputs:
%   dataset = smFRET dataset (struct)
%       .data
%           .donor
%           .acceptor
%           .time
%   order = filter order (default: 3)
%   window_length = filter window size (default: 11)
% output:
%   smoothed = denoised dataset (struct)
%       .data
%           .donor
%           .acceptor
%           .time

ds = dataset;

% check inputs
if ~exist('order', 'var')
    order = 3;
end

if ~exist('window_length', 'var')
    window_length = 11;
end

% denoise traces
for ii = 1:length(dataset.data)
    ds.data(ii) = smooth_trace(dataset.data(ii), order, window_length);
end
end