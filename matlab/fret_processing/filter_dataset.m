function ds = filter_dataset(dataset, correlation_cutoff, length_cutoff)
% filter_dataset: filter an smFRET data
%   ds = filter_dataset(dataset, correlation_cutoff, length_cutoff)
% inputs:
%   dataset = smFRET dataset (struct)
%       .data (1xN traces)
%   correlation_cutoff = maximum Pearson correlation (donor/acceptor)
%   length_cutoff = minimum trace length
% outputs:
%   ds = filtered dataset = smFRET dataset (struct)
%       .data (1xN traces)

% set the defaults
if ~exist('correlation_cutoff', 'var')
    correlation_cutoff = -0.5;
end

if ~exist('length_cutoff', 'var')
    length_cutoff = 50; % originally 50
end

for ii = length(dataset.data):-1:1
    % calculate correlation coefficient
    corr = correlation_coefficient(dataset.data(ii).acceptor, dataset.data(ii).donor);
    trace_length = length(dataset.data(ii).donor);
    
    % check if the trace is anti-correlated
    if corr > correlation_cutoff
        dataset.data(ii) = [];
        continue;
    end
    
    % check the trace length
    if trace_length < length_cutoff
        dataset.data(ii) = [];
        continue;
    end
end
ds = dataset;
end