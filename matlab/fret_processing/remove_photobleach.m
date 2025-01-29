function s = remove_photobleach(dataset)
% remove_photobleach: cut-off an smFRET trajectory at the point of
%                     photobleaching
% trace = remove_photobleach(trace)
% inputs:
%   dataset = smFRET dataset (struct)
%       .data = smFRET data (1xN struct)
%           .donor
%           .acceptor
%           .time
% outputs:
%   s = smFRET dataset (struct)
%       .data = smFRET data (1xN struct)
%           .donor
%           .acceptor
%           .time

for ii = 1:length(dataset.data)
    trace = dataset.data(ii);
    
    % update to use both donor and acceptor bleaching
    % find the photobleach index
    [idx,d] = external.ebfret.ebfret.analysis.photobleach_index(trace.acceptor(:));


    % create output
    corrected_trace = struct();

    corrected_trace.acceptor = trace.acceptor(1:idx);
    corrected_trace.donor = trace.donor(1:idx);
    corrected_trace.time = trace.time(1:idx);
    
    dataset.data(ii) = corrected_trace;
end
s = dataset;
end