function traces = extract_traces(path)
% extract_traces: extract HMM fitted FRET traces from a data file
%   traces = extract_traces(path)
% inputs:
%   path = filepath to data
% outputs:
%   traces = traces structure

% load the data file
data = readmatrix(path);

% create the traces struct
traces = struct();

% read and store each trace
num_traces = max(data(:,1));
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