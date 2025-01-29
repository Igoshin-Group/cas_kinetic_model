function dataset = load_binary_data(path)
% load_binary_data: load smFRET traces from binary file
% dataset = load_dataset_binary(path)
% inputs:
%   path = filepath (.traces file; binary)
% outputs:
%   dataset = loaded dataset (struct)

file_id = fopen(path, 'r');

if file_id == -1
    error('Invalid file.');
end

% load experiment metadata
trace_length = fread(file_id, 1, 'int32');
num_traces = fread(file_id, 1, 'int16');

% load trace raw data
data = fread(file_id, num_traces*trace_length, 'int16');

% process raw data
data = reshape(data, num_traces, trace_length);

fclose(file_id);

% extract donor and acceptor traces, create output structure
dataset = struct();
dataset.data = struct();
for ii = 1:(num_traces/2)
    dataset.data(ii).time = 0:0.1:(trace_length-1)*0.1;
    dataset.data(ii).donor = data(ii*2-1,:);
    dataset.data(ii).acceptor = data(ii*2,:);
end
end