function export_ascii_data(dataset, path)
% export_ascii_data: save smFRET data to ascii file
% export_ascii_data(dataset, path)
% inputs:
%   dataset = smFRET data (struct)
%   path = output file path (string)
%   format = ascii file format ('stacked', 'unstacked')
% outputs:
%   none

num_traces = length(dataset.data);

file_id = fopen(path, 'w');

% write data to file
for ii = 1:num_traces
    trace = dataset.data(ii);
    for jj = 1:length(trace.donor)
        fprintf(file_id, '%i\t%i\t%i\n', ii, trace.donor(jj), trace.acceptor(jj));
    end
end

fclose(file_id);
end