function data = combine_fret_and_high_throughput_data(data, fret_data)
% combine_fret_and_high_throughput_data: combine datasets
%   data = combine_fret_and_high_throughput_data(ht_data, fret_data)
% inputs:
%   data = high-throughput dataset (struct array)
%   fret_data = FRET dataset (struct array)
% outputs:
%   data = combined dataset

num_targets = length(data);
for ii=1:num_targets
    num_mismatches = data(ii).distal_mismatches;
    data(ii).fret = fret_data(num_mismatches+1);
end
end