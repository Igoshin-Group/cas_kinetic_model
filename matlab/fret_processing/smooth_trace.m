function smoothed = smooth_trace(trace, order, window_length)
% smooth_trace: denoise an smFRET trace using a SG filter
% smoothed = smooth_trace(trace, order, window_length)
% inputs:
%   trace = smFRET trace (struct)
%       .donor
%       .acceptor
%       .time
%   order = filter order (default: 3)
%   window_length = filter window size (default: 11)
% output:
%   smoothed = denoised trace (struct)
%       .donor
%       .acceptor
%       .time

if ~exist('order', 'var')
    order = 3;
end

if ~exist('window_length', 'var')
    window_length = 11;
end

filter = @(x) sgolayfilt(x, order, window_length);

% denoise the trace
smoothed = trace;
smoothed.donor = filter(trace.donor);
smoothed.acceptor = filter(trace.acceptor);
end