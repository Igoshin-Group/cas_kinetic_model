function fig = large_figure(width,height)
% large_figure: create a large figure

% create a figure with default settings
fig = figure('Units','inches');

% update figure size
fig.Position(1) = 0.01;
fig.Position(2) = 0.01;
fig.Position(3) = width;
fig.Position(4) = height;
end