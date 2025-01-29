function save_figure(fig, root, folder, name)
% save_figure: save a figure (both PNG and EPS) to a given folder/path
% inputs:
%   fig = figure handle
%   root = path to the project 'figures' directory
%   folder = folder name
%   name = file name (without extension)
% outputs:
%   none

% assemble the filepaths to the figure folder
folder_png = root + folder + "png/";
folder_eps = root + folder + "eps/";

% create the folders if necessary
[status, msg, msgID] = mkdir(folder_png);
[status, msg, msgID] = mkdir(folder_eps);

% assemble the full filepaths for the figure
path_png = folder_png + name + ".png";
path_eps = folder_eps + name + ".eps";

% save the PNG and EPS files (use default settings)
saveas(fig, path_png, "png");
saveas(fig, path_eps, "epsc"); % EPSC is required for color
end
