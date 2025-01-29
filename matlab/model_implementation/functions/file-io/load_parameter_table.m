function [table, scores] = load_parameter_table(folder)
% load_parameter_table: load the parameters for a set of fits
%   [table, scores]  = load_parameter_table(folder)
% inputs:
%   folder = folder path containing fitting run outputs
% outputs:
%   table = parameter table
%   scores = vector of fit scores

files = dir(folder);
num_files = length(files);
table = [];

for ii=1:num_files
    file = files(ii);
    
    % check if file is a folder
    if file.isdir
        continue;
    end
    
    % load parameters
    filepath = fullfile(folder, file.name);
    parameters = readmatrix(filepath);
    
    % store in parameter matrix
    table = [table; parameters];
end
scores = table(:,end);
table(:,end) = [];
end