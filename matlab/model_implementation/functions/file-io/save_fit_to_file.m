function save_fit_to_file(x, fx, filepath)

% write parameter set to file
file_id = fopen(filepath, 'w');

for ii=1:length(x)
    fprintf(file_id, '%4f,', x(ii));
end

fprintf(file_id, '%4f\n', fx);

% close the file
fclose(file_id);
end