function output_parameter_json(model, file)

% generate JSON from the model
if isMATLABReleaseOlderThan('R2021a')
    model_json = jsonencode(model);
else
    model_json = jsonencode(model, 'PrettyPrint', true);
end

% write to file
file_id = fopen(file, 'w');
fprintf(file_id,model_json);
fclose(file_id);


end