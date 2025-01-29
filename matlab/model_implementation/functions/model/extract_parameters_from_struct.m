function extract_parameters_from_struct(s)
% extract_parameters_from_struct: extract the parameters from a struct
% inputs:
% s = struct
% note: we assume that the structure contains only doubles

% get struct fields
fields = fieldnames(s);

% assign struct elements to caller workspace
for ii=1:length(fields)
    field = fields{ii};
    value = s.(field);
    assignin('caller',field,value);
end
end