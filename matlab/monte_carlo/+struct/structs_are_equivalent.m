function tf = structs_are_equivalent(a,b)
% structs_are_equivalent: check if two structs have the same fields
%   tf = structs_are_equivalent(a,b)
% inputs:
%   a, b = structs
% outputs:
%   tf = boolean

fields_a = fields(a);
fields_b = fields(b);

% compare fields
tf = isempty(setdiff(fields_a, fields_b)); 
end