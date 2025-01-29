function tf = struct_is_empty(s)
% struct_is_empty: check if a struct has no fields
%   tf = struct_is_empty(s)
% inputs:
%   s = struct
% outputs:
%   tf = boolean
tf = isempty(fields(s));
end