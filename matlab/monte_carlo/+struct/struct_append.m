function s = struct_append(s,a)
% struct_append: append a struct to a struct array
%   s = struct_append(s,a)
% inputs:
%   s = struct array
%   a = struct to be appended
% outputs:
%   s = struct array
if struct_is_empty(s)
    s = a;
else
    if structs_are_equivalent(s, a)
        s(end+1) = a;
    else
        error("Struct A is not compatible with struct array S\n");
    end
end
end