function [a,b] = swap_cleavage_rate(Pa,Pb)
% swap_cleavage_rate:
% inputs:
%   Pa, Pb = input parameter sets
% outputs:
%   a, b = swapped parameter sets

a = Pa;
b = Pb;

tmp = a;

a.kclv = b.kclv;
% a.fclv = b.fclv;

b.kclv = tmp.kclv;
% b.fclv = tmp.fclv;

end