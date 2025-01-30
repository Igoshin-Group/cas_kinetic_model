function [a,b] = swap_r_loop_transitions(Pa,Pb)
% swap_conformational_transitions:
% inputs:
%   Pa, Pb = input parameter sets
% outputs:
%   a, b = swapped parameter sets

% create output parameter sets
a = Pa;
b = Pb;

% swap conformational transitions
tmp = a;

a.k2 = b.k2;
a.k22 = b.k22;
a.k7 = b.k7;
a.k8 = b.k8;

% a.f2 = b.f2;
% a.f22 = b.f22;
% a.f7 = b.f7;
% a.f8 = b.f8;

b.k2 = tmp.k2;
b.k22 = tmp.k22;
b.k7 = tmp.k7;
b.k8 = tmp.k8;
% 
% b.f2 = tmp.f2;
% b.f22 = tmp.f22;
% b.f7 = tmp.f7;
% b.f8 = tmp.f8;
end