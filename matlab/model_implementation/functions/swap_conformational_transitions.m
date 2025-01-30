function [a,b] = swap_conformational_transitions(Pa,Pb)
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

a.k3 = b.k3;
a.k33 = b.k33;
a.k4 = b.k4;
a.k44 = b.k44;

% a.f3 = b.f3;
% a.f33 = b.f33;
% a.f4 = b.f4;
% a.f44 = b.f44;

a.k5 = b.k5;
a.k55 = b.k55;
a.k6 = b.k6;
a.k66 = b.k66;

b.k3 = tmp.k3;
b.k33 = tmp.k33;
b.k4 = tmp.k4;
b.k44 = tmp.k44;

% b.f3 = tmp.f3;
% b.f33 = tmp.f33;
% b.f4 = tmp.f4;
% b.f44 = tmp.f44;

b.k5 = tmp.k5;
b.k55 = tmp.k55;
b.k6 = tmp.k6;
b.k66 = tmp.k66;


% % test swapping k7 and k8
% a.k7 = b.k7;
% a.k8 = b.k8;
% a.k2 = b.k2;
% a.k22 = b.k22;
% a.k11a = b.k11a;
end