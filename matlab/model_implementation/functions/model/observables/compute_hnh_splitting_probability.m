function S = compute_hnh_splitting_probability(pset)
% compute_numerical_mfpt: compute the MFPT for a model numerically
%   T = compute_numerical_mfpt(pset)
% inputs:
%   pset = parameter set (struct)
% outputs:
%   T = mean first-passage time (MFPT)

% get the transition matrix
Mb = compute_backwards_transition_matrix(pset);

% compute the MFPT to the absorbing state
A = Mb;
b = zeros(size(Mb,1),1);
tol = 1e-15; % tolerance
T = lsqminnorm(A,b,tol,"warn");
disp(T);
T = T(1);
S=T;
end

function Mb = compute_backwards_transition_matrix(pset)
% compute_backwards_transition_matrix: compute Mb for a parameter set
%   Mb = compute_backwards_transition_matrix(pset)
% inputs:
%   pset = parameter set (struct)
% outputs:
%   Mb = transition matrix for backwards master equation

x = 100; % nM

% pre-allocate the transition matrix
num_states = 7; % number of discrete states
Mb = zeros(num_states);

% extract parameters
[k1a,k11a,k2,k22,k3,k33,k4,k44,...
          k5,k55,k6,k66,k7,k77,k8,k88,kclv,k1b,k11b,k9,k99,...
          f1a,f11a,f2,f22,f3,f33,f4,f44,...
          f5,f55,f6,f66,f7,f77,f8,f88,fclv,f1b,f11b,f9,f99]=extract_parameters(pset);

% set parameters by detailed balance
k77 = k7 * ((k33 * k22 * k5) / (k3 * k2 * k55));
f77 = f7 * ((f33 * f22 * f5) / (f2 * f3 * f55));

k88 = k8 * ((k44 * k33 * k22 * k5 * k6) / (k4 * k3 * k2 * k55 * k66));
f88 = f8 * ((f44 * f33 * f22 * f5 * f6) / (f4 * f3 * f2 * f55 * f66));

k9 = k99 * ((k1a * k5 * k1b) / (k11a * k55 * k11b));
f9 = f99 * ((f1a * f5 * f1b) / (f11a * f55 * f11b));

% initialize the transition matrix
Mb(2,1) = k11a*f11a;
Mb(8,1) = k11b*f11b;


Mb(1,2) = x*k1a*f1a;
Mb(3,2) = k22*f22;
Mb(4,2) = k55*f55;
Mb(4,2) = k55*f55;

Mb(2,3) = k2*f2;
Mb(5,3) = k33*f33;


Mb(2,4) = k5*f5;
Mb(5,4) = k77*f77;
Mb(6,4) = k66*f66;
Mb(8,4) = x*k9*f9;
Mb(4,4) = k55*f55 + k6*f6;

Mb(3,5) = k3*f3;
Mb(4,5) = k7*f7;
Mb(7,5) = k44*f44;
Mb(5,5) = k33*f33 + k4*f4;

Mb(4,6) = k6*f6;
Mb(7,6) = k88*f88;


Mb(5,7) = k4*f4;
Mb(6,7) = k8*f8;
% Mb(7,7) = kclv*fclv;

Mb(1,8) = k1b*f1b;
Mb(4,8) = k99*f99;

% compute diagonal elements
for ii=1:size(Mb,1)
    Mb(ii,ii) = -sum(Mb(ii,:));
end
end
