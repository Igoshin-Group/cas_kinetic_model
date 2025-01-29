function err = compute_numerical_cleavage_error(pset)

% get transition matrix
K = compute_backwards_transition_matrix(pset);
num_states = length(K);

% set source terms
source_r = zeros(num_states,1); source_r(7) = pset.kclv;
source_w = zeros(num_states,1); source_w(14) = pset.fclv .* pset.kclv;

zero = zeros(num_states,1);

% solve for splitting probabilities
pi_r = pinv(K)*(zero-source_r); pi_r = pi_r(1);
pi_w = pinv(K)*(zero-source_w); pi_w = pi_w(1);

% pi_r = K \ (zero - source_r); pi_r = pi_r(1);
% pi_w = K \ (zero - source_w); pi_w = pi_w(1);

% fprintf("PiR = %4.4f | PiW = %4.4f\n", pi_r, pi_w);

err = pi_w ./ pi_r;
end

function Mb = compute_backwards_transition_matrix(pset)
% compute_backwards_transition_matrix: compute Mb for a parameter set
%   Mb = compute_backwards_transition_matrix(pset)
% inputs:
%   pset = parameter set (struct)
% outputs:
%   Mb = transition matrix for backwards master equation

% pre-allocate the transition matrix
num_states = 14; % number of discrete states
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

k99 = k9 * ((k11a * k55 * k1b) / (k1a * k5 * k11b));
f99 = f9 * ((f11a * f55 * f1b) / (f1a * f5 * f11b));


% initialize the transition matrix
% right product branch
Mb(2,1) = k11a;
Mb(8,1) = k11b*f11b;
Mb(9,1) = k11a*f11a;

Mb(1,2) = k1a;
Mb(3,2) = k22;
Mb(4,2) = k55;

Mb(2,3) = k2;
Mb(5,3) = k33;

Mb(2,4) = k5;
Mb(5,4) = k77;
Mb(6,4) = k66;
Mb(8,4) = k9;

Mb(3,5) = k3;
Mb(4,5) = k7;
Mb(7,5) = k44;

Mb(4,6) = k6;
Mb(7,6) = k88;

Mb(5,7) = k4;
Mb(6,7) = k8;
Mb(7,7) = kclv;

Mb(1,8) = k1b*f1b;
Mb(4,8) = k99;
Mb(11,8) = k99*f99;

% wrong product branch
Mb(1,9) = k1a*f1a;
Mb(10,9) = k22*f22;
Mb(11,9) = k55*f55;

Mb(9,10) = k2*f2;
Mb(12,10) = k33*f33;

Mb(8,11) = k9*f9;
Mb(9,11) = k5*f5;
Mb(12,11) = k77*f77;
Mb(13,11) = k66*f66;

Mb(10,12) = k3*f3;
Mb(11,12) = k7*f7;
Mb(14,12) = k44*f44;

Mb(11,13) = k6*f6;
Mb(14,13) = k88*f88;

Mb(12,14) = k4*f4;
Mb(13,14) = k8*f8;
Mb(14,14) = kclv*fclv;

% compute diagonal elements
for ii=1:size(Mb,1)
    Mb(ii,ii) = -sum(Mb(ii,:));
end
end
