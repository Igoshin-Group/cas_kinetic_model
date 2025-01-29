function T = compute_numerical_mfpt(pset)
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
b = -ones(size(Mb,1),1);
tol = 1e-15; % tolerance
T = lsqminnorm(A,b,tol,"warn");
T = T(1);
end
