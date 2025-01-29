function Mf = create_forward_matrix_no_cleavage(model_parameters, x)
% compute_forward_matrix: compute Mf for a parameter set
%   Mf = create_forward_matrix(model_parameters)
% inputs:
%  model_parameters = parameter set (struct)
% outputs:
%   Mf = transition matrix for forward master equation

% pre-allocate the transition matrix
num_states = 8; % number of discrete states
Mf = zeros(num_states);

% extract parameters from parameter set
[k1a,k11a,k2,k22,k3,k33,k4,k44,...
          k5,k55,k6,k66,k7,k77,k8,k88,kclv,k1b,k11b,k9,k99,...
          f1a,f11a,f2,f22,f3,f33,f4,f44,...
          f5,f55,f6,f66,f7,f77,f8,f88,fclv,f1b,f11b,f9,f99]=extract_parameters(model_parameters);


% set parameters by detailed balance
k77 = k7 * ((k33 * k22 * k5) / (k3 * k2 * k55));
f77 = f7 * ((f33 * f22 * f5) / (f2 * f3 * f55));

k88 = k8 * ((k44 * k33 * k22 * k5 * k6) / (k4 * k3 * k2 * k55 * k66));
f88 = f8 * ((f44 * f33 * f22 * f5 * f6) / (f4 * f3 * f2 * f55 * f66));

k99 = k9 * ((k1b * k55 * k11a) / (k1a * k5 * k11b)); % * (x02/x01);
f99 = f9 * ((f1b * f55 * f11a) / (f1a * f5 * f11b));


% initialize the transition matrix
Mf(2,1) = x * f1a * k1a;
Mf(8,1) = f1b * k1b;

Mf(1,2) = f11a * k11a;
Mf(3,2) = f2 * k2;
Mf(4,2) = f5 * k5;

Mf(2,3) = f22 * k22;
Mf(5,3) = f3 * k3;

Mf(2,4) = f55 * k55;
Mf(5,4) = f7 * k7;
Mf(6,4) = f6 * k6;
Mf(8,4) = f99 * k99;

Mf(3,5) = f33 * k33;
Mf(4,5) = f77 * k77;
Mf(7,5) = f4 * k4;

Mf(4,6) = f66 * k66;
Mf(7,6) = f8 * k8;

Mf(5,7) = f44 * k44;
Mf(6,7) = f88 * k88;

Mf(1,8) = f11b * k11b;
Mf(4,8) = x * f9 * k9;

% compute diagonal elements
for ii=1:size(Mf,1)
    Mf(ii,ii) = -sum(Mf(:,ii));
end

% disp(Mf);
end