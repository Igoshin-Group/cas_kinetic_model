function off_pset = create_off_target_parameter_set(pset)
% create_off_target_parameter_set:
% inputs:
%   pset = parameter set structure (struct)
% outputs:
%   off_pset = parameter set for off-target substrate (struct)

% extract parameters
[k1a,k11a,k2,k22,k3,k33,k4,k44,...
          k5,k55,k6,k66,k7,k77,k8,k88,kclv,k1b,k11b,k9,k99,...
          f1a,f11a,f2,f22,f3,f33,f4,f44,...
          f5,f55,f6,f66,f7,f77,f8,f88,fclv,f1b,f11b,f9,f99]=extract_parameters(pset);


% create off-target parameter set
off_pset = struct();

% compute parameters based on discrimination factors and on-target values
off_pset.k1a = k1a * f1a;
off_pset.k11a = k11a * f11a;
off_pset.k1b = k1b * f1b;
off_pset.k11b = k11b * f11b;
off_pset.k2 = k2 * f2;
off_pset.k22 = k22 * f22;
off_pset.k3 = k3 * f3;
off_pset.k33 = k33 * f33;
off_pset.k4 = k4 * f4;
off_pset.k44 = k44 * f44;
off_pset.k5 = k5 * f5;
off_pset.k55 = k55 * f55;
off_pset.k6 = k6 * f6;
off_pset.k66 = k66 * f66;
off_pset.k7 = k7 * f7;
off_pset.k8 = k8 * f8;
off_pset.k9 = k9 * f9;
off_pset.kclv = kclv * fclv;

% calculate parameters set by detailed balance
k77 = k7 * ((k33 * k22 * k5) / (k3 * k2 * k55));
f77 = f7 * ((f33 * f22 * f5) / (f2 * f3 * f55));

k88 = k8 * ((k44 * k33 * k22 * k5 * k6) / (k4 * k3 * k2 * k55 * k66));
f88 = f8 * ((f44 * f33 * f22 * f5 * f6) / (f4 * f3 * f2 * f55 * f66));

k99 = k9 * ((k11a * k55 * k1b) / (k1a * k5 * k11b));
f99 = f9 * ((f11a * f55 * f1b) / (f1a * f5 * f11b));

off_pset.k77 = k77 * f77;
off_pset.k88 = k88 * f88;
off_pset.k99 = k99 * f99;

% set all discrimination factors to 1 (they are now already accounted for in the
% ki's )
off_pset.f1a = 1;
off_pset.f11a = 1;
off_pset.f1b = 1;
off_pset.f11b = 1;
off_pset.f2 = 1;
off_pset.f22 = 1;
off_pset.f3 = 1;
off_pset.f33 = 1;
off_pset.f4 = 1;
off_pset.f44 = 1;
off_pset.f5 = 1;
off_pset.f55 = 1;
off_pset.f6 = 1;
off_pset.f66 = 1;
off_pset.f7 = 1;
off_pset.f77 = 1;
off_pset.f8 = 1;
off_pset.f88 = 1;
off_pset.f9 = 1;
off_pset.f99 = 1;
off_pset.fclv = 1;
end