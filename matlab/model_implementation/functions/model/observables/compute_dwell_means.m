function mu = compute_dwell_means(pset)
% compute_dwell_means: compute the 1st-moment of state dwell time dist.
% inputs:
%   pset = parameter structure (struct)
% outputs:
%   mu = state dwell time means (3x1 vector)

% extract parameters
[k1a,k11a,k2,k22,k3,k33,k4,k44,...
          k5,k55,k6,k66,k7,k77,k8,k88,kclv,k1b,k11b,k9,k99,...
          f1a,f11a,f2,f22,f3,f33,f4,f44,...
          f5,f55,f6,f66,f7,f77,f8,f88,fclv,f1b,f11b,f9,f99]=extract_parameters(pset);

% pre-allocate output
mu = zeros(3,1);
% open state, low FRET
mu(1) = (f2.*k2+f22.*k22).*(f2.*f3.*k2.*k3+f22.*f5.*k22.*k5).^(-1);
% checkpoint state, intermediate FRET
mu(2) = (f22.*f33.*f5.*k22.*k33.*k5+f2.*f3.*f55.*k2.*k3.*k55).*(f2.*f3.* ...
  f55.*k2.*k3.*(f33.*k33+f4.*k4).*k55+f22.*f33.*f5.*k22.*k33.*k5.*( ...
  f55.*k55+f6.*k6)).^(-1);
% closed state, high FRET
mu(3) = (f22.*f33.*f44.*f5.*f6.*k22.*k33.*k44.*k5.*k6+f2.*f3.*f4.*f55.* ...
  f66.*k2.*k3.*k4.*k55.*k66).*(f2.*f3.*f4.*f44.*f55.*f66.*k2.*k3.* ...
  k4.*k44.*k55.*k66+f22.*f33.*f44.*f5.*f6.*f66.*k22.*k33.*k44.*k5.* ...
  k6.*k66).^(-1);
end
