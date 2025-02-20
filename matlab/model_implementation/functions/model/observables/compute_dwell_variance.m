function mu = compute_dwell_variance(pset)
% compute_dwell_variance: compute the 2nd-moment of state dwell time dist.
% inputs:
%   pset = parameter structure (struct)
% outputs:
%   mu = state dwell time variance (3x1 vector)

% extract parameters
[k1a,k11a,k2,k22,k3,k33,k4,k44,...
          k5,k55,k6,k66,k7,k77,k8,k88,kclv,k1b,k11b,k9,k99,...
          f1a,f11a,f2,f22,f3,f33,f4,f44,...
          f5,f55,f6,f66,f7,f77,f8,f88,fclv,f1b,f11b,f9,f99]=extract_parameters(pset);

% pre-allocate output
mu = zeros(3,1);
% open state, low FRET
mu(1) = (f2.*f3.*k2.*k3+f22.*f5.*k22.*k5).^(-2).*(f2.*f3.*k2.*k3+f5.*( ...
  f22.*k22+f3.*k3).*k5).^(-1).*(f2.*f3.*k2.*k3.*((f2.*k2+f22.*k22) ...
  .^2+2.*f22.*f3.*k22.*k3)+f5.*(f22.*k22.*(f2.*k2+f22.*k22).^2+f3.*( ...
  f2.*k2+(-1).*f22.*k22).^2.*k3).*k5+2.*f2.*f22.*f5.^2.*k2.*k22.* ...
  k5.^2);
% checkpoint state, intermediate FRET
mu(2) = (f2.*f3.*f55.*k2.*k3.*(f33.*k33+f4.*k4).*k55+f22.*f33.*f5.*k22.* ...
  k33.*k5.*(f55.*k55+f6.*k6)).^(-2).*(f22.*f33.*f5.*f7.*k22.*k33.* ...
  k5.*(f55.*k55+f6.*k6).*k7+f2.*f3.*f55.*k2.*k3.*(f33.*k33+f4.*k4).* ...
  k55.*(f55.*k55+f6.*k6+f7.*k7)).^(-1).*(f22.^3.*f33.^3.*f5.^3.*f7.* ...
  k22.^3.*k33.^3.*k5.^3.*(f55.*k55+f6.*k6).*k7+f2.^3.*f3.^3.* ...
  f55.^3.*k2.^3.*k3.^3.*(f33.*k33+f4.*k4).*k55.^3.*(f55.*k55+f6.*k6+ ...
  f7.*k7)+f2.*f22.^2.*f3.*f33.^2.*f5.^2.*f55.*k2.*k22.^2.*k3.* ...
  k33.^2.*k5.^2.*k55.*((f33.*k33+f4.*k4).*(f55.*k55+f6.*k6)+f7.*( ...
  f33.*k33+f4.*k4+2.*f55.*k55+2.*f6.*k6).*k7)+f2.^2.*f22.*f3.^2.* ...
  f33.*f5.*f55.^2.*k2.^2.*k22.*k3.^2.*k33.*k5.*k55.^2.*(2.*f33.^2.* ...
  k33.^2+2.*f4.^2.*k4.^2+(-2).*f4.*k4.*(f55.*k55+f6.*k6+(-1).*f7.* ...
  k7)+2.*f33.*k33.*(2.*f4.*k4+(-1).*f55.*k55+(-1).*f6.*k6+f7.*k7)+( ...
  f55.*k55+f6.*k6).*(2.*f55.*k55+2.*f6.*k6+f7.*k7)));
% closed state, high FRET
mu(3) = f44.^(-2).*f66.^(-2).*k44.^(-2).*(f2.*f3.*f4.*f55.*k2.*k3.*k4.* ...
  k55+f22.*f33.*f5.*f6.*k22.*k33.*k5.*k6).^(-2).*k66.^(-2).*(f22.* ...
  f33.*f5.*f6.*f8.*k22.*k33.*k5.*k6.*k8+f2.*f3.*f4.*f55.*k2.*k3.* ...
  k4.*k55.*(f66.*k66+f8.*k8)).^(-1).*(f22.^3.*f33.^3.*f44.^2.* ...
  f5.^3.*f6.^3.*f8.*k22.^3.*k33.^3.*k44.^2.*k5.^3.*k6.^3.*k8+f2.^3.* ...
  f3.^3.*f4.^3.*f55.^3.*f66.^2.*k2.^3.*k3.^3.*k4.^3.*k55.^3.* ...
  k66.^2.*(f66.*k66+f8.*k8)+f2.*f22.^2.*f3.*f33.^2.*f4.*f44.*f5.^2.* ...
  f55.*f6.^2.*k2.*k22.^2.*k3.*k33.^2.*k4.*k44.*k5.^2.*k55.*k6.^2.*( ...
  2.*f66.*f8.*k66.*k8+f44.*k44.*(f66.*k66+f8.*k8))+f2.^2.*f22.* ...
  f3.^2.*f33.*f4.^2.*f5.*f55.^2.*f6.*f66.*k2.^2.*k22.*k3.^2.*k33.* ...
  k4.^2.*k5.*k55.^2.*k6.*k66.*(2.*f44.^2.*k44.^2+2.*f44.*k44.*((-1) ...
  .*f66.*k66+f8.*k8)+f66.*k66.*(2.*f66.*k66+f8.*k8)));
end
