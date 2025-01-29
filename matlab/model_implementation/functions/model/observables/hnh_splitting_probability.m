function P21 = hnh_splitting_probability(model)
% extract parameters
[k1a,k11a,k2,k22,k3,k33,k4,k44,...
k5,k55,k6,k66,k7,k77,k8,k88,kclv,k1b,k11b,k9,k99,...
f1a,f11a,f2,f22,f3,f33,f4,f44,...
f5,f55,f6,f66,f7,f77,f8,f88,fclv,f1b,f11b,f9,f99]=extract_parameters(model);

P21 = f55.*k55.*(f2.*f3.*f55.*k2.*k3.*(f33.*k33+f4.*k4).*k55+f33.*f7.* ...
  k33.*(f2.*f3.*k2.*k3+f22.*f5.*k22.*k5).*k7).*(f2.*f3.*f55.*k2.* ...
  k3.*(f33.*k33+f4.*k4).*k55.*(f55.*k55+f6.*k6)+f7.*(f55.*(f2.*f3.* ...
  k2.*k3.*(f33.*k33+f4.*k4)+f22.*f33.*f5.*k22.*k33.*k5).*k55+f22.* ...
  f33.*f5.*f6.*k22.*k33.*k5.*k6).*k7).^(-1);
end