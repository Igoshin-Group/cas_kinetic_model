function Kd = compute_dissociation_constant(model)
equation = @(x) 0.5 - off_target_fraction_bound(model,10.^x);
options = optimset("FunValCheck","on");
Kd = 10 .^ fzero(equation,1,options);
end

function P = off_target_fraction_bound(model, x)

% extract parameters
[k1a,k11a,k2,k22,k3,k33,k4,k44,...
 k5,k55,k6,k66,k7,k8,kclv,k1b,k11b,k9,...
 f1a,f11a,f2,f22,f3,f33,f4,f44,...
 f5,f55,f6,f66,f7,f8,fclv,f1b,f11b,f9]=extract_parameters(model);

% compute bound fraction
P = f1a.*k11b.*k1a.*(f2.*f55.*f66.*k2.*(f33.*f44.*k33.*k44+f3.*k3.*( ...
  f4.*k4+f44.*k44)).*k55.*k66+f22.*f33.*f44.*k22.*k33.*k44.*(f55.* ...
  f66.*k55.*k66+f5.*k5.*(f6.*k6+f66.*k66))).*x.*(f22.*f33.*f44.* ...
  f55.*f66.*k11a.*(f11b.*k11b+f1b.*k1a).*k22.*k33.*k44.*k55.*k66+ ...
  f1a.*k11b.*k1a.*(f2.*f55.*f66.*k2.*(f33.*f44.*k33.*k44+f3.*k3.*( ...
  f4.*k4+f44.*k44)).*k55.*k66+f22.*f33.*f44.*k22.*k33.*k44.*(f55.* ...
  f66.*k55.*k66+f5.*k5.*(f6.*k6+f66.*k66))).*x).^(-1);  
end
