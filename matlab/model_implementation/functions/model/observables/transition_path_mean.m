function mu = transition_path_mean(pset)

k7 = pset.k7;
k33 = pset.k33;
k4 = pset.k4;
k55 = pset.k55;
k6 = pset.k6;
k22 = pset.k22;
k2 = pset.k2;
k3 = pset.k3;
k5 = pset.k5;

k77 = k7 * ((k33 * k22 * k5) / (k3 * k2 * k55));

mu = (k7+k77).^(-1).*((k33+k4).*(k55+k6+k7)+(k55+k6).*k77).^(-1).*( ...
  k77.*(k33.*k6+k4.*(k6+k7)+k6.*k77).^(-1).*(k33.^2.*k6+k4.^2.*(k6+ ...
  k7)+k4.*k7.*(k55+k6+k7)+k33.*k4.*(2.*k6+k7)+2.*k33.*k6.*k77+k4.*( ...
  2.*k6+k7).*k77+k6.*k77.*(k7+k77))+k7.*(k4.*(k55+k6+k7)+k6.*k77).^( ...
  -1).*(k6.*k77.*(k33+k55+k6+k7+k77)+k4.*((k55+k6+k7).^2+(k6+k7).* ...
  k77)));
end