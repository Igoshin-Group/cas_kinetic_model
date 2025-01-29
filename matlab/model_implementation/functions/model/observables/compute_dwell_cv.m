function cv = compute_dwell_cv(pset)

mu = compute_dwell_means(pset);
var = compute_dwell_variance(pset);
std = sqrt(var);

cv = std ./ mu;

end