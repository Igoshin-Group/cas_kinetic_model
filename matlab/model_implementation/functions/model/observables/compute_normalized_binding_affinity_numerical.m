function ndABA = compute_normalized_binding_affinity_numerical(model)
Kd_on = compute_dissociation_constant_numerical(model); % Kd for on-target (i.e.: fi = 1)
Kd_off = compute_dissociation_constant_numerical(model); % Kd for off-target (i.e.: ki,w = fi * ki,r)
ndABA = Kd_on ./ Kd_off;
end