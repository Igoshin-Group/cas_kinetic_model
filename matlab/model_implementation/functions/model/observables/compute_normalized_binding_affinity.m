function ndABA = compute_normalized_binding_affinity(model)
on_model = model;
off_model = create_off_target_parameter_set(model);
Kd_on = analytical_dissociation_constant(on_model);
Kd_off = analytical_dissociation_constant(off_model); % Kd for off-target (i.e.: ki,w = fi * ki,r)
ndABA = Kd_on ./ Kd_off;
end