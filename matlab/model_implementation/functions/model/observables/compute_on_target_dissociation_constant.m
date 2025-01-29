function Kd = compute_on_target_dissociation_constant(model)
equation = @(x) 0.5 - on_target_fraction_bound(model,10.^x);
options = optimset("FunValCheck","on");
Kd = 10 .^ fzero(equation,1,options);
end