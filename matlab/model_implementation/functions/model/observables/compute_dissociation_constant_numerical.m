function Kd = compute_dissociation_constant_numerical(model)
equation = @(x) 0.5 - compute_numerical_fraction_bound(model,10.^x);
options = optimset("FunValCheck","on");

% global observed_x
% global observed_P
%     
% observed_x = [];
% observed_P = [];

% max_Kd = 1e10;
% min_Kd = 1e-10;
% 
% left_end = equation(10.^min_Kd);
% right_end = equation(10.^max_Kd);
% 
% if (sign(left_end) ~= sign(right_end))
%     if (~isnan(left_end)) && (~isnan(right_end))
%         [log10_Kd, fx, exitflag] =  fzero(equation,log10([min_Kd, max_Kd]),options);
%     else
%         log10_Kd = log10(min_Kd);
%     end
% else
%     if sign(right_end) == 1
%         log10_Kd = log10(max_Kd);
%     else
%         log10_Kd = log10(min_Kd);
%     end
% end
% 
% Kd = 10 .^ log10_Kd;

% fprintf("Kd: %4.2f\n", Kd);
[log10_Kd, fx, exitflag] =  fzero(equation,1,options);
fprintf("Kd = %4.4f | F(x) = %4.4f | Exit Flag: %i\n", 10.^log10_Kd, fx, exitflag);
Kd = 10 .^ log10_Kd;
end
