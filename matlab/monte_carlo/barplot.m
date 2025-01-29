function barplot(x, y, err)
% Creates a bar plot from input data.
% function barplot(x, y, err)
% Inputs:
%   x = category labels or x-values (vector or array)
%   y = y-values (vector or array)
%   err = standard errors (vector or array)
% reference:
% https://www.mathworks.com/matlabcentral/answers/102220-how-do-i-place-errorbars-on-my-grouped-bar-graph-using-function-errorbar-in-matlab

hold on;

b = bar(y);

xticklabels(x);

ngroups = size(y,1);
nbars = size(y,2);

if ngroups > 1
    groupwidth = min(0.8, nbars/(nbars+1.5));
    for ii=1:nbars
        xt = (1:ngroups) - groupwidth/2 + (2*ii-1) * groupwidth / (2*nbars);
        er = errorbar(xt, y(:,ii), err(:,ii), 'k-');
        er.LineStyle = 'none';
    end 
else
    xt = 1:nbars;
    xticks(xt);
    er = errorbar(xt, y, err, 'k-');
    er.LineStyle = 'none';
end


hold off;
end
