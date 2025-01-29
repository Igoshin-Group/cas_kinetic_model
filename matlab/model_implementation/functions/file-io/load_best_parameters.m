function [x, fx] = load_best_parameters(folder)
x = load_parameter_table(folder);

scores = x(:,end);  % extract scores and delete from table
x(:,end) = [];

[fx, idx] = min(scores); % find best score

x = x(idx,:);

end