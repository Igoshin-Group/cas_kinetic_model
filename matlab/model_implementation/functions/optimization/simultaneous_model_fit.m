function score = simultaneous_model_fit(p, data, weights, model, alpha)
% model_fit: evaluates the fitness of an input parameter set
%   score = model_fit(p, n, y_obs)
% inputs:
%   p = model parameter vector
%   n = number of parameter sets / conditions
%   y_obs = vector of experimental data
%   y_norm = normalization vector (y_obs^2 with elements == 0 set to 1)
%   weights = weighting factors for each data point
%   model = function handle for model
%   alpha =
% outputs:
%   score = objective function score

% compute model predictions
[y, mfpt, psets_on, psets_off] = model(p);

% convert both to column vectors
y = y(:);
y_obs = data.y_obs(:);
y_norm = data.y_norm(:);
mfpt_obs = data.mfpt_obs;

% calculate error (normalized sum of squares)
objective_function_values = ((y - y_obs) .^ 2) ./ y_norm;
objective_function_values = weights .* objective_function_values;
score = sum(objective_function_values, 'omitnan');

% handle on-target parameter sets
for ii = 1:length(psets_on)
    score = score + qualitative_penalty(psets_on(ii), alpha); % add penalty for violating qualitative constraint
end

% % handle off-target parameter sets
% for ii = 1:length(psets_off)
%     score = score + qualitative_penalty(psets_off(ii), alpha); % add penalty for violating qualitative constraint
% end

% add penalty for cleavage MFPT
detection_limit = 6e4; % limit of detection based on Jones et al. experiment
for ii = 1:length(mfpt)
    if mfpt_obs(ii) > detection_limit
        mfpt_penalty = (heaviside(detection_limit - mfpt(ii)) * quadratic_penalty(mfpt_obs(ii), mfpt(ii)));
        score = score + mfpt_penalty;
    else
        mfpt_penalty = quadratic_penalty(mfpt_obs(ii), mfpt(ii));
        score = score + mfpt_penalty;
    end
end
end

%% Local Functions
function q = quadratic_penalty(a,b)
    q = ((b-a)./(a)).^2;
end

function H = heaviside(x)
    if x < 0
        H = 0.0;
    elseif x == 0
        H = 0.5;
    else
        H = 1.0;
    end
end

function penalty = qualitative_penalty(p, alpha)
% add qualitative penalty for initial fluxes J1 > alpha*J2
Z = @(p) p.k1b + p.k11b;    % partition function
P1 = @(p) p.k11b ./ Z(p);   % state S01 probability at t=0
P2 = @(p) 1 - P1(p);        % state S02 probability at t=0
J1 = @(p) p.k1a .* P1(p);   % state S01 to S11 flux
J2 = @(p) p.k9 .* P2(p);    % state S02 to S21 flux

g1 = 0.1; % quadratic penalty weight
g2 = 0.1; % count penalty weight

penalty = g1 .* heaviside(alpha*J2(p) - J1(p)) .* quadratic_penalty(J1(p), alpha*J2(p)) + ...
          g2 .* heaviside(alpha*J2(p) - J1(p));
end
