function s = objective_function_off_target_relaxed_hnh_assumptions_fixed(x, on_target_fit, on_data, off_data)
% objective_function: computes the objective function for a given parameter set
% inputs:
%   x = parameter vector (Nx1)
%   on_data = on-target substrate data (Mx1 struct array, M = number of Cas variants)
%   off_data = off-target substrate data (Mx1 struct array)

    C = 1; % [Cas9] (nM), necessary to compute binding affinities
    [~, energies] = generate_model_from_energies(on_target_fit, C);
    [model, energies] = generate_off_target_model_from_energies(x, energies, C);
    
    on_model_wt = model(1,1);
    on_model_hf1 = model(2,1);
    
    off_model_wt = model(1,2);
    off_model_hf1 = model(2,2);
    
    energies_wt = energies(1);
    energies_hf1 = energies(2);
    
    % =================================================================================================
    %   compute the base objective function (comparing fit to observed data)
    % =================================================================================================
    
    % compute the base objective function score for each variant on-target substrate
    s = 0;
    for ii = 1:2
        y_obs = [on_data(ii).dissociation_constant; ...
                 on_data(ii).fret.state_occupancies; ...
                 on_data(ii).fret.dwell_means; ...
                 on_data(ii).fret.dwell_cv];

        y_pred = [analytical_dissociation_constant(model(ii,1)); ...
                  compute_state_occupancies(model(ii,1)); ...
                  compute_dwell_means(model(ii,1)); ...
                  compute_dwell_cv(model(ii,1))];
              
        % compute the weights; weight = 1 if y_obs and y_pred differ by
        % less than one order-of-magnitude, and (log10(y_obs) - log10(y_pred))^2
        % otherwise (for each component of y_obs)
        weights = get_data_weights(y_pred, y_obs);              
              
        % compute the MSE for the current variant
        for jj = 1:length(y_obs)
            if isnan(y_obs(jj)); continue; end % skip NaN values to avoid blowing up the objective function            
            s = s + weights(jj) * ((y_pred(jj) - y_obs(jj))^2 / y_obs(jj)^2);
        end
    end
    
    % compute the base objective function score for each variant off-target substrate
    for ii = 1:2
        y_obs = [off_data(ii).normalized_binding; ...
                 off_data(ii).fret.state_occupancies; ...
                 off_data(ii).fret.dwell_means; ...
                 off_data(ii).fret.dwell_cv];

        y_pred = [compute_normalized_binding_affinity(model(ii,2)); ...
                  compute_state_occupancies(model(ii,2)); ...
                  compute_dwell_means(model(ii,2)); ...
                  compute_dwell_cv(model(ii,2))];
              
        % compute the MSE for the current variant
        for jj = 1:length(y_obs)
            if isnan(y_obs(jj)); continue; end % skip NaN values to avoid blowing up the objective function
            s = s + weights(jj) * ((y_pred(jj) - y_obs(jj))^2 / y_obs(jj)^2);
        end
    end
    
    % =================================================================================================
    %   enforce model constraints
    % =================================================================================================
    % enforce "shared" constraints (should be equivalent to prior model)
%     s = s + shared_constraint(energies_wt, energies_hf1, "E2", "E9", "E5R_dag", "E5W_dag");     % k5
%     s = s + shared_constraint(energies_wt, energies_hf1, "E5", "E11", "E5R_dag", "E5W_dag");    % k55
%     s = s + shared_constraint(energies_wt, energies_hf1, "E5", "E11", "E6R_dag", "E6W_dag");    % k6
%     s = s + shared_constraint(energies_wt, energies_hf1, "E6", "E13", "E6R_dag", "E6W_dag");    % k66
    
    % enforce "shared_discrimination" constraints (should be equivalent to prior model)
    s = s + shared_discrimination_constraint(off_model_wt, off_model_hf1, "f11a");
    
%     s = s + shared_discrimination_constraint(off_model_wt, off_model_hf1, "f3");
%     s = s + shared_discrimination_constraint(off_model_wt, off_model_hf1, "f33");   
%     s = s + shared_discrimination_constraint(off_model_wt, off_model_hf1, "f4");
%     s = s + shared_discrimination_constraint(off_model_wt, off_model_hf1, "f44");
    
    s = s + shared_discrimination_constraint(off_model_wt, off_model_hf1, "f2");
    s = s + shared_discrimination_constraint(off_model_wt, off_model_hf1, "f22");
    s = s + shared_discrimination_constraint(off_model_wt, off_model_hf1, "f7");
    s = s + shared_discrimination_constraint(off_model_wt, off_model_hf1, "f8");
    
    % enforce "fixed" constraints (should be equivalent to prior model)
    % why does the score decrease here?
    s = s + fixed_constraint(energies_wt, energies_hf1, "E1", "E1", "E1aR_dag", "E1aW_dag");    % k1a
    s = s + fixed_constraint(energies_wt, energies_hf1, "E8", "E8", "E9R_dag", "E9W_dag");      % k9
    s = s + fixed_constraint(energies_wt, energies_hf1, "E7", "E14", "EclvR_dag", "EclvW_dag"); % kclv
    
    % enforce "constant" constraints (should be equivalent to prior model)
    s = s + constant_constraint(on_model_wt, on_model_hf1, "k1b", 0.29);
    s = s + constant_constraint(on_model_wt, on_model_hf1, "k11b", 4.2);
    
    % =================================================================================================
    %   enforce other constraints
    % =================================================================================================
    % add a penalty for MFPTs greater than the limit of detection
    detection_threshold = 6E4; % limit of detection in Jones et al. NucleaSeq experiment
    for ii = 1:2 % on-target substrate
        mfpt_obs = on_data(ii).cleavage_mfpt;
        mfpt_pred = compute_off_target_mfpt(model(ii,1));
        if mfpt_obs > detection_threshold
            s = s + heaviside(detection_threshold - mfpt_pred) * quadratic_penalty(mfpt_obs, mfpt_pred);
        else
            s = s + quadratic_penalty(mfpt_obs, mfpt_pred);
        end
    end
    for ii = 1:2 % off-target substrate
        mfpt_obs = off_data(ii).cleavage_mfpt;
        mfpt_pred = compute_off_target_mfpt(model(ii,2));
        if mfpt_obs > detection_threshold
            s = s + heaviside(detection_threshold - mfpt_pred) * quadratic_penalty(mfpt_obs, mfpt_pred);
        else
            s = s + quadratic_penalty(mfpt_obs, mfpt_pred);
        end
    end    
    % add a penalty for relative binding fluxes (in low-FRET vs mid-FRET states)
    alpha = 10; % low-FRET to mid-FRET binding ratio
    for ii = 1:2
        s = s + qualitative_penalty(model(ii), alpha);
    end
    
    % =================================================================================================
    %   enforce soft upper bounds on parameters
    % =================================================================================================
    for ii = 1:2
    s = s + parameter_constraint(model(ii,2), "k1a", "f1a", 100);
    s = s + parameter_constraint(model(ii,2), "k11a", "f11a", 100);
    s = s + parameter_constraint(model(ii,2), "k2", "f2", 100);
    s = s + parameter_constraint(model(ii,2), "k22", "f22", 100);
    s = s + parameter_constraint(model(ii,2), "k3", "f3", 100);
    s = s + parameter_constraint(model(ii,2), "k33", "f33", 100);
    s = s + parameter_constraint(model(ii,2), "k4", "f4", 100);
    s = s + parameter_constraint(model(ii,2), "k44", "f44", 100);
    s = s + parameter_constraint(model(ii,2), "k5", "f5", 100);
    s = s + parameter_constraint(model(ii,2), "k55", "f55", 100);
    s = s + parameter_constraint(model(ii,2), "k6", "f6", 100);
    s = s + parameter_constraint(model(ii,2), "k66", "f66", 100);
    s = s + parameter_constraint(model(ii,2), "k7", "f7", 100);
    s = s + parameter_constraint(model(ii,2), "k77", "f77", 100);
    s = s + parameter_constraint(model(ii,2), "k8", "f8", 100);
    s = s + parameter_constraint(model(ii,2), "k88", "f88", 100);
    s = s + parameter_constraint(model(ii,2), "k9", "f9", 100);
    s = s + parameter_constraint(model(ii,2), "k99", "f99", 100);
    s = s + parameter_constraint(model(ii,2), "kclv", "fclv", 100);
    end
    
    % check for negative objective function value
    if s < 0; s = nan; end
    
    % set to high value if Inf or -Inf
    if isinf(s); s = 1e30; end
    
    % set to high value if NaN
    if isnan(s); s = 1e30; end
    
    % set to high value if complex
    if imag(s) ~= 0; s = 1e30; end
end

%% Local Functions
function s = substrate_constraint(energies, state_r, state_w, barrier_r, barrier_w)
    s = 0;
    s = s + ((energies.(state_r) - energies.(barrier_r)) - (energies.(state_w) - energies.(barrier_w)))^2;
end

function s = fixed_constraint(e1, e2, state_r, state_w, barrier_r, barrier_w)
% WT and HF1 should have the same rate
% On and Off-Target should have the same rate
    s = 0;
    s = s + ((e1.(state_r) - e1.(barrier_r)) - (e2.(state_r) - e2.(barrier_r)))^2;
    s = s + substrate_constraint(e1, state_r, state_w, barrier_r, barrier_w);
    s = s + substrate_constraint(e2, state_r, state_w, barrier_r, barrier_w);
end

function s = shared_constraint(e1, e2, state_r, state_w, barrier_r, barrier_w)
    s = 0;
    s = s + substrate_constraint(e1, state_r, state_w, barrier_r, barrier_w);
    s = s + substrate_constraint(e2, state_r, state_w, barrier_r, barrier_w);
end

function s = shared_discrimination_constraint(m1, m2, fi)
    s = 0;
    s = s + (m1.(fi) - m2.(fi))^2;
end

function s = constant_constraint(m1, m2, ki, value)
    s = 0;
    s = s + (m1.(ki) - value)^2;
    s = s + (m2.(ki) - value)^2;
end

function s = parameter_constraint(model, parameter, factor, threshold)
    rate = model.(parameter);
    off_rate = rate * model.(factor);
    alpha = 1000;
    
    s = alpha .* heaviside(rate - threshold) .* (rate - threshold).^2;
    s = s + (alpha .* heaviside(off_rate - threshold) .* (off_rate - threshold).^2);
end

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

function weights = get_data_weights(y, y_obs)
    weights = ones(size(y));
    for ii = 1:length(y)
        weights(ii) = data_weight(y(ii), y_obs(ii));
    end
end

function weight = data_weight(y, y_obs)
    alpha = 1;
    threshold = 1;
    difference = (log10(y) - log10(y_obs))^2; % get the order-of-magnitude difference
    if difference > threshold
        weight = alpha * difference;
    else
        weight = 1;
    end
end