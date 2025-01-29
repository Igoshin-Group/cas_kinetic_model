function [t,S] = SimulateModel2D(t0,tf,S0,model)
% SimulateModel2D: simulate the 2D model
%   [t,S] = SimulateModel2D(t0,tf,S0,model)
% inputs:
%   t0 = initial time (double)
%   tf = final time (double)
%   S0 = initial state (int)
%   model = model parameters (struct)
% outputs:
%   t = timepoints (vector)
%   S = states (vector)

num_reactions = 14;
num_states = 6;

% extract model parameters
extract_parameters_from_struct(model);

% initialize simulation
cur_t = 0;
cur_S = zeros(num_states,1); cur_S(S0) = 1;

t = [t0];
S = [S0];

propensities = zeros(num_reactions,1);
while cur_t <= tf

    % compute reaction propensities
    propensities(1) = k12 .* cur_S(1); % 1 -> 2
    propensities(2) = k21 .* cur_S(2); % 2 -> 1
    propensities(3) = k23 .* cur_S(2); % 2 -> 3
    propensities(4) = k32 .* cur_S(3); % 3 -> 2
    propensities(5) = k45 .* cur_S(4); % 4 -> 5
    propensities(6) = k54 .* cur_S(5); % 5 -> 4
    propensities(7) = k56 .* cur_S(5); % 5 -> 6
    propensities(8) = k65 .* cur_S(6); % 6 -> 5
    propensities(9) = k14 .* cur_S(1); % 1 -> 4
    propensities(10) = k41 .* cur_S(4); % 4 -> 1
    propensities(11) = k25 .* cur_S(2); % 2 -> 5
    propensities(12) = k52 .* cur_S(5); % 5 -> 2
    propensities(13) = k36 .* cur_S(3); % 3 -> 6
    propensities(14) = k63 .* cur_S(6); % 6 -> 3
    
    total_propensity = sum(propensities);

    % compute dt (time until the next reaction)
    r1 = rand();
    dt = (1 ./ total_propensity) .* log(1 ./ r1);

    % determine the next reaction
    r2 = rand();
    threshold = r2 .* total_propensity;
    
    reaction = 0;
    for ii = 1:num_reactions
        temp = 0;
        for jj = 1:ii
            temp = temp + propensities(jj);
        end
        if temp > threshold
            reaction = ii;
            break;
        end
    end
    
    % update the system state
    switch reaction
        case 1 % 1 -> 2
            cur_S(1) = cur_S(1) - 1;
            cur_S(2) = cur_S(2) + 1;
        case 2 % 2 -> 1
            cur_S(2) = cur_S(2) - 1;
            cur_S(1) = cur_S(1) + 1;            
        case 3 % 2 -> 3
            cur_S(2) = cur_S(2) - 1;
            cur_S(3) = cur_S(3) + 1;            
        case 4 % 3 -> 2
            cur_S(3) = cur_S(3) - 1;
            cur_S(2) = cur_S(2) + 1;
        case 5 % 4 -> 5
            cur_S(4) = cur_S(4) - 1;
            cur_S(5) = cur_S(5) + 1;            
        case 6 % 5 -> 4
            cur_S(5) = cur_S(5) - 1;
            cur_S(4) = cur_S(4) + 1;            
        case 7 % 5 -> 6
            cur_S(5) = cur_S(5) - 1;
            cur_S(6) = cur_S(6) + 1;
        case 8 % 6 -> 5
            cur_S(6) = cur_S(6) - 1;
            cur_S(5) = cur_S(5) + 1;            
        case 9 % 1 -> 4
            cur_S(1) = cur_S(1) - 1;
            cur_S(4) = cur_S(4) + 1;            
        case 10 % 4 -> 1
            cur_S(4) = cur_S(4) - 1;
            cur_S(1) = cur_S(1) + 1;            
        case 11 % 2 -> 5
            cur_S(2) = cur_S(2) - 1;
            cur_S(5) = cur_S(5) + 1;            
        case 12 % 5 -> 2
            cur_S(5) = cur_S(5) - 1;
            cur_S(2) = cur_S(2) + 1;            
        case 13 % 3 -> 6
            cur_S(3) = cur_S(3) - 1;
            cur_S(6) = cur_S(6) + 1;            
        case 14 % 6 -> 3
            cur_S(6) = cur_S(6) - 1;
            cur_S(3) = cur_S(3) + 1;            
        otherwise
            error("Invalid reaction: %i\n", reaction);
    end
    
    % update the simulation time
    cur_t = cur_t + dt;
    
    % save the updated system state and time
    state = 0;
    if (cur_S(1) || cur_S(4))
        state = 1;
    elseif (cur_S(2) || cur_S(5))
        state = 2;
    elseif (cur_S(3) || cur_S(6))
        state = 3;
    end
    S(end+1) = state; % add the index of the current state (which corresponds to state number)
    t(end+1) = cur_t;
end
end