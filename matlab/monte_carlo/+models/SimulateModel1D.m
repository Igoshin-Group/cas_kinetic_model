function [t,S] = SimulateModel1D(t0,tf,S0,model)
% SimulateModel1D: simulate the 1D model
%   [t,S] = SimulateModel1D(t0,tf,S0,model)
% inputs:
%   t0 = initial time (double)
%   tf = final time (double)
%   S0 = initial state (int)
%   model = model parameters (struct)
% outputs:
%   t = timepoints (vector)
%   S = states (vector)

num_reactions = 4;
num_states = 3;

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
    propensities(1) = k12 .* cur_S(1);
    propensities(2) = k21 .* cur_S(2);
    propensities(3) = k23 .* cur_S(2);
    propensities(4) = k32 .* cur_S(3);
    
    total_propensity = sum(propensities);

    % compute dt
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
        case 1
            cur_S(1) = cur_S(1) - 1;
            cur_S(2) = cur_S(2) + 1;
        case 2
            cur_S(2) = cur_S(2) - 1;
            cur_S(1) = cur_S(1) + 1;            
        case 3
            cur_S(2) = cur_S(2) - 1;
            cur_S(3) = cur_S(3) + 1;            
        case 4
            cur_S(3) = cur_S(3) - 1;
            cur_S(2) = cur_S(2) + 1;            
        otherwise
            error("Invalid reaction: %i\n", reaction);
    end
    
    % update the simulation time
    cur_t = cur_t + dt;
    
    % save the updated system state and time
    S(end+1) = find(cur_S); % add the index of the current state (which corresponds to state number)
    t(end+1) = cur_t;
end
end