function simulation = simulate_ensemble(t0, tf, S0, model, p, N, seed)
% simulate_ensemble: generate an N-member Monte-Carlo simulation
%   traces = SimulateEnsemble2D(t0, tf, S0, model, N)
% inputs:
%   t0, tf = initial and final time-points (double)
%   S0 = initial state (int)
%   model = model function (function handle)
%   p = model parameters (struct)
%   N = number of ensemble members (int)
%   seed = seed for random number generator (int)
% outputs:
%   simulation = struct
%       .seed = seed
%       .trajectories = struct array
%           .time = time-points (double)
%           .state = system states (int)

% output structure
simulation = struct();
simulation.seed = seed;

% set the RNG seed
rng(seed);

% simulate N trajectories
trajectories = struct();
for ii = 1:N
    [t, S] = model(t0,tf,S0,p);
    trajectory = struct("time",t,"state",S);
    trajectories = struct.struct_append(trajectories, trajectory);
end
simulation.trajectories = trajectories;
end