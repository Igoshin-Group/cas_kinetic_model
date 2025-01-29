function traces = SimulateEnsemble2D(t0, tf, S0, model, N, seed)
% SimulateEnsemble2D: simulate an N-member ensemble of traces
%   traces = SimulateEnsemble2D(t0, tf, S0, model, N)
% inputs:
%   t0, tf = initial and final time-points (double)
%   S0 = initial state (int)
%   model = model parameters (struct)
%   N = number of ensemble members (int)
%   seed = seed for random number generator (int)
% outputs:
%   traces = struct array
%       .time = time-points (double)
%       .state = system states (int)

% Set the seed
rng(seed);

% Simulate N trajectories
traces = struct();
for ii = 1:N
    [t, S] = models.SimulateModel2D(t0,tf,S0,model);
    trace = struct("time",t,"state",S);
    traces = struct.struct_append(traces, trace);
end
end