function [x, fx] = generate_fits(fit_function,options)
% generate_fits: generates model fits using particle swarm optimization
%   [fit, score] = generate_fits(fit_function, data_path, output_path, lb, ub, num_runs)
% inputs:
%   fit_function = function handle
%   data = fitting data (struct)
%   lb = lower bounds (vector)
%   ub = upper bounds (vector)
%   seed = RNG seed
% outputs:
%   x = best fit (vector)
%   fx = best fit score

% initialize default parallel pool for the cluster
parallel.defaultClusterProfile;
parallel_profile = parallel.defaultClusterProfile("local");

% get number of fitting parameters
n_vars = length(options.lb);

max_time = options.max_time;

switch options.algorithm
    case 'particleswarm'
        % hybrid optimizer settings - for refinement
        hybrid_opts = optimoptions('patternsearch');
        hybrid_opts.PlotFcn = {'psplotbestf'};

        % PSO optimizer settings
        opts = optimoptions('particleswarm', 'HybridFcn', {'patternsearch', hybrid_opts}); 

        opts.SwarmSize = 50 .* n_vars; % default is 10 particles per parameter
        opts.UseParallel = true;
        opts.FunctionTolerance = 1e-6;
        opts.MaxStallIterations = 20; % maximum iterations with relative change less than FunctionTolerance
        opts.MaxTime = max_time; % maximum running time (in seconds)
        opts.PlotFcn = "pswplotbestf"; % no plotting required when running from command line
        opts.MinNeighborsFraction = 0.1;
        opts.InertiaRange = [0.3, 1];
        opts.FunValCheck = 'on'; % check for NaN or Inf objective function values
        
        % initialize the initial swarm matrix
        initial_swarm_matrix = zeros(opts.SwarmSize, n_vars); % set all particles to the same point
        for row = 1:opts.SwarmSize
            initial_swarm_matrix(row,:) = options.initial_swarm_matrix;
        end
        opts.InitialSwarmMatrix = initial_swarm_matrix;

        % optimization problem structure
        problem.solver = 'particleswarm';
        problem.objective = @(p)fit_function(p);
        problem.nvars = n_vars;
        problem.lb = options.lb;
        problem.ub = options.ub;
        problem.options = opts;
        problem.rngstate = rng(options.seed);

        % run PSO optimizer
        [x,fx] = particleswarm(problem);
    case 'simulannealbnd'
        % hybrid optimizer settings - for refinement
        hybrid_opts = optimoptions('patternsearch');
        hybrid_opts.PlotFcn = {'psplotbestf'};

        % Simulated Annealing optimizer settings
        opts = optimoptions('simulannealbnd', 'HybridFcn', {'patternsearch', hybrid_opts}); 

        opts.InitialTemperature = 1000;
        opts.ReannealInterval = 100;
        opts.FunctionTolerance = 1e-6;
        opts.MaxStallIterations = 50*68; % maximum iterations with relative change less than FunctionTolerance
        opts.MaxTime = 4 .* 60. * 59; % maximum running time (in seconds); 4 minute margin
        opts.PlotFcn = {'saplotbestf','saplottemperature'}; % no plotting required when running from command line

        % optimization problem structure
        problem.solver = 'simulannealbnd';
        problem.objective = @(p)fit_function(p);
        problem.nvars = n_vars;
        problem.lb = options.lb;
        problem.ub = options.ub;
        problem.options = opts;
        problem.rngstate = rng(options.seed);
        problem.x0 = options.x0;
        % run PSO optimizer
        [x,fx] = simulannealbnd(problem);    
end
end

