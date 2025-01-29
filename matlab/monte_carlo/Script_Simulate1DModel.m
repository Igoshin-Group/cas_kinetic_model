clc;
clear;
close all;

% Load trace lengths
trace_lengths = readmatrix("./data/wt_0mm_trace_lengths.txt");
num_traces = length(trace_lengths);

% Simulation settings
t0 = 0;
tf = 100;
N = num_traces;
S0 = 1;
num_seeds = 25;
burnin = 800;
interpolate = false;
start_seed = 1231;
frequency = 10; % interpolation frequency, Hz

fprintf("Initial RNG Seed: %i\n", start_seed);
fprintf("Number of seeds: %i\n", num_seeds);
fprintf("Ensemble size: %i\n", N);
fprintf("Trace length range: [%i, %i]\n", min(trace_lengths), max(trace_lengths));

if interpolate
    fprintf("Interpolating traces!\n"); 
    fprintf("Interpolation frequency: %f\n", frequency);
    fprintf("Interpolation time-step: %f\n", 1./frequency);
end

% Simulate model
ensemble_paths = cell(num_seeds,2);
rev_ensemble_paths = cell(num_seeds, 1);

residence_times = cell(num_seeds, 3);
Z_statistic = cell(3,1);

model = @models.SimulateModel1D;
for seed = 0:num_seeds-1
    
    fprintf("Simulating seed (%i): ", seed+1);
    
    % Model parameters
    p = struct();

    % WT Parameters from Liu et al.
    p(1).k12 = 1;
    p(1).k21 = 3;
    p(1).k23 = 2.5;
    p(1).k32 = 1.2;
    
    % HF1 Parameters from Liu et al.
    p(2).k12 = 1;
    p(2).k21 = 0.63;
    p(2).k23 = 2.2;
    p(2).k32 = 0.2;

    for variant = 1:2
    simulation = models.simulate_ensemble(t0, tf+burnin, S0, model, p(variant), N, start_seed + seed);
    traces = simulation.trajectories;

    % Find transition paths
    paths = [];
    rev_paths = [];
    for ii = 1:N
        trace = traces(ii);
        trace_length = trace_lengths(ii);
        trace = trim_trace_to_length(trace, trace_length, burnin, interpolate, frequency);
        paths = [paths, find_transition_paths(trace)];
        rev_paths = [rev_paths, find_reverse_transition_paths(trace)];
        
        
        T = cell(3,1);
        Tn = fret.residence_times(trace, 3);
        for jj = 1:3
            T{jj} = [T{jj}, Tn{jj}];
        end        
    end
    
    
    for jj = 1:3
        residence_times{seed+1, jj} = T{jj};        
        Z = (std(T{jj}) / mean(T{jj}));
        Z_statistic{jj}(end+1) = Z;
    end
    
    ensemble_paths{seed+1,variant} = paths;
    rev_ensemble_paths{seed+1} = rev_paths;  
    fprintf("DONE\n");
    end
end

% % Draw transition path distribution
% figure("Renderer","painters");
% histogram(paths);
% xlabel("Transition Path Length (s)");
% ylabel("Count (#)");
% pbaspect([1,1,1]);
% saveas(gcf,"./1d_distribution.png","png");

% Compute transition path statistics
means = zeros(num_seeds,1);
stds = zeros(num_seeds,1);
cvs = zeros(num_seeds,2);

rev_means = zeros(num_seeds,1);
rev_stds = zeros(num_seeds,1);
rev_cvs = zeros(num_seeds,1);

transitions = zeros(num_seeds,1);
rev_transitions = zeros(num_seeds,1);

for variant = 1:2
    for seed = 1:num_seeds
    paths = ensemble_paths{seed,variant};
    rev_paths = rev_ensemble_paths{seed};
    
    means(seed) = mean(paths);
    stds(seed) = std(paths);
    
    rev_means(seed) = mean(rev_paths);
    rev_stds(seed) = std(rev_paths);
    
    cvs(seed,variant) = stds(seed) ./ means(seed);
    rev_cvs(seed) = rev_stds(seed) ./ rev_means(seed);
    
    transitions(seed) = length(paths);
    rev_transitions(seed) = length(rev_paths);
    
    end
end


mu_forw_cv = mean(cvs);
mu_back_cv = mean(rev_cvs);

fprintf("Forward: %f\nBackward: %f\n", mu_forw_cv, mu_back_cv);


%% Compute transition path statistics
figure("Renderer","painters");
hold on;
scatter(means, cvs);
yline(1,"k--","LineWidth",2);
xlabel("Mean");
ylabel("CV");
pbaspect([1,1,1]);
%saveas(gcf,"./cv_ensemble_without_interpolation.png","png");


% compute probability of observing CV <= 1 (rough)
cvs_round = round(cvs,2); % round to 2 digits to remove statistical noise
tf = cvs_round <= 1;
nans = isnan(cvs_round);
num_nans = sum(nans);
fprintf("Simulations with CV = NaN: %i\n", num_nans);

P = sum(tf) ./ (num_seeds - num_nans);
fprintf("Rough Probability: %f\n", P);

fprintf("Probability of observing 0 or 1 transitions: %f\n", num_nans ./ num_seeds);

figure();
hold on;
histogram(cvs);
xline(1,"r--","LineWidth",2);
xlabel("CV");
ylabel("Frequency (#)");
title("1D Model CV Distribution");
saveas(gcf,"./cv_distribution_without_interpolation.png","png");


figure();
hold on;
histogram(means);
xline(1,"r--","LineWidth",2);
xlabel("Mean (sec)");
ylabel("Frequency (#)");
title("1D Model Mean Distribution");
%saveas(gcf,"./mean_distribution_without_interpolation.png","png");

figure();
hold on;
histogram(stds);
xline(1,"r--","LineWidth",2);
xlabel("Standard Deviation (sec)");
ylabel("Frequency (#)");
title("1D Model STD Distribution");
%saveas(gcf,"./std_distribution_without_interpolation.png","png");


figure();
hold on;
histogram(transitions);
xline(1,"r--","LineWidth",2);
xlabel("Number of Transitions (#)");
ylabel("Frequency (#)");
title("1D Model Transition Count Distribution");
%saveas(gcf,"./transition_count_distribution_without_interpolation.png","png");


%% Draw CV bar plot and save simulated CV data
% Draw bar plot
figure("Renderer","painters");
x_labels = ["WT", "HF1"];
x = categorical(x_labels);
x = reordercats(x, x_labels);

data = zeros(2,1);
sem = zeros(2,1);
sim_ci = zeros(2,2);
for variant = 1:2
mu = mean(cvs(:,variant), "omitnan");
sigma = std(cvs(:,variant), "omitnan");
sem(variant) = sigma ./ sqrt(length(cvs(:,variant)));
data(variant) = mu;

ci = bootci(length(cvs(:,variant)), @(x) mean(x), cvs(:,variant));
sim_ci(:,variant) = ci;
end
hold on;
barplot(x, data', sem');
yline(1,"r--","LineWidth",2);
xlabel("Variant");
ylabel("CV");
pbaspect([1,1,1]);
% ylim([0,1.8]);
% xlim([0,2]);
% saveas(gcf,"./1d_statistics_without_interpolation.png","png");
% saveas(gcf,"./1d_statistics_without_interpolation.eps","epsc");

sim_1d_mu = data;
sim_1d_sem = sem;
sim_1d_ci = sim_ci;
save("./sim_1d_stats.mat", "sim_1d_mu", "sim_1d_sem", "sim_1d_ci");


%% Check residence time distribution

trajectories = simulation.trajectories;
T = cell(3,1);
for ii = 1:N
    trajectory = trajectories(ii);
    trace_length = trace_lengths(ii);
    trajectory = trim_trace_to_length(trajectory, trace_length, burnin, interpolate, frequency);
    Tn = fret.residence_times(trajectory, 3);
    for jj = 1:3
        T{jj} = [T{jj}, Tn{jj}];
    end
end

pvalues = zeros(3,1);
mu_data = zeros(3,1);
figure();
for ii = 1:3
    %figure("Name","Residence Times for State "+num2str(ii));
    data = T{ii};
    mu = mean(data);
    test_cdf = makedist("Exponential","mu",mu);
    rand_sample = exprnd(mu, length(data), 1);
    mu_data(ii) = mu;
    
    subplot(2,3,ii);
    probplot(test_cdf, data);
    xlabel("Residence Time (sec)");
    ylabel("Probability");
    title("Exponential Distribution");
    %legend("Reference", "Data");
    pbaspect([1,1,1]);
    

    
    subplot(2,3,ii+3);
    hold on;
    
    histogram(data, "DisplayStyle", "stairs");
    histogram(rand_sample, "DisplayStyle", "stairs");
    if ii == 3
    legend("Sim.", "Rand. Smp. from Exp(\lambda)");
    end
    xlabel("Residence Time (sec)");
    ylabel("Frequency (#)");
    pbaspect([1,1,1]);
    
    title("Histogram");

    [h,pval] = kstest(data, "CDF", test_cdf, "Alpha", 0.01);
    pvalues(ii) = pval;
    fprintf("P-Value: %f\n", pval);
    if ~h
        fprintf("We cannot reject the null hypothesis. State %i distribution is similar to the exponential distribution.\n", ii);
    else
        fprintf("Null hypothesis rejected. State %i distribution is different from the exponential distribution.\n", ii);
    end
    

end

% saveas(gca,"./1d_model_residence_time_exponential_without_interpolation.png","png");
mu_analytic = analytical_residence_means(p);

bar_data = [mu_data, mu_analytic];
figure();
bar(bar_data);
xlabel("State");
ylabel("Mean Residence Time (sec)");
legend("MC", "Analytical");
% saveas(gca,"./1d_residence_time_vs_analytical_without_interpolation.png","png");

%% Plot CDF
titles = ["Low-FRET","Mid-FRET","High-FRET"];
figure();
for ii = 1:3
    data = T{ii};
    
    mu = mean(data);
    test_cdf = makedist("Exponential","mu",mu);
    
    subplot(1,3,ii);
    hold on;
    histogram(data, "Normalization", "cdf");
    
    x = linspace(min(data), max(data));
    y = expcdf(x, mu);
    plot(x,y,"r-");
    ylabel("Cumulative Distribution");
    xlabel("Residence Time (sec)");
    title(titles(ii));
    pbaspect([1,1,1]);
    
    [h,pval] = kstest(data, "CDF", test_cdf, "Alpha", 0.01);
    pvalues(jj) = pval;
    fprintf("P-Value: %f\n", pval);
    if ~h
        fprintf("We cannot reject the null hypothesis. State %i distribution is similar to the exponential distribution.\n", jj);
    else
        fprintf("Null hypothesis rejected. State %i distribution is different from the exponential distribution.\n", jj);
    end    
end
% saveas(gca,"./1d_model_cum_dist.eps","epsc");

%% Check residence time distributions for all seeds
mu = zeros(3,1);
sem = zeros(3,1);
for ii = 1:3
    mu(ii) = mean(Z_statistic{ii}, "omitnan");
    sem(ii) = std(Z_statistic{ii}, "omitnan") / sqrt(num_seeds);
end

figure();
hold on;
titles = ["Low-FRET", "Mid-FRET", "High-FRET"];
barplot(titles,mu',sem');
yline(1,"r--");
xlabel("HNH Domain State");
ylabel("CV");
ylim([0,2]);
pbaspect([1,1,1]);
% saveas(gca,"./1d_model_residence_time_ensemble_stats_without_interpolation.eps","epsc");

%% Check residence time distribution
trajectories = simulation.trajectories;
T = cell(3,1);
titles = ["Low-FRET", "Mid-FRET", "High-FRET"];
for ii = 1:N
    trajectory = trajectories(ii);
    trace_length = trace_lengths(ii);
    trajectory = trim_trace_to_length(trajectory, trace_length, burnin, interpolate, frequency);
    Tn = fret.residence_times(trajectory, 3);
    for jj = 1:3
        T{jj} = [T{jj}, Tn{jj}];
    end
end

pvalues = zeros(3,1);
mu_data = zeros(3,1);
figure();
for ii = 1:3
    data = T{ii};
    mu = mean(data);
    test_cdf = makedist("Exponential","mu",mu);
    rand_sample = exprnd(mu, length(data), 1);
    mu_data(ii) = mu;
    
    subplot(1,3,ii);
    probplot(test_cdf, data);
    xlabel("Residence Time (sec)");
    
    if ii == 1
        ylabel("Probability");
    else
        ylabel("");
    end
    
    title(titles(ii));
    %legend("Reference", "Data");
    pbaspect([1,1,1]);
    


    [h,pval] = kstest(data, "CDF", test_cdf, "Alpha", 0.01);
    pvalues(ii) = pval;
    fprintf("P-Value: %f\n", pval);
    if ~h
        fprintf("We cannot reject the null hypothesis. State %i distribution is similar to the exponential distribution.\n", ii);
    else
        fprintf("Null hypothesis rejected. State %i distribution is different from the exponential distribution.\n", ii);
    end
    

end
% saveas(gca,"./1d_model_residence_time_exponential_probplot_without_interpolation.png","png");
% saveas(gca,"./1d_model_residence_time_exponential_probplot_without_interpolation.eps","epsc");




%% Local Functions
function paths = find_transition_paths(trace)
% 1 -> 2 -> 3 transitions
paths = [];

t = trace.time;
S = trace.state;
N = length(t);

t_start = 0;
in_transition = false;
for ii = 2:N
    cur_S = S(ii);
    prev_S = S(ii-1);
    
    if ((cur_S == 2) && (prev_S == 1)) % transition in progress
        t_start = t(ii);
        in_transition = true;
        continue;
    end
    
    if ((cur_S == 3) && (prev_S == 2) && in_transition) % transition completed
        paths(end+1) = t(ii) - t_start;
        in_transition = false;
        continue;
    end
end
end

function paths = find_reverse_transition_paths(trace)
% 3 -> 2 -> 1 transitions
paths = [];

t = trace.time;
S = trace.state;
N = length(t);

t_start = 0;
in_transition = false;
for ii = 2:N
    cur_S = S(ii);
    prev_S = S(ii-1);
    
    if ((cur_S == 2) && (prev_S == 3)) % transition in progress
        t_start = t(ii);
        in_transition = true;
        continue;
    end
    
    if ((cur_S == 1) && (prev_S == 2) && in_transition) % transition completed
        paths(end+1) = t(ii) - t_start;
        in_transition = false;
        continue;
    end
end
end

function [t,y] = interpolate_trajectory(t,y,frequency)
t0 = t(1);
tf = t(end);

tq = t0:frequency:tf;
yq = interp1(t,y,tq,"prev");

t = tq(:);
y = yq(:);
end


function mu = transition_mean(model)
k21 = model.k21;
k23 = model.k23;
mu = 1 ./ (k23);
end

function sigma = transition_standard_deviation(model)
k21 = model.k21;
k23 = model.k23;
sigma = sqrt(1 ./ (k23).^2);
end

function cv = transition_cv(model)
cv = transition_standard_deviation(model) ./ transition_mean(model);
end

function trace = trim_trace_to_length(trace, length, burnin, interpolate, frequency)
% input trace
time = trace.time;
state = trace.state;

% interpolate to 10Hz acquisition time
if interpolate
    [time,state] = interpolate_trajectory(time,state,1./frequency);
end

% get time-points
tf = (time >= burnin) & (time <= burnin+length);
time = time(tf);
state = state(tf);

% output trace
trace.time = time;
trace.state = state;
end


function mu = analytical_residence_means(model)
k12 = model.k12;
k21 = model.k21;
k23 = model.k23;
k32 = model.k32;

mu = zeros(3,1);
mu(1) = 1 ./ k12;
mu(2) = 1 ./ (k21 + k23);
mu(3) = 1 ./ k32;
end
