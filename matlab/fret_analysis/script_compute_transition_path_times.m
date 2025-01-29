clc;
clear;
close all;

%{
    This script is intended to look into whether the transition path time
    analysis put forward in Satija et al. 2020 (PNAS) can be used to show
    that the Cas9 HNH domain dynamics are mechanically uncoupled from the
    process of R-loop formation.
%}

clc;
clear;
close all;

%% Load the FRET data
fret_file_paths = ["../../data/fret/processed/wt/wt-0bp_traces.dat",...
                   "../../data/fret/processed/hf1/hf1-0bp_traces.dat"];

% load FRET datasets
num_states = 3; % number of discrete FRET states (see Chen et al. 2017)
dt = 0.1; % acquisition time-step, sec (see Chen et al. 2017)

fret_data = process_fret_data(fret_file_paths,num_states,dt);

%% Compute the HNH domain transition path times (state 1 to state 3)
path_times = cell(length(fret_data),1);
num_conditions = length(fret_data);
num_path_times = zeros(length(fret_data),1);
trace_lengths = cell(num_conditions,1);
for ii = 1:num_conditions
    cur_condition = fret_data(ii);
    num_traces = length(cur_condition.traces);
    condition_paths = [];
    for jj = 1:num_traces
        % get the current trace
        cur_trace = cur_condition.traces(jj);
        
        % compute trace length
        trace_length = length(cur_trace.state_mean);
        
        trace_lengths{ii}(end+1) = trace_length;
        
        % get the Viterbi state path
        state_path = get_fret_state_path(cur_trace);
        
        % count transition path lengths
        counter = 0;
        in_transition = false;
        
        % check if we start in state 2
        trash_first_path = false;
        if state_path(1) == 2
            trash_first_path = true;
        end
        
        current_paths = [];
        for kk = 1:length(state_path) - 2
            cur_state = state_path(kk);
            if cur_state == 1 && ~in_transition
                in_transition = true;
                counter = 0;
                continue;
            end
            
            if cur_state == 1 && in_transition
               counter = 0; % reset the counter because we have left the well 
               continue;
            end
            
            if cur_state ~= 3 && in_transition
                counter = counter + 1; % increment the counter
                continue;
            end
            
            if cur_state == 3 && in_transition
                in_transition = false;
                counter = counter + 1; % increment the counter
                current_paths(end+1) = counter*dt;
                counter = 0;
                continue;
            end
        end
        if length(current_paths) >= 2 && trash_first_path
            current_paths(1) = [];
        end
        condition_paths = [condition_paths, current_paths];
    end
    path_times{ii} = condition_paths;
    num_path_times(ii) = length(path_times{ii});
end

% compute statistics

cv_func = @(x) std(x) / mean(x);

mu_mean = zeros(size(path_times));
ci_mean = zeros(length(path_times),2);

mu_std = zeros(size(path_times));
ci_std = zeros(length(path_times),2);

mu_cv = zeros(size(path_times));
ci_cv = zeros(length(path_times),2);

nboot = 250;%10000;

for ii = 1:length(path_times)
    mu_mean(ii) = mean(path_times{ii});
    mu_std(ii) = std(path_times{ii});
    mu_cv(ii) = cv_func(path_times{ii});
    
    ci_mean(ii,:) = bootci(nboot,@(x)mean(x),path_times{ii});
    ci_std(ii,:) = bootci(nboot,@(x)std(x),path_times{ii});
    ci_cv(ii,:) = bootci(nboot,@(x)cv_func(x),path_times{ii});
end

%% Generate bar plot for transition path moments
figure("Renderer","painters");
subplot(1,2,1);
Xlabel = ["\mu", "\sigma"];
X = categorical(Xlabel);
X = reordercats(X, Xlabel);
data = [mu_mean, mu_std]';
b = bar(X, data);
hold on;

count = 1;
for ii = 1:2 % variant
    err_data = [ci_mean(ii,:); ci_std(ii,:)];
    xpos = [b.XEndPoints];
    %xpos = reshape(xpos, 2, 2);
    
    mu = [mu_mean(ii), mu_std(ii)];
    ci = err_data;
    
    for jj = 1:2 % statistic
        lower_bar = mu(jj) - ci(1,jj);
        upper_bar = ci(1,2) - mu(jj);

        errorbar(xpos(count), mu(jj),lower_bar,upper_bar,[],[],"ko");
        count = count + 1;
    end
end

%legend("WT", "HF1");
%legend("Location","north");
ylabel("Value (sec)");
pbaspect([1,1,1]);
ylim([0,2]);

subplot(1,2,2);

labels = {'WT','HF1'};
condition_labels = categorical(labels);
condition_labels = reordercats(condition_labels,labels);

b = plot_with_confidence_intervals(condition_labels, mu_cv, ci_cv);


% 
% data = [mu_cv]';
% Xlabel = ["CV"];
% X = categorical(Xlabel);
% X = reordercats(X, Xlabel);
% b = bar(X, data);
hold on;
yline(1, "r--");

ylabel("Value (unitless)");
pbaspect([1,1,1]);
%saveas(gcf,"../transition_path_moments.eps","epsc");

%% Generate plot for transition path distributions
variant = ["wt", "hf1"];
for ii = 1:2 % loop over variants
    figure("Renderer","painters");
    T = path_times{ii};
    histogram(T,"Normalization","probability","BinWidth",0.5);
    ylim([0,1]);
    xlim([0,10]);
    xlabel("Transition Path Time (sec)");
    ylabel("Probability");
    pbaspect([1,1,1]);
    %saveas(gcf,"../transition_path_distribution_"+variant(ii)+".eps","epsc");
end


%%
figure();
labels = {'WT','HF1'};
condition_labels = categorical(labels);
condition_labels = reordercats(condition_labels,labels);

b = plot_with_confidence_intervals(condition_labels, mu_cv, ci_cv);


% b = bar(condition_labels,mu_cv);
% hold on;
% 
% xpos = b.XData;
% 
% errorbar(xpos,mu_cv,abs(mu_cv - ci_cv(:,1)),abs(mu_cv - ci_cv(:,2)),[],[],"ko");

yline(1,"r--","LineWidth",2);
xlabel("Cas9 Variant");
ylabel("Coefficient of Variation, CV (unitless)");
title("CV of Transition Path Times");


% figure();
% b2 = plot_with_confidence_intervals(condition_labels, mu_mean, ci_mean);
% 
% xlabel("Cas9 Variant","FontSize",14);
% ylabel("Mean, \mu (sec)","FontSize", 14);
% title("Transition Path Means","FontSize",18);
% 
% figure();
% b3 = plot_with_confidence_intervals(condition_labels, mu_std, ci_std);
% 
% xlabel("Cas9 Variant","FontSize",14);
% ylabel("Standard Deviation, \sigma (sec)","FontSize", 14);
% title("Transition Path Standard Deviation","FontSize",18);

%{
    Current conclusion: the CVs suggest that the transition path
    for the HNH domain may be multi-dimensional (i.e.: no mechanical
    coupling between the R-loop and the HNH domain)
%}

fprintf("Num. State Paths Detected\n");
disp(num_path_times)

pbaspect([1,1,1]);

set(gcf,"renderer","Painters");
% %saveas(gca,"../transition_path_cv.eps");


%% Save WT trace lengths to file
fid = fopen("./wt_0mm_trace_lengths.txt", "w");
trace_lengths_wt = trace_lengths{1};

% correct trace lengths based on acquisition frequency
trace_lengths_wt = trace_lengths_wt .* dt;

for ii = 1:length(trace_lengths_wt)
    fprintf(fid, "%i\n", trace_lengths_wt(ii)); 
end

fclose(fid);

%% MC Simulation
% simulation settings
t_end = 5000;   % final time-point (sec)
N = 1;         % number of trajectories (#)

rng(10);

num_reactions = 14;
propensities = zeros(num_reactions,1); % pre-allocate for speed

num_states = 6;
initial_state = zeros(num_states,1);
initial_state(1) = 1;

trajectories = {};

% load parameters fopen("./model-02c-48-wt-4.json");
%fid = fopen("./model_02g_wt_best_fit.json");
raw = fread(fid,inf);
str = char(raw');
fclose(fid);
parameter_set = jsondecode(str);

% extract parameters
k2 = parameter_set.k2;
k22 = parameter_set.k22;
k3 = parameter_set.k3;
k33 = parameter_set.k33;
k4 = parameter_set.k4;
k44 = parameter_set.k44;
k5 = parameter_set.k5;
k55 = parameter_set.k55;
k6 = parameter_set.k6;
k66 = parameter_set.k66;
k7 = parameter_set.k7;
k8 = parameter_set.k8;

% compute parameters from detailed balance
k77 = k7 .* ((k33 .* k22 .* k5) ./ (k3 .* k2 .* k55));
k88 = k8 .* ((k44 .* k33 .* k22 .* k5 .* k6) ./ (k4 .* k3 .* k2 .* k55 .* k66));

% simulation loop
for ii = 1:N
    
    % reset for the next trajectory
    t = 0;
    state = initial_state;
    
    time = [t];
    state_trajectory = [state];
    
    % simulate the trajectory
    while (t <= t_end)
    
        % compute reaction propensities
        propensities(1) = k2 * state(1);
        propensities(2) = k22 * state(2);
        propensities(3) = k3 * state(2);
        propensities(4) = k33 * state(4);
        propensities(5) = k4 * state(4);
        propensities(6) = k44 * state(6);
        propensities(7) = k5 * state(1);
        propensities(8) = k55 * state(3);
        propensities(9) = k6 * state(3);
        propensities(10) = k66 * state(5);
        propensities(11) = k7 * state(3);
        propensities(12) = k77 * state(4);
        propensities(13) = k8 * state(5);
        propensities(14) = k88 * state(6);
        
        total_propensity = sum(propensities);
        
        % compute dt
        r1 = rand();
        dt = (1 / total_propensity) * log(1/r1);
        
        % determine the next reaction
        r2 = rand();
        threshold = r2 * total_propensity;
        
        reaction = 0;
        for jj = 1:length(propensities)
            temp = 0;
            for kk = 1:jj
                temp = temp + propensities(kk);
            end
            if temp > threshold
                reaction = jj;
                break;
            end
        end
        
        % update the system state
        switch reaction
            case 1
                state(1) = state(1) - 1;
                state(2) = state(2) + 1;
            case 2
                state(2) = state(2) - 1;
                state(1) = state(1) + 1;
            case 3
                state(2) = state(2) - 1;
                state(4) = state(4) + 1;
            case 4
                state(4) = state(4) - 1;
                state(2) = state(2) + 1;
            case 5
                state(4) = state(4) - 1;
                state(6) = state(6) + 1;
            case 6
                state(6) = state(6) - 1;
                state(4) = state(4) + 1;
            case 7
                state(1) = state(1) - 1;
                state(3) = state(3) + 1;
            case 8
                state(3) = state(3) - 1;
                state(1) = state(1) + 1;
            case 9
                state(3) = state(3) - 1;
                state(5) = state(5) + 1;
            case 10
                state(5) = state(5) - 1;
                state(3) = state(3) + 1;
            case 11
                state(3) = state(3) - 1;
                state(4) = state(4) + 1;
            case 12
                state(4) = state(4) - 1;
                state(3) = state(3) + 1;
            case 13
                state(5) = state(5) - 1;
                state(6) = state(6) + 1;
            case 14
                state(6) = state(6) - 1;
                state(5) = state(5) + 1;
            otherwise
                error("Invalid reaction: %i\n", reaction)
        end
        
        % update simulation time
        t = t + dt;
        
        % save current system state
        time = [time, t];
        state_trajectory = [state_trajectory, state];
    end
    
    trajectory.time = time;
    trajectory.state = state_trajectory;
    trajectories{end+1} = trajectory;
end

fret_state = determine_fret_state(trajectory);
sim_dwell_times = determine_state_dwell_times(trajectory);
distribution = determine_transition_path_distribution(trajectory);



dwell_mean = [mean(sim_dwell_times{1}), mean(sim_dwell_times{2}), mean(sim_dwell_times{3})];
dwell_variance = [var(sim_dwell_times{1}), var(sim_dwell_times{2}), var(sim_dwell_times{3})];

mu_x = mean(distribution);
std_x = std(distribution);
cv_x = std_x / mu_x;

fprintf("Mean: %f Std. Dev.: %f CV: %f\n", mu_x, std_x, cv_x);

%% Generate plot comparing MC and data
figure();
labels = {'WT'};
condition_labels = categorical(labels);
condition_labels = reordercats(condition_labels,labels);

cv_data = [mu_cv(1)];
std_data = [ci_cv(1,:)];

b = plot_with_confidence_intervals(condition_labels, cv_data, std_data);

hold on;

plot(1,cv_x,"r*","LineWidth",2,"DisplayName","Monte Carlo");

ylim([0.95,1.5]);

xlabel("Variant","FontSize",14);
ylabel("CV (\sigma/\mu)","FontSize",14);
title("Transition Path Distribution","FontSize",14);
legend("Data", "",  "Simulation");
%saveas(gca,"./model-02g-transition-path-simulation.png");

%% Generate plot for WT-statistics
figure();
labels = {'\mu', '\sigma', 'CV'};
condition_labels = categorical(labels);
condition_labels = reordercats(condition_labels,labels);

bar_data = [mu_mean(1), mu_std(1), mu_cv(1)]';
ci_data = [ci_mean(1,:); ci_std(1,:); ci_cv(1,:)];

b = plot_with_confidence_intervals(condition_labels, bar_data, ci_data);

xlabel("Statistic","FontSize",14);
ylabel("Value","FontSize",14);
title("WT On-Target Transition Path Statistics (Data)","FontSize",14);
ylim([0,5]);
%saveas(gcf,"./wt-on-target-transition-statistics-with-ci.png","png");

%% Local Functions
function b = plot_with_confidence_intervals(x,mu,ci)


b = bar(x,mu);
hold on;

xpos = b.XData;

lower_bar = mu - ci(:,1);
upper_bar = ci(:,2) - mu;

errorbar(xpos,mu,lower_bar,upper_bar,[],[],"ko");
end

function fret_state = determine_fret_state(trajectory)
    state_trajectory = trajectory.state;
    trajectory_length = length(state_trajectory);
    fret_state = zeros(trajectory_length, 1);
    for ii = 1:trajectory_length
        state = find(state_trajectory(:,ii) == 1);
        if ((state == 1) || (state == 2))
            fret_state(ii) = 1;
        elseif ((state == 3) || (state == 4))
            fret_state(ii) = 2;
        else
            fret_state(ii) = 3;
        end
    end
end

function dwell_times = determine_state_dwell_times(trajectory)
    time = trajectory.time;
    state_trajectory = trajectory.state;
    trajectory_length = length(state_trajectory);
    fret_state = determine_fret_state(trajectory);
    dwell_times = cell(3,1);
    prev_time = time(1);
    for ii = 2:trajectory_length
        state = fret_state(ii);
        prev_state = fret_state(ii-1);
        if state ~= prev_state
            dwell_times{prev_state}(end+1) = time(ii) - prev_time;
            prev_time = time(ii);
        end
    end
end

function distribution = determine_transition_path_distribution(trajectory)
    time = trajectory.time;
    state_trajectory = trajectory.state;
    trajectory_length = length(state_trajectory);
    fret_state = determine_fret_state(trajectory);
    distribution = [];
    prev_time = time(1);
    in_transition = false;
    first_transition_seen = false;
    
    % get forward transition paths
    for ii = 2:trajectory_length
        state = fret_state(ii);
        
        if ((state == 1) && ~in_transition)
            prev_time = time(ii);
            in_transition = true;
            continue;
        end
        
        if ((state == 1) && in_transition)
            prev_time = time(ii);
            continue;
        end
        
        if ((state == 2) && in_transition)
            continue;
        end
        
        if ((state == 3) && in_transition)
            in_transition = false;
            if (first_transition_seen)
                distribution(end+1) = time(ii) - prev_time;
            end
            
            first_transition_seen = true;
        end
    end
    
    % get reverse transition paths
    for ii = 2:trajectory_length
        state = fret_state(ii);
        
        if ((state == 3) && ~in_transition)
            prev_time = time(ii);
            in_transition = true;
            continue;
        end
        
        if ((state == 3) && in_transition)
            prev_time = time(ii);
            continue;
        end
        
        if ((state == 2) && in_transition)
            continue;
        end
        
        if ((state == 1) && in_transition)
            in_transition = false;
            if (first_transition_seen)
                distribution(end+1) = time(ii) - prev_time;
            end
            
            first_transition_seen = true;
        end
    end
end
