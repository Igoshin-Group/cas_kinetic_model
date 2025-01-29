clc;
clear;
close all;

%% Load the HF1-Cas9 FRET data
hf1_fret_file_paths = ["../../data/fret/processed/hf1/hf1-0bp_traces.dat",...
                      "../../data/fret/processed/hf1/hf1-1bp_traces.dat",...
                      "../../data/fret/processed/hf1/hf1-2bp_traces.dat",...
                      "../../data/fret/processed/hf1/hf1-3bp_traces.dat",...
                      "../../data/fret/processed/hf1/hf1-4bp_traces.dat"];
% load FRET datasets
num_states = 3; % number of discrete FRET states (see Chen et al. 2017)
dt = 0.1; % acquisition time-step (see Chen et al. 2017)

hf1_fret_data = process_fret_data(hf1_fret_file_paths,num_states,dt);

%% Check for evidence of 1 -> 3 and 3 -> 1 FRET transitions
num_conditions = length(hf1_fret_data);
transition_matrix = zeros(num_states+1); % matrix to hold transitions
trace_lengths = [];
total_transitions = 0;
for ii = 1:num_conditions
    cur_condition = hf1_fret_data(ii);
    num_traces = length(cur_condition.traces);
    for jj = 1:num_traces

        % get the current trace
        cur_trace = cur_condition.traces(jj);

        % compute trace length
        trace_length = length(cur_trace.state_mean);
        trace_lengths(end+1) = trace_length;
        
        % get the viterbi state path
        state_path = get_fret_state_path(cur_trace);
        
        % check for transitions
        for kk = 1:trace_length - 2
            % check if a transition has occured
            if state_path(kk) ~= state_path(kk+1)
                % transition has occured, increment transition matrix
                cur_state = state_path(kk);
                next_state = state_path(kk+1);
                transition_matrix(cur_state,next_state) = transition_matrix(cur_state,next_state) + 1;
                total_transitions = total_transitions + 1;
            end
        end
    end
end

% convert trace lengths to units of time
trace_lengths = trace_lengths .* dt;

% convert transition matrix to probability matrix
transition_matrix = transition_matrix ./ total_transitions;

%
% generate trace length histogram
%
figure();
histogram(trace_lengths,"BinMethod","auto");
xlabel("Trace Length (sec)","FontSize",14);
ylabel("Frequency (#)","FontSize",14);
title("HF1-Cas9 FRET Trace Length Distribution","FontSize",18);

%
% generate transition probability figure
%
figure();
x = 1:num_states+1;
y = 1:num_states+1;
colormap("copper");
pcolor(x,y,transition_matrix');

% add colorbar to plot
color_bar = colorbar();
color_bar.Label.String = "Probability";
color_bar.Label.FontSize = 14;

% add axis labels and title
xlabel("Initial FRET State (#)","FontSize",14);
ylabel("Final FRET State (#)","FontSize",14);
title("HF1-Cas9 Transition Probability Matrix","FontSize",18);

% set y-axis direction
set(gca,'YDir','reverse');

% fix axis tick labels
xtick_pos = [1.5,2.5,3.5];
ytick_pos = [1.5,2.5,3.5];
xticks(xtick_pos);
yticks(ytick_pos);
tick_labels = {'1','2','3'};
xticklabels(tick_labels);
yticklabels(tick_labels);

% write actually transition probability in each square
for ii = 1:length(xticks)
    for jj = 1:length(yticks)
        Pij = num2str(round(transition_matrix(ii,jj),3));
        xi = xtick_pos(ii);
        yi = ytick_pos(jj);

        text(xi,yi,Pij,'Color','white','FontSize',14);
    end
end

%% Check the HNH domain dwell time distribution and moments
dwell_times = cell(length(hf1_fret_data),num_states);
num_conditions = length(hf1_fret_data);
for ii = 1:num_conditions
    cur_condition = hf1_fret_data(ii);
    num_traces = length(cur_condition.traces);
    for jj = 1:num_traces

        % get the current trace
        cur_trace = cur_condition.traces(jj);

        % compute trace length
        trace_length = length(cur_trace.state_mean);
        trace_lengths(end+1) = trace_length;
        
        % get the viterbi state path
        state_path = get_fret_state_path(cur_trace);
        
        % determine fret dwell times
        counter = 0;
        current_dwells = cell(num_states,1);
        for kk=1:length(state_path)-2
            counter = counter + 1;
            % check if a transition has occured
            if state_path(kk) ~= state_path(kk+1)
                % transition has occured
                current_dwells{state_path(kk)}(end+1) = counter*dt; % convert to units of time
                counter = 0; % reset the counter
            end
        end
    
        % drop the first and last transition (we do not know start and end)
        for kk=1:num_states
            if length(current_dwells{kk}) > 2
                current_dwells{kk}(1) = []; % drop first
                current_dwells{kk}(end) = []; % drop last
            end
            % save current_dwells to dwell_times structure
            dwell_times{ii,kk} = [dwell_times{ii,kk}, current_dwells{kk}];
        end
    end
end

%
% plot residence time densities
%

line_colors = ["black","red","blue"];
line_styles = ["-","--","-."];
num_mismatches = [0,1,2,3,4];

states = {'State 1', 'State 2', 'State 3'};
means = zeros(num_conditions,num_states);
vars = zeros(num_conditions,num_states);
skews = zeros(num_conditions,num_states);

nboot = total_transitions;
means_ci = zeros(num_conditions,num_states,2);
vars_ci = zeros(num_conditions,num_states,2);
skews_ci = zeros(num_conditions,num_states,2);
for ii = 1:num_conditions
    figure();
    hold on;
    for jj = 1:num_states
        % get dwell times for current condition and state
        dwell_times_ij = dwell_times{ii,jj};

        % compute dwell time moments
        mean_ij = mean(dwell_times_ij);
        var_ij = var(dwell_times_ij);
        skew_ij = skewness(dwell_times_ij);
        boot_data = dwell_times_ij(~isnan(dwell_times_ij));

        means(ii,jj) = mean_ij;
        vars(ii,jj) = var_ij;
        skews(ii,jj) = skew_ij;
        
        % skip if empty (to avoid breaking bootci)
        if isempty(dwell_times_ij)
            fprintf("Skipping state: %s\n",states{jj});
            continue;
        end

        means_ci(ii,jj,:) = bootci(nboot,@mean,boot_data);
        vars_ci(ii,jj,:) = bootci(nboot,@var,boot_data);
        skews_ci(ii,jj,:) = bootci(nboot,@skewness,boot_data);
        
        % compute kernel density estimate and plot
        [f,xi] = ksdensity(dwell_times_ij,'Support',[0,inf],'BoundaryCorrection','reflection');
        plot(xi,f,"Color",line_colors(jj),"LineStyle","-","LineWidth",2);
    end
    % add axis labels and title
    title_label = join(["HF1-Cas9 HNH Residence Times, MM=",num2str(num_mismatches(ii))]);
    title(title_label,"FontSize",18);
    xlabel("Residence Time (sec)","FontSize",14);
    ylabel("Probability Density (units)","FontSize",14);
    % add legend
    legend("State 1", "State 2", "State 3");
    % set axis limits
    xlim([0,10]);
end

%
% plot residence time moments
%
figure();

width = 12;
height = 9;
set(gcf,'units','inches','position',[1,1,width,height])
subplot(1,3,1);
X = categorical({'0','1','2','3','4'});
X = reordercats(X,{'0','1','2','3','4'});
Y = means;
error_low = means_ci(:,:,1);
error_high = means_ci(:,:,2);

b = bar(X,Y);
hold on;
[num_groups,num_bars] = size(Y);

xpos = nan(num_groups, num_bars);
for ii = 1:num_bars
    xpos(:,ii) = b(ii).XEndPoints;
end

errorbar(xpos,Y,error_low,error_high,[],[],"k","LineStyle","none");



xlabel("Num. Mismatches (\#)","FontSize",14,"Interpreter","latex");
ylabel("Residence Time Mean ($\mu$)","FontSize",14,"Interpreter","latex");
title("Mean","FontSize",18,"Interpreter","latex");
lgd = legend("State 1","State 2", "State 3");
lgd.FontName = "Times";
pbaspect([1 1 1]);

subplot(1,3,2);
X = categorical({'0','1','2','3','4'});
X = reordercats(X,{'0','1','2','3','4'});
Y = vars;
error_low = vars_ci(:,:,1);
error_high = vars_ci(:,:,2);

b = bar(X,Y);
hold on;
[num_groups,num_bars] = size(Y);

xpos = nan(num_groups, num_bars);
for ii = 1:num_bars
    xpos(:,ii) = b(ii).XEndPoints;
end

errorbar(xpos,Y,error_low,error_high,[],[],"k","LineStyle","none");

xlabel("Num. Mismatches (\#)","FontSize",14,"Interpreter","latex");
ylabel("Residence Time Mean ($\sigma^2$)","FontSize",14,"Interpreter","latex");
title("Variance","FontSize",18,"Interpreter","latex");
lgd = legend("State 1","State 2", "State 3");
lgd.FontName = "Times";
pbaspect([1 1 1]);

subplot(1,3,3);

X = categorical({'0','1','2','3','4'});
X = reordercats(X,{'0','1','2','3','4'});
Y = skews;
error_low = skews_ci(:,:,1);
error_high = skews_ci(:,:,2);

b = bar(X,Y);
hold on;
[num_groups,num_bars] = size(Y);

xpos = nan(num_groups, num_bars);
for ii = 1:num_bars
    xpos(:,ii) = b(ii).XEndPoints;
end

errorbar(xpos,Y,error_low,error_high,[],[],"k","LineStyle","none");

xlabel("Num. Mismatches (\#)","FontSize",14,"Interpreter","latex");
ylabel("Residence Time Skewness ($\tilde{\mu}_3$)","FontSize",14,"Interpreter","latex");
title("Skewness","FontSize",18,"Interpreter","latex");
lgd=legend("State 1","State 2", "State 3");
lgd.FontName = "Times";
pbaspect([1 1 1]);

%% Check deviation from exponential distribution
CV = vars ./ means;
Z = CV - 1;

figure();
bar(X,Z);

title("Deviation from Exponential", "FontSize", 18);
xlabel("Mismatch (#)", "FontSize", 14);
ylabel("Deviation, Z", "FontSize", 14);