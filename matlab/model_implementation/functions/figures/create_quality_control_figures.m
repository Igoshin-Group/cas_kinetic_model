function create_quality_control_figures(psets,dataset,name)
create_dissociation_constant_figure(psets,dataset,name);
save_figure(gcf, "../figures/", "qc/", "model-on-target-dissociation-constant");
create_mfpt_figure(psets,dataset,name);
save_figure(gcf, "../figures/", "qc/", "model-on-target-mfpt");
create_fret_state_occupancy_figure(psets,dataset,name);
save_figure(gcf, "../figures/", "qc/", "model-on-target-fret-probability");
create_fret_dwell_mean_figure(psets,dataset,name);
save_figure(gcf, "../figures/", "qc/", "model-on-target-fret-mean");
create_fret_dwell_cv_figure(psets,dataset,name);
save_figure(gcf, "../figures/", "qc/", "model-on-target-fret-cv");
% create_fret_splitting_probability_figure(psets,dataset);
% save_figure(gcf, "../figures/", "qc/", "model-dissociation-constant");
% saveas(gcf,"../splitting_probability.png");
end


function create_dissociation_constant_figure(psets,dataset,name)

colors = ["#527722", "#A2142F", "#0072BD"];

num_substrates = length(psets); % required to be same length as variants
num_variants = length(psets);

dissociation_constant_predicted = nan(num_substrates,num_variants);
dissociation_constant_observed = nan(num_substrates,num_variants);

for ii = 1:num_variants
    % get the observed points
    dissociation_constant_observed(1,ii) = dataset(ii).dissociation_constant;
    % get the fitted points
    dissociation_constant_predicted(1,ii) = analytical_dissociation_constant(psets(ii));
end

% generate plot
figure("Name", name);
hold on;

s = scatter(dissociation_constant_observed, dissociation_constant_predicted,32,"filled","MarkerEdgeColor","k");

for ii = 1:num_variants
    s(ii).MarkerFaceColor = colors(ii);
end

% add reference line for perfect correlation
ref = fplot(@(x) 1.25.*x, [0, 1e2]);
ref.Color = 'k';
ref.LineStyle = '--';
ref.DisplayName = "Reference Line";
ref.HandleVisibility = "off";

ref = fplot(@(x) 0.75.*x, [0, 1e2]);
ref.Color = 'k';
ref.LineStyle = '--';
ref.DisplayName = "Reference Line";
ref.HandleVisibility = "off";

ref = fplot(@(x) x, [0, 1e2]);
ref.Color = 'k';
ref.LineStyle = '--';
ref.LineWidth = 2;
ref.DisplayName = "Reference Line";
ref.HandleVisibility = "off";
hold off;

pbaspect([1,1,1]);
% xlim([0,25]);
% ylim([0,25]);
set(gca,"XScale", "log", "YScale", "log");


legend("WT","HF1","Enh");
xlabel("Dissociation Constant - Observed (nM)", "FontSize", 14);
ylabel("Dissociation Constant - Fitted (nM)", "FontSize", 14);
title("Dissociation Constant Fit Quality", "FontSize", 18);
hold off;

end

function create_mfpt_figure(psets,dataset,name)
colors = ["#527722", "#A2142F", "#0072BD"];

num_substrates = length(psets); % required to be same length as variants
num_variants = length(psets);

cleavage_mfpt_predicted = nan(num_substrates,num_variants);
cleavage_mfpt_observed = nan(num_substrates,num_variants);

for ii = 1:num_variants
    % get the observed points
    cleavage_mfpt_observed(1,ii) = dataset(ii).cleavage_mfpt;
    % get the fitted points
    cleavage_mfpt_predicted(1,ii) = compute_off_target_mfpt(psets(ii));
end

% generate plot
figure("Name",name);
hold on;

s = scatter(cleavage_mfpt_observed, cleavage_mfpt_predicted,32,"filled","MarkerEdgeColor","k");

for ii = 1:num_variants
    s(ii).MarkerFaceColor = colors(ii);
end


% add reference line for perfect correlation
ref = fplot(@(x) 1.25.*x, [0, 1e2]);
ref.Color = 'k';
ref.LineStyle = '--';
ref.DisplayName = "Reference Line";
ref.HandleVisibility = "off";

ref = fplot(@(x) 0.75.*x, [0, 1e2]);
ref.Color = 'k';
ref.LineStyle = '--';
ref.DisplayName = "Reference Line";
ref.HandleVisibility = "off";

ref = fplot(@(x) x, [0, 1e2]);
ref.Color = 'k';
ref.LineStyle = '--';
ref.LineWidth = 2;
ref.DisplayName = "Reference Line";
ref.HandleVisibility = "off";
hold off;

pbaspect([1,1,1]);
% xlim([0,25]);
% ylim([0,25]);
set(gca,"XScale", "log", "YScale", "log");


legend("WT","HF1","Enh");
xlabel("Cleavage MFPT - Observed (sec)", "FontSize", 14);
ylabel("Cleavage MFPT - Fitted (sec)", "FontSize", 14);
title("Cleavage MFPT Fit Quality", "FontSize", 18);
hold off;
end

function create_fret_state_occupancy_figure(psets,dataset,name)
colors = ["#527722", "#A2142F", "#0072BD"];
variant = ["WT", "HF1", "Enh"];

num_variants = length(psets);
num_states = 3;

fit_data = zeros(num_variants,num_states);
obs_data = zeros(num_variants,num_states);

for ii = 1:num_variants
    obs_data(ii,:) = [dataset(ii).fret.state_occupancies];
    fit_data(ii,:) = [compute_state_occupancies(psets(ii))];
end

% generate plot
fig = large_figure(12,9);
fig.Name = name;

for ii = 1:num_states
    subplot(1,3,ii);
    hold on;
    for jj = 1:num_variants
        s = scatter(obs_data(jj,ii), fit_data(jj,ii),32,"filled","MarkerEdgeColor","k","MarkerFaceColor",colors(jj), "DisplayName", variant(jj));
    end

    ref = fplot(@(x) x, [0, 1e0]);
    ref.Color = 'k';
    ref.LineStyle = '--';
    ref.LineWidth = 2;
    ref.DisplayName = "Reference Line";
    ref.HandleVisibility = "off";
    hold off;
    
    pbaspect([1,1,1]);
    set(gca,"XScale", "log", "YScale", "log");
    legend();
    xlabel("Probability - Observed (unitless)", "FontSize", 14);
    ylabel("Probability - Fitted (unitless)", "FontSize", 14);
    title("State "+num2str(ii), "FontSize", 18);    
end
end

function create_fret_dwell_mean_figure(psets,dataset,name)
colors = ["#527722", "#A2142F", "#0072BD"];
variant = ["WT", "HF1", "Enh"];

num_variants = length(psets);
num_states = 3;

fit_data = zeros(num_variants,num_states);
obs_data = zeros(num_variants,num_states);

for ii = 1:num_variants
    obs_data(ii,:) = [dataset(ii).fret.dwell_means];
    fit_data(ii,:) = [compute_dwell_means(psets(ii))];
end

% generate plot
fig = large_figure(12,9);
fig.Name = name;

for ii = 1:num_states
    subplot(1,3,ii);
    hold on;
    for jj = 1:num_variants
        s = scatter(obs_data(jj,ii), fit_data(jj,ii),32,"filled","MarkerEdgeColor","k","MarkerFaceColor",colors(jj), "DisplayName", variant(jj));
    end

    ref = fplot(@(x) x, [0, 1e0]);
    ref.Color = 'k';
    ref.LineStyle = '--';
    ref.LineWidth = 2;
    ref.DisplayName = "Reference Line";
    ref.HandleVisibility = "off";
    hold off;
    
    pbaspect([1,1,1]);
    set(gca,"XScale", "log", "YScale", "log");
    legend();
    xlabel("Dwell Mean - Observed (unitless)", "FontSize", 14);
    ylabel("Dwell Mean - Fitted (unitless)", "FontSize", 14);
    title("State "+num2str(ii), "FontSize", 18);    
end
end

function create_fret_dwell_cv_figure(psets,dataset,name)
colors = ["#527722", "#A2142F", "#0072BD"];
variant = ["WT", "HF1", "Enh"];

num_variants = length(psets);
num_states = 3;

fit_data = zeros(num_variants,num_states);
obs_data = zeros(num_variants,num_states);

for ii = 1:num_variants
    obs_data(ii,:) = [dataset(ii).fret.dwell_cv];
    fit_data(ii,:) = [compute_dwell_cv(psets(ii))];
end

% generate plot
fig = large_figure(12,9);
fig.Name = name;

for ii = 1:num_states
    subplot(1,3,ii);
    hold on;
    for jj = 1:num_variants
        s = scatter(obs_data(jj,ii), fit_data(jj,ii),32,"filled","MarkerEdgeColor","k","MarkerFaceColor",colors(jj), "DisplayName", variant(jj));
    end

    ref = fplot(@(x) x, [0, 1e0]);
    ref.Color = 'k';
    ref.LineStyle = '--';
    ref.LineWidth = 2;
    ref.DisplayName = "Reference Line";
    ref.HandleVisibility = "off";
    hold off;
    
    pbaspect([1,1,1]);
    set(gca,"XScale", "log", "YScale", "log");
    legend();
    xlabel("Dwell Variance - Observed (unitless)", "FontSize", 14);
    ylabel("Dwell Variance - Fitted (unitless)", "FontSize", 14);
    title("State "+num2str(ii), "FontSize", 18);    
end
end

function create_fret_splitting_probability_figure(psets,dataset)

end