function visualize_off_target_substrate(model, data)

    % compute model observables
    Kd = analytical_dissociation_constant(model);
    ndABA = compute_normalized_binding_affinity(model);
    T = compute_off_target_mfpt(model);
    P = compute_state_occupancies(model);
    Dm = compute_dwell_means(model);
    Dv = compute_dwell_variance(model);

    % extract experimental data
    ndABA_obs = data.normalized_binding;
    T_obs = data.cleavage_mfpt;
    P_obs = [data.fret.state_occupancies];
    Dm_obs = [data.fret.dwell_means];
    Dv_obs = [data.fret.dwell_variance];
    
    % Create figure
    figure_name = "Substrate " + num2str(data.substrate_id) + " | " + ...
                  "Num. Distal Mismatches = " + num2str(data.distal_mismatches);
    f = figure("Name", figure_name, "Renderer", "painters", "Units", "inches", "NumberTitle", "off");
    
    % Resize figure
    width = 12; % inches
    height = 10; % inches
    pos = f.Position;
    x0 = pos(1);
    y0 = pos(2);
    x1 = width;
    y1 = height;
    pos = [x0, y0, x1, y1];
    f.Position = pos;
    
    % Plot: Normalized Binding Affinity
    subplot(4,4,[1,5]);
    data = [ndABA_obs; ndABA];
    x = categorical({'Obs.', 'Fit.'});
    x = reordercats(x,{'Obs.', 'Fit.'});
    bar(x,data);
    ylabel("\DeltaABA (unitless)")
    title("Binding");
    
    % Plot: Cleavage MFPT
    subplot(4,4,[2,6]);
    data = [T_obs; T];
    x = categorical({'Obs.', 'Fit.'});
    x = reordercats(x,{'Obs.', 'Fit.'});
    bar(x,data);
    ylabel("\tau (sec)")
    title("MFPT");
    
    % Plot: HNH State Occupancies
    subplot(4,4,[3,4,7,8]);
    data = [P_obs, P];
    x = categorical({'State 1', 'State 2', 'State 3'});
    x = reordercats(x,{'State 1', 'State 2', 'State 3'});
    bar(x,data);
    ylabel("P (unitless)")
    legend("Obs","Fit");
    title("State Occupancy");
    
    % Plot: HNH Dwell Means
    subplot(4,4,[9,10,13,14]);
    data = [Dm_obs, Dm];
    x = categorical({'State 1', 'State 2', 'State 3'});
    x = reordercats(x,{'State 1', 'State 2', 'State 3'});
    bar(x,data);
    ylabel("\mu (sec)")
    legend("Obs","Fit");
    title("Dwell Mean");
    
    % Plot: HNH Dwell Variance
    subplot(4,4,[11,12,15,16]);
    data = [Dv_obs, Dv];
    x = categorical({'State 1', 'State 2', 'State 3'});
    x = reordercats(x,{'State 1', 'State 2', 'State 3'});
    bar(x,data);
    ylabel("\sigma (sec^{2})")
    legend("Obs","Fit");
    title("Dwell Variance");
end