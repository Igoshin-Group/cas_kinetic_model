function model = adjust_model_energies(model, energies, C)
% adjust the rate constants of a model, keeping the fi's intact

    % extract on-target model state free-energies
    E1 = 0;
    E2 = energies.E2;
    E3 = energies.E3;
    E4 = energies.E4;
    E5 = energies.E5;
    E6 = energies.E6;
    E7 = energies.E7;
    E8 = energies.E8;

%     % extract off-target model state free-energies
%     E9 = energies.E9;
%     E10 = energies.E10;
%     E11 = energies.E11;
%     E12 = energies.E12;
%     E13 = energies.E13;
%     E14 = energies.E14;
    
    % extract on-target model transition state free-energies
    E1aR_dag = energies.E1aR_dag;
    E2R_dag = energies.E2R_dag;
    E5R_dag = energies.E5R_dag;
    E3R_dag = energies.E3R_dag;
    E7R_dag = energies.E7R_dag;
    E6R_dag = energies.E6R_dag;
    E9R_dag = energies.E9R_dag;
    E4R_dag = energies.E4R_dag;
    E8R_dag = energies.E8R_dag;
    EclvR_dag = energies.EclvR_dag;
    E1bR_dag = energies.E1bR_dag;
    
%     % extract off-target model transition state free-energies
%     E1aW_dag = energies.E1aW_dag;
%     E2W_dag = energies.E2W_dag;
%     E5W_dag = energies.E5W_dag;
%     E3W_dag = energies.E3W_dag;
%     E7W_dag = energies.E7W_dag;
%     E6W_dag = energies.E6W_dag;
%     E9W_dag = energies.E9W_dag;
%     E4W_dag = energies.E4W_dag;
%     E8W_dag = energies.E8W_dag;
%     EclvW_dag = energies.EclvW_dag;

    % other thermodynamic terms
    mu = 0; % Cas9 chemical potential, set to 0 for consistency with prior fit
    
    % Setting K*T = 1 implies energies are in units of K*T (for simplicity)
    k = 1; % Boltzmann constant
    T = 1; % system temperature
    
    w = 1; % frequency factor

    B = 1 / (k * T); % Boltzmann factor
    
    % compute rate constants via transition state theory
    % non-specific rates:
    k1b = w * exp(B * (E1 - E1bR_dag));
    k11b = w * exp(B * (E8 - E1bR_dag));
    
    % compute new on-target rates:
    k1a_R = w * exp(B * (E1 - E1aR_dag)) * exp(mu + log(C));
    k11a_R = w * exp(B * (E2 - E1aR_dag));

    k2_R = w * exp(B * (E2 - E2R_dag));
    k22_R = w * exp(B * (E3 - E2R_dag));

    k3_R = w * exp(B * (E3 - E3R_dag));
    k33_R = w * exp(B * (E5 - E3R_dag));

    k4_R = w * exp(B * (E5 - E4R_dag));
    k44_R = w * exp(B * (E7 - E4R_dag));

    k5_R = w * exp(B * (E2 - E5R_dag));
    k55_R = w * exp(B * (E4 - E5R_dag));

    k6_R = w * exp(B * (E4 - E6R_dag));
    k66_R = w * exp(B * (E6 - E6R_dag));

    k7_R = w * exp(B * (E4 - E7R_dag));
    k77_R = w * exp(B * (E5 - E7R_dag));

    k8_R = w * exp(B * (E6 - E8R_dag));
    k88_R = w * exp(B * (E7 - E8R_dag));

    k9_R = w * exp(B * (E8 - E9R_dag)) * exp(mu + log(C));
    k99_R = w * exp(B * (E4 - E9R_dag));

    kclv_R = w * exp(B * (E7 - EclvR_dag));
    
%     % off-target rates:
%     k1a_W = w * exp(B * (E1 - E1aW_dag)) * exp(mu + log(C));
%     k11a_W = w * exp(B * (E9 - E1aW_dag));
% 
%     k2_W = w * exp(B * (E9 - E2W_dag));
%     k22_W = w * exp(B * (E10 - E2W_dag));
% 
%     k3_W = w * exp(B * (E10 - E3W_dag));
%     k33_W = w * exp(B * (E12 - E3W_dag));
% 
%     k4_W = w * exp(B * (E12 - E4W_dag));
%     k44_W = w * exp(B * (E14 - E4W_dag));
% 
%     k5_W = w * exp(B * (E9 - E5W_dag));
%     k55_W = w * exp(B * (E11 - E5W_dag));
% 
%     k6_W = w * exp(B * (E11 - E6W_dag));
%     k66_W = w * exp(B * (E13 - E6W_dag));
% 
%     k7_W = w * exp(B * (E11 - E7W_dag));
%     k77_W = w * exp(B * (E12 - E7W_dag));
% 
%     k8_W = w * exp(B * (E13 - E8W_dag));
%     k88_W = w * exp(B * (E14 - E8W_dag));
% 
%     k9_W = w * exp(B * (E8 - E9W_dag)) * exp(mu + log(C));
%     k99_W = w * exp(B * (E11 - E9W_dag));
% 
%     kclv_W = w * exp(B * (E14 - EclvW_dag));
    
    % update model on-target rate constants
    model.k1a = k1a_R;
    model.k11a = k11a_R;
    model.k2 = k2_R;
    model.k22 = k22_R;
    model.k3 = k3_R;
    model.k33 = k33_R;
    model.k4 = k4_R;
    model.k44 = k44_R;
    model.k5 = k5_R;
    model.k55 = k55_R;
    model.k6 = k6_R;
    model.k66 = k66_R;
    model.k7 = k7_R;
    model.k77 = k77_R;
    model.k8 = k8_R;
    model.k88 = k88_R;
    model.k9 = k9_R;
    model.k99 = k99_R;
    model.kclv = kclv_R;
    
%     % update discrimination factors
%     model.f1a = k1a_W / k1a_R;
%     model.f11a = k11a_W / k11a_R;
%     model.f2 = k2_W / k2_R;
%     model.f22 = k22_W / k22_R;
%     model.f3 = k3_W / k3_R;
%     model.f33 = k33_W / k33_R;
%     model.f4 = k4_W / k4_R;
%     model.f44 = k44_W / k44_R;
%     model.f5 = k5_W / k5_R;
%     model.f55 = k55_W / k55_R;
%     model.f6 = k6_W / k6_R;
%     model.f66 = k66_W / k66_R;
%     model.f7 = k7_W / k7_R;
%     model.f77 = k77_W / k77_R;
%     model.f8 = k8_W / k8_R;
%     model.f88 = k88_W / k88_R;
%     model.f9 = k9_W / k9_R;
%     model.f99 = k99_W / k99_R;
%     model.fclv = kclv_W / kclv_R;
    
    % non-specific parameters
    model.k1b = k1b;
    model.k11b = k11b;
    model.f1b = 1;
    model.f11b = 1;    
end