function [model, energies] = generate_model_from_energies(x, C)
    % re-shape x to extract the parameters for each variant
    x = reshape(x,[],2); % x should have shape (34,2)
      
    % compute model parameter sets
    [off_model_wt, energies_wt] = helper_function(x(:,1), C);
    [off_model_hf1, energies_hf1] = helper_function(x(:,2), C);
    
    on_model_wt = set_on_target_discrimination(off_model_wt);
    on_model_hf1 = set_on_target_discrimination(off_model_hf1);

    % package model and energies for output
    model(1,1) = on_model_wt;
    model(2,1) = on_model_hf1;
    model(1,2) = off_model_wt;
    model(2,2) = off_model_hf1;
    
    energies(1) = energies_wt;
    energies(2) = energies_hf1;
end


%% Local Functions
function [model, energies] = helper_function(x, C)

% x = vector of energy terms (in units of kT)
% C = [Cas9] (nM)
% Vector: [E2=1, E3=2, E4=3, E5=4, E6=5, E7, E8, E9, E10, E11, E12, E13, E14=13, ...
%          E1_2=14, E2_3=15, E2_4=16, E4_5=17, E3_5=18, E4_6=19, E5_7=20, E6_7=21, E1_8=22, E4_8=23, E7p=24, ...
%          E1_9=25, E9_10=26, E9_11=27, E10_12=28, E11_12=29, E11_13=30, E12_14=31, E13_14=32, E8_11=33, E14p=34]
    % state free-energies
    E1 = 0;     % we can take E1 to be the reference energy
    E2 = x(1);
    E3 = x(2);
    E4 = x(3);
    E5 = x(4);
    E6 = x(5);
    E7 = x(6);
    E8 = x(7);
    E9 = x(8);
    E10 = x(9);
    E11 = x(10);
    E12 = x(11);
    E13 = x(12);
    E14 = x(13);

    % transition state free-energies
    % E1_2=14, E2_3=15, E2_4=16, E4_5=17, E3_5=18, E4_6=19, E5_7=20, E6_7=21, E1_8=22, E4_8=23, E7p=24, ...
    %          E1_9=25, E9_10=26, E9_11=27, E10_12=28, E11_12=29, E11_13=30, E12_14=31, E13_14=32, E8_11=33, E14p=34]
    E1aR_dag = x(14);
    E2R_dag = x(15);
    E5R_dag = x(16);
    E3R_dag = x(17);
    E7R_dag = x(18);
    E6R_dag = x(19);
    E9R_dag = x(20);
    E4R_dag = x(21);
    E8R_dag = x(22);
    EclvR_dag = x(23);
    E1bR_dag = x(24);
    E1aW_dag = x(25);
    E9W_dag = x(26);
    E2W_dag = x(27);
    E5W_dag = x(28);
    E3W_dag = x(29);
    E7W_dag = x(30);
    E6W_dag = x(31);
    E4W_dag = x(32);
    E8W_dag = x(33);
    EclvW_dag = x(34);

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

    % on-target rates:
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

    % off-target rates:
    k1a_W = w * exp(B * (E1 - E1aW_dag)) * exp(mu + log(C));
    k11a_W = w * exp(B * (E9 - E1aW_dag));

    k2_W = w * exp(B * (E9 - E2W_dag));
    k22_W = w * exp(B * (E10 - E2W_dag));

    k3_W = w * exp(B * (E10 - E3W_dag));
    k33_W = w * exp(B * (E12 - E3W_dag));

    k4_W = w * exp(B * (E12 - E4W_dag));
    k44_W = w * exp(B * (E14 - E4W_dag));

    k5_W = w * exp(B * (E9 - E5W_dag));
    k55_W = w * exp(B * (E11 - E5W_dag));

    k6_W = w * exp(B * (E11 - E6W_dag));
    k66_W = w * exp(B * (E13 - E6W_dag));

    k7_W = w * exp(B * (E11 - E7W_dag));
    k77_W = w * exp(B * (E12 - E7W_dag));

    k8_W = w * exp(B * (E13 - E8W_dag));
    k88_W = w * exp(B * (E14 - E8W_dag));

    k9_W = w * exp(B * (E8 - E9W_dag)) * exp(mu + log(C));
    k99_W = w * exp(B * (E11 - E9W_dag));

    kclv_W = w * exp(B * (E14 - EclvW_dag));


    % =================================================================================================
    %   create the parameters struct
    % =================================================================================================
    model = struct();

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

    model.f1a = k1a_W / k1a_R;
    model.f11a = k11a_W / k11a_R;
    model.f2 = k2_W / k2_R;
    model.f22 = k22_W / k22_R;
    model.f3 = k3_W / k3_R;
    model.f33 = k33_W / k33_R;
    model.f4 = k4_W / k4_R;
    model.f44 = k44_W / k44_R;
    model.f5 = k5_W / k5_R;
    model.f55 = k55_W / k55_R;
    model.f6 = k6_W / k6_R;
    model.f66 = k66_W / k66_R;
    model.f7 = k7_W / k7_R;
    model.f77 = k77_W / k77_R;
    model.f8 = k8_W / k8_R;
    model.f88 = k88_W / k88_R;
    model.f9 = k9_W / k9_R;
    model.f99 = k99_W / k99_R;
    model.fclv = kclv_W / kclv_R;

    % non-specific parameters
    model.k1b = k1b;
    model.k11b = k11b;
    model.f1b = 1;
    model.f11b = 1;
    
    % =================================================================================================
    %   create the energies struct
    % =================================================================================================
    energies.E1 = 0;
    energies.E2 = E2;
    energies.E3 = E3;
    energies.E4 = E4;
    energies.E5 = E5;
    energies.E6 = E6;
    energies.E7 = E7;
    energies.E8 = E8;
    energies.E9 = E9;
    energies.E10 = E10;
    energies.E11 = E11;
    energies.E12 = E12;
    energies.E13 = E13;
    energies.E14 = E14;
    
    energies.E1aR_dag = E1aR_dag;
    energies.E2R_dag = E2R_dag;
    energies.E5R_dag = E5R_dag;
    energies.E3R_dag = E3R_dag;
    energies.E7R_dag = E7R_dag;
    energies.E6R_dag = E6R_dag;
    energies.E9R_dag = E9R_dag;
    energies.E4R_dag = E4R_dag;
    energies.E8R_dag = E8R_dag;
    energies.EclvR_dag = EclvR_dag;
    energies.E1bR_dag = E1bR_dag;
    energies.E1aW_dag = E1aW_dag;
    energies.E9W_dag = E9W_dag;
    energies.E2W_dag = E2W_dag;
    energies.E5W_dag = E5W_dag;
    energies.E3W_dag = E3W_dag;
    energies.E7W_dag = E7W_dag;
    energies.E6W_dag = E6W_dag;
    energies.E4W_dag = E4W_dag;
    energies.E8W_dag = E8W_dag;
    energies.EclvW_dag = EclvW_dag;
end

function model = set_on_target_discrimination(model)
    model.f1a = 1;
    model.f11a = 1;
    model.f2 = 1;
    model.f22 = 1;
    model.f3 = 1;
    model.f33 = 1;
    model.f4 = 1;
    model.f44 = 1;
    model.f5 = 1;
    model.f55 = 1;
    model.f6 = 1;
    model.f66 = 1;
    model.f7 = 1;
    model.f77 = 1;
    model.f8 = 1;
    model.f88 = 1;
    model.f9 = 1;
    model.f99 = 1;
    model.fclv = 1;
end