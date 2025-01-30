function x = create_off_target_initial_condition(energies)
% fill in x for each Cas variant
x = zeros(16,2);
for ii = 1:length(energies)
    % microstate free-energies
    x(1,ii) = energies(ii).E9;
    x(2,ii) = energies(ii).E10;
    x(3,ii) = energies(ii).E11;
    x(4,ii) = energies(ii).E12;
    x(5,ii) = energies(ii).E13;
    x(6,ii) = energies(ii).E14;
    % transition state free-energies
    x(7,ii) = energies(ii).E1aW_dag;
    x(8,ii) = energies(ii).E9W_dag;
    x(9,ii) = energies(ii).E2W_dag;
    x(10,ii) = energies(ii).E5W_dag;
    x(11,ii) = energies(ii).E3W_dag;
    x(12,ii) = energies(ii).E7W_dag;
    x(13,ii) = energies(ii).E6W_dag;
    x(14,ii) = energies(ii).E4W_dag;
    x(15,ii) = energies(ii).E8W_dag;
    x(16,ii) = energies(ii).EclvW_dag;
end

% flatten x
x = reshape(x, 32, []);
end