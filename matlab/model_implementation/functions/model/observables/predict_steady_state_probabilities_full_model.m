function Pss = predict_steady_state_probabilities_full_model(model, x)
    Mf = create_forward_matrix_full_model(model, x);
    P0 = zeros(14,1);
    P0(1) = 1;
    end_t = 60*10;
    Pss = expm(end_t*Mf)*P0;
end