function Pss = predict_steady_state_probabilities(model, x)
    Mf = create_forward_matrix_no_cleavage(model, x);
    P0 = zeros(8,1);
    P0(1) = 1;
    end_t = 60*10;
    Pss = expm(end_t*Mf)*P0;
end