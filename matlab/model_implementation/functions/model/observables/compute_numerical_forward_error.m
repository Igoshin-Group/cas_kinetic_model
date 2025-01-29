function err = compute_numerical_forward_error(pset)

conc = 10; % concentration (nM)
Pss = predict_steady_state_probabilities_full_model(pset, conc);

right = Pss(7);
wrong = Pss(14);
err = wrong ./ right;
end