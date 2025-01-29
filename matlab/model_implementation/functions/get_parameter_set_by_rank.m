function [fit, score] = get_parameter_set_by_rank(fits, scores, rank)
    [~,indices] = sort(scores);
    idx = indices(rank);
    fit = fits(idx,:);
    score = scores(idx);
end