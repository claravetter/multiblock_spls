%% DP function for standardizing
% function to standardize X according to method specified in IN and then
% applying the betas of X onto Y
function [TRAIN_s, TEST_s] = dp_standardize_comb(TRAIN, TEST, IN)

[TRAIN_s,IN] = dp_standardize(TRAIN, IN);

[TEST_s,~] = dp_standardize(TEST, IN);

end