function [TRAIN_imp, TEST_imp] = cv_impute_comb(TRAIN, TEST)
%IN_imp.method = 'euclidean';
%IN_imp.k = 5; 
[TRAIN_imp, IN_imp] = nk_PerfImputeObj(TRAIN);
IN_imp.X = TRAIN_imp;
[TEST_imp, ~] = nk_PerfImputeObj(TEST, IN_imp);
end