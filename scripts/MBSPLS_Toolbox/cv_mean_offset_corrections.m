function [TRAIN_c, TEST_c] = cv_mean_offset_corrections(TRAIN, TEST, TRAIN_covariates, TEST_covariates)

IN = struct;
IN.sTrInd = TRAIN_covariates(:,1);

if size(TRAIN_covariates,2)>1
    IN.dTrInd = TRAIN_covariates(:,2);
end

[TRAIN_c, IN] = nk_PerfRemMeanDiffObj(TRAIN, IN);

IN.sTsInd = TEST_covariates(:,1);

if size(TRAIN_covariates,2)>1
    IN.dTsInd = TEST_covariates(:,2);
end

[TEST_c, ~] = nk_PerfRemMeanDiffObj(TEST, IN);

end