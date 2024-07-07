%% function for correction keep in and hold out data for covariates
% regress out covariate effects from TRAIN and then apply this correction to
% TEST, with betas from TRAIN and covariates from TEST

function [TRAIN_c, TEST_c] = dp_corrections_multi(TRAIN, TEST, TRAIN_covariates, TEST_covariates)

switch corr_option
    case 1 % use partial correlations to correct data
        IN = struct;
        IN.TrCovars = TRAIN_covariates;
        [TRAIN_c, IN] = nk_PartialCorrelationsObj(TRAIN, IN);
        IN.TsCovars = TEST_covariates;
        [TEST_c, ~] = nk_PartialCorrelationsObj(TEST, IN);
        
    case 2 % use PCA correction
        IN.S = TRAIN;
        IN.G = TRAIN_covariates;
        [TRAIN_c, IN, ~] = nk_PerfAdjForCovarsUsingPCAObj(TRAIN, IN);
        
        IN.S = TRAIN;
        IN.G = TEST_covariates;
        [TEST_c, ~, ~] = nk_PerfAdjForCovarsUsingPCAObj(TEST, IN);
end

end