%% DP script for all multiple comparison techniques

function [c_pvalues, h, extra] = dp_multiple_testing(test_name, pvalues, alpha, plotting)

switch test_name
    case 'Bonferroni'
        [c_pvalues, ~, h, extra] = fwer_bonf(pvalues, alpha, plotting);
    case 'Sidak'
        [c_pvalues, ~, h, extra] = fwer_sidak(pvalues, alpha, plotting);
    case 'Holm_Bonferroni'
        [c_pvalues, ~, h, extra] = fwer_holmbonf(pvalues, alpha, plotting);
    case 'Benjamini_Hochberg'
        [c_pvalues, ~, h, extra] = fdr_BH(pvalues, alpha, plotting);
    case 'Benjamini_Yekutieli'
        [c_pvalues, ~, h, extra] = fdr_BY(pvalues, alpha, correction, plotting);
    case 'Storey'
        [qvalues, ~, h, extra] = fdr_storey(pvalues, alpha, plotting);
        c_pvalues = qvalues;
        disp('Be aware: in Storey FDR correction, cpvalues are in fact qvalues');
    case 'Fisher'
        [c_pvalues, ~, h, extra] = mt_fisher(pvalues, alpha, plotting);
end

end