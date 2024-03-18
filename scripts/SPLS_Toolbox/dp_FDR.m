%% DP function for FDR correction
function [pvalue_FDR] = dp_FDR(pvalues, FDRvalue)

if ~exist('FDRvalue','var')
    FDRvalue = 0.05;
end

pvalues_sorted = sort(pvalues);
m = numel(pvalues);
k = 1:m;
ratio = size(pvalues,1)/size(pvalues,2);

if ratio >= 1
    FDRthreshold = ((k*FDRvalue)/m)';
else
    FDRthreshold = (k*FDRvalue)/m;
end

decision = pvalues_sorted <= FDRthreshold;
kmax = find(decision, 1, 'last');
pvalue_FDR = pvalues_sorted(kmax);

if isempty(pvalue_FDR)
    pvalue_FDR = 0;
end

end
