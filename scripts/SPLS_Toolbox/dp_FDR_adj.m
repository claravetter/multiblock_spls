%% DP function for FDR correction
function [pvalues_adj, pvalue_FDR] = dp_FDR_adj(pvalues, FDRvalue)

pvalues_adj = pvalues;

temp = pvalues(:);

try temp = temp(:);
catch
end

temp_save = temp;
temp(isnan(temp))=[];

if ~exist('FDRvalue','var')
    FDRvalue = 0.05;
end

pvalues_sorted = sort(temp);
m = numel(temp);
k = 1:m;
ratio = size(temp,1)/size(temp,2);

if ratio >= 1
    FDRthreshold = ((k*FDRvalue)/m)';
else
    FDRthreshold = (k*FDRvalue)/m;
end

decision = pvalues_sorted < FDRthreshold;
kmax = find(decision, 1, 'last');
pvalue_FDR = FDRthreshold(kmax);

if isempty(pvalue_FDR)
    pvalue_FDR = FDRvalue/m;
end

pvalues_adj_temp = temp.*(FDRvalue/pvalue_FDR);
temp_save(~isnan(temp_save))=pvalues_adj_temp;

pvalues_adj(:) = temp_save;

pvalues_adj(pvalues_adj>1)=rescale(pvalues_adj(pvalues_adj>1), 0.9, 1);

end
