%% DP function for training/retraining and merging training results

function values_merged = dp_trainmerge_single(values, type_merge, dim, weights)

switch type_merge
    case 'mean'
        values_merged = mean(values, dim);
    case 'median'
        values_merged = median(values, dim);
    case 'weighted_mean'
        values_merged = wmean(values, weights, dim);
    case 'best'
        values_merged = max(values, dim);
end

end