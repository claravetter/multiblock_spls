%% DP function to compute AUC and test whether one value belongs to a distribution

function [p_val, h] = dp_auc_testing(IN)

dist = IN.dist;
val = IN.val;
if isfield(IN, 'bin_width')
    bin_width = IN.bin_width;
else
    bin_width = 0.001;
end

if isfield(IN, 'testing_precision')
    testing_precision = IN.testing_precision;
else
    testing_precision = 'lenient';
end

[N, edges] = histcounts(dist, 'BinWidth', bin_width);
x = 0:bin_width:((size(edges,2)-2)*bin_width);
y = N;

[~, index] = min(abs(x-val));
switch testing_precision
    case 'lenient'
        p_val = sum(cumtrapz(x(index:end), y(index:end)))/sum(cumtrapz(x,y));
    case 'conservative'
        p_val = trapz(x(index:end), y(index:end))/trapz(x,y);
end

if p_val<0.05
    h=true;
else
    h=false;
end

end