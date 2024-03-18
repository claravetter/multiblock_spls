%% DP function to quickly do fitlm and extract relevant features

function [result] = dp_fitlm(X, Y)

for v=1:size(Y,2)
    mdl = fitlm(X, Y(:,v));
    p(v,1) = mdl.coefTest;
    rsquared(v,1) = mdl.Rsquared.Adjusted;
    output.post_hoc_mdl.test.rsquared = [output.post_hoc_mdl.test.rsquared, rsquared];
    output.post_hoc_mdl.test.p = [output.post_hoc_mdl.test.p, p];
    output.post_hoc_mdl.test.mdl = {output.post_hoc_mdl.test.mdl, mdl};
end

[output.post_hoc_mdl.test.p, ~] = dp_FDR_adj(output.post_hoc_mdl.test.p);


end