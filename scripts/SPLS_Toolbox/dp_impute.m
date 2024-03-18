%% DP function for imputation using nk_PerfImputeObj
function data_imputed = dp_impute(data_tobeimputed, method)

IN.method = method;
IN.X = data_tobeimputed;
IN.Y = data_tobeimputed;
[data_imputed, ~] = nk_PerfImputeObj(data_tobeimputed, IN);

end