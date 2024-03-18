%% function for permutation testing

function dp_ICV_bootstrap_csv(i, analysis_folder)

m_setup = matfile([analysis_folder '/bootstrap_setup.mat']);

if m_setup.selection_train == 1
    m_data = matfile([analysis_folder, '/bootstrap_partition_fold.mat']);
    m_opt = matfile([analysis_folder '/bootstrap_opt.mat']);
end

% 1) retrain on single test splits within
% folds, then merge, 2) retrain on all inner folds
% separately, then merge with mean or median, 3)
% retrain on entirety of inner folds, 4) use already
% existing u and v from inner folds without retraining

X = m_data.train_data_x;
Y = m_data.train_data_y;
covars = m_data.train_covariates;
labels = m_data.train_DiagNames;

cs_method_bootstrap = m_data.cs_method;
correction_target = m_setup.correction_target;

[~,bootsam] = bootstrp(m_setup.size_sets_bootstrap,[],1:size(Y,1));

% perform procrustean transformation to minimize rotation effects of
% permutated y matrix, if V_opt available
RHO_boot=[]; u_boot=[]; v_boot=[];
for ii=1:size(bootsam,2)
    X_boot = X(bootsam(:,ii),:);
    Y_boot = Y(bootsam(:,ii),:);
    covars_boot = covars(bootsam(:,ii),:);
    labels_boot = labels(bootsam(:,ii),:);
    
    cs_method_bootstrap.subgroup_train = matches(labels_boot, cs_method_bootstrap.correction_subgroup);
    cs_method_bootstrap.subgroup_test = matches(labels_boot, cs_method_bootstrap.correction_subgroup);
    [OUT_x, OUT_y] = dp_master_correctscale(X_boot, Y_boot, covars_boot, cs_method_bootstrap, correction_target);
    cu = m_opt.cu_opt;
    cv = m_opt.cv_opt;
    
    if ~islogical(m_opt.V_opt)
        [RHO_boot(1,ii), u_boot(:,ii), v_boot(:,ii), ~, ~, ~] = dp_spls_full(OUT_x,OUT_y,OUT_x, OUT_y, cu, cv, m_setup.correlation_method, m_opt.V_opt);
    else
        [RHO_boot(1,ii), u_boot(:,ii), v_boot(:,ii), ~, ~, ~] = dp_spls_full(OUT_x,OUT_y,OUT_x, OUT_y, cu, cv, m_setup.correlation_method);
    end
end

writematrix(RHO_boot,[analysis_folder, '/RHO_results_', i, '.csv'],'Delimiter','tab')
writematrix(u_boot,[analysis_folder, '/u_results_', i, '.csv'],'Delimiter','tab')
writematrix(v_boot,[analysis_folder, '/v_results_', i, '.csv'],'Delimiter','tab')

end