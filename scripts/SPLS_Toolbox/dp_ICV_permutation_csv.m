%% function for permutation testing

function dp_ICV_permutation_csv(i, analysis_folder)

m_setup = matfile([analysis_folder '/permutation_setup.mat']);

if m_setup.selection_train == 1
    m_data = matfile([analysis_folder, '/permutation_partition_fold.mat']);
    m_opt = matfile([analysis_folder '/permutation_opt.mat']);
elseif m_setup.selection_train == 2
    w = m_setup.perm_coding(2, (find(m_setup.perm_coding(1,:)==str2double(i))));
    m_data = matfile([analysis_folder, '/permutation_partition_fold_', num2str(w), '.mat']);
    m_opt = matfile([analysis_folder '/permutation_opt.mat']); 
end

% 1) retrain on single test splits within
% folds, then merge, 2) retrain on all inner folds
% separately, then merge with mean or median, 3)
% retrain on entirety of inner folds, 4) use already
% existing u and v from inner folds without retraining

IN_x.train = m_data.train_data_x;
IN_x.test = m_data.test_data_x;
COV.train = m_data.train_covariates;
COV.test = m_data.test_covariates;
IN_y.train = m_data.train_data_y;
IN_y.test = m_data.test_data_y;

cs_method_permutation = m_data.cs_method;

for pp=1:str2double(i)
    permmat = nk_PermInd2(m_setup.size_sets_permutation, m_data.train_Diag);
end

RHO_collection_ICV = nan(m_setup.size_sets_permutation,1);

for ii=1:m_setup.size_sets_permutation
    % perform procrustean transformation to minimize rotation effects of
    % permutated y matrix, if V_opt available
    IN_y.train = IN_y.train(permmat(ii,:),:);
    
    [OUT_x, OUT_y] = dp_master_correctscale(IN_x, IN_y, COV, cs_method_permutation, m_setup.correction_target);
    
    if ~islogical(m_opt.V_opt)
        RHO_collection_ICV(ii,1) = dp_spls_slim(OUT_x.train,OUT_y.train,OUT_x.test, OUT_y.test, m_opt.cu_opt, m_opt.cv_opt, m_setup.correlation_method, m_opt.V_opt);
    else
        RHO_collection_ICV(ii,1) = dp_spls_slim(OUT_x.train,OUT_y.train,OUT_x.test, OUT_y.test, m_opt.cu_opt, m_opt.cv_opt, m_setup.correlation_method);
    end
end

writematrix(RHO_collection_ICV,[analysis_folder, '/RHO_results_', i, '.csv'],'Delimiter','tab')

end