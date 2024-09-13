%% function for permutation testing

function cv_mbspls_ICV_permutation_csv(i, analysis_folder)

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

IN_matrices.train = m_data.train_data_matrices;
IN_matrices.test = m_data.test_data_matrices;

train_covariates = m_data.train_covariates; 
test_covariates = m_data.test_covariates; 
for num_m=1:size(IN_matrices.train,2)
    COV{num_m}.train = train_covariates{num_m};
    COV{num_m}.test = test_covariates{num_m};
end

cs_method_permutation = m_data.cs_method;

% permut all but last matrix
% for num_m=1:length(IN_matrices.train)-1
%     for pp=1:str2double(i)
%         permmat = nk_PermInd(m_setup.size_sets_permutation, m_data.train_Diag); % CV: nk_PermInd if NeuroMiner_1.3
%     end
%     permmats{num_m} = permmat; 
% end

% permute all matrices
for num_m=1:length(IN_matrices.train)
    for pp=1:str2double(i)
        permmat = nk_PermInd(m_setup.size_sets_permutation, m_data.train_Diag); 
    end
    permmats{num_m} = permmat; 
end
clear permmat

RHO_collection_ICV = nan(m_setup.size_sets_permutation,1);

for ii=1:m_setup.size_sets_permutation
    % perform procrustean transformation to minimize rotation effects of
    % permutated y matrix, if V_opt available
    for num_m=1:length(IN_matrices.train)-1
        IN_matrices.train{num_m+1} = IN_matrices.train{num_m+1}(permmats{num_m}(ii,:),:);%CV: does it matter which matrix is being permuted? Only one? 
    % CV 10.4.2024: alle anderen Matrizen außer die erste müssen permutiert werden --> permmat{} pro Matrix  
    end

    [OUT_matrices] = cv_master_correctscale(IN_matrices, COV, cs_method_permutation, m_setup.correction_target);
    
    if ~islogical(m_opt.Vs_opt) && m_setup.procrustes_flag
        RHO_collection_ICV(ii,1) = cv_mbspls_wrapper(OUT_matrices.train, ...
            OUT_matrices.test, ...
            m_opt.c_weights_opt, ...
            m_setup.mbspls_params, ...
            m_setup.correlation_method, ...
            m_setup.matrix_norm, ...
            m_opt.Vs_opt); % slim
    else
        RHO_collection_ICV(ii,1) = cv_mbspls_wrapper(OUT_matrices.train, ...
            OUT_matrices.test, ...
            m_opt.c_weights_opt, ...
            m_setup.mbspls_params, ...
            m_setup.correlation_method, ...
            m_setup.matrix_norm, ...
            []); % slim
    end
end

writematrix(RHO_collection_ICV,[analysis_folder, '/RHO_results_', num2str(i), '.csv'],'Delimiter','tab')

end