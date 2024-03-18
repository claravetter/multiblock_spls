%% function for permutation testing

function dp_ICV_permutation(i, analysis_folder)

m_setup = matfile([analysis_folder '/permutation_setup.mat']);
% load([analysis_folder '/permutation_setup.mat'],'selection_train', 'correlation_method', 'size_sets_permutation', 'perm_coding', 'scaling_method');

if m_setup.selection_train == 1
    m_data = matfile([analysis_folder, '/permutation_partition_fold.mat']);
    m_opt = matfile([analysis_folder '/permutation_opt.mat']);
    %test={'IN_x','IN_y','OCV_train_Diag','COV','cu_opt','cv_opt','V_opt'};
elseif m_setup.selection_train == 2
    w = m_setup.perm_coding(2, (find(m_setup.perm_coding(1,:)==str2double(i))));
    m_data = matfile([analysis_folder, '/permutation_partition_fold_', num2str(w), '.mat']);
    m_opt = matfile([analysis_folder '/permutation_opt.mat']); %,'cu_opt','cv_opt','V_opt');
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

for pp=1:str2double(i)
    permmat = nk_PermInd2(m_setup.size_sets_permutation, m_data.train_Diag);
end

% RHO_collection = nan(m_setup.size_sets_permutation,1);
% u_collection = cell(m_setup.size_sets_permutation,size(IN_x.train,2));
% v_collection = cell(m_setup.size_sets_permutation,size(IN_y.train,2));

% intervall_fill = m_setup.perm_coding(3,str2double(i)):(m_setup.perm_coding(3,str2double(i))+m_setup.size_sets_permutation-1);

RHO_collection_ICV = nan(m_setup.size_sets_permutation,1);
u_collection_ICV = nan(m_setup.size_sets_permutation,size(m_data.train_data_x,2));
v_collection_ICV = nan(m_setup.size_sets_permutation,size(m_data.train_data_y,2));

for ii=1:m_setup.size_sets_permutation
    % perform procrustean transformation to minimize rotation effects of
    % permutated y matrix, if V_opt available
    IN_y.train = IN_y.train(permmat(ii,:),:);
    
    OUT_x = dp_correctscale(IN_x,COV,m_setup.scaling_method);
    OUT_y = dp_correctscale(IN_y,COV,m_setup.scaling_method);
    
    if ~islogical(m_opt.V_opt)
        [RHO_collection_ICV(ii,1), u_collection_ICV(ii,:), v_collection_ICV(ii,:), ~, ~, ~] = dp_spls_full(OUT_x.train,OUT_y.train,OUT_x.test, OUT_y.test, m_opt.cu_opt, m_opt.cv_opt, m_setup.correlation_method, m_opt.V_opt);
    else
        [RHO_collection_ICV(ii,1), u_collection_ICV(ii,:), v_collection_ICV(ii,:), ~, ~, ~] = dp_spls_full(OUT_x.train,OUT_y.train,OUT_x.test, OUT_y.test, m_opt.cu_opt, m_opt.cv_opt, m_setup.correlation_method);
    end
    
end

save([analysis_folder, '/RHO_results_', i, '.mat'], 'RHO_collection_ICV', 'u_collection_ICV', 'v_collection_ICV');

% errorcount=1;
% while errorcount>0
%     try m_coll = matfile([analysis_folder, '/RHO_results.mat'],'Writable',true);
%         m_coll.RHO_collection(intervall_fill,:) = RHO_collection_ICV;
%         m_coll.u_collection(intervall_fill,1) = u_collection_ICV;
%         m_coll.v_collection(intervall_fill,1) = v_collection_ICV;
%         errorcount=0;
%     catch ME
%         errorcount=1;
%         pause(1)
%     end
% end

end