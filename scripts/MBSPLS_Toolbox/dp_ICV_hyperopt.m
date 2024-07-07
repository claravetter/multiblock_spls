%% new function for permutation testing

function dp_ICV_hyperopt(i, analysis_folder)

m = matfile([analysis_folder '/hyperopt_partition.mat']);

extract_target = (str2double(i)-1)*m.size_sets_hyperopt+1;
try
    cu_cv_extract = m.cu_cv_combination(extract_target:(extract_target+m.size_sets_hyperopt-1),:);
catch
    cu_cv_extract = m.cu_cv_combination(extract_target:end,:);
end

% RHO_collection_ICV = nan(size(cu_cv_extract,1),size(m.cv_inner_TrainInd,1)*size(m.cv_inner_TrainInd,2));
u_collection_ICV = cell(size(cu_cv_extract,1),1);
v_collection_ICV = cell(size(cu_cv_extract,1),1);
train_data_x = m.train_data_x;
train_data_y = m.train_data_y;
train_covariates = m.train_covariates;

for ii=1:size(cu_cv_extract,1)    
    cu = cu_cv_extract(ii,1);
    cv = cu_cv_extract(ii,2);
    nn=1;
    for ib=1:size(m.cv_inner_TrainInd,1)
        for k=1:size(m.cv_inner_TrainInd,2)
            
            IN_x.train = train_data_x(cell2mat(m.cv_inner_TrainInd(ib,k)),:);
            IN_x.test = train_data_x(cell2mat(m.cv_inner_TestInd(ib,k)),:);
            
            IN_y.train = train_data_y(cell2mat(m.cv_inner_TrainInd(ib,k)),:);
            IN_y.test = train_data_y(cell2mat(m.cv_inner_TestInd(ib,k)),:);

            COV.train = train_covariates(cell2mat(m.cv_inner_TrainInd(ib,k)),:);
            COV.test = train_covariates(cell2mat(m.cv_inner_TestInd(ib,k)),:);
            
            OUT_x = dp_correctscale(IN_x,COV,m.scaling_method);
            OUT_y = dp_correctscale(IN_y,COV,m.scaling_method);
            
            [RHO_collection_ICV(ii,nn), u_collection_ICV{ii,1}(:,nn), v_collection_ICV{ii,1}(:,nn), ~, ~, ~] = dp_spls_full(OUT_x.train,OUT_y.train,OUT_x.test,OUT_y.test, cu, cv, m.correlation_method);
            nn=nn+1;
        end
    end
end

% errorcount=1;
% while errorcount>0
%     try  m_coll = matfile([analysis_folder, '/RHO_results.mat'],'Writable',true);
%         m_coll.RHO_collection(extract_target:(extract_target+m.size_sets_hyperopt-1),:) = RHO_collection_ICV;
%         m_coll.u_collection(extract_target:(extract_target+m.size_sets_hyperopt-1),1) = u_collection_ICV;
%         m_coll.v_collection(extract_target:(extract_target+m.size_sets_hyperopt-1),1) = v_collection_ICV;
%         errorcount=0;
%     catch ME
%         errorcount=1;
%         pause(1)
%     end
% end

save([analysis_folder, '/RHO_results_', i, '.mat'], 'RHO_collection_ICV', 'u_collection_ICV', 'v_collection_ICV');

end