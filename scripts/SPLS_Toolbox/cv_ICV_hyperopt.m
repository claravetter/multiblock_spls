%% new function for permutation testing

function cv_ICV_hyperopt(i, analysis_folder)

m = matfile([analysis_folder '/hyperopt_partition.mat']);

extract_target = (i-1)*m.size_sets_hyperopt+1;
try
    weights_extract = m.weights_combination(extract_target:(extract_target+m.size_sets_hyperopt-1),:);
catch
    weights_extract = m.weights_combination(extract_target:end,:);
end

% RHO_collection_ICV = nan(size(weights_extract,1),size(m.cv_inner_TrainInd,1)*size(m.cv_inner_TrainInd,2));
weights_collection_ICV = cell(size(weights_extract,1),1);
%v_collection_ICV = cell(size(weights_extract,1),1);
train_data_matrices = m.train_data_matrices;
%train_data_y = m.train_data_y;
train_covariates = m.train_covariates;

for ii=1:size(weights_extract,1)    
    
    weights = weights_extract(ii,:);
    %cv = weights_extract(ii,2);
    nn=1;
    for ib=1:size(m.cv_inner_TrainInd,1)
        for k=1:size(m.cv_inner_TrainInd,2)
            
            for mat_ind = 1:length(train_data_matrices)
                IN.train{mat_ind} = train_data_matrices{mat_ind}(cell2mat(m.cv_inner_TrainInd(ib,k)),:);
                IN.test{mat_ind} = train_data_matrices{mat_ind}(cell2mat(m.cv_inner_TestInd(ib,k)),:);
            end
            
            %IN_y.train = train_data_y(cell2mat(m.cv_inner_TrainInd(ib,k)),:);
            %IN_y.test = train_data_y(cell2mat(m.cv_inner_TestInd(ib,k)),:);

            COV.train = train_covariates(cell2mat(m.cv_inner_TrainInd(ib,k)),:);
            COV.test = train_covariates(cell2mat(m.cv_inner_TestInd(ib,k)),:);
            
            if isfield(m, 'scaling_method')
                OUT = dp_correctscale(IN,COV,m.scaling_method);
            %OUT_y = dp_correctscale(IN_y,COV,m.scaling_method);
            else
                OUT = IN;
            end
            
            [RHO_collection_ICV(ii,nn), weights_collection_ICV{ii,:}(:,nn), ~, ~] = cv_gspls_full(OUT.train, OUT.test, weights, m.correlation_method);
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

save([analysis_folder, '/RHO_results_', num2str(i), '.mat'], 'RHO_collection_ICV', 'weights_collection_ICV');

end