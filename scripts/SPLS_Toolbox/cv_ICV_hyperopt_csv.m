%% new function for permutation testing

function cv_ICV_hyperopt_csv(i,analysis_folder)

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

train_DiagNames = m.train_DiagNames;
cs_method_hyperopt = m.cs_method;

for ii=1:size(weights_extract,1)
    weights = weights_extract(ii,:);

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

            DiagNames.train = train_DiagNames(cell2mat(m.cv_inner_TrainInd(ib,k)),:);
            DiagNames.test = train_DiagNames(cell2mat(m.cv_inner_TestInd(ib,k)),:);

            if ~isempty(cs_method_hyperopt.correction_subgroup)
                cs_method_hyperopt.subgroup_train = contains(DiagNames.train, cs_method_hyperopt.correction_subgroup);
                cs_method_hyperopt.subgroup_test = contains(DiagNames.test, cs_method_hyperopt.correction_subgroup);
            else
                cs_method_hyperopt.subgroup_train = [];
                cs_method_hyperopt.subgroup_test = [];
            end

            OUT = cv_master_correctscale(IN, COV, cs_method_hyperopt, m.correction_target);

           
            RHO_collection_ICV(ii,nn) = cv_gspls_slim(OUT.train,OUT.test, weights, m.correlation_method);
        
            nn=nn+1;
        end
    end
end

writematrix(RHO_collection_ICV,[analysis_folder, '/RHO_results_', num2str(i), '.csv'],'Delimiter','tab')

end