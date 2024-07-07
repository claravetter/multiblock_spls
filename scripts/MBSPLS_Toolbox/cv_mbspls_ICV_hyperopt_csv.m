%% new function for permutation testing

function cv_mbspls_ICV_hyperopt_csv(i,analysis_folder)

m = matfile(fullfile(analysis_folder, 'hyperopt_partition.mat'));

if isnumeric(i)
    extract_target = (i-1)*m.size_sets_hyperopt+1;
else
    extract_target = (str2double(i)-1)*m.size_sets_hyperopt+1;
end

disp(['i = ', i, "; extract_target = ", extract_target]);
try
    %disp(extract_target:(extract_target+m.size_sets_hyperopt-1));
    weights_extract = m.weights_combination(extract_target:(extract_target+m.size_sets_hyperopt-1),:);
catch
    %disp(extract_target);
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
                COV{mat_ind}.train = train_covariates{mat_ind}(cell2mat(m.cv_inner_TrainInd(ib,k)),:);
                COV{mat_ind}.test = train_covariates{mat_ind}(cell2mat(m.cv_inner_TestInd(ib,k)),:);

            end

            %IN_y.train = train_data_y(cell2mat(m.cv_inner_TrainInd(ib,k)),:);
            %IN_y.test = train_data_y(cell2mat(m.cv_inner_TestInd(ib,k)),:);

           
            DiagNames.train = train_DiagNames(cell2mat(m.cv_inner_TrainInd(ib,k)),:);
            DiagNames.test = train_DiagNames(cell2mat(m.cv_inner_TestInd(ib,k)),:);

            for num_m=1:size(COV,2)
                if ~isempty(cs_method_hyperopt{num_m}.correction_subgroup)
                    cs_method_hyperopt{num_m}.subgroup_train = contains(DiagNames.train, cs_method_hyperopt{num_m}.correction_subgroup);
                    cs_method_hyperopt{num_m}.subgroup_test = contains(DiagNames.test, cs_method_hyperopt{num_m}.correction_subgroup);
                else
                    cs_method_hyperopt{num_m}.subgroup_train = [];
                    cs_method_hyperopt{num_m}.subgroup_test = [];
                end
            end

            OUT = cv_master_correctscale(IN, COV, cs_method_hyperopt, m.correction_target);

           
            RHO_collection_ICV(ii,nn) = cv_mbspls_slim(OUT.train,OUT.test, weights, m.mbspls_params, m.correlation_method, [], m.matrix_norm);
        
            nn=nn+1;
        end
    end
end

writematrix(RHO_collection_ICV,[analysis_folder, '/RHO_results_', num2str(i), '.csv'],'Delimiter','tab')

end