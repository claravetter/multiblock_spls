%% new function for permutation testing

function dp_ICV_hyperopt_csv(i, analysis_folder)

m = matfile([analysis_folder '/hyperopt_partition.mat']);

extract_target = (str2double(i)-1)*m.size_sets_hyperopt+1;
try
    cu_cv_extract = m.cu_cv_combination(extract_target:(extract_target+m.size_sets_hyperopt-1),:);
catch
    cu_cv_extract = m.cu_cv_combination(extract_target:end,:);
end

RHO_collection_ICV = nan(size(cu_cv_extract,1),size(m.cv_inner_TrainInd,1)*size(m.cv_inner_TrainInd,2));
train_data_x = m.train_data_x;
train_data_y = m.train_data_y;
train_covariates = m.train_covariates;
train_DiagNames = m.train_DiagNames;
cs_method_hyperopt = m.cs_method;

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
            
            DiagNames.train = train_DiagNames(cell2mat(m.cv_inner_TrainInd(ib,k)),:);
            DiagNames.test = train_DiagNames(cell2mat(m.cv_inner_TestInd(ib,k)),:);
            
            if ~isempty(cs_method_hyperopt.correction_subgroup)
                cs_method_hyperopt.subgroup_train = contains(DiagNames.train, cs_method_hyperopt.correction_subgroup);
                cs_method_hyperopt.subgroup_test = contains(DiagNames.test, cs_method_hyperopt.correction_subgroup);
            else
                cs_method_hyperopt.subgroup_train = [];
                cs_method_hyperopt.subgroup_test = [];
            end
            
            [OUT_x, OUT_y] = dp_master_correctscale(IN_x, IN_y, COV, cs_method_hyperopt, m.correction_target);

            RHO_collection_ICV(ii,nn) = dp_spls_slim(OUT_x.train,OUT_y.train,OUT_x.test,OUT_y.test, cu, cv, m.correlation_method);
            nn=nn+1;
        end
    end
end

writematrix(RHO_collection_ICV,[analysis_folder, '/RHO_results_', i, '.csv'],'Delimiter','tab')

end