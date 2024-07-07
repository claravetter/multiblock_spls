%% new function for permutation testing

function dp_hyperopt_ICV(i, size_sets_str, analysis_folder)

load([analysis_folder '/keep_in_partition.mat']);

size_sets = str2double(size_sets_str);

for ii=1:size_sets
    if exist([analysis_folder '/cu_' i '_' num2str(ii) '.txt'],'file')
        RHO_collection = nan(size(cv_inner.TestInd,2),ii);
        cu = dp_txtscan([analysis_folder '/cu_', i, '_', num2str(ii), '.txt'], '%f');
        cv = dp_txtscan([analysis_folder '/cv_', i, '_', num2str(ii), '.txt'], '%f');
        
        if cu > sqrt(size(keep_in_data_x,2))
            cu = sqrt(size(keep_in_data_x,2));
        end
        
        if cv > sqrt(size(keep_in_data_y,2))
            cv = sqrt(size(keep_in_data_y,2));
        end
        
        for k=1:size(cv_inner.TestInd,2)
            test_data_x = keep_in_data_x(cv_inner.TestInd{k},:);
            test_data_y = keep_in_data_y(cv_inner.TestInd{k},:);
            training_data_x = keep_in_data_x(cv_inner.TrainInd{k},:);
            training_data_y = keep_in_data_y(cv_inner.TrainInd{k},:);
            RHO_collection(k,ii) = dp_k_split(training_data_x,training_data_y,test_data_x, test_data_y, cu, cv, correlation_method);
        end
    end
end

try
    save([analysis_folder, '/RHO_', i, '.mat'], 'RHO_collection');
end

end