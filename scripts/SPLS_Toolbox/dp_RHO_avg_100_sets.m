%% new function for permutation testing

function dp_RHO_avg_100_sets(i, size_sets_str, analysis_folder)

load([analysis_folder '/keep_in_partition.mat']);
dp_txt_write(analysis_folder, ['init_' i],'initialized','%s \n');

size_sets = str2double(size_sets_str);
RHO_avg_collection = nan(size_sets,1);

for ii=1:size_sets
    if exist([analysis_folder '/cu_' i '_' num2str(ii) '.txt'],'file')
        RHO_collection = nan(K,1);
        cu = dp_txtscan([analysis_folder '/cu_', i, '_', num2str(ii), '.txt'], '%f');
        cv = dp_txtscan([analysis_folder '/cv_', i, '_', num2str(ii), '.txt'], '%f');
        
        if cu > sqrt(size(keep_in_data_x,2))
            cu = sqrt(size(keep_in_data_x,2));
        end
        
        if cv > sqrt(size(keep_in_data_y,2))
            cv = sqrt(size(keep_in_data_y,2));
        end
        
        % set up inner fold CV partitions, balanced for sites
        
        repeat=true; K_I = 10;
        while repeat
            try cv_inner_test = nk_CVpartition2_dp(1, K_I, keep_in_labels)
                repeat=false;
            catch
                repeat = true;
                K_I=K_I-1;
            end
        end
        
        nn=1;
        for iii=1:((round(K/K_I))+1)
            
            cv_inner = nk_CVpartition2_dp(1, K_I, keep_in_labels);
            
            for k=1:K_I
                test_data_x = keep_in_data_x(cv_inner.TestInd{k},:);
                test_data_y = keep_in_data_y(cv_inner.TestInd{k},:);
                training_data_x = keep_in_data_x(cv_inner.TrainInd{k},:);
                training_data_y = keep_in_data_y(cv_inner.TrainInd{k},:);
                RHO_collection(nn) = dp_k_split(training_data_x,training_data_y,test_data_x, test_data_y, cu, cv);
                nn=nn+1;
            end
            
        end
        
        RHO_avg_collection(ii,1) = mean(RHO_collection);
        %     end
        FID = fopen([analysis_folder, '/init_' i '.txt'], 'a');
        fprintf(FID, '%d \n', ii);
        fclose(FID);
    end
    
end

RHO_avg_collection(isnan(RHO_avg_collection))=[];

dp_txt_write(analysis_folder, ['RHO_' i], RHO_avg_collection, '%.4f\n');

end