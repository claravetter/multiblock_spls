%% DP function to compute NaN values for RHO_b

function [RHO_b_collection_completed] = dp_perm_nan(RHO_b_collection, spls_standalone_path, permutation_folder, utilities_folder, mem_total, max_sim_jobs, queue_name)
   
nan_log = isnan(RHO_b_collection);
    
if sum(nan_log)>0
    if ~exist([permutation_folder, '/addition'],'dir')
        mkdir([permutation_folder, '/addition']);   
        permutation_addition_folder=[permutation_folder, '/addition']; 
    else
        permutation_addition_folder=[permutation_folder, '/addition'];    
    end

    copyfile([permutation_folder '/opt_param.mat'], [permutation_addition_folder '/opt_param.mat']);
    load([permutation_addition_folder '/opt_param.mat']);

    while sum(nan_log)>0

        for i=1:sum(nan_log)
            temp = randperm(size(keep_in_data_y,1))';
            FID = fopen([permutation_addition_folder '/perm_' num2str(i) '.txt'],'w');
            fprintf(FID,'%d\n',temp);
            fclose(FID);
        end

        % write a bash script for hyperparameter optimization, which can later be called
        permutation_bash = dp_bash_para_new(spls_standalone_path, permutation_addition_folder, utilities_folder, 'permutation', mem_total, sum(nan_log), max_sim_jobs, queue_name);

        % submit the bash script for parallel computing of RHO_b values to the queue
        cd(permutation_addition_folder);
        system(['qsub ' permutation_bash]);

        % collect all files with RHO_b values, open them and store the RHO values in a vector
        filecount_add = size(dir([permutation_addition_folder '/RHO_b_*.txt']),1);
        RHO_b_collection_addition = nan(sum(nan_log),1);

        while filecount_add < sum(nan_log)
            filecount_add = size(dir([permutation_addition_folder '/RHO_b_*.txt']),1);    
        end 

        for i=1:size(RHO_b_collection_addition,1)
            if exist([permutation_addition_folder '/RHO_b_',num2str(i),'.txt'],'file')
                RHO_b_collection_addition(i,1) = dp_txtscan([permutation_addition_folder '/RHO_b_' num2str(i) '.txt'], '%f');            
                delete([permutation_addition_folder '/RHO_b_',num2str(i),'.txt']);
            else
                RHO_b_collection_addition(i,1) = NaN;
            end
        end

        RHO_b_collection(nan_log) = RHO_b_collection_addition;
        nan_log = isnan(RHO_b_collection);

    end

    RHO_b_collection_completed = RHO_b_collection;

else
    
    RHO_b_collection_completed = RHO_b_collection;
    
end

end
 