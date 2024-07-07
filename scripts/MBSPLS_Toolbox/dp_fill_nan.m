%% DP function to compute NaN values for RHO_b

function [RHO_collection_completed] = dp_fill_nan(RHO_collection, spls_standalone_path, analysis_folder, type_analysis, mem_total, max_sim_jobs, queue_name, tp, cu, cv)
   
nan_log = isnan(RHO_collection);

switch type_analysis
    case 'hyperopt'
        search_prefix = 'RHO_';
    case 'permutation'
        search_prefix = 'RHO_b_';
        load([analysis_folder '/opt_param.mat']);
        for i=1:sum(nan_log)
            temp = randperm(size(keep_in_data_y,1))';
            dp_txt_write(permutation_folder, ['perm_' num2str(i)], temp, '%d\n');
        end
end

while sum(nan_log)>0

    % write a bash script for hyperparameter optimization, which can later be called
    bash_file = dp_bash_parallel(spls_standalone_path, analysis_folder, type_analysis, mem_total, max_sim_jobs, queue_name, sum(nan_log), tp, cu, cv);

    % submit the bash script for parallel computing of RHO_b values to the queue
    cd(analysis_folder);
    system(['qsub ' bash_file]);

    % collect all files with RHO_b values, open them and store the RHO values in a vector
    filecount_add = size(dir([addition_folder '/' search_prefix '*.txt']),1);
    RHO_collection_addition = nan(sum(nan_log),1);

    while filecount_add < sum(nan_log)
        filecount_add = size(dir([addition_folder '/' search_prefix '*.txt']),1);    
    end

    for i=1:size(RHO_collection_addition,1)
        if exist([addition_folder '/' search_prefix, num2str(i),'.txt'],'file')
            RHO_collection_addition(i,1) = dp_txtscan([addition_folder '/' search_prefix, num2str(i),'.txt'], '%f');            
            delete([addition_folder '/' search_prefix, num2str(i),'.txt']);
        else
            RHO_collection_addition(i,1) = NaN;
        end
    end

    RHO_collection(nan_log) = RHO_collection_addition;
    nan_log = isnan(RHO_collection);

end

RHO_collection_completed = RHO_collection;

end
 