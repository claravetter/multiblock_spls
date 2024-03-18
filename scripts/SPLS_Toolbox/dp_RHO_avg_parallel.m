%% function for RHO_avg BO

function [RHO_avg] = dp_RHO_avg_parallel(spls_standalone_path, analysis_folder, type_analysis, mem_total, max_sim_jobs, queue_name, K, tp, cu, cv) 

RHO_avg_bash = dp_bash_parallel(spls_standalone_path, analysis_folder, type_analysis, mem_total, max_sim_jobs, queue_name, K, tp, cu, cv);
    
cd(analysis_folder);
system(['qsub ' RHO_avg_bash]);

filecount = size(dir([analysis_folder '/RHO_*']),1);
RHO_collection = nan(K,1);

% see hyperparameter optimization
while filecount < K
    filecount = size(dir([analysis_folder '/RHO_*.txt']),1);    
end

for i=1:size(RHO_collection,1)
    if exist([analysis_folder '/RHO_',num2str(i),'.txt'],'file')    
        RHO_collection(i,1) = dp_txtscan([permutation_folder '/RHO_' num2str(i) '.txt'], '%f');                    
        delete([analysis_folder '/RHO_',num2str(i),'.txt']);       
    else                
        RHO_collection(i,1) = NaN;       
    end    
end

RHO_collection_completed = dp_fill_nan(RHO_collection, spls_standalone_path, analysis_folder, type_analysis, mem_total, max_sim_jobs, queue_name, tp, cu, cv);
RHO_avg = mean(RHO_collection_completed);

end