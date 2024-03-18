%% function for RHO parallel computing

function [RHO_output] = dp_RHO_parallel(spls_standalone_path, analysis_folder, type_analysis, mem_total, max_sim_jobs, queue_name, total_jobs, tp_size_sets, cu, cv) 

switch type_analysis
    case 'hyperopt'
        search_prefix = 'RHO_';
        RHO_bash = dp_bash_parallel(spls_standalone_path, analysis_folder, type_analysis, mem_total, max_sim_jobs, queue_name, total_jobs, tp_size_sets, cu, cv);
        RHO_collection = zeros(total_jobs,1);
    case 'permutation'
        search_prefix = 'RHO_b_';
        RHO_bash = dp_bash_parallel(spls_standalone_path, analysis_folder, type_analysis, mem_total, max_sim_jobs, queue_name, total_jobs, tp_size_sets); 
        RHO_collection = zeros(total_jobs*tp_size_sets,1);
end
   
cd(analysis_folder);
system(['qsub ' RHO_bash]);

% see hyperparameter optimization

while sum(sum(RHO_collection == 0))>0;
    for i=1:size(RHO_collection,1)
        if exist([analysis_folder '/' search_prefix, num2str(i),'.txt'],'file')
            try
                RHO_collection(i,1) = dp_txtscan([analysis_folder '/' search_prefix, num2str(i),'.txt'], '%f');                    
                delete([analysis_folder '/' search_prefix, num2str(i),'.txt']);  
            catch ME
                disp(ME.identifier);
                RHO_collection(i,1) = NaN;
            end
        end  
    end
end


switch type_analysis
    case 'hyperopt'
%         RHO_collection_completed = dp_fill_nan(RHO_collection, spls_standalone_path, analysis_folder, type_analysis, mem_total, max_sim_jobs, queue_name, tp_size_sets, cu, cv);
        RHO_output = mean(RHO_collection);
    case 'permutation'
%         RHO_collection_completed = dp_fill_nan(RHO_collection, spls_standalone_path, analysis_folder, type_analysis, mem_total, max_sim_jobs, queue_name);
        RHO_output = RHO_collection;
end

end