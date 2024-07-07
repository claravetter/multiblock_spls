%% function for one hyperoptimization step

function [RHO_avg] = dp_hyperopt_single(folders, settings, cu, cv)
K=100; tp=10;
no_jobs = K;

hyperopt_bash = dp_bash_para_new(folders.spls_standalone_path, folders.hyperopt_folder, folders.utilities_folder, 'hyperopt', settings.mem_total, settings.no_jobs, settings.max_sim_jobs, settings.queue_name, cu, cv);

% access the temp_RHO folder and initialize the bash script for
% parallel computing of the hyperparameter optimization
cd(hyperopt_folder);
system(['qsub ' hyperopt_bash]);

% set a variable which counts the files in the RHO folder
filecount = size(dir([hyperopt_folder '/RHO_*.txt']),1);
RHO_collection = nan(size(cu_cv_combination,2),1);

% use a while loop which recounts all the files in the RHO folder
% until all computations from the parallelized hyperparameter
% optimization step are completed and all files are saved in that
% folder, the while loop will be exited as soon as all files are
% saved in the folder

while filecount < size(cu_cv_combination,2)
    filecount = size(dir([hyperopt_folder '/RHO_avg_*.txt']),1);
end

% load all RHO values which were saved in mat files into a
% RHO_collection vector

for i=1:size(cu_cv_combination,2)
    if exist([hyperopt_folder '/RHO_avg_' num2str(i) '.txt'],'file') 
        RHO_avg_collection(i,1) = dp_txtscan([hyperopt_folder '/RHO_avg_' num2str(i) '.txt'], '%f');        
        delete([hyperopt_folder '/RHO_avg_' num2str(i) '.txt']);         
    else        
        RHO_avg_collection(i,1) = NaN;        
    end   
end

        
for k=1:K
    RHO_collection(k) = dp_k_split(tp, keep_in_data_x, keep_in_data_y, cu, cv);
end

nan_log = isnan(RHO_collection);

while sum(nan_log)>0
    for ii=1:size(RHO_collection,1)
        if isnan(RHO_collection(ii))
            RHO_collection(ii) = dp_k_split(tp, keep_in_data_x, keep_in_data_y, cu, cv);
        end
    end
    nan_log = isnan(RHO_collection);
end

RHO_avg = mean(RHO_collection);

end