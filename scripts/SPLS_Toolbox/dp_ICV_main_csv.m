%% function for RHO parallel computing

function RHO_ICV = dp_ICV_main_csv(spls_standalone_path, queue_name, analysis_folder, type_analysis, total_jobs, parallel_jobs, mem_request)

RHO_bash = dp_ICV_bash_job(spls_standalone_path, queue_name, analysis_folder, type_analysis, total_jobs, parallel_jobs, mem_request);

cd(analysis_folder);
system(['qsub ' RHO_bash]);

% see hyperparameter optimization
target = [analysis_folder '/RHO_results_*.csv'];
mydir = size(dir(target),1);
RHO_collection = [];

while mydir<total_jobs
    mydir = size(dir(target),1);
end

pause(30);

for i=1:total_jobs
    path_file = [analysis_folder '/RHO_results_' num2str(i),'.csv'];
    RHO_collection_temp = readmatrix(path_file);
    RHO_collection = [RHO_collection; RHO_collection_temp];
    delete(path_file);
end

RHO_ICV = RHO_collection;

end