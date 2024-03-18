%% function for RHO parallel computing

function RHO_ICV = dp_ICV_main_txt(spls_standalone_path, queue_name, analysis_folder, type_analysis, total_jobs, parallel_jobs, mem_request, time_limit)

RHO_bash = dp_ICV_bash_job(spls_standalone_path, queue_name, analysis_folder, type_analysis, total_jobs, parallel_jobs, mem_request);

cd(analysis_folder);
system(['qsub ' RHO_bash]);

% see hyperparameter optimization
target = [analysis_folder '/RHO_results_*.csv'];
mydir = size(dir(target),1);
RHO_collection = [];

tic
while mydir<total_jobs
    mydir = size(dir(target),1);
    elapsedTime = toc;
    if elapsedTime > (time_limit*60)
        for i=1:total_jobs
            if ~exist([analysis_folder '/RHO_results_' num2str(i) '.csv'])
                switch type_analysis
                    case 'hyperopt'
                        dp_ICV_hyperopt(num2str(i), analysis_folder);
                    case 'permutation'
                        dp_ICV_permutation(num2str(i), analysis_folder);
                end
            end
        end
        mydir = size(dir(target),1);
    end
end

for i=1:total_jobs
    path_file = [analysis_folder '/RHO_results_' num2str(i),'.csv'];
    RHO_collection_temp = readmatrix(path_file);
    RHO_collection = [RHO_collection; RHO_collection_temp];
    delete(path_file);
end

RHO_ICV = RHO_collection;

end