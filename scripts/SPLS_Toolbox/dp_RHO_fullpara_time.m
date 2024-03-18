%% function for RHO parallel computing

function [RHO_output, success] = dp_RHO_fullpara_time(spls_standalone_path, analysis_folder, type_analysis, mem_total, max_sim_jobs, queue_name, total_jobs, size_sets, time_limit)

RHO_bash = dp_bash_parallel(spls_standalone_path, analysis_folder, type_analysis, mem_total, max_sim_jobs, queue_name, total_jobs, size_sets);

cd(analysis_folder);
system(['qsub ' RHO_bash]);

% see hyperparameter optimization
mydir = size(dir([analysis_folder '/RHO_*.txt']),1);
RHO_collection = [];

tic
while mydir<total_jobs
    mydir = size(dir([analysis_folder '/RHO_*.txt']),1);
    elapsedTime = toc;
    if elapsedTime > (time_limit*60)
        for i=1:total_jobs
            if ~exist([analysis_folder '/RHO_' num2str(i) '.txt'])
                switch type_analysis
                    case 'hyperopt'
                        dp_RHO_avg_10_sets(num2str(i), num2str(size_sets), analysis_folder);
                    case 'permutation'
                        dp_RHO_b_100_sets(num2str(i), num2str(size_sets), analysis_folder);
                end
            end
        end
        mydir = size(dir([analysis_folder '/RHO_*.txt']),1);
    end
end



for i=1:total_jobs
    RHO_collection_temp = [];
    RHO_collection_temp = dp_txtscan([analysis_folder '/RHO_' num2str(i),'.txt'], '%f\n');
    delete([analysis_folder '/RHO_' num2str(i),'.txt']);
    RHO_collection = [RHO_collection; RHO_collection_temp];
end

RHO_output = RHO_collection;


if sum(isnan(RHO_output))>(0.1*size(RHO_output,1))
    success=false;
else
    success=true;
end

end