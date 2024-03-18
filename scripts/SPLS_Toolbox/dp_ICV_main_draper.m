%% function for RHO parallel computing

function [RHO_ICV, u_ICV, v_ICV] = dp_ICV_main_draper(spls_standalone_path, queue_name, analysis_folder, type_analysis, total_jobs, time_limit)

RHO_bash = dp_ICV_bash_job(spls_standalone_path, queue_name, analysis_folder, type_analysis, total_jobs);

cd(analysis_folder);
system(['qsub ' RHO_bash]);

% see hyperparameter optimization
read_in = 0;
RHO_collection = [];
u_collection = [];
v_collection = [];

tic
while read_in<total_jobs
    elapsedTime = toc;
    if elapsedTime > (time_limit*60)
        for i=1:total_jobs
            if ~exist([analysis_folder '/RHO_results_' num2str(i) '.mat'])
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
    
    pause('on');
    pause(10);
    
    for i=1:total_jobs
        if exist([analysis_folder '/RHO_results_' num2str(i) '.mat'], 'file')
            path_mat = [analysis_folder '/RHO_results_' num2str(i),'.mat'];
            m = matfile(path_mat);
            RHO_collection_temp = m.RHO_collection_ICV;
            u_collection_temp = m.u_collection_ICV;
            v_collection_temp = m.v_collection_ICV;
            
            RHO_collection = [RHO_collection; RHO_collection_temp];
            u_collection = [u_collection; u_collection_temp];
            v_collection = [v_collection; v_collection_temp];
            
            temp = load(path_mat);
            temp_names = fieldnames(temp);
            
            delete(path_mat);
            
            read_in = read_in+1;
        end
    end
    
end
RHO_ICV = RHO_collection;
u_ICV = u_collection;
v_ICV = v_collection;
% epsilon_ICV = [];
% omega_ICV = [];

end