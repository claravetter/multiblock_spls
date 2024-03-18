%% function for RHO parallel computing

function [RHO_ICV, u_ICV, v_ICV] = dp_ICV_main(spls_standalone_path, queue_name, analysis_folder, type_analysis, total_jobs, parallel_jobs, mem_request, time_limit)

RHO_bash = dp_ICV_bash_job(spls_standalone_path, queue_name, analysis_folder, type_analysis, total_jobs, parallel_jobs, mem_request);

cd(analysis_folder);
system(['qsub ' RHO_bash]);

% see hyperparameter optimization
target = [analysis_folder '/RHO_results_*.mat'];
mydir = size(dir(target),1);
RHO_collection = [];
u_collection = [];
v_collection = [];

tic
while mydir<total_jobs
    mydir = size(dir(target),1);
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
end

pause('on');
pause(10);

for i=1:total_jobs
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
    
%     RHO_collection_temp = temp.(temp_names{~cellfun(@isempty,strfind(temp_names, 'RHO'))});
%     RHO_collection = [RHO_collection; RHO_collection_temp];
%     
%     u_collection_temp = temp.(temp_names{~cellfun(@isempty,strfind(temp_names, 'u_collection'))});
%     v_collection_temp = temp.(temp_names{~cellfun(@isempty,strfind(temp_names, 'v_collection'))});
%     
%     u_collection = [u_collection; u_collection_temp];
%     v_collection = [v_collection; v_collection_temp];
    
%     try epsilon_collection_temp = temp.(temp_names{contains(temp_names, 'epsilon_collection')});
%         omega_collection_temp = temp.(temp_names{contains(temp_names, 'omega_collection')});
%         
%         epsilon_collection = [epsilon_collection; epsilon_collection_temp];
%         omega_collection = [omega_collection; omega_collection_temp];
%     catch
%         epsilon_collection=[];omega_collection=[];
%     end
    delete(path_mat);
end

RHO_ICV = RHO_collection;
u_ICV = u_collection;
v_ICV = v_collection;
% epsilon_ICV = [];
% omega_ICV = [];

end