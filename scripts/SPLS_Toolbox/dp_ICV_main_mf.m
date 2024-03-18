%% function for sending and collecting job array results

function [RHO_ICV, u_ICV, v_ICV] = dp_ICV_main_mf(spls_standalone_path, analysis_folder, type_analysis, total_jobs)

RHO_bash = dp_ICV_bash_job(spls_standalone_path, analysis_folder, type_analysis, total_jobs);

cd(analysis_folder);
system(['qsub ' RHO_bash]);

% m = matfile([analysis_folder '/RHO_results.mat']);
% while any(cellfun(@isempty,m.u_collection))
% end
% 
% RHO_ICV = m.RHO_collection;
% u_ICV = m.u_collection;
% v_ICV = m.v_collection;

end