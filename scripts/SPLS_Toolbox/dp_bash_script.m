
%% create main bash file to submit to queue
addpath /volume/DP_FEF/ScrFun/ScriptsRepository/SPLS_Toolbox;
analysis_folder = '/volume/HCStress/Analysis/20-Apr-2018';
project_name = 'DP_CISS_parallel';
queue_name = 'mitnvp1';
full_path_script = '/volume/DP_FEF/ScrFun/ScriptsRepository';
script_name = 'DP_SPLS_para_FEF.m';
main_job_bash = dp_bash(analysis_folder, project_name, queue_name, full_path_script, script_name);