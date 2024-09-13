% compilation

%% first add required code repositories 
clear all
addpath(genpath('/volume/projects/CV_gs_PLS/ScrFun/multiblock_spls/scripts/MBSPLS_Toolbox'))
addpath(genpath('/opt/PrecisionCodeRep/SPLS_Toolbox/SPLS_Toolbox'))
addpath(genpath('/opt/NM/NeuroMiner_Current'))

matlab_version = ['R', version('-release')]; 

%% change
compiled_mbspls_directory = '/volume/projects/CV_gs_PLS/ScrFun/Compilations/';
compiled_mbspls_name = ['MBSPLS_DEV_Aug2024_bootsample_' matlab_version];
compiled_mbspls_version = 'v1';
mbspls_directory = '/volume/projects/CV_gs_PLS/ScrFun/multiblock_spls/scripts/MBSPLS_Toolbox';
mbspls_main_script = fullfile(mbspls_directory, 'cv_mbspls_main.m');
mbspls_hyperopt_script = fullfile(mbspls_directory, 'cv_mbspls_ICV_hyperopt_csv.m');
mbspls_permutation_script = fullfile(mbspls_directory, 'cv_mbspls_ICV_permutation_csv.m');
mbspls_bootstrap_script = fullfile(mbspls_directory, 'cv_mbspls_ICV_bootstrap_csv.m');

%% if hyperopt, etc are recompiled, they need to be added to the SLURM directory too 
% good to know: when recompiling the SLURM main script, the previously added folders (hyperopt etc) are
% not being removed
compile_main = 1;
compile_hyperopt = 1; 
compile_permutation = 1;
compile_bootstrap = 1; 

%% add compilation notes 
log_file = fullfile(compiled_mbspls_directory, 'compilation.log');
[~, username] = system('whoami');
username = strtrim(username); 
user_note = 'bootstrap debugging'; % change
current_time = datestr(now, 'yyyy-mm-dd HH:MM:SS'); 
log_entry = sprintf('\n\n%s - %s\nCompiled by %s\nChanges: %s \nMain: %d, Hyperopt: %d, Permutation: %d, Bootstrap: %d\nPath: %s', ...
    current_time, compiled_mbspls_name, username, user_note, compile_main, compile_hyperopt, compile_permutation, compile_bootstrap, fullfile(compiled_mbspls_directory, compiled_mbspls_name));

fid = fopen(log_file, 'a');
if fid == -1
    error('Could not open log file for writing.');
end
fprintf(fid, '%s', log_entry);
fclose(fid);

%% compile main mbspls run script
if compile_main
    eval(['mcc -m', ...
        ' -o ', compiled_mbspls_name, ...
        ' -W main:', compiled_mbspls_name, '_', compiled_mbspls_version, '_', matlab_version, ...
        ' -T link:exe -d ', fullfile(compiled_mbspls_directory, compiled_mbspls_name, 'for_testing'), ...
        ' -v ', mbspls_main_script]);


end

%% compile helper modules (hyperopt, bootstrap, permutation)
if compile_hyperopt
    eval(['mcc -m', ...
        ' -o hyperopt_mbspls', ...
        ' -W main:', compiled_mbspls_name, ...
        ' -T link:exe -d ', fullfile(compiled_mbspls_directory, compiled_mbspls_name, 'hyperopt_mbspls', 'for_testing'), ...
        ' -v ', mbspls_hyperopt_script]);
end

if compile_permutation
    eval(['mcc -m', ...
        ' -o permutation_mbspls', ...
        ' -W main:', compiled_mbspls_name, ...
        ' -T link:exe -d ', fullfile(compiled_mbspls_directory, compiled_mbspls_name, 'permutation_mbspls', 'for_testing'), ...
        ' -v ', mbspls_permutation_script]);
end

if compile_bootstrap
    eval(['mcc -m', ...
        ' -o bootstrap_mbspls', ...
        ' -W main:', compiled_mbspls_name, ...
        ' -T link:exe -d ', fullfile(compiled_mbspls_directory, compiled_mbspls_name, 'bootstrap_mbspls', 'for_testing'), ...
        ' -v ', mbspls_bootstrap_script]);
end

%%
disp(['scp -r ', fullfile(compiled_mbspls_directory, compiled_mbspls_name), ' clvetter@srvcorem2:/data/core-psy-archive/projects/CV_mbspls/mbspls/'])