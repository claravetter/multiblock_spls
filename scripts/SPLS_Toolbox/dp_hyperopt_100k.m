%% function for one hyperoptimization step

function [RHO_avg] = dp_hyperopt_100k(folders, settings, K, tp, cu, cv)

hyperopt_folder = folders.analysis_folder;
hyperopt_bash = dp_bash_100k(folders, settings, K, tp, cu, cv);

% access the temp_RHO folder and initialize the bash script for
% parallel computing of the hyperparameter optimization
cd(hyperopt_folder);
system(['qsub ' hyperopt_bash]);

% set a variable which counts the files in the RHO folder
filecount = size(dir([hyperopt_folder '/RHO_*.txt']),1);
RHO_collection = nan(K,1);

% use a while loop which recounts all the files in the RHO folder
% until all computations from the parallelized hyperparameter
% optimization step are completed and all files are saved in that
% folder, the while loop will be exited as soon as all files are
% saved in the folder

while filecount < K
    filecount = size(dir([hyperopt_folder '/RHO_*.txt']),1);
end

% load all RHO values which were saved in mat files into a
% RHO_collection vector

for i=1:K
    if exist([hyperopt_folder '/RHO_' num2str(i) '.txt'],'file') 
        RHO_collection(i,1) = dp_txtscan([hyperopt_folder '/RHO_' num2str(i) '.txt'], '%f');        
        delete([hyperopt_folder '/RHO_' num2str(i) '.txt']);         
    else        
        RHO_collection(i,1) = NaN;        
    end   
end

RHO_avg = mean(RHO_collection);

end