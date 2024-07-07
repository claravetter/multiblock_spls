%% fill in missing RHO values


filecount = size(dir([analysis_folder '/RHO_avg_*']),1);
RHO_collection = nan(size(cu_cv_combination,2),1);  
dir_RHO = dir([analysis_folder '/RHO_avg_*']);
    
for i=1:size(RHO_collection,1)
    if exist(['RHO_avg_' num2str(i) '.mat'],'file')
        load([analysis_folder '/RHO_avg_',num2str(i),'.mat']);
        RHO_collection(i,1) = RHO_avg;
%         delete([analysis_folder '/RHO_avg_',num2str(i),'.mat']); 
    else
        RHO_collection(i,1) = NaN;
    end
end

missing_rho = isnan(RHO_collection);

temp = 1:1600;
ID_missing = temp(missing_rho)';
dlmwrite('missing_ID.txt',ID_missing)

   
cd(analysis_folder) %'/temp_RHO']);
system('qsub /volume/DP_FEF/ScrFun/ScriptsRepository/SPLS_Toolbox/Hyperparameter_optimization/dp_rho_missing.sh');