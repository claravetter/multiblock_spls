filecount = size(dir([hyperparameter_folder '/RHO_avg_*']),1);

for i=1:filecount
           
    if exist([hyperparameter_folder '/RHO_avg_' num2str(i) '.mat'])
        load([hyperparameter_folder '/RHO_avg_',num2str(i),'.mat']);
        RHO_avg_collection(i,1) = RHO_avg;
    else
        RHO_avg_collection(i,1) = NaN;
    end    
end
    

filecount = size(dir([permutation_folder '/RHO_b_*']),1);

for i=1:filecount
    if exist([permutation_folder '/RHO_b_' num2str(i) '.mat'])
        load([permutation_folder '/RHO_b_',num2str(i),'.mat']);
        RHO_b_collection(i,1) = RHO_b;
    else
        RHO_b_collection(i,1) = NaN;
    end
end