%% new function for cu/cv combination and 100 splits

function dp_RHO_avg_k_AR_DP(i, hyperopt_folder, utilities_folder)

load([hyperopt_folder '/keep_in_partition.mat']);

%% start of the function
RHO_avg=[]; cu=[]; cv=[]; tp=[]; K=[];

fileID = fopen([hyperopt_folder '/cu_' num2str(i) '.txt']);
cu = fscanf(fileID,'%f');
fclose(fileID);

fileID = fopen([hyperopt_folder '/cv_' num2str(i) '.txt']);
cv = fscanf(fileID,'%f');
fclose(fileID);

fileID = fopen([utilities_folder '/tp.txt']);
tp = fscanf(fileID,'%d');
fclose(fileID);

fileID = fopen([utilities_folder '/K.txt']);
K = fscanf(fileID,'%d');
fclose(fileID);

RHO_collection = nan(K,1);

if cu > sqrt(size(keep_in_data_x,2))
    cu = sqrt(size(keep_in_data_x,2));
end

if cv > sqrt(size(keep_in_data_y,2))
    cv = sqrt(size(keep_in_data_y,2));
end
    
for k=1:K
    RHO_collection(k) = dp_k_split(tp, keep_in_data_x, keep_in_data_y, cu, cv);
end

nan_log = isnan(RHO_collection);
while sum(nan_log)>0
    for ii=1:size(RHO_collection,1)
        if isnan(RHO_collection(ii))
            RHO_collection(ii) = dp_k_split(tp, keep_in_data_x, keep_in_data_y, cu, cv);
        end
    end
    nan_log = isnan(RHO_collection);
end

RHO_avg = mean(RHO_collection);

FID_RHO = fopen(['RHO_avg_' num2str(i) '.txt'],'w');
fprintf(FID_RHO, '%.4f', RHO_avg);
fclose(FID_RHO);

end
