%% new function for cu/cv combination and 100 splits

function dp_RHO_avg_100k(i, cu_str, cv_str, tp_str, hyperopt_folder)

cu = str2double(cu_str); cv = str2double(cv_str); tp = str2double(tp_str);
load([hyperopt_folder '/keep_in_partition.mat']);
RHO = dp_k_split(tp, keep_in_data_x, keep_in_data_y, cu, cv);

FID_RHO = fopen([hyperopt_folder '/RHO_' i '.txt'],'w');
fprintf(FID_RHO, '%.4f', RHO);
fclose(FID_RHO);

end
