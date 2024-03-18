%% new function for permutation testing

function dp_RHO_b_100_sets(i, size_sets_str, analysis_folder)

load([analysis_folder '/opt_param.mat']);
% dp_txt_write(analysis_folder, ['init_' i],'initialized','%s \n');

size_sets = str2double(size_sets_str);
RHO_b_collection = nan(size_sets,1);

for pp=1:str2double(i)
    permmat = nk_PermInd2(size_sets, keep_in_Diag);
end

for ii=1:size_sets
    
    % perform procrustean transformation to minimize rotation effects of
    % permutated y matrix
    perm_data_y = keep_in_data_y(permmat(ii,:),:);
    
    [u_b, v_b, ~]=dp_spls_resample(keep_in_data_x, perm_data_y, cu_opt, cv_opt, V_opt);
    
    % compute the absolute correlation between the hold_out data
    % and the permuted u_b and v_b
    RHO_b_collection(ii,1) = abs(corr(hold_out_data_x*u_b,hold_out_data_y*v_b, 'Type', correlation_method));
%     FID = fopen([analysis_folder, '/init_' i '.txt'], 'a');
%     fprintf(FID, '%d \n', ii);
%     fclose(FID);
end

% save the calculated RHO_b value
try
    dp_txt_write(analysis_folder, ['RHO_' i], RHO_b_collection, '%.4f\n');
catch ME
    ME.message
    save([analysis_folder, '/RHO_', i, '.mat'], 'RHO_b_collection');
end

end