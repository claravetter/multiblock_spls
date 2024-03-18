%% new function for permutation testing

function dp_perm_testing(i, analysis_folder)

load([analysis_folder '/opt_param.mat']);

perm_rows = dp_txtscan([analysis_folder '/perm_' i '.txt'], '%f');
perm_data_y = keep_in_data_y(perm_rows,:);
       
[u_b, v_b, ~]=spls_suppressed_display(keep_in_data_x, perm_data_y, cu_opt, cv_opt);

% compute the absolute correlation between the hold_out data
% and the permuted u_b and v_b
RHO_b = abs(corr(hold_out_data_x*u_b,hold_out_data_y*v_b));

% save the calculated RHO_b value
FID_RHO = fopen([analysis_folder '/RHO_b_' i '.txt'],'w');
fprintf(FID_RHO, '%.4f', RHO_b);
fclose(FID_RHO);

delete([analysis_folder '/perm_' i '.txt']);

end