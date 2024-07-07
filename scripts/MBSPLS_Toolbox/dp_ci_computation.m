%% function for post analysis CI computation

function [output_CI_u, output_CI_v] = dp_ci_computation(opt_datamatrix, final_LV, number_bootstraps)

%% compute CI for item weights using bootstrapping of 10 opt_parameter iterations

% first find inverted weight vectors in turn them around in order to
% avoid problems with mean computation during bootstrapping
inv=struct; inv.to_extract = 'v';
[opt_datamatrix, inv] = dp_inverse_array(opt_datamatrix, inv);
inv.to_extract = 'u';
[opt_datamatrix, ~] = dp_inverse_array(opt_datamatrix, inv);

disp('checkpoint inversion');

% set features to zero, which are also zero in the final LV
temp = opt_datamatrix(:,4:5);
log_zero = {final_LV{4}==0, final_LV{5}==0};
for i=1:size(temp,1)
    temp{i,1}(log_zero{1}) = 0;
    temp{i,2}(log_zero{2}) = 0;
end

opt_datamatrix(:,4:5)=temp;

disp('checkpoint zero');

% vector u
[ci_collection_temp, ~] = dp_bootstrap_opt(opt_datamatrix, number_bootstraps, 'u');
output_CI_u =[abs(ci_collection_temp(:,1)-final_LV{4}),abs(ci_collection_temp(:,2)-final_LV{4})];

disp('checkpoint bootstrap u');

% vector v
[ci_collection_temp, ~] = dp_bootstrap_opt(opt_datamatrix, number_bootstraps, 'v');
output_CI_v =[abs(ci_collection_temp(:,1)-final_LV{5}),abs(ci_collection_temp(:,2)-final_LV{5})];

disp('checkpoint bootstrap v');

end