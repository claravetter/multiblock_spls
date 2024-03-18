%% function to find best opt_param combination

function [best_p, index_best_p] = dp_find(matrix, matrix_names)

opt_RHO = strcmp(matrix_names, 'RHO');
opt_p = strcmp(matrix_names, 'p');
% remove parameter combinations with nans as RHO_opt
to_rmv=isnan([matrix{:,opt_RHO}]');
matrix(to_rmv,:)=[];
matrix_temp = matrix;

%change RHO values to absolute values
matrix_temp(:, opt_RHO) = num2cell(abs(cell2mat(matrix_temp(:, opt_RHO))));

%create a logical indexing the lowest p-value(s)
p_min = min([matrix_temp{:,opt_p}]);
I_pmin = [matrix_temp{:,opt_p}] == p_min;
temp1 = matrix_temp(I_pmin,:); % placeholder to make syntax easier

if size(temp1,1) > 1 % if there is more than one p-value that is equal to the lowest p-value, ie there are two absolute minimal p-values, then continue in the loop
    RHO_max = max([temp1{:,opt_RHO}]);
    I_RHO_max = [temp1{:,opt_RHO}] == RHO_max; % logical indexing the highest RHO value
    if sum(I_RHO_max)>1
        I_RHO_max_new = find(I_RHO_max==1,1);
        best = temp1(I_RHO_max_new,:); 
    else
        best = temp1(I_RHO_max,:);% take the row with the lowest p-value and the highest RHO value
    end
else
    best = temp1; % if there is just one absolute minimal p value, then take this row
end

index_best_p = best{1};
best_p = matrix(index_best_p, :);

end