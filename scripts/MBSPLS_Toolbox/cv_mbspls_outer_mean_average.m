function preliminary_result_mean_average = cv_mbspls_outer_mean_average(resultfile, analysisfolder)

%resultfile = [analysisfolder, 'final_results/preliminary_result.mat'];
load(resultfile)
%%
lv_name = ['LV_' num2str(input.n_LV)];
temp = output.opt_parameters.(lv_name);
temp_cat = [];
for ii=1:size(temp,1)
    temp_cat=cat(1, temp_cat, temp{ii,3});
end
%%
[numRows, numCols] = size(temp_cat);

% Initialize a cell array to store mean vectors
mean_vectors = cell(1, numCols);

% Iterate over columns
for col = 1:numCols
    % Extract the column data
    col_data = temp_cat(:, col);

    % Concatenate data vertically
    concatenated_data = cat(2, col_data{:});

    % Compute mean along rows
    mean_vectors{col} = mean(concatenated_data, 2);
end

%%
output.final_parameters{input.n_LV, 3} = mean_vectors; 
output.final_parameters{input.n_LV,5} = mean([temp{:,5}],2);

%%

preliminary_result_mean_average = [analysisfolder '/preliminary_result_mean_LV' num2str(input.n_LV) '.mat'];
save(preliminary_result_mean_average, 'input', 'setup', 'output' )

end