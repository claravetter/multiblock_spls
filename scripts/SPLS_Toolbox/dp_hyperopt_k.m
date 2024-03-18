%% new function for cu/cv combination and 100 splits

function [RHO_avg] = dp_hyperopt_k(output_folder, input_folder, cu_cv_file, keepin_X_split, columns_X_file, keepin_Y_split, columns_Y_file)

input_names = {cu_cv_file keepin_X_split keepin_Y_split};

for i=1:numel(input_names)
    formatSpec = '%f';
    fileID = fopen([input_folder input_names{i}],'r'); 
    switch input_names{i}
        case keepin_X_split
            fileID_columns = fopen([input_folder columns_X_file],'r'); 
            columns = fscanf(fileID_columns,formatSpec);
            size=[columns Inf];
            keep_in_data_x = [fscanf(fileID, formatSpec, size)]';
        case keepin_Y_split
            fileID_columns = fopen([input_folder columns_Y_file],'r'); 
            columns = fscanf(fileID_columns,formatSpec);
            size=[columns Inf];
            keep_in_data_y = [fscanf(fileID, formatSpec, size)]';
        case cu_cv_file
            cu_cv = [fscanf(fileID, formatSpec)]';    
    end
end

% create an inner loop with k splits and a tp test perecentage set within
% each split

% RHO_avg is a vector where all average RHO values of the different cu/cv
% combinations are stored, preallocated for speed       

% separate the keep_in data into training and test data according to the
% chosen test percentage tp
[test_data_x, training_data_x, test_data_y, training_data_y] = dp_partition_holdout(tp, keep_in_data_x, keep_in_data_y);

%perform SPLS on the training data using the current cu/cv combination
[u, v, ~] = spls(training_data_x,training_data_y,cu_index,cv_index); 

%compute the correlation between the projections of the training and test
%matrices onto the SPLS latent space spanned by the weight vectors
RHO = abs(corr(test_data_x*u,test_data_y*v)); 

dlmwrite([output_folder cu_cv_file '_RHO.txt'],M,'delimiter','\t');

end