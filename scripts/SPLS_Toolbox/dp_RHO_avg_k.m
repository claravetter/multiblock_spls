%% new function for cu/cv combination and 100 splits

function [RHO_avg] = dp_RHO_avg_k(output_folder, input_folder, cu_cv_file, K_file, keepin_X_split, columns_X_file, keepin_Y_split, columns_Y_file)

input_names = {cu_cv_file K_file keepin_X_split keepin_Y_split};

for i=1:numel(input_names)
    formatSpec = '%f';
    fileID = fopen([input_folder input_names{i}],'r'); 
    switch input_names{i}
        case keepin_X_split
            fileID_columns = fopen([input_folder columns_X_file],'r'); 
            columns = fscanf(fileID_columns,formatSpec);
            size=[columns Inf];
            keep_in_data_x = fscanf(fileID, formatSpec, size)';
        case keepin_Y_split
            fileID_columns = fopen([input_folder columns_Y_file],'r'); 
            columns = fscanf(fileID_columns,formatSpec);
            size=[columns Inf];
            keep_in_data_y = fscanf(fileID, formatSpec, size)';
        case cu_cv_file
            cu_cv = fscanf(fileID, formatSpec)';  
        case K_file
            K = fscanf(fileID, formatSpec);  
    end
end

RHO_collection = zeros(K,1);

for k=1:K
            
    % separate the keep_in data into training and test data according to
    % the chosen test percentage tp
    [test_data_x, training_data_x, test_data_y, training_data_y] = dp_partition_holdout(tp, keep_in_data_x, keep_in_data_y);
        
    %perform SPLS on the training data using the current cu/cv combination
    [u, v, ~] = spls(training_data_x,training_data_y,cu_cv_combination(index_cu,i),cu_cv_combination(index_cv,i)); 
                    
    %compute the correlation between the projections of the training and
    %test matrices onto the SPLS latent space spanned by the weight vectors
    RHO = abs(corr(test_data_x*u,test_data_y*v)); 
                
    %store the results  
    RHO_collection(k) = RHO;
                
end


% computing the mean value of all RHO calculated for the k iterations for
% one specific cu/cv combination
            
RHO_avg = mean(cell2mat(RHO_collection(:,index_RHO_k))); 

M = [cu_cv RHO_avg];
dlmwrite([output_folder cu_cv_file '_RHO.txt'],M,'delimiter','\t');

end
