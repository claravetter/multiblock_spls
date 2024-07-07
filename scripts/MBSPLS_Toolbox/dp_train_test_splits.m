 %% DP K-Splits
function [train_test_matrix] = dp_train_test_splits(output_folder, X, Y, K)

% addpath /volume/DP_FEF/ScrFun/ScriptsRepository/
% 
% input_names = {hp_file K_file X_file Y_file};
% 
% for i=1:numel(input_names)
%     formatSpec = '%f';
%     fileID = fopen([input_folder input_names{i}],'r'); 
%     switch input_names{i}
%         case hp_file
%             hp = fscanf(fileID,formatSpec);
%         case K_file
%             K = fscanf(fileID,formatSpec);
%         case X_file
%             fileID_columns = fopen([input_folder columns_X_file],'r'); 
%             columns = fscanf(fileID_columns,formatSpec);
%             size=[columns Inf];
%             X = [fscanf(fileID, formatSpec, size)]';
%         case Y_file
%             fileID_columns = fopen([input_folder columns_Y_file],'r'); 
%             columns = fscanf(fileID_columns,formatSpec);
%             size=[columns Inf];
%             Y = [fscanf(fileID, formatSpec, size)]';
%     end
% end


temp = num2cell(zeros(K,4));
temp_columns = zeros(K,4);
for k=1:K
    % separate the keep_in data into training and test data according
    % to the chosen test percentage tp
    [test_X, train_X, test_Y, train_Y] = dp_partition_holdout(hp, X, Y);
    temp(k,:) = {test_X, train_X, test_Y, train_Y};
    temp_columns(k,:) = [size(test_X,2), size(train_X,2), size(test_Y,2), size(train_Y,2)];
end

train_test_matrix = temp;
train_test_matrix_names = {'holdout_X' 'keepin_X' 'holdout_Y' 'keepin_Y'};

for i=1:size(temp,1)
    for ii=1:size(temp,2)
        M = temp{i,ii};
        dlmwrite([output_folder train_test_matrix_names{ii} '_split_' num2str(i) '.txt'],M,'delimiter','\t');
        C = temp_columns(i,ii);
        dlmwrite([output_folder 'columns_' train_test_matrix_names{ii} '_split_' num2str(i) '.txt'],C); 
    end
end

end
