%% DP K-Splits
function [holdout_X, keepin_X, holdout_Y, keepin_Y] = dp_k__holdout_splits(output_folder, input_folder, hp_file, K_file, X_file, Y_file)

addpath /volume/DP_FEF/ScrFun/ScriptsRepository/

input_names = {hp_file K_file X_file Y_file};

for i=1:numel(input_names)
    formatSpec = '%f';
    fileID = fopen([input_folder input_names{i}],'r'); 
    switch input_names{i}
        case hp_file
            hp = fscanf(fileID,formatSpec);
        case K_file
            K = fscanf(fileID,formatSpec);
        case X_file
            fileID_columns = fopen([input_folder columns_X_file],'r'); 
            columns = fscanf(fileID_columns,formatSpec);
            size=[columns Inf];
            X = [fscanf(fileID, formatSpec, size)]';
        case Y_file
            fileID_columns = fopen([input_folder columns_Y_file],'r'); 
            columns = fscanf(fileID_columns,formatSpec);
            size=[columns Inf];
            Y = [fscanf(fileID, formatSpec, size)]';
    end
end


temp = num2cell(zeros(K,4));
temp_columns = zeros(K,4);
for k=1:K
    % separate the keep_in data into training and test data according
    % to the chosen test percentage tp
    [holdout_X, keepin_X, holdout_Y, keepin_Y] = dp_partition_holdout(hp, X, Y);
    temp(k,:) = {holdout_X, keepin_X, holdout_Y, keepin_Y};
    temp_columns(k,:) = [size(holdout_X,2), size(keepin_X,2), size(holdout_Y,2), size(keepin_Y,2)];
end

output_names = {'holdout_X' 'keepin_X' 'holdout_Y' 'keepin_Y'};

for i=1:size(temp,1)
    for ii=1:size(temp,2)
        M = temp{i,ii};
        dlmwrite([output_folder output_names{ii} '_split_' num2str(i) '.txt'],M,'delimiter','\t');
        C = temp_columns(i,ii);
        dlmwrite([output_folder 'columns_' output_names{ii} '_split_' num2str(i) '.txt'],C); 
    end
end
