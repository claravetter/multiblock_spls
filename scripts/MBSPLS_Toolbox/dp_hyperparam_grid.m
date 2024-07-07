%% DP cu cv grid

function [cu_cv_combination] = dp_hyperparam_grid(output_folder, X, Y, range_points)

% input_names = {X_file Y_file range_points};
% 
% for i=1:numel(input_names)
%     formatSpec = '%f';
%     fileID = fopen([input_folder input_names{i}],'r'); 
%     switch input_names{i}
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

% cu_range and cv_range define the range for the grid search. The grid
% search is performed along 40 points from 1 to sqrt(number of variables)
% as proposed in Monteiro et al. 2016.
cu_range = linspace(1,sqrt(size(X,2)),range_points);
cv_range = linspace(1,sqrt(size(Y,2)),range_points);

% compile a matrix with separate row for all possible cu and cv
% combinations by taking cu and repeating every single element 40 times and
% then takin cv and repeating the entire vector 40 times
cu_cv_combination = [repelem(cu_range,numel(cu_range));repmat(cv_range,1,numel(cv_range))]; % this matrix contains all the possible combinations between cu and cv

for i=1:size(cu_cv_combination,2)
    M = [cu_cv_combination(:,i)]';
    dlmwrite([output_folder 'cu_cv_' num2str(i) '.txt'],M,'delimiter','\t');
end

end


