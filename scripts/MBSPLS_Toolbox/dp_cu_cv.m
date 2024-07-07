%% DP function for cu and cv setup

function [cu_cv_combination, hyperopt_sets] = dp_cu_cv(X, Y, grid_x, grid_y, size_sets_hyperopt)

% % check if previous cu/cv files exist and clear the hyperopt folder
% cu_dir = dir([hyperopt_folder, '/cu*']);
% cv_dir = dir([hyperopt_folder, '/cv*']);
% if size(cu_dir,1)>0
%     for i=1:size(cu_dir,1)
%         delete([hyperopt_folder, '/' cu_dir(i).name]);
%     end
% end
% 
% if size(cv_dir,1)>0
%     for i=1:size(cv_dir,1)
%         delete([hyperopt_folder, '/' cv_dir(i).name]);
%     end
% end

% cu_range and cv_range define the range for the grid search. The grid
% search is performed along 20 points from 1 to sqrt(number of variables)
% as proposed in Monteiro et al. 2016.

% create the original 100-point grid between 1 and sqrt(size(x,2))
cu_range_temp_1 = linspace(1,sqrt(size(X,2)),20);
cu_range_temp = linspace(cu_range_temp_1(2),cu_range_temp_1(end),100);
cv_range_temp_1 = linspace(1,sqrt(size(Y,2)),20);
cv_range_temp = linspace(cv_range_temp_1(2),cv_range_temp_1(end),100);

% check if there are new start and end points for the grid, if not, then
% use defaults
if ~isfield(grid_x, 'start')
    grid_x.start = 1;
end

if ~isfield(grid_x, 'end')
    grid_x.end = 0;
end

if ~isfield(grid_x, 'density')
    grid_x.density = 20;
end

if ~isfield(grid_y, 'start')
    grid_y.start = 1;
end

if ~isfield(grid_y, 'end')
    grid_y.end = 0;
end

if ~isfield(grid_y, 'density')
    grid_y.density = 20;
end

% apply start and end points to generic grid
cu_range_points = linspace(cu_range_temp(grid_x.start),cu_range_temp(end-grid_x.end),grid_x.density);
cv_range_points = linspace(cv_range_temp(grid_y.start),cv_range_temp(end-grid_y.end),grid_y.density);

% compile a matrix with separate row for all possible cu and cv
% combinations by taking cu and repeating every single element X times and
% then takin cv and repeating the entire vector X times
cu_cv_combination = zeros(size(cu_range_points,2)*size(cv_range_points,2),2);
nn=1;
for i=1:size(cu_range_points,2)
    for ii=1:size(cv_range_points,2)
        cu_cv_combination(nn,:) = [cu_range_points(i), cv_range_points(ii)];
        nn=nn+1;
    end
end

% set up cu and cv for hyperparameter optimization
rest_hyperopt = mod(size(cu_cv_combination,1),size_sets_hyperopt);
if rest_hyperopt>0
    hyperopt_sets = ((size(cu_cv_combination,1) - rest_hyperopt)/size_sets_hyperopt)+1;
else
    hyperopt_sets = ((size(cu_cv_combination,1) - rest_hyperopt)/size_sets_hyperopt);
end

% nn=1;
% for i=1:hyperopt_sets
%     for ii=1:size_sets_hyperopt
%         try
%             temp1 = cu_cv_combination(nn,1);
%             dp_txt_write(hyperopt_folder, ['cu_', num2str(i), '_', num2str(ii)], temp1, '%.4f');
%             temp2 = cu_cv_combination(nn,2);
%             dp_txt_write(hyperopt_folder, ['cv_', num2str(i), '_', num2str(ii)], temp2, '%.4f');
%             nn=nn+1;
%         catch
%             break
%         end
%     end
% end

end