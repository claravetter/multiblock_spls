%% DP function for cu and cv setup

function [weights_combination, hyperopt_sets] = cv_weights_ext(matrices, grid, size_sets_hyperopt)

% cu_range and cv_range define the range for the grid search. The grid
% search is performed along 20 points from 1 to sqrt(number of variables)
% as proposed in Monteiro et al. 2016.

% create the original 100-point grid between 1 and sqrt(size(x,2))
range_temps = cellfun(@create_range, matrices, UniformOutput=false);

% check if there are new start and end points for the grid, if not, then
% use defaults

for i=1:length(matrices)
    if ~isfield(grid, 'start')
        grid(i).start = 1;
    end
    if ~isfield(grid, 'end')
        grid.end = 0;
    end
    if ~isfield(grid, 'density')
        grid.density = 20;
    end
    range_temp = range_temps{i};
    range_points{i} = linspace(range_temp(grid(i).start),range_temp(end-grid(i).end),grid(i).density);

end

% compile a matrix with separate row for all possible cu and cv
% combinations by taking cu and repeating every single element X times and
% then takin cv and repeating the entire vector X times
weights_combination = allcomb(range_points); % needs R2023b




% set up cu and cv for hyperparameter optimization
rest_hyperopt = mod(size(weights_combination,1),size_sets_hyperopt);
if rest_hyperopt>0
    hyperopt_sets = ((size(weights_combination,1) - rest_hyperopt)/size_sets_hyperopt)+1;
else
    hyperopt_sets = ((size(weights_combination,1) - rest_hyperopt)/size_sets_hyperopt);
end

end

function range = create_range(matrix)
range_temp_1 = linspace(1,sqrt(size(matrix,2)),20);
range = linspace(range_temp_1(2), range_temp_1(end),100);
end

