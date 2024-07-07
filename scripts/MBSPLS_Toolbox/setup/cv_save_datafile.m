function input = cv_save_datafile(input, setup, custom_path)

if exist('custom_path', 'var') && ~isempty(custom_path)
    analysis_folder = custom_path;
else
    analysis_folder = setup.analysis_folder;
end
modalities_str = strjoin(input.Xs_names, '_');

if isnumeric(input.matrix_norm)
    matrix_norm_str = num2str(input.matrix_norm);
else 
    matrix_norm_str = input.matrix_norm;
end

name_part1 = ['CV_mbspls_', ...
    setup.cluster, '_', ...
    input.project_name, '_', ...
    modalities_str, '_', ...
    num2str(input.outer_folds), 'x', num2str(input.inner_folds), '_', ...
    num2str(input.permutation_testing), 'perm_', ...
    num2str(input.bootstrap_testing), 'boot_', ...
    matrix_norm_str, '_matrixnorm'];

switch input.optimization_strategy
    case 'randomized_search'
        input.name = [name_part1, '_', input.optimization_strategy, '_', num2str(input.randomized_search_params.randomized_search_iterations), '_iterations'];
    case 'grid_search'
        density_str = strjoin(cellfun(@num2str, input.density, 'UniformOutput', false), '_');
        input.name = [name_part1, '_', input.optimization_strategy, '_', density_str, '_densities'];
end

input.datafile = fullfile(analysis_folder, [setup.date, '_', input.name, '_datafile.mat']); % Path for storing datafile containing input and setup

mkdir(analysis_folder); % Creates analysis folder
save([input.datafile], 'setup', 'input'); % Saves datafile in analysis folder

end