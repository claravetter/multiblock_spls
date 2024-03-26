%% DP function to set up parameters

function [input, setup, matrices, X, Y, B, K, W, OB, IB, size_sets_permutation, size_sets_bootstrap, correlation_method, cs_method, selection_train, selection_retrain, correction_target] = dp_setup_parameters(input, setup)

% TO DO: check if paths are specified instead of matrices, cell fun 
if isfield(input, 'matrices')
    matrices = cellfun(@load_matrices,input.matrices, UniformOutput=false);
    X = []; 
    Y = [];
else 
    if isfield(input, 'X')
        try temp = load(input.X);
            % CV: what does this do?
            field = fieldnames(temp);
            X = temp.(field{1});
            clear('temp');
        catch
            X = input.X;
        end
    end

    if isfield(input, 'Y')
        try temp = load(input.Y);
            field = fieldnames(temp);
            Y = temp.(field{1});
            clear('temp');
        catch
            Y = input.Y;
        end
    end
    matrices = [];
end


if ~isfield(input, 'correct_limit')
    input.correct_limit=1;
end

if ~isfield(setup, 'max_sim_jobs')
    setup.max_sim_jobs = 40;
end

if ~isfield(input, 'permutation_testing')
    input.permutation_testing = 5000;
end
B = input.permutation_testing;

if ~isfield(input, 'inner_folds')
    input.inner_folds = 10;
end
K = input.inner_folds;

if ~isfield(input, 'outer_folds')
    input.outer_folds = 10;
end
W = input.outer_folds;

if ~isfield(input, 'outer_permutations')
    input.outer_permutations = 1;
end
OB = input.outer_permutations;

if ~isfield(input, 'inner_permutations')
    input.inner_permutations = 1;
end
IB = input.inner_permutations;

if ~isfield(input, 'validation_train')
    input.validation_train = 1;
end

if ~isfield(input, 'size_sets_permutation')
    input.size_sets_permutation = round(input.permutation_testing/setup.max_sim_jobs);
end
size_sets_permutation = input.size_sets_permutation;

if ~isfield(input, 'permutation_testing_precision')
    input.permutation_testing_precision = 'lenient';
end

if ~isfield(input, 'size_sets_bootstrap')
    input.size_sets_bootstrap = round(input.bootstrap_testing/setup.max_sim_jobs);
end
size_sets_bootstrap = input.size_sets_bootstrap;

if ~isfield(input, 'statistical_testing')
    input.statistical_testing = 1;
end

if ~isfield(input, 'correlation_method')
    input.correlation_method = 'Spearman';
end
correlation_method = input.correlation_method;

if ~isfield(input, 'cs_method')
    input.cs_method.method = 'mean-centering';
    input.cs_method.correction_subgroup = '';
end
cs_method = input.cs_method;

if ~isfield(input, 'selection_train')
    input.selection_train = 1;
end
selection_train = input.selection_train;

if ~isfield(input, 'selection_retrain')
    input.selection_retrain = 1;
end
selection_retrain = input.selection_retrain;

if ~isfield(input, 'type_correction')
    input.type_correction = 'correct';
end

if ~isfield(input, 'merge_train')
    input.merge_train = 'median';
end

if ~isfield(input, 'merge_retrain')
    input.merge_retrain = 'median';
end

if ~isfield(input, 'correct_limit')
    input.correct_limit = 1;
end

if ~isfield(input, 'final_merge')
    input.final_merge.type = 'best';
end

if ~isfield(input, 'correction_target')
    input.correction_target = 3;
end
correction_target = input.correction_target;

if ~isfield(input, 'coun_ts_limit')
    input.coun_ts_limit = 1;
end

if ~isfield(input, 'alpha_value')
    input.alpha_value = 0.05;
end

if ~isfield(input, 'grid_static') && ~isfield(input, 'grid_dynamic')
    input.grid_static = 20;
end

if isfield(input, 'grid_static')
    input.grid_dynamic.onset    = 1;
    input.grid_dynamic.LV_1.x   = struct('start', 1, 'end', 0, 'density', input.grid_static);
    input.grid_dynamic.LV_1.y   = struct('start', 1, 'end', 0, 'density', input.grid_static);
end

if ~isfield(input, 'additional_NCV')
    input.additional_NCV = false;
end

end


function matrix = load_matrices(matrix)
try temp = load(matrix);
    field = fieldnames(temp);
    matrix = temp.(field{1});
    clear('temp');
catch
    matrix = matrix;
end
end