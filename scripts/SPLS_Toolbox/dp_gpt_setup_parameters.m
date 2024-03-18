% DP function to setup parameters, GPT improved

function [input, setup, X, Y, B, K, W, OB, IB, size_sets_permutation, size_sets_bootstrap, correlation_method, cs_method, selection_train, selection_retrain, correction_target] = dp_gpt_setup_parameters(input, setup)

X = get_field(input, 'X');
Y = get_field(input, 'Y');

input.correct_limit = get_default(input, 'correct_limit', 1);
setup.max_sim_jobs = get_default(setup, 'max_sim_jobs', 40);
input.permutation_testing = get_default(input, 'permutation_testing', 5000);
B = input.permutation_testing;
input.inner_folds = get_default(input, 'inner_folds', 10);
K = input.inner_folds;
input.outer_folds = get_default(input, 'outer_folds', 10);
W = input.outer_folds;
input.outer_permutations = get_default(input, 'outer_permutations', 1);
OB = input.outer_permutations;
input.inner_permutations = get_default(input, 'inner_permutations', 1);
IB = input.inner_permutations;
input.validation_train = get_default(input, 'validation_train', 1);
input.size_sets_permutation = get_default(input, 'size_sets_permutation', round(input.permutation_testing/setup.max_sim_jobs));
size_sets_permutation = input.size_sets_permutation;
input.permutation_testing_precision = get_default(input, 'permutation_testing_precision', 'lenient');
input.size_sets_bootstrap = get_default(input, 'size_sets_bootstrap', round(input.bootstrap_testing/setup.max_sim_jobs));
size_sets_bootstrap = input.size_sets_bootstrap;
input.statistical_testing = get_default(input, 'statistical_testing', 1);
correlation_method = get_default(input, 'correlation_method', 'Spearman');
cs_method = get_default(input, 'cs_method', struct('method', 'mean-centering', 'correction_subgroup', ''));
selection_train = get_default(input, 'selection_train', 1);
selection_retrain = get_default(input, 'selection_retrain', 1);
input.type_correction = get_default(input, 'type_correction', 'correct');
input.merge_train = get_default(input, 'merge_train', 'median');
input.merge_retrain = get_default(input, 'merge_retrain', 'median');
input.correct_limit = get_default(input, 'correct_limit', 1);
input.final_merge.type = get_default(input.final_merge, 'type', 'best');
input.correction_target = get_default(input, 'correction_target', 3);
correction_target = input.correction_target;
input.coun_ts_limit = get_default(input, 'coun_ts_limit', 1);
input.alpha_value = get_default(input, 'alpha_value', 0.05);

if ~isfield(input, 'grid_static') && ~isfield(input, 'grid_dynamic')
    input.grid_static = 20;
end

if isfield(input, 'grid_static')
    input.grid_dynamic.onset    = 1;
    input.grid_dynamic.LV_1.x   = struct('start', 1, 'end', 0, 'density', input.grid_static);
    input.grid_dynamic.LV_1.y   = struct('start', 1, 'end', 0, 'density', input.grid_static);
end

input.additional_NCV = get_default(input, 'additional_NCV', false);

end

function value = get_field(s, field)
    value = [];
    if isfield(s, field)
        try
            temp = load(s.(field));
            value = temp.(field);
            clear('temp');
        catch
            value = s.(field);
        end
    end
end

function value = get_default(s, field, default)
if isfield(s, field)
    value = s.(field);
else
    value = default;
end
end
