function cv_mbspls_outer_mean_perm(analysisfolder, n_LV)

resultfile = [analysisfolder, 'final_results/preliminary_result_mean.mat'];
load(resultfile)

[~, ~, matrices,  ~, ~, ~, ~, ~, ~, size_sets_bootstrap, correlation_method, cs_method, selection_train, selection_retrain, correction_target, matrix_norm, search_params, mbspls_params] = cv_mbspls_setup_parameters(input, setup);


input.inner_folds = size(input.CV.cv_inner_indices{1,1}.TrainInd,2);
input.outer_folds = size(input.CV.cv_outer_indices.TrainInd,2);
input.inner_permutations = size(input.CV.cv_inner_indices{1,1}.TrainInd,1);
input.outer_permutations = size(input.CV.cv_outer_indices.TrainInd,1);
K = input.inner_folds;
W = input.outer_folds;
IB = input.inner_permutations;
OB = input.outer_permutations;


% if isfield(input, 'procrustes_flag')
procrustes_flag = input.procrustes_flag;

weights_log = false;

weights{1,num_matrices} = [];
weights_Vs = [];

weights_V = [];
weights_RHO=[];


matrices_ana{num_matrices} = [];
matrices_val{num_matrices} = [];
for num_m=1:num_matrices
    matrix = matrices{num_m};
    matrices_val{num_m} = matrix(validation_set,:);
    matrices_ana{num_m} = matrix(analysis_set,:);
end


Diag_val = input.Diag(validation_set,:);
DiagNames_val = input.DiagNames(validation_set,:);
sites_val = input.sites(validation_set,:);


Diag_ana = input.Diag(analysis_set,:);
DiagNames_ana = input.DiagNames(analysis_set,:);
sites_ana = input.sites(analysis_set,:);

for num_m= 1:size(input.Xs,2)
    try
        Covars_ana{num_m} = input.covariates{num_m}(analysis_set,:);
        Covars_val{num_m} = input.covariates{num_m}(validation_set,:);

    catch

        Covars_ana{num_m} = nan(size(input.Diag(analysis_set,:),1),1);
        Covars_val{num_m} = nan(size(input.Diag(validation_set,:),1),1);
    end
end

if n_LV>input.correct_limit
    input.corrected_log(n_LV) = false;
elseif strcmp(input.type_correction, 'uncorrected')
    input.corrected_log(n_LV) = false;
else
    input.corrected_log(n_LV) = true;
end

if ~input.corrected_log(n_LV) % CV: what does this do?
    for num_m= 1:size(input.Xs,2)

        Covars_ana{num_m} = nan(size(input.Diag(analysis_set,:),1),1);
        Covars_val{num_m} = nan(size(input.Diag(validation_set,:),1),1);

    end

    correction_target = ones(size(input.Xs,2),1);

end


for num_m=1:num_matrices
                    test_covariates{num_m} = Covars_ana{num_m}(output.CV.cv_outer_indices.TestInd{ob,w},:);
                    train_covariates{num_m} = Covars_ana{num_m}(output.CV.cv_outer_indices.TrainInd{ob,w},:);
                end


                wrapper_test = @(x) get_split_matrix(x, output.CV.cv_outer_indices.TestInd{ob,w}) ;
                wrapper_train = @(x) get_split_matrix(x, output.CV.cv_outer_indices.TrainInd{ob,w}) ;
                train_data_matrices = cellfun(wrapper_train, matrices_ana, UniformOutput=false);

                test_data_matrices = cellfun(wrapper_test, matrices_ana, UniformOutput=false);

                train_Diag = Diag_ana(output.CV.cv_outer_indices.TrainInd{ob,w},:);

                % save the optimized parameters and the permutation matrix
                save([permutation_folder '/permutation_partition_fold.mat'],...
                    'train_data_matrices', 'test_data_matrices', 'train_covariates', 'test_covariates',...
                    'train_Diag', 'train_DiagNames', 'test_DiagNames', 'cs_method', '-v7.3');

                save([permutation_folder '/permutation_opt.mat'],...
                    'c_weights_opt','Vs_opt','-v7.3');

datafile = [analysisfolder, 'final_results/preliminary_result_mean.mat'];
save(datafile, 'input', 'setup', 'output')

end