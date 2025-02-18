function cv_mbspls_outer_mean_bootstrap(analysisfolder, n_LV)

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

rest_boot = mod(input.bootstrap_testing, size_sets_bootstrap);
if rest_boot>0
    disp(['Please choose a number of bootstrap sets, which can be divided by ', num2str(setup.max_sim_jobs)]);
end
boot_sets = (input.bootstrap_testing - rest_boot)/size_sets_bootstrap;

% if isfield(input, 'procrustes_flag')
procrustes_flag = input.procrustes_flag;

save([bootstrap_folder '/bootstrap_setup.mat'],'selection_train', 'mbspls_params', 'correlation_method', 'matrix_norm', 'size_sets_bootstrap', 'correction_target',  'procrustes_flag','-v7.3');


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

for ob=1:OB
    for w=1:W
      

        w_opt = output.opt_parameters.(lv_name){w, matches(output.parameters_names, 'w')}; 

        for num_m=1:size(input.Xs,2)
            train_covariates{num_m} = Covars_ana{num_m}(output.CV.cv_outer_indices.TrainInd{ob,w_opt},:);
            %train_data_matrices{num_m} = matrices_ana{num_m}(output.CV.cv_outer_indices.TrainInd{ob,w_opt},:);
        end
        train_DiagNames = DiagNames_ana(output.CV.cv_outer_indices.TrainInd{ob,w_opt},:);

        c_weights_opt = output.final_parameters{n_LV,matches(output.parameters_names, 'c_weights')};
        Vs_opt = output.final_parameters{n_LV,matches(output.parameters_names, 'Vs_opt')};

        save([bootstrap_folder '/bootstrap_opt.mat'], 'c_weights_opt','Vs_opt','-v7.3');

        wrapper_train_opt = @(x) get_split_matrix(x, output.CV.cv_outer_indices.TrainInd{ob,w_opt}) ;
        train_data_matrices = cellfun(wrapper_train_opt, matrices_ana, UniformOutput=false);


        % save the optimized parameters and the permutation matrix
        save([bootstrap_folder '/bootstrap_partition_fold.mat'],...
            'train_data_matrices', 'train_covariates', 'train_DiagNames', 'cs_method', '-v7.3');


        [RHO_boot, weights_boot] = cv_mbspls_ICV_main_job_mult(bootstrap_folder, 'bootstrap', boot_sets, setup, num_matrices);


        lv_name = ['LV_', num2str(n_LV)];

        
        % Optimizsation criterion
        RHO_analysis = output.opt_parameters.(lv_name){w, matches(output.opt_parameters_names, 'RHO')};
        RHO_mean = mean(RHO_boot);
        RHO_SE = std(RHO_boot)/(sqrt(input.bootstrap_testing));
        ci_RHO = [RHO_mean - 1.96 * RHO_SE, RHO_mean + 1.96 * RHO_SE];
        bs_ratio_RHO = RHO_analysis/RHO_SE;
        output.bootstrap_results.(lv_name).ci_RHO = ci_RHO;
        output.bootstrap_results.(lv_name){1,w}{1,w}.bs_ratio_RHO = bs_ratio_RHO;
        output.bootstrap_results.(lv_name){1,w}.RHO_boot_size = size(RHO_boot);


        % weights
        for num_m=1:num_matrices
            weights_analysis{num_m} = output.final_parameters{n_LV, matches(output.opt_parameters_names, 'weights')}{num_m};

            weights_mean{num_m} = mean(weights_boot{num_m},2);
            weights_SE{num_m} = std(weights_boot{num_m},0,2)/(sqrt(input.bootstrap_testing));

            bs_ratio_weights{num_m} = weights_analysis{num_m}./weights_SE{num_m};
            bs_ratio_weights{num_m}(isnan(bs_ratio_weights{num_m})) = 0;
            bs_ratio_weights{num_m}(bs_ratio_weights{num_m} == Inf) = 0;
            bs_ratio_weights{num_m}(bs_ratio_weights{num_m} == -Inf) = 0; % CV added; correct?
            log_bs_weights{num_m} = abs(bs_ratio_weights{num_m})<=1.96;

            ci_weights{num_m} = [weights_mean{num_m} - 1.96 * weights_SE{num_m}, weights_mean{num_m} + 1.96 * weights_SE{num_m}];
            log_ci_weights{num_m} = ((sum(ci_weights{num_m}>0, 2) == 2) + (sum(ci_weights{num_m}<0, 2) == 2)) == 0;

            output.bootstrap_results.(lv_name){1,w}.ci_weights{num_m} = ci_weights{num_m};
            output.bootstrap_results.(lv_name){1,w}.bs_ratio_weights{num_m} = bs_ratio_weights{num_m};
            output.bootstrap_results.(lv_name){1,w}.weights_boot_size{num_m} = size(weights_boot{num_m});
            output.bootstrap_results.(lv_name){1,w}.log_bs_weights{num_m} = log_bs_weights{num_m};
            output.bootstrap_results.(lv_name){1,w}.sum_bs_weights{num_m} = sum(log_bs_weights{num_m});
            output.bootstrap_results.(lv_name){1,w}.log_ci_weights{num_m} = log_ci_weights{num_m};
            output.bootstrap_results.(lv_name){1,w}.sum_ci_weights{num_m} = sum(log_ci_weights{num_m});


        end

    end


end

datafile = [analysisfolder, 'final_results/preliminary_result_mean.mat'];
save(datafile, 'input', 'setup', 'output')

end