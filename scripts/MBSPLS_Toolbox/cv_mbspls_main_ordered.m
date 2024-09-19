%% DP SPLS standalone function
function cv_mbspls_main(datafile)

%% initialize analysis folders
% load the datafile that contains the data and the analysis information
load(datafile);
%cv_get_analysis_info(input); % TO DO: print setup and input to log file

if ~exist('n_LV', 'var') % i.e., analysis is fresh, and not one to be continued and same directory should be used
    nn=1;
    analysis_name = input.name;
    analysis_folder = [setup.analysis_folder  '/' analysis_name];
    if ~exist(analysis_folder,'dir')
        mkdir(analysis_folder);
    else % CV: could this be done more elegantly/ is it necessary? in some cases it may be better to overwrite folder
        while nn < 100
            analysis_name = [input.name '_' num2str(nn)];
            analysis_folder = [setup.analysis_folder  '/' analysis_name];
            if ~exist(analysis_folder,'dir')
                mkdir(analysis_folder)
                nn=100;
            else
                nn=nn+1;
            end
            %setup.analysis_name = analysis_name;
            %setup.analysis_folder = analysis_folder;
        end
    end
    [permutation_folder, hyperopt_folder, bootstrap_folder, detailed_results, final_results] = dp_create_folders(setup.scratch_space, analysis_folder, analysis_name);

else
%     analysis_folder = setup.analysis_folder;
%     analysis_name = setup.analysis_name;
%     permutation_folder = aux.permutation_folder;
%     hyperopt_folder = aux.hyperopt_folder;
%     bootstrap_folder = aux.boostrap_folder;
%     n_LV = aux.n_LV;
%     count_ns = aux.count_ns;
%     matrices_ana = aux.matrices_ana;
%     matrices_val = aux.matrices_val;

end

%% Prepare analysis folders
% sets the subfolders for the different analysis steps
% temporary files for hyperopt and permutation will be saved to the common scratch space "/volume/mitnvp1_scratch"


cd(final_results);

%% Set parameters for SPLS analysis
[input, setup, matrices,  B, K, W, OB, IB, size_sets_permutation, size_sets_bootstrap, correlation_method, cs_method, selection_train, selection_retrain, correction_target, matrix_norm, search_params, mbspls_params] = cv_mbspls_setup_parameters(input, setup);

% define column names for matrices so that you can access them later by
% indexing

% depending on how many matrices are being used, the number of columns
% differs

num_matrices = length(matrices);
% arbitrary amount of matrices
%weights_colnames = cellstr("c_" + (1:num_matrices));
% hyp_colnames = cellstr("hyp_" + (1:num_matrices));
%output.parameters_names = {'w', weights_colnames{:}, hyp_colnames{:}, 'success', 'RHO', 'p', 'epsilon', 'omega', 'V_opt'};
output.parameters_names = {'w', 'c_weights', 'weights', 'success', 'RHO', 'p', 'lVs', 'Vs_opt'};


% indices for later use
opt_weights = strcmp(output.parameters_names,'weights');
opt_p = strcmp(output.parameters_names,'p');
opt_c_weights = strcmp(output.parameters_names,'c_weights');
opt_RHO = strcmp(output.parameters_names,'RHO');

output.opt_parameters_names = output.parameters_names;

% CV: what is V_opt?


% preallocate placeholder matrices for higher speed
% output.final_parameters = num2cell(nan(10, numel(output.parameters_names))); % cell array to store final parameters for each LV iteration


%% Set up analysis and validation sets, including covariates

% n_LV counts the iteration through the LV loops, so that the results for each
% LV can be stored in separate rows
if ~exist('n_LV', 'var')
    n_LV = 1;
end

% check if a validation set is set up => extract it before starting the
% analysis
val_labels=[];
if isfield(input, 'validation_set')
    if ~islogical(input.validation_set)
        validation_folds = round(100/input.validation_set);
        switch input.val_stratification % take out validation set depending on 1) diagnosis, 2) sites, 3) both
            case 1
                val_labels = input.Diag;
            case 2
                for i=1:size(input.sites,1)
                    val_labels(i,1) = find(input.sites(i,:));
                end
            case 3
                for i=1:size(input.sites,1)
                    val_labels(i,1) = input.Diag(i,1)*10 + find(input.sites(i,:));
                end
        end
        output.validation_indices = nk_CVpartition2(1, validation_folds, val_labels);
        analysis_set = sort(output.validation_indices.TrainInd{1,1});
        validation_set = sort(output.validation_indices.TestInd{1,1});
    else
        analysis_set = ones(size(input.Diag,1),1)>0;
        validation_set = ~analysis_set;
    end

else
    analysis_set = ones(size(input.Diag,1),1)>0;
    validation_set = ~analysis_set;
end
%%
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

if input.framework == 4
    input.framework = 1;
    input.validation_set = 50;
end

IN.type = input.framework;

if IN.type ~= 3
    IN.labels = Diag_ana;
else
    IN.labels = sites_ana;
    W = size(sites_ana,2);
    if input.additional_NCV
        IN.sublabels = Diag_ana;
    else
        IN.sublabels = [];
        K = W-1;
    end
end

IN.OB = OB; IN.IB = IB; IN.OF = W; IN.IF = K;

if ~isfield(input, 'CV')
    output.CV = dp_setup_framework(IN);
else 
    output.CV = input.CV; 
end

if input.save_CV
    CV = output.CV; 
    save(fullfile(analysis_folder, 'CVstruct.mat'), 'CV');
    clear CV;
end
%save permutation setup files; CV:this causes a bug if max_sim_jobs = 1
% 
if setup.max_sim_jobs > 1
    temp_pc = repelem(1:W, ceil(setup.num_jobs_perm/W)); %repelem(1:W, ceil(setup.max_sim_jobs/W));
    temp_pc(randi(setup.num_jobs_perm,1,size(temp_pc,2)-setup.num_jobs_perm))=[]; % temp_pc(randi(setup.max_sim_jobs,1,size(temp_pc,2)-setup.max_sim_jobs))=[];
    perm_coding = [[1:setup.num_jobs_perm;temp_pc]; 1:(B/setup.num_jobs_perm):(B-1)]; % perm_coding = [[1:setup.max_sim_jobs;temp_pc]; 1:(B/setup.max_sim_jobs):(B-1)];
    rest_perm = mod(B,size_sets_permutation);
    if rest_perm>0
        disp(['Please choose a number of permutation sets, which can be divided by ', num2str(setup.max_sim_jobs)]);
    end
    perm_sets = (B - rest_perm)/size_sets_permutation;
else
    perm_coding =  [1;1;1];
    perm_sets = 1; % also if 'sequential'
end




rest_boot = mod(input.bootstrap_testing, size_sets_bootstrap);
if rest_boot>0
    disp(['Please choose a number of bootstrap sets, which can be divided by ', num2str(setup.max_sim_jobs)]);
end
boot_sets = (input.bootstrap_testing - rest_boot)/size_sets_bootstrap;

% if isfield(input, 'procrustes_flag')
procrustes_flag = input.procrustes_flag; 
% else
%     procrustes_flag = true;
% end
save([permutation_folder '/permutation_setup.mat'],'selection_train', 'mbspls_params','correlation_method', 'matrix_norm', 'size_sets_permutation', 'perm_coding', 'selection_retrain', 'correction_target', 'procrustes_flag', '-v7.3');
save([bootstrap_folder '/bootstrap_setup.mat'],'selection_train', 'mbspls_params', 'correlation_method', 'matrix_norm', 'size_sets_bootstrap', 'correction_target',  'procrustes_flag','-v7.3');

% set count for not significant LVs to 0, as soon as one LV is not
% significant, the LV loop stops
count_ns = 0;

if strcmp(input.merge_retrain, 'weighted_mean')
    weights_log = true;
else
    weights_log = false;

    weights{1,num_matrices} = [];
    weights_Vs = [];

    weights_V = [];
    weights_RHO=[];
end

%% Latent Variable loop begins: setup of grid and cu/cv combination for analysis

% here starts the outer loop for each single LV, the next LV is only
% computed if the previous LV was significant (with FDR correction)


% CV 04.08.2024: should previous mbspls analysis be continued? In this case, input,
% output, preliminary results should be loaded instead of datafile; or what
% other fields do we need? What other fields change in loop? 
% count_ns, matrices_ana, matrices_val
if count_ns > 0
    
    if ~islogical(input.validation_set)


    else


    end
end


while count_ns<input.coun_ts_limit && (input.max_n_LVs == -1 || n_LV < input.max_n_LVs) % if max_n_LV = -1, no maximum set
    % repeats the entire process W times to obtain W different final p values,
    % which go into the omnibus hypothesis

    tic % track time elapsed for one LV 
    
    % matrix to store the optimal parameters of the w-th loop
    output.opt_parameters.(['LV_', num2str(n_LV)]) = num2cell(nan(W,numel(output.parameters_names)));


    if (strcmp(input.optimization_strategy, 'grid_search') && any(n_LV==input.grid_dynamic.onset)) || strcmp(input.optimization_strategy, 'randomized_search')


        switch input.optimization_strategy
            case 'grid_search'
                %grid = input.grid_dynamic.LVs;
                temp_density = 1;
                for num_m=1:num_matrices
                    temp_density = search_params(num_m).density*temp_density;
                end
                input.size_sets_hyperopt = round(temp_density/setup.num_jobs_hyperopt); %round(temp_density/setup.max_sim_jobs);
            case 'randomized_search'
                input.size_sets_hyperopt = round(input.randomized_search_params.randomized_search_iterations/setup.num_jobs_hyperopt); % round(input.randomized_search_params.randomized_search_iterations/setup.max_sim_jobs);

        end


        size_sets_hyperopt = input.size_sets_hyperopt;
        if size_sets_hyperopt==0
            size_sets_hyperopt=1;
        end

        
        [weights_combination, hyperopt_sets] = cv_weights_ext(input.optimization_strategy, matrices_ana, search_params, size_sets_hyperopt);
        
        if strcmp(input.optimization_strategy, 'randomized_search')
            input.randomized_search_params.weights_combination = weights_combination; 
        end
%         if strcmp(setup.cluster, 'sequential')
%             hyperopt_sets = 1; 
%         end
        weights_ICV_collection = cell(size(weights_combination,1),1);

    end

    RHO_ICV_collection=[];
    RHO_opt_collection=[];



    lVs_opt_collection{num_matrices}=[];
    Vs_opt_collection=[];
    weights{num_matrices} = []; weights_Vs=[];



    weights_RHO=[];


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

        correction_target = ones(size(input.matrices,2),1);

    end

    %% Outer folds begin: hyperopt

    for w=1:W

        % remove hp% of the data randomly and keep it in a hold-out
        % dataset, the rest is called the keep_in data for training/testing
        % the algorithm
        RHO_ICV_collection_temp=[];


        ICV_collection_temp{num_matrices}={};

        for ob=1:OB

            wrapper_test = @(x) get_split_matrix(x, output.CV.cv_outer_indices.TestInd{ob,w}) ;
            wrapper_train = @(x) get_split_matrix(x, output.CV.cv_outer_indices.TrainInd{ob,w}) ;

            test_data_matrices = cellfun(wrapper_test, matrices_ana, UniformOutput=false);
            train_data_matrices = cellfun(wrapper_train, matrices_ana, UniformOutput=false);


            test_Diag = Diag_ana(output.CV.cv_outer_indices.TestInd{ob,w},:);
            test_DiagNames = DiagNames_ana(output.CV.cv_outer_indices.TestInd{ob,w},:);
            for num_m=1:size(input.Xs,2)
                test_covariates{num_m} = Covars_ana{num_m}(output.CV.cv_outer_indices.TestInd{ob,w},:);
            end

            train_Diag = Diag_ana(output.CV.cv_outer_indices.TrainInd{ob,w},:);
            train_DiagNames = DiagNames_ana(output.CV.cv_outer_indices.TrainInd{ob,w},:);
            for num_m=1:size(input.Xs,2)
                train_covariates{num_m} = Covars_ana{num_m}(output.CV.cv_outer_indices.TrainInd{ob,w},:);
            end
            cv_inner_TrainInd = output.CV.cv_inner_indices{ob,w}.TrainInd;
            cv_inner_TestInd = output.CV.cv_inner_indices{ob,w}.TestInd;
            cv_outer = output.CV.cv_outer_indices;


            save([hyperopt_folder '/hyperopt_partition.mat'],...
                'train_data_matrices', 'train_covariates', 'cs_method',...
                'train_DiagNames', 'test_DiagNames',...
                'cv_inner_TrainInd', 'cv_inner_TestInd','mbspls_params', 'correlation_method', 'matrix_norm', ...
                'weights_combination', 'size_sets_hyperopt', 'correction_target', '-v7.3');

            if selection_train == 2 && ob == OB

                save([permutation_folder '/permutation_partition_fold_', num2str(w), '.mat'],...
                    'train_data_matrices', 'train_covariates', 'test_covariates',...
                    "test_covariates", 'test_data_matrices', 'train_Diag', 'train_DiagNames', 'test_DiagNames', 'cs_method', '-v7.3');
            end


            %% Hyperparameter optimization starts
            % this performs training and testing within the inner training
            % and testing folds of one OCV => according to W, K
            RHO_ICV_collection_temp = cv_mbspls_ICV_main_job_mult(hyperopt_folder, 'hyperopt', hyperopt_sets, setup, num_matrices);

            RHO_ICV_collection = cat(2, RHO_ICV_collection, RHO_ICV_collection_temp);
            %end

        end

        switch selection_train % 1=within one OCV fold, 2=across all OCV folds
            case 1 % within one OCV fold


                save([hyperopt_folder, '/checkpoint_hyperopt.mat'], 'RHO_ICV_collection', 'weights_combination');



                switch input.merge_train
                    case 'mean'
                        f_hyperopt = @(x)mean(x,2);
                    case 'median'
                        f_hyperopt = @(x)median(x,2);
                    case 'best'
                        f_hyperopt = @(x)max(x,2);
                end


                for icv=1:size(weights_combination,1)
                    RHO_selection(icv,1) = f_hyperopt(RHO_ICV_collection(icv,:));
                end


                [~,I_max] = max(RHO_selection);

                if sum(I_max)>1

                end
                c_weights_opt = weights_combination(I_max,:);


                % retrain the inner folds of this one outer fold with
                % cu_max and cv_max => options: along all inner folds
                % or with all inner folds combined
                % optqions: 1) retrain on single test splits within
                % folds, then merge, 2) retrain on all inner folds
                % separately, then merge with mean or median, 3)
                % retrain on entirety of inner folds, 4) use already
                % existing u and v from inner folds without retraining

                for ob=1:OB


                    for num_m=1:num_matrices
                        IN_matrices.train{num_m} = matrices_ana{num_m}(output.CV.cv_outer_indices.TrainInd{ob,w},:);
                        IN_matrices.test{num_m} = matrices_ana{num_m}(output.CV.cv_outer_indices.TestInd{ob,w},:);
                        COV{num_m}.test = Covars_ana{num_m}(output.CV.cv_outer_indices.TestInd{ob,w},:);
                        COV{num_m}.train = Covars_ana{num_m}(output.CV.cv_outer_indices.TrainInd{ob,w},:);
                    end

                    DiagNames_train = DiagNames_ana(output.CV.cv_outer_indices.TrainInd{ob,w},:);
                    DiagNames_test = DiagNames_ana(output.CV.cv_outer_indices.TestInd{ob,w},:);

                    for num_m=1:num_matrices
                        if ~isempty(cs_method{num_m}.correction_subgroup)
                            cs_method{num_m}.subgroup_train = contains(DiagNames_train, cs_method{num_m}.correction_subgroup);
                            cs_method{num_m}.subgroup_test = contains(DiagNames_test, cs_method{num_m}.correction_subgroup);
                        else
                            cs_method{num_m}.subgroup_train = [];
                            cs_method{num_m}.subgroup_test = [];
                        end
                    end


                    [OUT_matrices] = cv_master_correctscale(IN_matrices, COV, cs_method, correction_target);


                    switch selection_retrain

                        case 1 % retrain on entirety of folds

                            [RHO_opt_collection(1,ob), weights_opt_collection(:,ob), Vs_opt_collection_temp, lVs_opt_temp] = cv_mbspls_wrapper(OUT_matrices.train, OUT_matrices.test, c_weights_opt, mbspls_params, correlation_method, matrix_norm);
                            Vs_opt_collection = cat(3, Vs_opt_collection, Vs_opt_collection_temp);



                            if weights_log
                                for num_m=1:num_matrices
                                    weights{num_m} = RHO_opt_collection(1,ob)*ones(size(weights_opt_collection(:,ob),1),1)';
                                end
                                weights_Vs_1 = RHO_opt_collection(1,ob)*ones(size(Vs_opt_collection_temp));
                                weights_Vs = cat(3, weights_Vs, weights_Vs_1);
                            end
                            save([hyperopt_folder, '/weights_check.mat'], 'weights', 'weights_Vs');
                            clear V_opt_collection_temp;

                            for num_m=1:num_matrices
                                lVs_opt_collection{num_m} = cat(1, lVs_opt_collection{num_m}, lVs_opt_temp(:,num_m));
                            end

                        case 2 % use already existing u and v from inner folds without retraining

                            weights_collection = weight_selection(I_max,:)';
                            Vs_opt_collection = false;
                            test_data = cellfun(@(x) (x.test), OUT_matrices);
                            [RHO_opt_collection(1,ob), lVs_opt_temp, weights_opt_collection] = cv_mbspls_projection(test_data, weights_opt_collection, correlation_method, matrix_norm);

                            for num_m=1:num_matrices
                                lVs_opt_collection{num_m} = cat(1, lVs_opt_collection{num_m}, lVs_opt_temp{num_m});
                            end



                    end

                end

                switch input.merge_retrain
                    case 'best'
                        [RHO_opt,I_opt] = max(RHO_opt_collection);


                        for num_m=1:num_matrices
                            weights_opt{num_m} = weights_opt_collection{num_m}(:, I_opt);
                            for num_m2=num_m+1:num_matrices
                                Vs_opt{num_m, num_m2} = Vs_opt_collection{num_m, num_m2}(:, :, I_opt);
                            end
                        end


                    otherwise
                        RHO_opt = f_hyperopt(RHO_opt_collection);

                        Vs_opt = dp_trainmerge_single(Vs_opt_collection, input.merge_retrain, 3, weights_Vs);
                        for num_m= 1:num_matrices
                            weights_opt{num_m} = dp_trainmerge_single(weights_opt_collection{num_m}, input.merge_retrain, 2, weights{num_m});
                            for num_m=1:num_matrices
                                for num_m2=num_m+1:num_matrices
                                    Vs_opt{num_m, num_m2} = dp_trainmerge_single(Vs_opt_collection{num_m, num_m2}, input.merge_retrain, 3, weights_Vs);
                                end
                            end
                        end
                        save([hyperopt_folder '/opt_vars_check.mat'],...
                            'RHO_opt_collection','weights_opt_collection', 'RHO_opt', 'weights_opt', 'Vs_opt', 'RHO_ICV_collection', '-v7.3');


                end


                %% Permutation Testing within fold

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

                % choose whether permutation testing shall be conducted right
                % within the current fold or whether permutation testing shall
                % be done across all outer folds

                RHO_b_collection = cv_mbspls_ICV_main_job_mult(permutation_folder, 'permutation', perm_sets, setup,num_matrices);

                % calculate the p value, how the optimized model performs
                % agains the permuted model: 2 options => counting or AUC
                % testing
                switch input.statistical_testing
                    case 1 % counting
                        RHO_count_b = sum(RHO_b_collection > RHO_opt);
                        p = (RHO_count_b+1)/(B+1);
                    case 2 % AUC testing
                        IN.dist = RHO_b_collection;
                        IN.val = RHO_opt;
                        IN.testing_precision = input.permutation_testing_precision;
                        [p, ~] = dp_auc_testing(IN);
                end

                % store all parameters of the w-th wide loop into a matrix
                success_ICV=1;

                output.opt_parameters.(['LV_', num2str(n_LV)])(w,:) = {w c_weights_opt weights_opt success_ICV RHO_opt p lVs_opt_collection Vs_opt};

                weights_ICV_collection = cell(size(weights_combination,1),1);
                Vs_opt_collection=[];
                weights{1,num_matrices} = []; weights_Vs=[];

                save([detailed_results '/preliminary_results.mat'],'input', 'setup', 'output');

                % clear collections for next W fold
                RHO_ICV_collection=[];


                RHO_opt_collection=[];


                lVs_opt_collection{num_matrices}=[];


                RHO_selection = [];

        end
    end

    %% Hyperopt selection across all OCV folds

    if selection_train == 2 % select optimal hyperparameter combination across all folds and all permutations

        %         input.final_merge.p_value_adjust = 'off';
        input.final_merge.significant_only = 'off';
        input.final_merge.majority_vote = 'off';
        input.final_merge.type = 'best';

        switch input.merge_train % defines how the performance values (RHO) of all inner folds will be collected
            case 'mean'
                f_hyperopt = @(x)mean(x,2);
            case 'median'
                f_hyperopt = @(x)median(x,2);
            case 'best'
                f_hyperopt = @(x)max(x,2);
        end


        for icv=1:size(weights_combination,1)
            RHO_selection(icv,1) = f_hyperopt(RHO_ICV_collection(icv,:));
        end



        [~,I_max] = max(RHO_selection);

        c_weights_opt = weights_combination(I_max,:);

        Vs_opt_collection=[];


        % retrain the inner folds of this one outer fold with
        % cu_max and cv_max => options: along all inner folds
        % or with all inner folds combined
        cv_outer = output.CV.cv_outer_indices;


        nn=1;
        for w=1:W
            for ob=1:OB

                for num_m=1:size(input.Xs,2)
                    COV{num_m}.test = Covars_ana{num_m}(output.CV.cv_outer_indices.TestInd{ob,w},:);
                    COV{num_m}.train = Covars_ana{num_m}(output.CV.cv_outer_indices.TrainInd{ob,w},:);
                end

                IN_matrices.train = matrices_ana(output.CV.cv_outer_indices.TrainInd{ob,w},:);
                IN_matrices.test = matrices_ana(output.CV.cv_outer_indices.TestInd{ob,w},:);

                [OUT_matrices] = cv_master_correctscale(IN_matrices, COV, cs_method, correction_target);



                switch selection_retrain
                    case 1 % retrain on entirety of outer folds

                        [RHO_opt_collection(1,nn), weights_opt_collection(:,nn), Vs_opt_collection_temp, lVs_opt_temp] = cv_mbspls_wrapper(OUT_matrices.train, OUT_matrices.test, c_weights_opt, correlation_method, matrix_norm);
                        Vs_opt_collection = cat(3, Vs_opt_collection, Vs_opt_collection_temp);



                        if weights_log
                            if input.framework == 3
                                ratio = size(output.CV.cv_outer_indices.TestInd{1,w},1)/size(cat(1,output.CV.cv_outer_indices.TestInd{1,:}),1);
                            else
                                ratio = RHO_opt_collection(1,nn);
                            end
                            weights_RHO(1,nn) = ratio;

                            for num_m=1:num_matrices
                                weights{num_m}(:,nn) = ratio*ones(size(weights_opt_collection(:,nn),1),1)';
                            end
                            weights_V_1 = ratio*ones(size(Vs_opt_collection_temp));
                            weights_Vs = cat(3, weights_Vs, weights_Vs_1);

                        end

                        clear V_opt_collection_temp;
                        clear Vs_opt_collection_temp;


                        for num_m=1:num_matrices
                            lVs_opt_collection{num_m} = cat(1, lVs_opt_collection{num_m}, lVs_opt_temp{num_m});
                        end


                    case 2 % do not retrain, use already existing u and v



                        weights_opt_collection = weights_selection(I_max,:)';
                        Vs_opt_collection = false;

                        [RHO_opt_collection(1,nn), lVs_opt_temp, weights_opt_collection] = cv_mbspls_projection(OUT_matrices.test, weights_opt_collection, correlation_method, matrix_norm);


                        for num_m=1:num_matrices
                            lVs_opt_collection{num_m} = cat(1, lVs_opt_collection{num_m}, lVs_opt_temp{num_m});
                        end





                end
                nn=nn+1;
            end
        end

        if weights_log


            save([hyperopt_folder, '/weights_check.mat'], 'ratio', 'weights_RHO', 'weights', 'weights_Vs');

        end

        switch input.merge_retrain % defines how the results on the outer folds after retraining on entire inner folds will be condensed
            case 'best'
                [RHO_opt,I_opt] = max(RHO_opt_collection);


                weights_opt = weights_opt_collection(:, I_opt);
                Vs_opt = Vs_opt_collection(:, :, I_opt);


            otherwise
                RHO_opt = dp_trainmerge_single(RHO_opt_collection, input.merge_retrain, 2, weights_RHO);


                weights_opt = dp_trainmerge_single(u_opt_collection, input.merge_retrain, 2, weights);
                Vs_opt = Vs_opt_collection(:, :, I_opt);

        end


        save([detailed_results '/opt_vars_LV', num2str(n_LV), '.mat'],...
            'RHO_opt_collection','weights_opt_collection', ...
            'RHO_opt', 'weights_opt', 'RHO_ICV_collection');

        % save the optimized parameters and the permutation matrix
        save([permutation_folder '/permutation_opt.mat'], 'c_weights_opt', 'Vs_opt','-v7.3');

        RHO_b_collection = cv_mbspls_ICV_main_job_mult(permutation_folder, 'permutation', perm_sets, num_matrices);





        % calculate the p value to test whether the optimized model is significantly
        % better than the permutation model
        switch input.statistical_testing
            case 1 % counting
                RHO_count_b = sum(RHO_b_collection > RHO_opt);
                p = (RHO_count_b+1)/(B+1);
            case 2 % AUC testing
                IN.dist = RHO_b_collection;
                IN.val = RHO_opt;
                IN.testing_precision = input.permutation_testing_precision;
                [p, ~] = dp_auc_testing(IN);
        end
        success_ICV=true;

        % store all parameters of the w-th wide loop into a matrix


        output.opt_parameters.(['LV_', num2str(n_LV)]) = {1 c_weights_opt weights_opt success_ICV RHO_opt p lVs_opt_collection, Vs_opt};
        save([detailed_results '/preliminary_results.mat'],'input', 'setup', 'output');

        % clear collections for next W fold
        RHO_ICV_collection=[]; weights_ICV_collection=cell(size(weights_combination,1),1);
        RHO_opt_collection=[]; lVs_opt_collection{num_matrices}=[]; Vs_opt_collection=[];
        RHO_selection = []; weights_selection = [];
        weights{num_matrices}=[]; weights_Vs=[]; weights_RHO=[];

    end

    %% Final merging of folds (if applicable)

    % calculate FDR-corrected p value
    %     switch input.final_merge.mult_test
    %         case 'FDR'
    %             output.pvalue_FDR(n_LV,1) = dp_FDR([output.opt_parameters.(['LV_', num2str(n_LV)]){:,opt_p}], input.alpha_value);
    %         case 'SFisher'
    output.pvalue_FDR(n_LV,1) = 0.05;
    %     end

    % use different options to find optimal solution,
    % 1) collapse all folds using mean/median with/without majority vote
    % 2) collapse only significant folds using mean/median with/without majority vote
    % 3) take only best fold

    % set up function for collapsing/selecting folds => mean, median,
    % weighted mean, best(which automatically ends all further selection)
    temp  = output.opt_parameters.(['LV_', num2str(n_LV)]);
    switch input.final_merge.type
        case 'best'
            input.final_merge.significant_only = 'off';
            input.final_merge.majority_vote = 'off';
            log_p_min = [temp{:,opt_p}]==min([temp{:,opt_p}]);
            temp = temp(log_p_min,:);
            if sum(log_p_min)>1
                log_RHO_max = [temp{:,opt_RHO}]==max([temp{:,opt_RHO}]);
                % if still several solutions are equal, choose the sparsest
                temp = temp(log_RHO_max,:);
                if sum(log_RHO_max)>1
                    % Define a function that converts a cell array to a numeric array and computes the mean
                    %computeMean = @(x) mean(x(:));
                    %log_min_sparsity = arrayfun(computeMean, [temp{:,opt_weights}])==min(arrayfun(computeMean, [temp{:,opt_weights}]));
                    %temp = temp(log_min_sparsity,:);
                    %if sum(log_min_sparsity)>1  % then take first one as all are same model / same sparsities
                        temp = temp(1,:);
                    %end
                end
            end
            output.final_parameters(n_LV,:)=temp;
            p_threshold = output.pvalue_FDR(n_LV,1);
        otherwise
            p_threshold = 0.05;
    end

    % generate log for significant folds
    switch input.final_merge.significant_only
        case 'on'
            log_sig = [temp{:,opt_p}] <= output.pvalue_FDR(n_LV,1);
            if sum(log_sig)==0
                log_sig = [temp{:,opt_p}] <= 2;
                count_ns = count_ns+1;
            end
        case 'off'
            log_sig = [temp{:,opt_p}] <= 2;
    end

    temp = temp(log_sig,:);

    % generate logs for u and v vectors for majority voting on sparsity of
    % each single value
    switch input.final_merge.majority_vote
        case 'on'
            for i=1:size(temp,1)



                for num_m=1:num_matrices
                    weights_temp = temp{i, opt_weights};
                    log_weights{num_m, :} = weights_temp{num_m}~=0;
                    log_fm_weights{num_m} = sum(log_weights{num_m}, 1)./size(log_weights{num_m},1)<0.5;
                end

            end

        case 'off'



            log_fm_weights{num_m}= temp{1,opt_weights}{num_m}>100;


    end

    for i=1:size(temp,1)

        for num_m=1:num_matrices
            temp{i, opt_weights}{num_m}(log_fm_weights{num_m})= 0;
        end

    end

    % merge the folds according to the previous selection
    switch input.final_merge.type
        case 'mean'
            for i=1:size(temp,2)
                try output.final_parameters{n_LV,i}=mean([temp{:,i}],2);
                catch
                    temp_cat = [];
                    for ii=1:size(temp,1)
                        temp_cat=cat(1, temp_cat, temp{ii,i});
                    end
                    output.final_parameters{n_LV,i}=temp_cat;
                end
            end
        case 'median'
            for i=1:size(temp,2)
                try output.final_parameters{n_LV,i}=median([temp{:,i}],2);
                catch
                    temp_cat = [];
                    for ii=1:size(temp,1)
                        temp_cat=cat(1, temp_cat, temp{ii,i});
                    end
                    output.final_parameters{n_LV,i}=temp_cat;
                end
            end
        case 'weighted_mean'
            weights = 1-[temp{:,opt_p}];
            if size(temp,1)==1
                output.final_parameters(n_LV,:) = temp;
            else
                for i=1:size(temp,2)
                    try temp_wm = [temp{:,i}];
                    catch
                        temp_wm = []; temp_wm_cat=[];
                        for ii=1:size(temp,1)
                            temp_wm_cat=cat(1, temp_wm_cat, temp{ii,i});
                        end
                    end
                    if isempty(temp_wm)
                        output.final_parameters{n_LV,i}=temp_wm_cat;
                    elseif isvector(temp_wm)
                        output.final_parameters{n_LV,i}=wmean(temp_wm, weights);
                    elseif ismatrix(temp_wm)
                        temp_weights=[];
                        for iii=1:size(weights,2)
                            temp_weights(:,iii) = ones(size(temp_wm,1),1)*weights(iii);
                        end
                        output.final_parameters{n_LV,i}=wmean(temp_wm, temp_weights,2);
                    end
                end
                output.final_parameters{n_LV,1}=n_LV;
            end
    end

    [c_pvalues, ~, ~] = dp_multiple_testing(input.final_merge.mult_test, [output.opt_parameters.(['LV_', num2str(n_LV)]){:,opt_p}], 0.05, false);

    for ppp=1:size(output.opt_parameters.(['LV_', num2str(n_LV)]),1)
        output.opt_parameters.(['LV_', num2str(n_LV)]){ppp,opt_p} = c_pvalues(ppp);
    end

    p_LV = min(c_pvalues);

    output.final_parameters{n_LV,opt_p} = p_LV;

    %% Bootstrapping for chosen LV weight vectors
    w_opt = output.final_parameters{n_LV,matches(output.parameters_names, 'w')}; %CV: what is w

    for num_m=1:size(input.Xs,2)
        train_covariates{num_m} = Covars_ana{num_m}(output.CV.cv_outer_indices.TrainInd{ob,w_opt},:);
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
    
    % CV 20.08.2024: check bootstrap stability (delete later)
    save([detailed_results '/bootstrap_results_' lv_name '.mat'], 'RHO_boot', 'weights_boot');

    % Optimizsation criterion
    RHO_analysis = output.final_parameters{n_LV, matches(output.opt_parameters_names, 'RHO')};
    RHO_mean = mean(RHO_boot);
    RHO_SE = std(RHO_boot)/(sqrt(input.bootstrap_testing));
    ci_RHO = [RHO_mean - 1.96 * RHO_SE, RHO_mean + 1.96 * RHO_SE];
    bs_ratio_RHO = RHO_analysis/RHO_SE;
    output.bootstrap_results.(lv_name).ci_RHO = ci_RHO;
    output.bootstrap_results.(lv_name).bs_ratio_RHO = bs_ratio_RHO;
    output.bootstrap_results.(lv_name).RHO_boot_size = size(RHO_boot);





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

        output.bootstrap_results.(lv_name).ci_weights{num_m} = ci_weights{num_m};
        output.bootstrap_results.(lv_name).bs_ratio_weights{num_m} = bs_ratio_weights{num_m};
        output.bootstrap_results.(lv_name).weights_boot_size{num_m} = size(weights_boot{num_m});
        output.bootstrap_results.(lv_name).log_bs_weights{num_m} = log_bs_weights{num_m};
        output.bootstrap_results.(lv_name).sum_bs_weights{num_m} = sum(log_bs_weights{num_m});
        output.bootstrap_results.(lv_name).log_ci_weights{num_m} = log_ci_weights{num_m};
        output.bootstrap_results.(lv_name).sum_ci_weights{num_m} = sum(log_ci_weights{num_m});

    end

    %% check if LV is significant and apply it to validation set (if applicable)

    if ~islogical(input.validation_set)

        % apply the final parameters to the validations sets

        IN_matrices.train = matrices_ana;
        IN_matrices.test = matrices_val;

        COV.train = Covars_ana;
        COV.test = Covars_val;

        cs_method_val = cs_method;
        if ~isempty(cs_method_val.correction_subgroup)
            cs_method_val.subgroup_train = contains(DiagNames_ana, cs_method_val.correction_subgroup);
            cs_method_val.subgroup_test = contains(DiagNames_val, cs_method_val.correction_subgroup);
        else
            cs_method_val.subgroup_train = [];
            cs_method_val.subgroup_test = [];
        end

        temp = output.final_parameters(n_LV,:);


        [OUT_matrices] = cv_master_correctscale(IN_matrices, COV, cs_method_val, correction_target);
        for num_m=1:num_matrices
            weights_opt{num_m} = temp{1, opt_weights{num_m}};
            c_weights_opt{num_m} = temp{1, opt_c_weights{num_m}};
        end


        switch input.validation_train
            case 1 % retrain

                [RHO_opt, weights_opt, ~, lVs_val] = cv_mbspls_wrapper(OUT_matrices.train,OUT_matrices.test, c_weights_opt, correlation_method, matrix_norm);

                train_data_matrices = OUT_matrices.train;
                test_data_matrices = OU_matrices.test;


                % save the optimized parameters and the permutation matrix
                test_covariates = COV.test;
                train_covariates = COV.train;

                train_Diag = Diag_ana;

                selection_train=input.selection_train; %1

               % save([permutation_folder '/permutation_setup.mat'],'selection_train', 'mbspls_params','correlation_method', 'matrix_norm', 'size_sets_permutation', 'perm_coding', 'correction_target', 'procrustes_flag', '-v7.3');
               % save([bootstrap_folder '/bootstrap_setup.mat'],'selection_train', 'mbspls_params','correlation_method', 'matrix_norm', 'size_sets_bootstrap', 'selection_retrain', 'correction_target', 'bs_method', 'procrustes_flag', '-v7.3');

                save([permutation_folder '/permutation_partition_fold.mat'],...
                    'train_data_matrices', 'test_data_matrices', 'train_covariates', 'test_covariates',...
                    'train_Diag', 'cs_method_val', '-v7.3');

                RHO_b_collection = cv_mbspls_ICV_main_job_mult(permutation_folder, 'permutation', perm_sets, num_matrices);


            case 2 % apply optimal u and v


                [RHO_opt, lVs_val, weights_opt] = cv_mbspls_projection(OUT_matrices.test, weights_opt, correlation_method, matrix_norm);

                for pp=1:size(weights_b_collection,1)
                    [RHO_b_collection(pp,1), ~, ~, ~, ~] = cv_mbspls_projection(OUT_matrices.test, weights_b_collection(pp,:)', correlation_method, matrix_norm);
                end




        end

        switch input.statistical_testing
            case 1 % counting
                RHO_count_b = sum(RHO_b_collection > RHO_opt);
                p = (RHO_count_b+1)/(B+1);
            case 2 % AUC testing
                IN.dist = RHO_b_collection;
                IN.val = RHO_opt;
                IN.testing_precision = input.permutation_testing_precision;
                [p, ~] = dp_auc_testing(IN);
            case 3 % regular correlation testing
                % CV: Change
                [RHO_opt, p, epsilon_val, omega_val, u_opt, v_opt] = dp_projection_ext(OUT_x.test, OUT_y.test, u_opt, v_opt, correlation_method);
        end

        success_val=true;


        output.validation_results(n_LV,:) = {n_LV c_weights_opt weights_opt success_val RHO_opt p lVs_val};


        if p >= 0.05
            count_ns = count_ns+1;
        end
    else
        if p_LV > p_threshold
            count_ns = count_ns+1;
        end

    end

    %% Matrix Deflation => Projection Deflation
    % This method removes the covariance explained by u and v from X
    % and Y by projecting the data matrices onto the space spanned by the
    % corresponding weight vector, and subtracting this from the data:
    % choose whether X and Y should be deflated separately or together


    if ~islogical(input.validation_set)
        for num_m=1:num_matrices
            weights{num_m} = output.validation_results{n_LV,opt_weights{num_m}}; % weight vector for matrix num_m
            [matrices_val{num_m}] = proj_def_single(matrices_val{num_m}, weights{num_m});
        end
        
    else
        for num_m=1:num_matrices
            weights{num_m} = output.final_parameters{n_LV,opt_weights}; % weight vector for matrix num_m
            matrices_ana{num_m} = proj_def_single(matrices_ana{num_m}, weights{1}{num_m});
        end

    end

    disp('checkpoint deflation');
    toc
    
    n_LV = n_LV+1;
   
    if isfield(setup, 'saveall') && setup.saveall
        save([final_results, '/preliminary_result.mat']);
    else 
        save([final_results, '/preliminary_result.mat'], 'input', 'output', 'setup');
    end

end

% after the LVs are computed, clean up the final parameters matrix by removing empty rows
% output.final_parameters(n_LV:end,:) = [];

save([final_results, '/result.mat'], 'input', 'output', 'setup');
delete([final_results, '/preliminary_result.mat']);

rmdir(hyperopt_folder, 's');
rmdir(permutation_folder,'s');
rmdir(bootstrap_folder,'s');

cd(setup.analysis_folder)
end

function matrix = get_split_matrix(matrix,inds)
matrix = matrix(inds,:);
end
