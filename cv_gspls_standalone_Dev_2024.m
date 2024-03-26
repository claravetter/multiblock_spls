%% DP SPLS standalone function

function cv_gspls_standalone_Dev_2024(datafile)

%% initialize analysis folders
% load the datafile that contains the data and the analysis information
load(datafile, 'input', 'setup');

% this sets up the directory where the results of the current analysis are
% saved
nn=1;
analysis_name = input.name;
analysis_folder = [setup.analysis_folder  '/' analysis_name];
if ~exist(analysis_folder,'dir')
    mkdir(analysis_folder);
else
    while nn < 100
        analysis_name = [input.name '_' num2str(nn)];
        analysis_folder = [setup.analysis_folder  '/' analysis_name];
        if ~exist(analysis_folder,'dir')
            mkdir(analysis_folder)
            nn=100;
        else
            nn=nn+1;
        end
    end
end

%% Prepare analysis folders
% sets the subfolders for the different analysis steps
% temporary files for hyperopt and permutation will be saved to the common scratch space "/volume/mitnvp1_scratch"

[permutation_folder, hyperopt_folder, bootstrap_folder, detailed_results, final_results] = dp_create_folders(setup.scratch_space, analysis_folder, analysis_name);

cd(final_results);

%% Set parameters for SPLS analysis
[input, setup, X, Y, B, K, W, OB, IB, size_sets_permutation, size_sets_bootstrap, correlation_method, cs_method, selection_train, selection_retrain, correction_target] = dp_setup_parameters(input, setup);

% define column names for matrices so that you can access them later by
% indexing
output.parameters_names = {'w', 'cu', 'cv', 'u', 'v', 'success', 'RHO', 'p', 'epsilon', 'omega', 'V_opt'}; % names of parameters for optimal cu/cv combination of the w-th loop
output.opt_parameters_names = {'w', 'cu', 'cv', 'u', 'v', 'success', 'RHO', 'p', 'epsilon', 'omega', 'V_opt'};

% preallocate placeholder matrices for higher speed
% output.final_parameters = num2cell(nan(10, numel(output.parameters_names))); % cell array to store final parameters for each LV iteration

% indices for later use
opt_u = strcmp(output.parameters_names,'u');
opt_v = strcmp(output.parameters_names,'v');
opt_p = strcmp(output.parameters_names,'p');
opt_cu = strcmp(output.parameters_names,'cu');
opt_cv = strcmp(output.parameters_names,'cv');
opt_RHO = strcmp(output.parameters_names,'RHO');

%% Set up analysis and validation sets, including covariates

% ff counts the iteration through the LV loops, so that the results for each
% LV can be stored in separate rows
ff = 1;

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

X_val = X(validation_set,:);
Y_val = Y(validation_set,:);
Diag_val = input.Diag(validation_set,:);
DiagNames_val = input.DiagNames(validation_set,:);
sites_val = input.sites(validation_set,:);

X_ana = X(analysis_set,:);
Y_ana = Y(analysis_set,:);
Diag_ana = input.Diag(analysis_set,:);
DiagNames_ana = input.DiagNames(analysis_set,:);
sites_ana = input.sites(analysis_set,:);

try Covars_ana = input.covariates(analysis_set,:);
    Covars_val = input.covariates(validation_set,:);
catch
    Covars_ana = nan(size(input.Diag(analysis_set,:),1),1);
    Covars_val = nan(size(input.Diag(validation_set,:),1),1);
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
output.CV = dp_setup_framework(IN);

% save permutation setup files
temp_pc = repelem(1:W, ceil(setup.max_sim_jobs/W));
temp_pc(randi(setup.max_sim_jobs,1,size(temp_pc,2)-setup.max_sim_jobs))=[];
perm_coding = [[1:setup.max_sim_jobs;temp_pc]; 1:(B/setup.max_sim_jobs):(B-1)];
rest_perm = mod(B,size_sets_permutation);
if rest_perm>0
    disp(['Please choose a number of permutation sets, which can be divided by ', num2str(setup.max_sim_jobs)]);
end
perm_sets = (B - rest_perm)/size_sets_permutation;

rest_boot = mod(input.bootstrap_testing, size_sets_bootstrap);
if rest_boot>0
    disp(['Please choose a number of bootstrap sets, which can be divided by ', num2str(setup.max_sim_jobs)]);
end
boot_sets = (input.bootstrap_testing - rest_boot)/size_sets_bootstrap;

save([permutation_folder '/permutation_setup.mat'],'selection_train', 'correlation_method', 'size_sets_permutation', 'perm_coding', 'selection_retrain', 'correction_target', '-v7.3');
save([bootstrap_folder '/bootstrap_setup.mat'],'selection_train', 'correlation_method', 'size_sets_bootstrap', 'correction_target', '-v7.3');

% set count for not significant LVs to 0, as soon as one LV is not
% significant, the LV loop stops
count_ns = 0;

if strcmp(input.merge_retrain, 'weighted_mean')
    weights_log = true;
else
    weights_log = false;
    weights_u = [];
    weights_v = [];
    weights_V = [];
    weights_RHO=[];
end

%% Latent Variable loop begins: setup of grid and cu/cv combination for analysis

% here starts the outer loop for each single LV, the next LV is only
% computed if the previous LV was significant (with FDR correction)
while count_ns<input.coun_ts_limit
    
    % repeats the entire process W times to obtain W different final p values,
    % which go into the omnibus hypothesis
    
    % matrix to store the optimal parameters of the w-th loop
    output.opt_parameters.(['LV_', num2str(ff)]) = num2cell(nan(W,numel(output.parameters_names)));
    
    if any(ff==input.grid_dynamic.onset)
        grid_x = input.grid_dynamic.(['LV_', num2str(ff)]).x;
        grid_y = input.grid_dynamic.(['LV_', num2str(ff)]).y;
        input.size_sets_hyperopt = round((grid_x.density*grid_y.density)/setup.max_sim_jobs);
        size_sets_hyperopt = input.size_sets_hyperopt;
        if size_sets_hyperopt==0
            size_sets_hyperopt=1;
        end
        [cu_cv_combination, hyperopt_sets] = dp_cu_cv_ext(X_ana, Y_ana, grid_x, grid_y, size_sets_hyperopt);
    end
    
    RHO_ICV_collection=[]; u_ICV_collection=cell(size(cu_cv_combination,1),1); v_ICV_collection=cell(size(cu_cv_combination,1),1);
    RHO_opt_collection=[]; epsilon_opt_collection=[]; omega_opt_collection=[]; V_opt_collection=[];
    weights_u = []; weights_v = []; weights_V=[]; weights_RHO=[];
    
    if ff>input.correct_limit
        input.corrected_log(ff) = false;
    elseif strcmp(input.type_correction, 'uncorrected')
        input.corrected_log(ff) = false;
    else
        input.corrected_log(ff) = true;
    end
    
    if ~input.corrected_log(ff)
        Covars_ana = nan(size(input.Diag(analysis_set,:),1),1);
        Covars_val = nan(size(input.Diag(validation_set,:),1),1);
        correction_target = 3;
    end
    
    %% Outer folds begin: hyperopt
    
    for w=1:W
        
        % remove hp% of the data randomly and keep it in a hold-out
        % dataset, the rest is called the keep_in data for training/testing
        % the algorithm
        RHO_ICV_collection_temp=[]; u_ICV_collection_temp={}; v_ICV_collection_temp={};
        
        for ob=1:OB
            test_data_x = X_ana(output.CV.cv_outer_indices.TestInd{ob,w},:);
            test_data_y = Y_ana(output.CV.cv_outer_indices.TestInd{ob,w},:);
            test_Diag = Diag_ana(output.CV.cv_outer_indices.TestInd{ob,w},:);
            test_DiagNames = DiagNames_ana(output.CV.cv_outer_indices.TestInd{ob,w},:);
            test_covariates = Covars_ana(output.CV.cv_outer_indices.TestInd{ob,w},:);
            
            train_data_x = X_ana(output.CV.cv_outer_indices.TrainInd{ob,w},:);
            train_data_y = Y_ana(output.CV.cv_outer_indices.TrainInd{ob,w},:);
            train_Diag = Diag_ana(output.CV.cv_outer_indices.TrainInd{ob,w},:);
            train_DiagNames = DiagNames_ana(output.CV.cv_outer_indices.TrainInd{ob,w},:);
            train_covariates = Covars_ana(output.CV.cv_outer_indices.TrainInd{ob,w},:);
            
            cv_inner_TrainInd = output.CV.cv_inner_indices{ob,w}.TrainInd;
            cv_inner_TestInd = output.CV.cv_inner_indices{ob,w}.TestInd;
            cv_outer = output.CV.cv_outer_indices;
            
            save([hyperopt_folder '/hyperopt_partition.mat'],...
                'train_data_x', 'train_data_y', 'train_covariates', 'cs_method',...
                'train_DiagNames', 'test_DiagNames',...
                'cv_inner_TrainInd', 'cv_inner_TestInd','correlation_method',...
                'cu_cv_combination', 'size_sets_hyperopt', 'correction_target', '-v7.3');
            
            if selection_train == 2 && ob == OB
                
                save([permutation_folder '/permutation_partition_fold_', num2str(w), '.mat'],...
                    'train_data_x', 'test_data_x', 'train_covariates', 'test_covariates',...
                    'train_data_y', 'test_data_y', 'train_Diag', 'train_DiagNames', 'test_DiagNames', 'cs_method', '-v7.3');
            end
            
            %% Hyperparameter optimization starts
            % this performs training and testing within the inner training
            % and testing folds of one OCV => according to W, K
            RHO_ICV_collection_temp = dp_ICV_main_csv_mult(setup.spls_standalone_path, setup.queue_name, hyperopt_folder, 'hyperopt', hyperopt_sets, setup.parallel_jobs, setup.mem_request, setup.matlab_version, setup.compilation_subpath, setup.cache_path);
            RHO_ICV_collection = cat(2, RHO_ICV_collection, RHO_ICV_collection_temp);
            
        end
        
        switch selection_train % 1=within one OCV fold, 2=across all OCV folds
            case 1 % within one OCV fold
                save([hyperopt_folder, '/checkpoint_hyperopt.mat'], 'RHO_ICV_collection', 'cu_cv_combination');
                
                switch input.merge_train
                    case 'mean'
                        f_hyperopt = @(x)mean(x,2);
                    case 'median'
                        f_hyperopt = @(x)median(x,2);
                    case 'best'
                        f_hyperopt = @(x)max(x,2);
                end
                
                for icv=1:size(cu_cv_combination,1)
                    RHO_selection(icv,1) = f_hyperopt(RHO_ICV_collection(icv,:));
                end
                
                [~,I_max] = max(RHO_selection);
                cu_opt = cu_cv_combination(I_max, 1);
                cv_opt = cu_cv_combination(I_max, 2);
                
                % retrain the inner folds of this one outer fold with
                % cu_max and cv_max => options: along all inner folds
                % or with all inner folds combined
                % options: 1) retrain on single test splits within
                % folds, then merge, 2) retrain on all inner folds
                % separately, then merge with mean or median, 3)
                % retrain on entirety of inner folds, 4) use already
                % existing u and v from inner folds without retraining
                
                for ob=1:OB
                    
                    COV.test = Covars_ana(output.CV.cv_outer_indices.TestInd{ob,w},:);
                    COV.train = Covars_ana(output.CV.cv_outer_indices.TrainInd{ob,w},:);
                    
                    IN_x.train = X_ana(output.CV.cv_outer_indices.TrainInd{ob,w},:);
                    IN_x.test = X_ana(output.CV.cv_outer_indices.TestInd{ob,w},:);
                    
                    IN_y.train = Y_ana(output.CV.cv_outer_indices.TrainInd{ob,w},:);
                    IN_y.test = Y_ana(output.CV.cv_outer_indices.TestInd{ob,w},:);
                    
                    DiagNames_train = DiagNames_ana(output.CV.cv_outer_indices.TrainInd{ob,w},:);
                    DiagNames_test = DiagNames_ana(output.CV.cv_outer_indices.TestInd{ob,w},:);
                    
                    if ~isempty(cs_method.correction_subgroup)
                        cs_method.subgroup_train = contains(DiagNames_train, cs_method.correction_subgroup);
                        cs_method.subgroup_test = contains(DiagNames_test, cs_method.correction_subgroup);
                    else
                        cs_method.subgroup_train = [];
                        cs_method.subgroup_test = [];
                    end
                    
                    [OUT_x, OUT_y] = dp_master_correctscale(IN_x, IN_y, COV, cs_method, correction_target);
                    
                    switch selection_retrain
                        
                        case 1 % retrain on entirety of folds
                            [RHO_opt_collection(1,ob), u_opt_collection(:,ob), v_opt_collection(:,ob), V_opt_collection_temp, epsilon_opt_temp, omega_opt_temp] = dp_spls_full(OUT_x.train,OUT_y.train,OUT_x.test, OUT_y.test, cu_opt, cv_opt, correlation_method);
                            V_opt_collection = cat(3, V_opt_collection, V_opt_collection_temp);
                            
                            if weights_log
                                weights_u(:,ob) = RHO_opt_collection(1,ob)*ones(size(u_opt_collection(:,ob),1),1)';
                                weights_v(:,ob) = RHO_opt_collection(1,ob)*ones(size(v_opt_collection(:,ob),1),1)';
                                weights_V_1 = RHO_opt_collection(1,ob)*ones(size(V_opt_collection_temp));
                                weights_V = cat(3, weights_V, weights_V_1);
                            end
                            
                            save([hyperopt_folder, '/weights_check.mat'], 'weights_u', 'weights_v', 'weights_V');
                            clear V_opt_collection_temp;
                            
                            epsilon_opt_collection = cat(1, epsilon_opt_collection, epsilon_opt_temp);
                            omega_opt_collection = cat(1, omega_opt_collection, omega_opt_temp);
                            
                        case 2 % use already existing u and v from inner folds without retraining
                            u_opt_collection = u_selection(I_max,:)';
                            v_opt_collection = v_selection(I_max,:)';
                            V_opt_collection = false;
                            
                            [RHO_opt_collection(1,ob), epsilon_opt_temp, omega_opt_temp, u_opt_collection, v_opt_collection] = dp_projection(OUT_x.test, OUT_y.test, u_opt_collection, v_opt_collection, correlation_method);
                            
                            epsilon_opt_collection = cat(1, epsilon_opt_collection, epsilon_opt_temp);
                            omega_opt_collection = cat(1, omega_opt_collection, omega_opt_temp);
                            
                    end
                    
                end
                
                switch input.merge_retrain
                    case 'best'
                        [RHO_opt,I_opt] = max(RHO_opt_collection);
                        u_opt = u_opt_collection(:, I_opt);
                        v_opt = v_opt_collection(:, I_opt);
                        V_opt = V_opt_collection(:, :, I_opt);
                    otherwise
                        RHO_opt = f_hyperopt(RHO_opt_collection);
                        u_opt = dp_trainmerge_single(u_opt_collection, input.merge_retrain, 2, weights_u);
                        v_opt = dp_trainmerge_single(v_opt_collection, input.merge_retrain, 2, weights_v);
                        V_opt = dp_trainmerge_single(V_opt_collection, input.merge_retrain, 3, weights_V);
                end
                
                save([hyperopt_folder '/opt_vars_check.mat'],...
                    'RHO_opt_collection','u_opt_collection','v_opt_collection','V_opt_collection', 'RHO_opt', 'u_opt', 'v_opt', 'V_opt', 'RHO_ICV_collection', '-v7.3');
                
                %% Permutation Testing within fold
                
                test_covariates = Covars_ana(output.CV.cv_outer_indices.TestInd{ob,w},:);
                train_covariates = Covars_ana(output.CV.cv_outer_indices.TrainInd{ob,w},:);
                
                train_data_x = X_ana(output.CV.cv_outer_indices.TrainInd{ob,w},:);
                test_data_x = X_ana(output.CV.cv_outer_indices.TestInd{ob,w},:);
                train_data_y = Y_ana(output.CV.cv_outer_indices.TrainInd{ob,w},:);
                test_data_y = Y_ana(output.CV.cv_outer_indices.TestInd{ob,w},:);
                
                train_Diag = Diag_ana(output.CV.cv_outer_indices.TrainInd{ob,w},:);
                
                % save the optimized parameters and the permutation matrix
                save([permutation_folder '/permutation_partition_fold.mat'],...
                    'train_data_x', 'test_data_x', 'train_covariates', 'test_covariates',...
                    'train_data_y', 'test_data_y', 'train_Diag', 'train_DiagNames', 'test_DiagNames', 'cs_method', '-v7.3');
                
                save([permutation_folder '/permutation_opt.mat'],...
                    'cu_opt','cv_opt','V_opt','-v7.3');
                
                % choose whether permutation testing shall be conducted right
                % within the current fold or whether permutation testing shall
                % be done across all outer folds
                RHO_b_collection = dp_ICV_main_csv_mult(setup.spls_standalone_path, setup.queue_name, permutation_folder, 'permutation', perm_sets, setup.parallel_jobs, setup.mem_request, setup.matlab_version, setup.compilation_subpath, setup.cache_path);
                
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
                output.opt_parameters.(['LV_', num2str(ff)])(w,:) = {w cu_opt cv_opt u_opt v_opt success_ICV RHO_opt p epsilon_opt_collection omega_opt_collection, V_opt};
                save([detailed_results '/preliminary_results.mat'],'input', 'setup', 'output');
                
                % clear collections for next W fold
                RHO_ICV_collection=[]; u_ICV_collection=cell(size(cu_cv_combination,1),1); v_ICV_collection=cell(size(cu_cv_combination,1),1);
                RHO_opt_collection=[]; epsilon_opt_collection=[]; omega_opt_collection=[]; V_opt_collection=[];
                RHO_selection = []; u_selection = []; v_selection = [];
                weights_u=[]; weights_v=[]; weights_V=[];
                
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
        
        for icv=1:size(cu_cv_combination,1)
            RHO_selection(icv,1) = f_hyperopt(RHO_ICV_collection(icv,:));
        end
        
        [~,I_max] = max(RHO_selection);
        cu_opt = cu_cv_combination(I_max, 1);
        cv_opt = cu_cv_combination(I_max, 2);
        
        % retrain the inner folds of this one outer fold with
        % cu_max and cv_max => options: along all inner folds
        % or with all inner folds combined
        cv_outer = output.CV.cv_outer_indices;
        V_opt_collection=[];
        nn=1;
        for w=1:W
            for ob=1:OB
                
                COV.test = Covars_ana(output.CV.cv_outer_indices.TestInd{ob,w},:);
                COV.train = Covars_ana(output.CV.cv_outer_indices.TrainInd{ob,w},:);
                
                IN_x.train = X_ana(output.CV.cv_outer_indices.TrainInd{ob,w},:);
                IN_x.test = X_ana(output.CV.cv_outer_indices.TestInd{ob,w},:);
                
                IN_y.train = Y_ana(output.CV.cv_outer_indices.TrainInd{ob,w},:);
                IN_y.test = Y_ana(output.CV.cv_outer_indices.TestInd{ob,w},:);
                
                [OUT_x, OUT_y] = dp_master_correctscale(IN_x, IN_y, COV, cs_method, correction_target);
                
                switch selection_retrain
                    case 1 % retrain on entirety of outer folds
                        [RHO_opt_collection(1,nn), u_opt_collection(:,nn), v_opt_collection(:,nn), V_opt_collection_temp, epsilon_opt_temp, omega_opt_temp] = dp_spls_full(OUT_x.train,OUT_y.train,OUT_x.test, OUT_y.test, cu_opt, cv_opt, correlation_method);
                        V_opt_collection = cat(3, V_opt_collection, V_opt_collection_temp);
                        
                        if weights_log
                            if input.framework == 3
                                ratio = size(output.CV.cv_outer_indices.TestInd{1,w},1)/size(cat(1,output.CV.cv_outer_indices.TestInd{1,:}),1);
                            else
                                ratio = RHO_opt_collection(1,nn);
                            end
                            weights_RHO(1,nn) = ratio;
                            weights_u(:,nn) = ratio*ones(size(u_opt_collection(:,nn),1),1)';
                            weights_v(:,nn) = ratio*ones(size(v_opt_collection(:,nn),1),1)';
                            weights_V_1 = ratio*ones(size(V_opt_collection_temp));
                            weights_V = cat(3, weights_V, weights_V_1);
                        end
                        
                        clear V_opt_collection_temp;
                        
                        epsilon_opt_collection = cat(1, epsilon_opt_collection, epsilon_opt_temp);
                        omega_opt_collection = cat(1, omega_opt_collection, omega_opt_temp);
                        
                    case 2 % do not retrain, use already existing u and v
                        u_opt_collection = u_selection(I_max,:)';
                        v_opt_collection = v_selection(I_max,:)';
                        V_opt_collection = false;
                        
                        [RHO_opt_collection(1,nn), epsilon_opt_temp, omega_opt_temp, u_opt_collection, v_opt_collection] = dp_projection(OUT_x.test, OUT_y.test, u_opt_collection, v_opt_collection, correlation_method);
                        
                        epsilon_opt_collection = cat(1, epsilon_opt_collection, epsilon_opt_temp);
                        omega_opt_collection = cat(1, omega_opt_collection, omega_opt_temp);
                        
                end
                nn=nn+1;
            end
        end
        
        if weights_log
            save([hyperopt_folder, '/weights_check.mat'], 'ratio', 'weights_RHO', 'weights_u', 'weights_v', 'weights_V');
        end
        
        switch input.merge_retrain % defines how the results on the outer folds after retraining on entire inner folds will be condensed
            case 'best'
                [RHO_opt,I_opt] = max(RHO_opt_collection);
                u_opt = u_opt_collection(:, I_opt);
                v_opt = v_opt_collection(:, I_opt);
                V_opt = V_opt_collection(:, :, I_opt);
            otherwise
                RHO_opt = dp_trainmerge_single(RHO_opt_collection, input.merge_retrain, 2, weights_RHO);
                u_opt = dp_trainmerge_single(u_opt_collection, input.merge_retrain, 2, weights_u);
                v_opt = dp_trainmerge_single(v_opt_collection, input.merge_retrain, 2, weights_v);
                V_opt = dp_trainmerge_single(V_opt_collection, input.merge_retrain, 3, weights_V);
        end
        save([detailed_results '/opt_vars_LV', num2str(ff), '.mat'],...
            'RHO_opt_collection','u_opt_collection','v_opt_collection',...
            'RHO_opt', 'u_opt', 'v_opt', 'RHO_ICV_collection');
        
        % save the optimized parameters and the permutation matrix
        save([permutation_folder '/permutation_opt.mat'], 'cu_opt','cv_opt','V_opt','-v7.3');
        
        RHO_b_collection = dp_ICV_main_csv_mult(setup.spls_standalone_path, setup.queue_name, permutation_folder, 'permutation', perm_sets, setup.parallel_jobs, setup.mem_request, setup.matlab_version, setup.compilation_subpath, setup.cache_path);
        
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
        output.opt_parameters.(['LV_', num2str(ff)]) = {1 cu_opt cv_opt u_opt v_opt success_ICV RHO_opt p epsilon_opt_collection omega_opt_collection, V_opt};
        save([detailed_results '/preliminary_results.mat'],'input', 'setup', 'output');
        
        % clear collections for next W fold
        RHO_ICV_collection=[]; u_ICV_collection=cell(size(cu_cv_combination,1),1); v_ICV_collection=cell(size(cu_cv_combination,1),1);
        RHO_opt_collection=[]; epsilon_opt_collection=[]; omega_opt_collection=[]; V_opt_collection=[];
        RHO_selection = []; u_selection = []; v_selection = [];
        weights_u=[]; weights_v=[]; weights_V=[]; weights_RHO=[];
    end
    
    %% Final merging of folds (if applicable)
    
    % calculate FDR-corrected p value
    %     switch input.final_merge.mult_test
    %         case 'FDR'
    %             output.pvalue_FDR(ff,1) = dp_FDR([output.opt_parameters.(['LV_', num2str(ff)]){:,opt_p}], input.alpha_value);
    %         case 'SFisher'
    output.pvalue_FDR(ff,1) = 0.05;
    %     end
    
    % use different options to find optimal solution,
    % 1) collapse all folds using mean/median with/without majority vote
    % 2) collapse only significant folds using mean/median with/without majority vote
    % 3) take only best fold
    
    % set up function for collapsing/selecting folds => mean, median,
    % weighted mean, best(which automatically ends all further selection)
    temp = output.opt_parameters.(['LV_', num2str(ff)]);
    switch input.final_merge.type
        case 'best'
            input.final_merge.significant_only = 'off';
            input.final_merge.majority_vote = 'off';
            log_p_min = [temp{:,opt_p}]==min([temp{:,opt_p}]);
            if sum(log_p_min)>1
                temp = temp(log_p_min,:);
                log_RHO_max = [temp{:,opt_RHO}]==max([temp{:,opt_RHO}]);
                temp = temp(log_RHO_max,:);
            else
                temp = temp(log_p_min,:);
            end
            output.final_parameters(ff,:)=temp;
            p_threshold = output.pvalue_FDR(ff,1);
        otherwise
            p_threshold = 0.05;
    end
    
    % generate log for significant folds
    switch input.final_merge.significant_only
        case 'on'
            log_sig = [temp{:,opt_p}] <= output.pvalue_FDR(ff,1);
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
                u_temp=temp{i,opt_u};
                log_u(i,:)=u_temp~=0;
                v_temp=temp{i,opt_v};
                log_v(i,:)=v_temp~=0;
            end
            log_fm_u = sum(log_u,1)./size(log_u,1)<0.5;
            log_fm_v = sum(log_v,1)./size(log_v,1)<0.5;
        case 'off'
            log_fm_u = temp{1,opt_u} > 100;
            log_fm_v = temp{1,opt_v} > 100;
    end
    
    for i=1:size(temp,1)
        temp{i,opt_u}(log_fm_u) = 0;
        temp{i,opt_v}(log_fm_v) = 0;
    end
    
    % merge the folds according to the previous selection
    switch input.final_merge.type
        case 'mean'
            for i=1:size(temp,2)
                try output.final_parameters{ff,i}=mean([temp{:,i}],2);
                catch
                    temp_cat = [];
                    for ii=1:size(temp,1)
                        temp_cat=cat(1, temp_cat, temp{ii,i});
                    end
                    output.final_parameters{ff,i}=temp_cat;
                end
            end
        case 'median'
            for i=1:size(temp,2)
                try output.final_parameters{ff,i}=median([temp{:,i}],2);
                catch
                    temp_cat = [];
                    for ii=1:size(temp,1)
                        temp_cat=cat(1, temp_cat, temp{ii,i});
                    end
                    output.final_parameters{ff,i}=temp_cat;
                end
            end
        case 'weighted_mean'
            weights = 1-[temp{:,opt_p}];
            if size(temp,1)==1
                output.final_parameters(ff,:) = temp;
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
                        output.final_parameters{ff,i}=temp_wm_cat;
                    elseif isvector(temp_wm)
                        output.final_parameters{ff,i}=wmean(temp_wm, weights);
                    elseif ismatrix(temp_wm)
                        temp_weights=[];
                        for iii=1:size(weights,2)
                            temp_weights(:,iii) = ones(size(temp_wm,1),1)*weights(iii);
                        end
                        output.final_parameters{ff,i}=wmean(temp_wm, temp_weights,2);
                    end
                end
                output.final_parameters{ff,1}=ff;
            end
    end
    
    [c_pvalues, ~, ~] = dp_multiple_testing(input.final_merge.mult_test, [output.opt_parameters.(['LV_', num2str(ff)]){:,opt_p}], 0.05, false);
    
    for ppp=1:size(output.opt_parameters.(['LV_', num2str(ff)]),1)
        output.opt_parameters.(['LV_', num2str(ff)]){ppp,opt_p} = c_pvalues(ppp);
    end
    
    p_LV = min(c_pvalues);
    
    output.final_parameters{ff,opt_p} = p_LV;
    
    %% Bootstrapping for chosen LV u and v vectors
    w_opt = output.final_parameters{ff,matches(output.parameters_names, 'w')};
    cu_opt = output.final_parameters{ff,matches(output.parameters_names, 'cu')};
    cv_opt = output.final_parameters{ff,matches(output.parameters_names, 'cv')};
    V_opt = output.final_parameters{ff,matches(output.parameters_names, 'V_opt')};
    
    save([bootstrap_folder '/bootstrap_opt.mat'], 'cu_opt','cv_opt','V_opt','-v7.3');
    
    train_data_x = X_ana(output.CV.cv_outer_indices.TrainInd{ob,w_opt},:);
    train_data_y = Y_ana(output.CV.cv_outer_indices.TrainInd{ob,w_opt},:);
    train_covariates = Covars_ana(output.CV.cv_outer_indices.TrainInd{ob,w_opt},:);
    
    train_DiagNames = DiagNames_ana(output.CV.cv_outer_indices.TrainInd{ob,w_opt},:);
    
    % save the optimized parameters and the permutation matrix
    save([bootstrap_folder '/bootstrap_partition_fold.mat'],...
        'train_data_x', 'train_covariates', 'train_data_y', 'train_DiagNames', 'cs_method', '-v7.3');
    
    [RHO_boot, u_boot, v_boot] = dp_ICV_main_csv_mult(setup.spls_standalone_path, setup.queue_name, bootstrap_folder, 'bootstrap', boot_sets, setup.parallel_jobs, setup.mem_request, setup.matlab_version, setup.compilation_subpath, setup.cache_path);
    
    lv_name = ['LV_', num2str(ff)];

    % RHO
    RHO_analysis = output.final_parameters{ff, matches(output.opt_parameters_names, 'RHO')};
    RHO_mean = mean(RHO_boot);
    RHO_SE = std(RHO_boot)/(sqrt(input.bootstrap_testing));
    ci_RHO = [RHO_mean - 1.96 * RHO_SE, RHO_mean + 1.96 * RHO_SE];
    bs_ratio_RHO = RHO_analysis/RHO_SE;
    output.bootstrap_results.(lv_name).ci_RHO = ci_RHO;
    output.bootstrap_results.(lv_name).bs_ratio_RHO = bs_ratio_RHO;
    output.bootstrap_results.(lv_name).RHO_boot_size = size(RHO_boot);
    
    % u
    u_analysis = output.final_parameters{ff, matches(output.opt_parameters_names, 'u')};
    u_mean = mean(u_boot,2);
    u_SE = std(u_boot,0,2)/(sqrt(input.bootstrap_testing));
    
    bs_ratio_u = u_analysis./u_SE;
    bs_ratio_u(isnan(bs_ratio_u)) = 0;
    bs_ratio_u(bs_ratio_u == Inf) = 0;
    log_bs_u = abs(bs_ratio_u)<=2;
    
    ci_u = [u_mean - 1.96 * u_SE, u_mean + 1.96 * u_SE];
    log_ci_u = ((sum(ci_u>0, 2) == 2) + (sum(ci_u<0, 2) == 2)) == 0;

    output.bootstrap_results.(lv_name).ci_u = ci_u;
    output.bootstrap_results.(lv_name).bs_ratio_u = bs_ratio_u;
    output.bootstrap_results.(lv_name).u_boot_size = size(u_boot);
    output.bootstrap_results.(lv_name).log_bs_u = log_bs_u;
    output.bootstrap_results.(lv_name).sum_bs_u = sum(log_bs_u);
    output.bootstrap_results.(lv_name).log_ci_u = log_ci_u;
    output.bootstrap_results.(lv_name).sum_ci_u = sum(log_ci_u);
    
    % v    
    v_analysis = output.final_parameters{ff, matches(output.opt_parameters_names, 'v')};
    v_mean = mean(v_boot,2);
    v_SE = std(v_boot,0,2)/(sqrt(input.bootstrap_testing));
    
    bs_ratio_v = v_analysis./v_SE;
    bs_ratio_v(isnan(bs_ratio_v)) = 0;
    bs_ratio_v(bs_ratio_v == Inf) = 0;
    log_bs_v = abs(bs_ratio_v)<=2;
    
    ci_v = [v_mean - 1.96 * v_SE, v_mean + 1.96 * v_SE];
    log_ci_v = ((sum(ci_v>0, 2) == 2) + (sum(ci_v<0, 2) == 2)) == 0;

    output.bootstrap_results.(lv_name).ci_v = ci_v;
    output.bootstrap_results.(lv_name).bs_ratio_v = bs_ratio_v;
    output.bootstrap_results.(lv_name).v_boot_size = size(v_boot);
    output.bootstrap_results.(lv_name).log_bs_v = log_bs_v;
    output.bootstrap_results.(lv_name).sum_bs_v = sum(log_bs_v);
    output.bootstrap_results.(lv_name).log_ci_v = log_ci_v;
    output.bootstrap_results.(lv_name).sum_ci_v = sum(log_ci_v);
    
    % testing between BSR and CI
    output.testing.(lv_name).u_comparison = sum(log_bs_u == log_ci_u)/size(log_bs_u,1);
    output.testing.(lv_name).v_comparison = sum(log_bs_v == log_ci_v)/size(log_bs_v,1);
    
    %% check if LV is significant and apply it to validation set (if applicable)
    
    if ~islogical(input.validation_set)
        
        % apply the final parameters to the validations sets
        IN_x.train = X_ana;
        IN_x.test = X_val;
        IN_y.train = Y_ana;
        IN_y.test = Y_val;
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
        
        [OUT_x, OUT_y] = dp_master_correctscale(IN_x, IN_y, COV, cs_method_val, correction_target);
        
        temp = output.final_parameters(ff,:);
        cu_opt = temp{1,opt_cu};
        cv_opt = temp{1,opt_cv};
        u_opt = temp{1,opt_u};
        v_opt = temp{1,opt_v};
        
        switch input.validation_train
            case 1 % retrain
                [RHO_opt, u_opt, v_opt, ~, epsilon_val, omega_val] = dp_spls_full(OUT_x.train,OUT_y.train,OUT_x.test, OUT_y.test, cu_opt, cv_opt, correlation_method);
                
                % save the optimized parameters and the permutation matrix
                test_covariates = COV.test;
                train_covariates = COV.train;
                train_data_x = OUT_x.train;
                test_data_x = OUT_x.test;
                train_data_y = OUT_y.train;
                test_data_y = OUT_y.test;
                train_Diag = Diag_ana;
                
                selection_train=1;
                save([permutation_folder '/permutation_setup.mat'],'selection_train', 'correlation_method', 'size_sets_permutation', 'perm_coding', 'correction_target', '-v7.3');
                save([bootstrap_folder '/bootstrap_setup.mat'],'selection_train', 'correlation_method', 'size_sets_bootstrap', 'selection_retrain', 'correction_target', 'bs_method', '-v7.3');
                
                save([permutation_folder '/permutation_partition_fold.mat'],...
                    'train_data_x', 'test_data_x', 'train_covariates', 'test_covariates',...
                    'train_data_y', 'test_data_y', 'train_Diag', 'cs_method_val', '-v7.3');
                
                RHO_b_collection = dp_ICV_main_csv_mult(setup.spls_standalone_path, setup.queue_name, permutation_folder, 'permutation', perm_sets, setup.parallel_jobs, setup.mem_request, setup.matlab_version, setup.compilation_subpath, setup.cache_path);
                
                selection_train=input.selection_train;
                save([permutation_folder '/permutation_setup.mat'],'selection_train', 'correlation_method', 'size_sets_permutation', 'perm_coding', 'correction_target', '-v7.3');
                save([bootstrap_folder '/bootstrap_setup.mat'],'selection_train', 'correlation_method', 'size_sets_bootstrap', 'selection_retrain', 'correction_target', 'bs_method', '-v7.3');
                
            case 2 % apply optimal u and v
                [RHO_opt, epsilon_val, omega_val, u_opt, v_opt] = dp_projection(OUT_x.test, OUT_y.test, u_opt, v_opt, correlation_method);
                for pp=1:size(u_b_collection,1)
                    [RHO_b_collection(pp,1), ~, ~, ~, ~] = dp_projection(OUT_x.test, OUT_y.test, u_b_collection(pp,:)', v_b_collection(pp,:)', correlation_method);
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
                [RHO_opt, p, epsilon_val, omega_val, u_opt, v_opt] = dp_projection_ext(OUT_x.test, OUT_y.test, u_opt, v_opt, correlation_method);
        end
        
        success_val=true;
        output.validation_results(ff,:) = {ff cu_opt cv_opt u_opt v_opt success_val RHO_opt p epsilon_val omega_val};
        
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
        u = output.validation_results{ff,opt_u}; % weight vector for X
        v = output.validation_results{ff,opt_v}; % weight vector for Y
        [X_val,Y_val] = proj_def(X_val, Y_val, u, v);
    else
        u = output.final_parameters{ff,opt_u}; % weight vector for X
        v = output.final_parameters{ff,opt_v}; % weight vector for Y
    end
    
    [X_ana,Y_ana] = proj_def(X_ana, Y_ana, u, v);
    
    disp('checkpoint deflation');
    
    save([final_results, '/preliminary_result.mat'], 'input', 'output', 'setup');
    
    ff = ff+1;
    
end

% after the LVs are computed, clean up the final parameters matrix by removing empty rows
% output.final_parameters(ff:end,:) = [];

save([final_results, '/result.mat'], 'input', 'output', 'setup');
delete([final_results, '/preliminary_result.mat']);

rmdir(hyperopt_folder, 's');
rmdir(permutation_folder,'s');
rmdir(bootstrap_folder,'s');
end