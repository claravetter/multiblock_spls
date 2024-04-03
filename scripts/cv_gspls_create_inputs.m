%% generalized sPLS Analysis: 
% based on: volume/projects/RU_DP_immune/ScrFun/SPLS_data_prep_immune_07_2022.m
% Date: 18.03.2024

%% load PLEXUS test data
clear all
load('/Users/claravetter/local/Projects/mSPLS/20240131_testrun_realdata_mri_plexus_clinical/mri_plexus_clinical_matrices.mat')


%% Assign data for analysis in SPLS toolbox and create data and setup structure
% Parameters for IT infrastructure
setup.spls_standalone_path  = '/Users/claravetter/local/Code/multiblock_spls/'; % Path of the SPLS Toolboxsetup.date                  = date; % automatic date
setup.date                  = date; % automatic date
setup.analysis_folder       = ['/Users/claravetter/local/Projects/mSPLS/', setup.date]; % analysis folder is created automatically using the current date
setup.queue_name            = 'all.q'; % Choose queue for master and slave jobs, choose between psy0cf20 and mitnvp1-2 and all.q
setup.email                 = 'clara.vetter@med.uni-muenchen.de'; % your email address
setup.max_sim_jobs          = 1; % Define how many parallel jobs are created
setup.parallel_jobs         = 20; % Define how many jobs run in parallel at the same time (soft threshold)
setup.mem_request           = 5; % Memory request for master and slave jobs, in GB, maybe decrease to 2 or 3
setup.matlab_version        = 'R2022b'; % Define the runtime engine, currently its R2020b
setup.cache_path            = '/Users/claravetter/local/Projects/mSPLS/temp/'; % Path for output text files during hyperopt, permutation, bootstrapping => generally same as scratch space
setup.scratch_space         = '/Users/claravetter/local/Projects/mSPLS/nmcache/clvetter'; % Path for temporary file storage (hyperopt, permutation, bootstrapping) during analysis, please insert your own folder in the scratch space
setup.compilation_subpath   = 'for_testing'; % default

% Accessory data input
input.final_ID              = num2cell(1:267); % ID of your subjects (cell array)
%input.data_complete         = data_complete; % Optional, not needed for analysis: I store a large table of the entire sample in my files, mostly for post hoc analyses, 
input.type_correction       = 'uncorrected'; % Define whether you want to correct for covariates, choose between: correct, uncorrected

% TO DO: input should be an array with logical indices which matrices
% should be corrected 
input.correction_target     = 1; % Define whether you want to remove covariate effects from 1) X, 2) Y or 3) both matrices

input.sites                 = ones(267,1); % Dummy coded vector for sites, if only one site, then enter a column vector of ones
%input.sites_names           = input.data_complete.foranalysis.sites.Properties.VariableNames; % Names of the sites
%input.covariates            = [covariate]; % [input.sites]; % Matrix with covariates (double)
%input.covariates_names      = [covariate_name]; % Names of covariates (cell array)
input.Diag                  = [ones(150,1); ones(117,1)+1]; % input.data_complete.foranalysis.basic{PSN_MRI_final, 'DiagNum'}; % Column vector with diagnoses coded via numbers, i.e., [1 3 2 3 4 ]
input.DiagNames             = num2cell(input.Diag); % input.data_complete.foranalysis.basic{input.Y_final.Properties.RowNames, 'Labels'}; % Column cell array with diagnoses/labels, i.e., {'HC', 'ROD', 'CHR', 'HC', 'ROP'}

% define X and Y
% TO DO: change to matrices (cell array)
input.matrices              = matrices; % 1st data matrix, usually for MRI/biological data, if applicable. If MRI data, then either put in the path to a Matlab file, containing only one variable with vectorized MRI data, or put in vectorized MRI data itself, otherwise just put in the matrix (double format)
input.matrix_names               = []; % define names of features in X, if MRI data, or no names applicable, leave empty
%input.Y                     = Y_final.Variables; % 2nd data matrix, usually for behavioral/phenotypical data (double format) 
%input.Y_names               = Y_final.Properties.VariableNames; % define names of features in X, if no names applicable, leave empty
%input.subscales             = input.Y_names; % Optional, only for post hoc visualization needed

% Define ML framework
input.framework             = 1; % Cross-validation setup: 1 = nested cross-validation, 2 = random hold-out splits, 3 = LOSOCV, 4 = random split-half
input.outer_folds           = 5; % Applicable only for nested cross-validation and Random Hold-Out Splits: Define Outer folds CV2
input.inner_folds           = 10; % Applicable only for nested cross-validation and Random Hold-Out Splits: Define Inner folds CV1
input.permutation_testing   = 5000; % Number of permutations for significance testing of each LV, default: 5000
input.bootstrap_testing     = 500; % Number of bootstrap samples to measure Confidence intervals and bootstrap ratios for feature weights within LV: default 500 (100 also possible)
input.correlation_method    = 'Spearman'; % Define which correlation method is used to compute correlation between latent scores of X and Y (used for significance testing of LV): default 'Spearman', also possible 'Pearson'
input.cs_method.method      = 'mean-centering'; % Scaling of features, default: mean-centering, also possible 'min_max' (scaling from 0 to 1) => preferred scaling is mean-centering!
input.cs_method.correction_subgroup = 'HC'; % Define whether you want to correct the covariates based on the betas of a subgroup, or simply across all individuals => for subgroup-based correction use the label, i.e., 'HC' or 'ROD, etc. Otherwise leave as empty string: ''.
input.coun_ts_limit         = 1; % Define after how many non-significant LVs the algorithm should stop, default: 1 (means that as soon as one LV is not significant, the operation ends)
input.outer_permutations    = 1; % Define number of permutations in the CV2 folds, default: 1 (Toolbox is so far not optimized for permutations on folds, also, permutating the folds would severely increase computation time and is therefore not recommended
input.inner_permutations    = 1; % Define number of permutations in the CV1 folds, default: 1 (Toolbox is so far not optimized for permutations on folds, also, permutating the folds would severely increase computation time and is therefore not recommended

% training setup within ML framework
input.selection_train       = 1; % 1) Define how the RHO values between X and Y are collected across the cross-validation structure, default: 1, possible options: 1) within one CV2 fold, 2) across all CV2 folds (option 2 not recommended)
input.selection_retrain     = 1; % 1) Define whether you want to pool data from all CV1 folds and retrain the model on these before applying on CV2 testing fold, default: 1, possible options: 1) retrain on all CV1 folds, 2) no retraining, use already existing model
input.merge_train           = 'median'; % Define how the RHO values are collected, default: median, possible options: mean, median
input.merge_retrain         = 'best'; % Define how the best hyperparameters will be chosen on the CV1 level, default: 'best' (winner takes all, the best performing CV1 hyperparameters will be chosen), possible options: mean, median, weighted_mean, best => mean, median and weighted mean lead to a merging of all CV1 models

% Additional validation
input.validation_set        = false; % Define whether you want to hold out a validation set, default: false, possible options: false or number as percentage of the whole sample, i.e., 25, 50, etc. => 50 means 50% means random split half
input.val_stratification    = 1; % If applicable, define how you want to extract the validation set, options: 1) diagnosis, 2) sites, 3) both
input.validation_train      = 1; % If applicable, define how you want to test performance of the model on the validation set, options: 1) Retrain optimal model on permutations of the all samples, except for validation set, 2) use already computed permuted performances from the CV structure, default: 1

% Significance testing
input.alpha_value           = 0.05; % Define overall threshold for significance
input.final_merge.type      = 'best'; % Define how the final LV model will be chosen on the CV2 level, default: 'best' (winner takes all, the best performing CV2 model will be the new LV), possible options: mean, median, weighted_mean, best => mean, median and weighted mean lead to a merging of all CV2 models => the feature weights of all CV2 models are merged via mean, median or weighted mean (based on RHO values)
input.final_merge.mult_test = 'Benjamini_Hochberg'; % Define how correction for multiple testing across CV2 folds is done, default: 'Benjamini_Hochberg', possible options: Bonferroni, Sidak, Holm_Bonferroni, Benjamini_Hochberg, Benjamini_Yekutieli, Storey, Fisher
input.final_merge.significant_only  = 'on'; % Only applicable if input.final_merge.type is not set to best! Defines type of CV2 fold merging: options: 'on' use only significant folds for merging, 'off' use all folds for merging
input.final_merge.majority_vote     = 'on'; % Only applicable if input.final_merge.type is not set to best! options: 'on' use majority voting across folds to determine whether a value in u or v should be zero or non-zero, 'off' no majority vote, merging is done for all features, irrespective of whether in the majority of folds the feature was zero
input.correct_limit         = 1; % Define in which iteration of the process covariate correction should be done, default: 1 (means that covariate correction is done before computing the first LV, then no more correction)
input.statistical_testing   = 2; % Define how the P value is computed during permutation testing: 1) Counting method (number of instance where permuted models outperformed opimized model/number of permutations), 2) AUC method (permuted RHO values are used to compute AUC for optimal RHO value => option 2 usually gives slightly lower P values

% Hyperopt grid search specifics: Define in which LV iteration you want to
% use which grid, default: stable grid across all iterations
input.grid_dynamic.onset    = 1; % Choose the marks for grid applications, default: 1 (means that one grid is defined at first iteration and then not changed at all in later iterations)


% TO DO: grid needs to be defined based on number of matrices; better to
% add this as cell also 
input.density               = num2cell([50, 50, 50]); % should this be the same for all matrices or per matrix? 
input.grid_dynamic.LVs   = cellfun(@create_grid, input.density); 
%input.grid_dynamic.LV_1.x   = struct('start', 1, 'end', 0, 'density', input.density); 
%input.grid_dynamic.LV_1.y   = struct('start', 1, 'end', 0, 'density', input.density);
% Define grid for hyperparameter search for the X matrix (cu parameter) =>
% 'start': 1 means start is at value 1, 10 means it starts at the lower 10
% percentile of the grid, etc. 
% 'end' defines the upper limit of the hyperparameter search, 0 means all
% the way to the end, 10 means to stop at the upper 10 percentile, etc.,
% 'density' defines the number of data points which are tested during the
% grid, i.e., 20 means that between start and end point, 20 equidistant
% values are tested for the hyperparameter

%% Create analysis datafile

input.name = ['CV_gspls_test_', num2str(input.outer_folds), 'x', num2str(input.inner_folds), '_', num2str(input.permutation_testing), 'perm_', num2str(input.bootstrap_testing), 'boot_', num2str(50), 'density']; 

input.datafile = [setup.analysis_folder, '/' setup.date, '_', input.name, '_datafile.mat']; % Path for storing datafile containing input and setup

mkdir(setup.analysis_folder); % Creates analysis folder
save([input.datafile], 'setup', 'input'); % Saves datafile in analysis folder

%try rmdir('/volume/mitnvp1_scratch/CW_Med/.mcrCache9.9', 's'); % Cleans your scratch space, be careful: only do this when you have nothing running in the queue
%catch
%    disp('no deletion');
%end

%dp_bash_main_job_slim_addedruntime_mult(input.datafile); % Start the analysis


%%
function grid = create_grid(density)
grid = struct('start', 1, 'end', 0, 'density', density); 
end
