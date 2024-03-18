%% sPLS Analysis: Medication - GMV
% based on: volume/projects/RU_DP_immune/ScrFun/SPLS_data_prep_immune_07_2022.m
% Date: 06.02.2023

%% Add Paths
clear all
addpath(genpath('/opt/NM/NeuroMiner_1.1/NeuroMiner_1.1/'));
addpath('/opt/SPM/spm12_v6685_cat12_r1207/');
addpath(genpath('/volume/DP_FEF/ScrFun/ScriptsRepository'));
addpath(genpath('/opt/PRONIASoftware/Developpment/DataAllCron/PRONIA_JuliaCode'));
addpath(genpath('/volume/projects/CV_StructCovNetw/'))
cd('/volume/projects/CV_StructCovNetw/')

%% Define Input Options
input.timepoint = 'T1'; % 'T0' || 'T1'
input.dosage_type = 'CD'; % 'CD' || 'DD'
input.selected_studygroups = {'ROP', 'CHR', 'ROD', 'HC'}; % {'ROP', 'CHR', 'ROD', 'HC'} || {'ROP', 'CHR', 'ROD'}
input.approach = 'THLD'; % 'THLD' || 'EQ'; 
if contains(input.approach, 'EQ')
    input.approach_eqtype = 'OLA'; % 'OLA' || 'CPZE'
end 
input.density = 40; 

%% Load Matrix containing MRI and Behavioral Data 
matrix_path = ['/volume/projects/CW_Med/Data/BehavioralData/MED_DATA_4sPLS_', input.approach, '.mat']; 
load(matrix_path)
clear S
if contains(input.approach, 'EQ')
    data_table = data_table.(input.approach_eqtype); 
end 

PSN_for_MRI = data_table.PSN;
BOG_all = data_table.(['BOGEN_ID_sMRI_', input.timepoint]); 
PSN_BOG_all = [data_table.PSN, data_table.(['BOGEN_ID_sMRI_', input.timepoint])]; 

% PSN_cellstr = num2str(data_table.PSN); cellstr(string(PSN_all)); %cell of char

%% MRI extraction, excluding the missing MRI data
MainDataFolder = '/volume/data/PRONIA/MRI/20-Dec-2017/Data/';

% Paths_xml = {}; % initialize path file
Paths = {}; % initialize the path file

PipelineName = 'MRI_sMRI_cat12_spm12_v6685_cat12_r1207_ar21112017';
image_modality1 = 's8rmwp1';
image_modality2 = 'mwp1';

TIV=[];
Mask_MRI=true;

Paths_names = {image_modality1, image_modality2, 'report_xml', 'report_mat'};

for i = 1:size(PSN_BOG_all,1)
    %     Paths{i,1} = PSN_BOG_all{i,1};
    Paths{i,1} = [MainDataFolder,PSN_BOG_all{i,2},'/',PipelineName,'/mri/' image_modality1 ,PSN_BOG_all{i,2},'_MRI_sMRI_',num2str(PSN_BOG_all{i,1}),'.nii'];
    Paths{i,2} = [MainDataFolder,PSN_BOG_all{i,2},'/',PipelineName,'/mri/' image_modality2 ,PSN_BOG_all{i,2},'_MRI_sMRI_',num2str(PSN_BOG_all{i,1}),'.nii'];
    Paths{i,3} = [MainDataFolder,PSN_BOG_all{i,2},'/',PipelineName,'/report/cat_',PSN_BOG_all{i,2},'_MRI_sMRI_',num2str(PSN_BOG_all{i,1}),'.xml'];
    Paths{i,4} = [MainDataFolder,PSN_BOG_all{i,2},'/',PipelineName,'/report/cat_',PSN_BOG_all{i,2},'_MRI_sMRI_',num2str(PSN_BOG_all{i,1}),'.mat'];
end

for i = 1:size(Paths,1)
    if exist(Paths{i,1})
        if ~exist(Paths{i,2})
            try gunzip([Paths{i,2}, '.gz']);
            catch
                Mask_MRI(i,1)=false;
                disp([Paths{i,2} ' cannot be unzipped.']);
            end
        end
        try load(Paths{i,4});
            TIV(i,1) = S.subjectmeasures.vol_TIV;
            Mask_MRI(i,1)=true;
        catch
            Mask_MRI(i,1)=false;
            disp([Paths{i,4} ' cannot be loaded.']);
        end
    elseif exist([Paths{i,1},'.gz'])
        try gunzip([Paths{i,1},'.gz']);
            if ~exist(Paths{i,2})
                gunzip([Paths{i,2}, '.gz']);
            end
            try load(Paths{i,4});
                TIV(i,1) = S.subjectmeasures.vol_TIV;
                Mask_MRI(i,1)=true;
            catch
                Mask_MRI(i,1)=false;
                disp([Paths{i,4} ' cannot be loaded.']);
            end
        catch
            Mask_MRI(i,1)=false;
            disp([Paths{i,1} ' cannot be unzipped.']);
        end
    else
        Mask_MRI(i,1) = false;
        disp([Paths{i,1} ' does not exist']);
    end
end
clear i 
mri_table = array2table(Paths, 'VariableNames', Paths_names, 'RowNames', PSN_for_MRI);
mri_table.Mask_MRI = Mask_MRI; clear Mask_MRI Paths_names

BOG_new=[];
for i=1:size(BOG_all,1)
    try BOG_new(i,1) = str2num(BOG_all{i,1});
    catch
        BOG_new(i,1)=NaN;
    end
end
mri_table.BOG = BOG_new; clear BOG_new BOG_all

% get data quality measures
% qm_names = {['Cat12_ICR_sMRI_', input.timepoint], ['Cat12_IQR_sMRI_', input.timepoint], ...
%     ['Cat12_NCR_sMRI_', input.timepoint]}; temp=[];
qm_names = {['Cat12_IQR_sMRI_', input.timepoint]}; 
temp=[];
for q=1:size(qm_names,2)
    for i=1:size(data_table.(qm_names{q}),1)
        temp(i,1) = data_table.(qm_names{q})(i);
    end
    mri_table.(qm_names{q}) = temp;
end

mri_table.TIV = TIV;
clear PSN_BOG_all image_modality1 image_modality2 PipelineName S temp qm_names TIV matrix_path Paths 
%% Cat12 homogeneity check
% Implemented here: /volume/projects/CW_Med/ScrFun/SCRIPT_CW_MEDICATION_2023_02_17.m
% Mahalanobis distance

PSN_final = data_table.PSN; 
clear temp1 files_not_found

%% Diagnosis
data_table_names = data_table.Properties.VariableNames; 
Diag_names = data_table_names(contains(data_table_names, 'Studygroup_', 'IgnoreCase', true));
Diagfull = table2array(data_table(:,contains(data_table_names, Diag_names))); % double
basic_table = array2table(Diagfull, 'VariableNames', Diag_names, 'RowNames', PSN_final); 

Labels = cell(size(data_table,1),1);DiagNum=[]; % cell
for i=1:size(Diagfull,1)
    temp = strfind(Diag_names(1, Diagfull(i,:)>0),'_');
    Labels{i,1} = Diag_names{1, Diagfull(i,:)>0}((temp{1}+1):end);
    DiagNum(i,1) = find(Diagfull(i,:)); % double
end
basic_table.Labels = Labels; clear Labels
basic_table.DiagNum = DiagNum; clear DiagNum Diagfull DiagNames temp temp1

%% Create a matrix with Age & Sex as covariates

% AGE
temp = data_table.(['AGE_', input.timepoint]); age_data=[];
for i=1:size(temp,1)
    age_data(i,1) = temp(i);
end
basic_table.age = age_data; 

% SEX
temp = data_table.SEX;
sex_data = zeros(size(age_data,1),2);
for i=1:size(temp,1)
    try sex_data(i,temp(i)) = 1;
    catch
        sex_data(i,:) = NaN;
    end
end

sex_names = {'male_sex', 'female_sex'};
for i=1:size(sex_names,2)
    basic_table.(sex_names{i}) = sex_data(:,i);
end
clear sex_data sex_names age_data
% Check original script for additional covariates (CTQ, PANSS, ....) 

% PSN_MRI_complete = PSN_all; % cell

%% Site
sites_names = unique(data_table.INSTITUTE_SHORTNAME)';sites_data=[];
for i=1:size(data_table,1)
    sites_data(i,:) = double(contains(sites_names, data_table.INSTITUTE_SHORTNAME{i}));
end
sites_table  = array2table(sites_data, 'VariableNames', sites_names, 'RowNames', PSN_final);
clear sites_daa sties_names

%% Y data
% Choose studygroups
% PSN_studygroup = PSN_MRI_complete(ismember(data_table.Studygroup, input.selected_studygroups)); % cell
% PSN_studygroup = PSN_MRI_complete(ismember(data_table.Studygroup, input.selected_studygroups)); % cell

selected = {'age', 'sex', 'studygroup'}; % 'age', 'sex', 'studygroup', 'bmi'
% basic_names: cell
basic_names = basic_table.Properties.VariableNames(contains(basic_table.Properties.VariableNames, selected, 'IgnoreCase', true));
% basic_table: table
basic_data = basic_table(PSN_final, basic_names); clear basic_names

% Medication Data 
for m = 1:numel(mednames.medlist_abbr)
    placeholder = [mednames.medlist_abbr{m}, '_', input.timepoint, '_', input.dosage_type];
    med_data(:,m) = data_table.(placeholder); clear placeholder
end
med_data = array2table(cell2mat(med_data), 'VariableNames', mednames.medlist_abbr);

% MRI quality measures
qm_names = ['Cat12_IQR_sMRI_', input.timepoint];
qm_data = mri_table(PSN_final, qm_names);

% Assemble Y_temp
Y_raw = [med_data, basic_data, qm_data];

% Check for outliers in Y data
% Y_outlier = Y_raw(:, ~contains(Y_raw.Properties.VariableNames, 'Labels'));
% TF = isoutlier(Y_outlier.Variables);

%& ADAPT: set imputation threshold
% input.threshold_imputation  = 4/5;
% 
% % take only subjects that have complete data above threshold
% log_complete = sum(~isnan(Y_raw.Variables),2)/size(Y_raw.Variables,2)>input.threshold_imputation;
% PSN_complete = Y_raw.Properties.RowNames(log_complete);
% Y_temp_r = Y_raw(PSN_complete,:);
% 
% % take only features that have complete data above threshold
% check_complete_var = sum(~isnan(Y_temp_r.Variables),1)/size(Y_temp_r.Variables,1);
% 
% % update Y_names with pruned variables
% Y_temp_v = Y_temp_r(:, check_complete_var >= input.threshold_imputation);
% PSN_for_MRI = Y_temp_v.Properties.RowNames;

%% Save data sets in struct
data_complete.foranalysis.basic = basic_table;
data_complete.foranalysis.sites = sites_table;
data_complete.foranalysis.mri = mri_table;
data_complete.foranalysis.medication = med_data;
clear basic_table sites_table mri_table med_data

%% generate GM mask from all available MRIs
% PSN_MRI_final = PSN_cellchar;

mask_name = 'Med_GM_Mask';
GM_mask{1}.spm.tools.masking{1}.makeavg.innames = data_complete.foranalysis.mri.s8rmwp1(PSN_final);
GM_mask{1}.spm.tools.masking{1}.makeavg.avgexpr = 'mean(X)';
GM_mask{1}.spm.tools.masking{1}.makeavg.outname = ['CW_', mask_name, '_', num2str(size(PSN_final,1)), '_', input.timepoint, '_average.nii'];
GM_mask{1}.spm.tools.masking{1}.makeavg.outdir = {'/volume/CW_Med/Data/MRIData/Masks/'};
GM_mask{2}.spm.tools.masking{1}.optthr.inname(1) = cfg_dep('Make Average: Average Image', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
GM_mask{2}.spm.tools.masking{1}.optthr.optfunc = '@opt_thr_antimode';
GM_mask{2}.spm.tools.masking{1}.optthr.outname = ['CW_', mask_name, '_', num2str(size(PSN_final,1)),'_', input.timepoint, '_average_optthr.nii'];
GM_mask{2}.spm.tools.masking{1}.optthr.outdir = {'/volume/CW_Med/Data/MRIData/Masks/'};

spm_jobman('run', GM_mask);
% spm_jobman('interactive', GM_mask);

input.GM_mask_path = [GM_mask{2}.spm.tools.masking{1}.optthr.outdir{1}, '/', GM_mask{2}.spm.tools.masking{1}.optthr.outname];

%% Get vectorized brain data
NM_images = data_complete.foranalysis.mri.s8rmwp1(PSN_final);
NM_TIV = data_complete.foranalysis.mri.TIV(PSN_final);
% nm
% ADAPT
input.NM_structure = ['/volume/projects/CW_Med/Data/MRIData/NM_GMV_', input.timepoint, '_n', num2str(height(PSN_final)), '.mat'];
load(input.NM_structure);
MRI_for_analysis = NM.Y{1,1};
input.MRI_path = [input.NM_structure(1:(strfind(input.NM_structure,'.mat')-1)), '_X.mat'];
save(input.MRI_path, 'MRI_for_analysis');

%% Impute missing Y Data
% ADD UPDATED
[Y_final] = impute_nan(Y_raw, data_table);
input.Y_final = Y_final; clear Y_raw

%% Covariate: Software Version
covariate = data_table.(['dicom_sMRI_SoftwareVersion_sMRI_', input.timepoint, '_num']); 
covariate_name = {'software'};

%% Assign data for analysis in SPLS toolbox and create data and setup structure
% Parameters for IT infrastructure
setup.spls_standalone_path  = '/volume/DP_FEF/ScrFun/ScriptsRepository/SPLS_Toolbox_Dev_2023'; % Path of the SPLS Toolboxsetup.date                  = date; % automatic date
setup.date                  = date; % automatic date
setup.analysis_folder       = ['/data/core-psy/data/pronia/projects/core_spls_test_cv/Analysis/', setup.date]; % analysis folder is created automatically using the current date
setup.queue_name            = 'all.q'; % Choose queue for master and slave jobs, choose between psy0cf20 and mitnvp1-2 and all.q
setup.email                 = 'clara.weyer@med.uni-muenchen.de'; % your email address
setup.max_sim_jobs          = 40; % Define how many parallel jobs are created
setup.parallel_jobs         = 20; % Define how many jobs run in parallel at the same time (soft threshold)
setup.mem_request           = 5; % Memory request for master and slave jobs, in GB, maybe decrease to 2 or 3
setup.matlab_version        = 'R2022a'; % Define the runtime engine, currently its R2020b
setup.cache_path            = '/data/core-psy/data/pronia/opt/temp/clvetter'; % Path for output text files during hyperopt, permutation, bootstrapping => generally same as scratch space
setup.scratch_space         = '/data/core-psy/data/pronia/opt/temp/nmcache/clvetter'; % Path for temporary file storage (hyperopt, permutation, bootstrapping) during analysis, please insert your own folder in the scratch space
setup.compilation_subpath   = 'for_testing'; % default

% Accessory data input
input.final_ID              = input.Y_final.Properties.RowNames; % ID of your subjects (cell array)
input.data_complete         = data_complete; % Optional, not needed for analysis: I store a large table of the entire sample in my files, mostly for post hoc analyses, 
input.type_correction       = 'correct'; % Define whether you want to correct for covariates, choose between: correct, uncorrected
input.correction_target     = 1; % Define whether you want to remove covariate effects from 1) X, 2) Y or 3) both matrices
input.sites                 = input.data_complete.foranalysis.sites{PSN_final,:}; % Dummy coded vector for sites, if only one site, then enter a column vector of ones
input.sites_names           = input.data_complete.foranalysis.sites.Properties.VariableNames; % Names of the sites
input.covariates            = [covariate]; % [input.sites]; % Matrix with covariates (double)
input.covariates_names      = [covariate_name]; % Names of covariates (cell array)
input.Diag                  = input.data_complete.foranalysis.basic.DiagNum; % input.data_complete.foranalysis.basic{PSN_MRI_final, 'DiagNum'}; % Column vector with diagnoses coded via numbers, i.e., [1 3 2 3 4 ]
input.DiagNames             = input.data_complete.foranalysis.basic.Labels; % input.data_complete.foranalysis.basic{input.Y_final.Properties.RowNames, 'Labels'}; % Column cell array with diagnoses/labels, i.e., {'HC', 'ROD', 'CHR', 'HC', 'ROP'}

% define X and Y
input.X                     = input.MRI_path; % 1st data matrix, usually for MRI/biological data, if applicable. If MRI data, then either put in the path to a Matlab file, containing only one variable with vectorized MRI data, or put in vectorized MRI data itself, otherwise just put in the matrix (double format)
input.X_names               = []; % define names of features in X, if MRI data, or no names applicable, leave empty
input.Y                     = Y_final.Variables; % 2nd data matrix, usually for behavioral/phenotypical data (double format) 
input.Y_names               = Y_final.Properties.VariableNames; % define names of features in X, if no names applicable, leave empty
input.subscales             = input.Y_names; % Optional, only for post hoc visualization needed

% Define ML framework
input.framework             = 1; % Cross-validation setup: 1 = nested cross-validation, 2 = random hold-out splits, 3 = LOSOCV, 4 = random split-half
input.outer_folds           = 10; % Applicable only for nested cross-validation and Random Hold-Out Splits: Define Outer folds CV2
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
input.grid_dynamic.LV_1.x   = struct('start', 1, 'end', 0, 'density', input.density); 
input.grid_dynamic.LV_1.y   = struct('start', 1, 'end', 0, 'density', input.density);
% Define grid for hyperparameter search for the X matrix (cu parameter) =>
% 'start': 1 means start is at value 1, 10 means it starts at the lower 10
% percentile of the grid, etc. 
% 'end' defines the upper limit of the hyperparameter search, 0 means all
% the way to the end, 10 means to stop at the upper 10 percentile, etc.,
% 'density' defines the number of data points which are tested during the
% grid, i.e., 20 means that between start and end point, 20 equidistant
% values are tested for the hyperparameter

%% Create analysis datafile

input.name = ['CW_MED_GMV_', input.timepoint, '_', input.approach, '_', num2str(input.outer_folds), 'x', num2str(input.inner_folds), '_', num2str(input.permutation_testing), 'perm_', num2str(input.bootstrap_testing), 'boot_', num2str(input.density), 'density']; %'Immune_only_678_IQRadd_HCcorr_33_noval_min10_2020_1000AUC_1000boot';

% input.name = 'Immune_only_678_IQRadd_HCcorr_55_noval_min10_4040_5000AUC_500boot'; % Choose a filename

input.datafile = [setup.analysis_folder, '/' setup.date, '_', input.name, '_datafile.mat']; % Path for storing datafile containing input and setup

% mkdir(setup.analysis_folder); % Creates analysis folder
save([ setup.date, '_', input.name, '_datafile.mat'], 'setup', 'input'); % Saves datafile in analysis folder

%try rmdir('/volume/mitnvp1_scratch/CW_Med/.mcrCache9.9', 's'); % Cleans your scratch space, be careful: only do this when you have nothing running in the queue
%catch
%    disp('no deletion');
%end

%dp_bash_main_job_slim_addedruntime_mult(input.datafile); % Start the analysis
