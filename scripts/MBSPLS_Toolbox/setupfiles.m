%% testing

spls_standalone_path = '/volume/DP_FEF/ScrFun/ScriptsRepository/SPLS_Toolbox_standalone_BO_5k_sets_FDR';
analysis_folder = '/volume/HCStress/Analysis/12-Jul-2018/DP_CISS_RSA_allgroups_correctedforsites_1';
queue_name = 'psy0cf20';
variables_for_analysis = '/volume/HCStress/Analysis/06-Jul-2018/DP_CISS_RSA_allgroups_T0_80P_ASincluded.mat';
covariates_file = '/volume/HCStress/Analysis/06-Jul-2018/DP_CISS_RSA_allgroups_T0_80P_sites.mat';

dp_spls_standalone_BO_5k_sets_FDR(spls_standalone_path, analysis_folder, queue_name_perm, variables_for_analysis)