%% DP updated script for SPLS visualization 
% now including:
% - Automated results collection
% - bootstrapping and plotting of bootstrapped results
% - adaptive columns for legends, so that it always fits to the number of
% weights
% - selection whether visualizations should be saved to results path or to
% overall folder

directory = '/volume/projects/CV_gs_PLS/Analysis/PRONIA_sites/';
folder_spec = 'PRONIA_sites'; % keyword that is used in all your SPLS analysis, it must be found in all folders that you want to assess
specific_analysis_dates = {'25-Jun-2024_CORE'}; % you can specify the analysis dates (as created by the SPLS toolbox (i.e., 23-Feb-2023), otherwise leave empty as []
% results_collection = dp_find_spls_results_files(directory, folder_spec, specific_analysis_dates);%{'/volume/projects/DP_FEF/FKBP5/Analysis/24-Nov-2023/FKBP5_361_cellcovHC_LSOVNCV55_4040_100BS_5000count/final_results/result.mat'};%
results_collection = {'/volume/projects/CV_gs_PLS/Analysis/PRONIA_sites/25-Jun-2024_CORE/CV_mbspls_SLURM_PRONIA_sites_rsmri_muse_demog_10x10_5000perm_500boot_fro_matrixnorm_randomized_search_1000_iterations/preliminary_result.mat'};%{'/volume/projects/DP_FEF/FKBP5/Analysis/24-Nov-2023/FKBP5_361_cellcovHC_LSOVNCV55_4040_100BS_5000count/final_results/result.mat'};%

save_mode = 'results_path'; % where to save the visualizations: results_path => saves it in the same folder as the result.mat file, overall_folder => saves it in an overall collection folder in the Analysis main directory

for i=1:size(results_collection,1)
    load(results_collection{i});
    k_temp = strfind(results_collection{i},'/','ForceCellOutput',true);
    switch save_mode
        case 'results_path'
            folder_temp = results_collection{i}(1:max(k_temp{1}));
        case 'overall_folder'
            folder_temp = [directory, '/Visualizations/', setup.date, '_', input.name, '/'];
    end
    vis_folder = [folder_temp, 'vis'];
    mkdir(vis_folder);
    selection_vectors = {1:3};
    plot_all = true;
    additional_plotting = false;
    to_plot_collection = [];
    val_log = false;

    dp_vis_general(vis_folder, results_collection{i}, selection_vectors, plot_all, additional_plotting, to_plot_collection, val_log);
    
    weights_collection =[];lv_names=[];
    for r=1:size(output.final_parameters,1)
        for n_mat = 1:size(input.Xs,2)
        weights_collection{n_mat}(r,:) = output.final_parameters{r, matches(output.parameters_names, 'weights')}{n_mat};
        end
        lv_names{r,1} = ['LV_', num2str(r)];
        
    end
    u_table = array2table(u_collection, 'RowNames', lv_names, 'VariableNames', input.X_names);
    v_table = array2table(v_collection, 'RowNames', lv_names, 'VariableNames', input.Y_names);
    writetable(u_table, [vis_folder, '/LV_results.xlsx'], 'WriteRowNames', true, 'Sheet', 'u_vector');
    writetable(v_table, [vis_folder, '/LV_results.xlsx'], 'WriteRowNames', true, 'Sheet', 'v_vector');
    
    switch input.type_correction
        case 'correct'
            correct_log = true;
        case 'uncorrected'
            correct_log = false;
    end

    latent_scores_table = dp_get_latent_scores(results_collection{i}, 'all', correct_log);
    writetable(latent_scores_table, [vis_folder, '/LV_results.xlsx'], 'WriteRowNames', true, 'Sheet', 'latent_scores_table');

    boot_options = {'CI', 'BS'};
    for ii=1:size(boot_options,2)
        boot_results_file = dp_bootstrap_pruning(results_collection{i}, boot_options{ii});
        load(boot_results_file);
        boot_vis_folder = [folder_temp, 'boot_vis'];
        mkdir(boot_vis_folder);
        dp_vis_general(boot_vis_folder, boot_results_file, selection_vectors, plot_all, additional_plotting, to_plot_collection, val_log);

        u_collection = []; v_collection = []; lv_names=[];
        for r=1:size(output.final_parameters,1)
            u_collection(r,:) = output.final_parameters{r, matches(output.parameters_names, 'u')};
            v_collection(r,:) = output.final_parameters{r, matches(output.parameters_names, 'v')};
            lv_names{r,1} = ['LV_', num2str(r)];
        end
        u_table = array2table(u_collection, 'RowNames', lv_names, 'VariableNames', input.X_names);
        v_table = array2table(v_collection, 'RowNames', lv_names, 'VariableNames', input.Y_names);
        writetable(u_table, [boot_vis_folder, '/LV_results.xlsx'], 'WriteRowNames', true, 'Sheet', 'u_vector');
        writetable(v_table, [boot_vis_folder, '/LV_results.xlsx'], 'WriteRowNames', true, 'Sheet', 'v_vector');

    boot_latent_scores_table = dp_get_latent_scores(boot_results_file, 'all', correct_log);
    writetable(boot_latent_scores_table, [boot_vis_folder, '/LV_results.xlsx'], 'WriteRowNames', true, 'Sheet', ['latent_scores_table', '_', boot_options{ii}]);

    end
end

