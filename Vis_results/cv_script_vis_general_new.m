%% DP updated script for SPLS visualization 
% now including:
% - Automated results collection
% - bootstrapping and plotting of bootstrapped results
% - adaptive columns for legends, so that it always fits to the number of
% weights
% - selection whether visualizations should be saved to results path or to
% overall folder

addpath(genpath('/volume/projects/CV_gs_PLS/ScrFun/multiblock_spls'))



directory = '/volume/projects/CV_gs_PLS/Analysis/PRONIA_sites/';
folder_spec = 'PLEXUS_sites'; % keyword that is used in all your SPLS analysis, it must be found in all folders that you want to assess
specific_analysis_dates = {'25-Jun-2024-CORE'}; % you can specify the analysis dates (as created by the SPLS toolbox (i.e., 23-Feb-2023), otherwise leave empty as []
results_collection = {'/volume/projects/CV_gs_PLS/Analysis/PLEXUS_new_segmentation/16-Jul-2024/CV_mbspls_SLURM_PLEXUS_new_segmentation_plexus_mri_clinical_5x5_5000perm_500boot_fro_matrixnorm_grid_search_10_10_10_densities_2/final_results/result.mat'}%, ...
 %results_collection ={'/volume/projects/CV_gs_PLS/Analysis/PRONIA_sites/25-Jun-2024_CORE/CV_mbspls_SLURM_PRONIA_sites_rsmri_muse_10x10_5000perm_500boot_0_matrixnorm_randomized_search_1000_iterations/final_results/preliminary_result.mat'};

save_mode = 'results_path'; % where to save the visualizations: results_path => saves it in the same folder as the result.mat file, overall_folder => saves it in an overall collection folder in the Analysis main directory

for i=1:size(results_collection,2)
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
    selection_vectors = 1:size(input.Xs,2);
    plot_all = false;
    additional_plotting = false;
    to_plot_collection = [];
    val_log = false;

    %cv_vis_general(vis_folder, results_collection{i}, selection_vectors, plot_all, additional_plotting, to_plot_collection, val_log);
    
    weights_collection = []; [];lv_names=[];
    for r=1:size(output.final_parameters,1)
        for num_m =1:size(input.Xs,2)
            weights_collection{num_m}(r,:) = output.final_parameters{r, matches(output.parameters_names, 'weights')}{num_m};
            lv_names{r,1} = ['LV_', num2str(r)];
            weights_table = array2table(weights_collection{num_m}, 'RowNames', lv_names, 'VariableNames', input.Xs_feature_names{num_m});
            writetable(weights_table, [vis_folder, '/LV_results.xlsx'], 'WriteRowNames', true, 'Sheet', ['weights_vector_', num_m]);
 
        end
           

    end
        
    switch input.type_correction
        case 'correct'
            correct_log = true;
        case 'uncorrect'
            correct_log = false;
    end

    latent_scores_table = cv_get_latent_scores(results_collection{i}, 'all', correct_log, flip_flag);
    %writetable(latent_scores_table, [vis_folder, '/LV_results.xlsx'], 'WriteRowNames', true, 'Sheet', 'latent_scores_table');
%%
    boot_options = {'CI', 'BS'};
    for ii=1:size(boot_options,2)
        boot_results_file = cv_mbspls_bootstrap_pruning(results_collection{i}, boot_options{ii});
        load(boot_results_file);
        boot_vis_folder = [folder_temp, 'boot_vis_', boot_options{ii}];
        mkdir(boot_vis_folder);
        
        %cv_vis_general(boot_vis_folder, boot_results_file, selection_vectors, plot_all, additional_plotting, to_plot_collection, val_log);

        weights_collection = []; [];lv_names=[];
        for r=1:size(output.final_parameters,1)
            for num_m =1:size(input.Xs,2)
                weights_collection{num_m}(r,:) = output.final_parameters{r, matches(output.parameters_names, 'weights')}{num_m};
                lv_names{r,1} = ['LV_', num2str(r)];
                weights_table = array2table(weights_collection{num_m}, 'RowNames', lv_names, 'VariableNames', input.Xs_feature_names{num_m});
                writetable(weights_table, [boot_vis_folder, '/LV_results_' boot_options{ii} '.xlsx'], 'WriteRowNames', true, 'Sheet', ['weights_vector_', num2str(num_m)]);
                writetable(weights_table, [boot_vis_folder, '/LV_results_' boot_options{ii} '_weights_vector_', num2str(num_m), '.csv'], 'WriteRowNames', true);

            end
            

        end

    %boot_latent_scores_table = cv_get_latent_scores(boot_results_file, 'all', correct_log);
    %writetable(boot_latent_scores_table, [boot_vis_folder, '/LV_results.xlsx'], 'WriteRowNames', true, 'Sheet', 'latent_scores_table');

                writetable(weights_table, [boot_vis_folder, '/LV_results_latent_scores_table.csv'], 'WriteRowNames', true);



    end
end

