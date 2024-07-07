%% correlation analyses between phenotypic features and GMV
addpath(genpath('/volume/projects/DP_FEF/ScrFun/ScriptsRepository/'));
analysis_folder = '/volume/projects/ST_Trauma_MDD/Analysis/30-May-2022/MDD_singleitems_633_IQRadd_NCV55_noval_min10_4040_5000AUC_boot/final_results/Phen_voxels_corr/voxels_y/'; %Output folder

% load results
load('/volume/projects/ST_Trauma_MDD/Analysis/30-May-2022/MDD_singleitems_633_IQRadd_NCV55_noval_min10_4040_5000AUC_boot/final_results/result_BS_final_vis_final_vis.mat');
% load cortical/cerebellar atlases
temp = load('/volume/HCStress/Data/MRI/Atlases/Brainnetome_Atlas/brainnetome_3mm_633_MDD_Trauma_NM_X.mat');
fields=fieldnames(temp);
brain_atlas = temp.(fields{1});

temp = load('/volume/HCStress/Data/MRI/Atlases/Cerebellum-MNIflirt-MRICroN/Cerebellum-MNIflirt_633_MDD_Trauma_X.mat');
fields=fieldnames(temp);
cerebellum_atlas = temp.(fields{1});

temp = load('/volume/HCStress/Data/MRI/Atlases/Brainnetome_Atlas/BNA_table.mat');
fields=fieldnames(temp);
BNA_table = temp.(fields{1});

temp = load(input.X);
fields=fieldnames(temp);
X = temp.(fields{1});
Y = input.Y;

% choose whether you want to use brain regions or just single voxels for
% correlation
correlation_setup = 1; % 1 = single voxel with highest salience, 2 = brain regions with highest percentage

% choose whether you want deflation steps after first LV
deflation_setup = 2; % 1 = use deflation, 2 = do not use deflation

% choose whether rescaling shall be activated
correct_setup = 2; % 1 = use correcting, 2 = do not use correcting
cs_method.correction_subgroup = input.cs_method.correction_subgroup;
if ~isempty(cs_method.correction_subgroup)
    cs_method.subgroup_train = contains(input.DiagNames, cs_method.correction_subgroup);
    cs_method.subgroup_test = contains(input.DiagNames, cs_method.correction_subgroup);
else
    cs_method.subgroup_train = [];
    cs_method.subgroup_test = [];
end

% choose scaling
scaling_setup = 2; % 1 = use rescaling, 2 = do not use rescaling
cs_method.method      = 'min_max'; % Scaling of features, default: mean-centering, also possible 'min_max' (scaling from 0 to 1) => preferred scaling is mean-centering!

% define phenotypic items
variables_selected =      {'age', 'BDI2_20', 'CTQ_09', 'female_sex', 'SPI_A_A2_1_1_red', 'male_sex'};

% define brain regions (used only when correlation_setup is set at 2 =
% brain regions, if correlation_setup is set at 1, this will become
% arbitrary since the loop finds the voxels itself and locates it in the
% right region and hemisphere
regions_selected =        {'MTG', 'Hipp', 'MTG', 'pSTS','Amyg', 'BG'};
hemispheres_selected =    {'right', 'right', 'right', 'right', 'left', 'left'};

% Use nicer colors
nice_blue = [0 0.4470 0.7410];
nice_red = [0.8500 0.3250 0.0980];

log_u = matches(output.parameters_names, 'u');
log_v = matches(output.parameters_names, 'v');

for i=1:(size(output.final_parameters,1)-1)
    
    switch correlation_setup
        case 1
            u = output.final_parameters{i,log_u};
            [max_value, index_region_extract] = max(abs(u));
            region_extract = zeros(size(u))';
            region_extract(index_region_extract) = 1;
            
            region_number_temp = brain_atlas(index_region_extract);
            
            if region_number_temp == 0
                region_number_temp = cerebellum_atlas(index_region_extract);
                
            else
                hemispheres_names = {'left', 'right'};
                [row_temp, column_temp] = find(BNA_table{:, hemispheres_names} == region_number_temp);
                regions_selected{1,i} = BNA_table{row_temp, 'regions'}{1};
                hemispheres_selected{1,i} = hemispheres_names{column_temp};
            end
            
            
        case 2
            region_indices = BNA_table{contains(BNA_table.regions, regions_selected{i}, 'IgnoreCase', false), hemispheres_selected{i}};
            region_extract = ismember(round(brain_atlas), region_indices);
            
            
    end
    
    if i==1
        switch correct_setup % perform correction and scaling
            case 1 % with site correction
                COV = input.covariates;
            case 2 % without site correction
                COV = nan(size(input.covariates,1),1);
        end
    else
        COV = nan(size(input.covariates,1),1);
    end
    
    switch scaling_setup % perform scaling
        case 1 % with scaling
            X = dp_correctscale(X,COV, cs_method);
            Y = dp_correctscale(Y,COV, cs_method);
    end
    
    switch deflation_setup
        case 1
            if i > 1
                % deflation step
                u = output.final_parameters{i-1, 4};
                v = output.final_parameters{i-1, 5};
                [X,Y] = proj_def(X, Y, u, v);
                switch scaling_setup % perform scaling
                    case 1 % with scaling
                        X = dp_correctscale(X,COV, cs_method);
                        Y = dp_correctscale(Y,COV, cs_method);
                end
            end
            
    end
    
    brain_volumes = X*region_extract';
    log_y = matches(input.Y_names, variables_selected{i});
    x = brain_volumes;
    y = Y(:, log_y);
    if dp_binary_check(y)
        [p(i), h, stats] = ranksum(x(y==0),x(y==1));
        boxplot(x,y);
        title(['LV', num2str(i), ', Ranksum = ', num2str(stats.zval), ', P value = ', num2str(p(i))]); % add third line
    else
        [rho(i), p(i)] = corr(x,y, 'type', 'Spearman');
        if rho(i) < 0
            scatter(y,x, 'filled', 'blue');
        else
            scatter(y,x, 'filled', 'red');
        end
        lsline
        title(['LV', num2str(i), ', Spearman''s Rho = ', num2str(rho(i)), ', P value = ', num2str(p(i))]); % add third line
    end
    
    fontsize = 8;
    switch correlation_setup
        case 1
            ylabel(['voxel in ', num2str(index_region_extract), ' ', regions_selected{i}, ' ', hemispheres_selected{i}], 'FontSize', fontsize);
        case 2
            ylabel([regions_selected{i}, ' ', hemispheres_selected{i}], 'FontSize', fontsize);
    end
    xlabel(strrep(variables_selected{i}, '_', ' '), 'FontSize', fontsize);
    set(gcf,'Position', get(0,'Screensize'));
    set(gcf,'PaperPositionMode','auto')
    print(gcf, [analysis_folder, '/LV_' num2str(i), '_GMV_phenotype_corr'], '-dpng', '-r0');
    saveas(gcf, [analysis_folder, '/LV_' num2str(i), '_GMV_phenotype_corr.fig']);
    saveas(gcf,[analysis_folder, '/LV_' num2str(i), '_GMV_phenotype_corr'],'epsc');
    close all
end