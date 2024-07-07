%% function to visualize SPLS output

function dp_visualize_data_Devmod(IN)

load(IN.results_path);

switch IN.overall_analysis
    case 'Stress'
        overall_folder = '/volume/HCStress/Analysis/Stress';
    case 'Resilience'
        overall_folder = '/volume/HCStress/Analysis/Resilience';
    case 'Schizotypy'
        overall_folder = '/volume/MU_Pronia_PLS_Schizotypy/Analysis/SPLS/Schizotypy';
    case 'Mimics'
        overall_folder = '/volume/DP_Mimics/Analysis/Mimics';
    case 'immune'
        overall_folder = '/volume/RU_DP_immune/Analysis/immune';
    case 'MDD_Trauma'
        overall_folder = '/volume/projects/ST_Trauma_MDD/Analysis/MDD_Trauma';
    case 'FF'
        overall_folder = '/volume/projects/DP_JL_Forensics/Analysis/FF';
    case 'COPE'
        overall_folder = '/volume/projects/ES_become/ES_COPE/Analysis/SPLS/COPE';
end

folder_name = [setup.date, '_', input.name];
if contains(IN.results_path, 'final_vis')
    collection_folder = [overall_folder, '/', folder_name, '/final_vis'];
    val_log = false;
elseif contains(IN.results_path, 'validation_vis')
    collection_folder = [overall_folder, '/', folder_name, '/validation_vis'];
    val_log = true;
end

mkdir(collection_folder);

if ~isfield(IN, 'specific')
    all_jobs = true;
    IN.specific = 'empty';
else
    all_jobs = false;
end

marker_size = 14;
font_size = 14;
LV_x = 'clinical features';
LV_y = 'weights';
LS_epsilon = 'brain score';
LS_omega = 'behavior score';

%% get min and max voxels
output.min_max_table = dp_min_max_voxels(output.final_parameters);

%% compute sociodemographic and clinical outcome correlations
if any(strcmp(IN.specific, 'sociodemographic') | all_jobs)
    [input, output, setup] = dp_sociodemographic_2020_full(IN);
end

if any(strcmp(IN.specific, 'post_hoc') | all_jobs)
    p_folder = [collection_folder, '/post_hoc_correlation'];
    mkdir(p_folder);
%     load('/volume/RU_DP_immune/Analysis/24-Apr-2021/Immune_only_678_IQRadd_NCV55_noval_min10_2020_5000AUC_Dev/final_results/result_final_vis.mat');
    group_select = 'all'; % all OR fold
    correct = 1; % 1=yes, 2=no
    output.RHO_all=[];
    epsilon=[]; omega=[]; latent_scores=[];latent_scores_names=[];
    switch group_select
        case 'all'
            temp = load(input.MRI_path);
            name_temp = fieldnames(temp);
            X = temp.(name_temp{1});
            Y = input.Y;
            
            for i=1:(size(output.final_parameters,1)-1)
                
                IN_x.train = X;
                IN_x.test = X;
                IN_y.train = Y;
                IN_y.test = Y;
                
                switch correct
                    
                    case 1 % with site correction
                        if i==1
                            COV.train = input.sites;
                            COV.test = input.sites;
                        else
                            COV.train = nan(size(input.sites,1),1);
                            COV.test = nan(size(input.sites,1),1);
                            input.correction_target = 3;
                        end
                    case 2 % without site correction
                        COV.train = nan(size(input.sites,1),1);
                        COV.test = nan(size(input.sites,1),1);
                        input.correction_target = 3;
                end
                
                [OUT_x, OUT_y] = dp_master_correctscale(IN_x, IN_y, COV, input.scaling_method, input.correction_target);
                
                u = output.final_parameters{i,4};
                v = output.final_parameters{i,5};
                
                [output.RHO_all(i), epsilon, omega, u, v] = dp_projection(OUT_x.train, OUT_y.train, u, v, input.correlation_method);
                
                [X,Y] = proj_def(X, Y, u, v);
                latent_scores = [latent_scores, omega, epsilon];
                latent_scores_names = [latent_scores_names, {['omega_LV', num2str(i)], ['epsilon_LV', num2str(i)]}];
            end
            
            output.latent_scores_table = array2table(latent_scores, 'RowNames', input.final_ID, 'VariableNames', latent_scores_names);
            
        case 'fold'
            
            brain_scores_LV3 = epsilon{3};
            brain_scores_LV4 = epsilon{4};
            output.latent_scores_table = array2table([brain_scores_LV3, brain_scores_LV4], 'RowNames', input.final_ID, 'VariableNames', {'epsilon_LV3', 'epsilon_LV4'});
    end
    
    correlation_names=[];
    
    if any(contains(IN.SD_selection, 'GAF'))
        try temp_names = input.data_complete.foranalysis.additional.Properties.VariableNames(contains(input.data_complete.foranalysis.additional.Properties.VariableNames, 'TGF', 'IgnoreCase', false));
            input.data_complete.foranalysis.additional(:, temp_names)=[];
        end
        correlation_names = [correlation_names, input.data_complete.foranalysis.additional.Properties.VariableNames(contains(input.data_complete.foranalysis.additional.Properties.VariableNames, {'GAF', 'GF'}, 'IgnoreCase', false))];
    end
    
    if any(contains(IN.SD_selection, 'BDI'))
        correlation_names = [correlation_names, input.data_complete.foranalysis.additional.Properties.VariableNames(contains(input.data_complete.foranalysis.additional.Properties.VariableNames, 'BDI', 'IgnoreCase', false))];
    end
    
    if any(contains(IN.SD_selection, 'NEO'))
        correlation_names = [correlation_names, input.data_complete.foranalysis.additional.Properties.VariableNames(contains(input.data_complete.foranalysis.additional.Properties.VariableNames, 'NEO', 'IgnoreCase', false))];
    end
    
    if any(contains(IN.SD_selection, 'WHO'))
        correlation_names = [correlation_names, input.data_complete.foranalysis.additional.Properties.VariableNames(contains(input.data_complete.foranalysis.additional.Properties.VariableNames, 'WHO', 'IgnoreCase', false))];
    end
    
    if any(contains(IN.SD_selection, 'PANSS'))
        correlation_names = [correlation_names, input.data_complete.foranalysis.additional.Properties.VariableNames(contains(input.data_complete.foranalysis.additional.Properties.VariableNames, 'PANSS', 'IgnoreCase', false))];
    end
    
    if any(contains(IN.SD_selection, 'COGDIS'))
        load('/volume/RU_DP_immune/Data/Immune_megapaper_request_Disc/Immune_megapaper_request_Disc/DATA/18-Dec-2020/Immune_megapaper_request_Disc_Data_all_18-Dec-2020.mat');
        temp_data_table = data_table_all;
        
        load('/volume/RU_DP_immune/Data/Immune_megapaper_request_Repl/Immune_megapaper_request_Repl/DATA/18-Dec-2020/Immune_megapaper_request_Repl_Data_all_18-Dec-2020.mat');
        temp_names = temp_data_table.Properties.VariableNames(matches(temp_data_table.Properties.VariableNames, data_table_all.Properties.VariableNames));
        data_table_all = [temp_data_table(:, temp_names); data_table_all(:, temp_names)];
        cogdis_items_names = {'SPI_A_COGDIS_B1_1', 'SPI_A_COGDIS_C2_1', 'SPI_A_COGDIS_C3_1', 'SPI_A_COGDIS_C4_1', 'SPI_A_COGDIS_C5_1',...
            'SPI_A_COGDIS_D3_1', 'SPI_A_COGDIS_D4_1', 'SPI_A_COGDIS_O3_1', 'SPI_A_COGDIS_O7_1'};
        cogdis_names = data_table_all.Properties.VariableNames(contains(data_table_all.Properties.VariableNames, cogdis_items_names, 'IgnoreCase', true) & contains(data_table_all.Properties.VariableNames, 'red', 'IgnoreCase', true));
        
        input.data_complete.foranalysis.cogdis_data = data_table_all(:, cogdis_names);
        input.data_complete.foranalysis.cogdis_data.cogdis_fulfilled = sum(input.data_complete.foranalysis.cogdis_data.Variables>=3,2)>=2;
        input.data_complete.foranalysis.cogdis_data.Properties.RowNames = input.data_complete.foranalysis.basic.Properties.RowNames;
        correlation_names = [correlation_names, 'cogdis_fulfilled'];
    end
    
    if any(contains(IN.SD_selection, 'neurocognition'))
        % get extraction data
        load('/volume/RU_DP_immune/Data/Immune_megapaper_request_Disc/Immune_megapaper_request_Disc/DATA/18-Dec-2020/Immune_megapaper_request_Disc_Data_all_18-Dec-2020.mat');
        temp_data_table = data_table_all;
        
        load('/volume/RU_DP_immune/Data/Immune_megapaper_request_Repl/Immune_megapaper_request_Repl/DATA/18-Dec-2020/Immune_megapaper_request_Repl_Data_all_18-Dec-2020.mat');
        temp_names = temp_data_table.Properties.VariableNames(matches(temp_data_table.Properties.VariableNames, data_table_all.Properties.VariableNames));
        data_table_all = [temp_data_table(:, temp_names); data_table_all(:, temp_names)];
        
        NC = rs_neurocognition(data_table_all);
        input.data_complete.foranalysis.neurocognition=table;
        for i=1:size(input.data_complete.foranalysis.basic.Properties.RowNames,1)
            log_temp = str2num(cell2mat(data_table_all.PSN)) == str2num(input.data_complete.foranalysis.basic.Properties.RowNames{i,1});
            if sum(log_temp)~=0
                input.data_complete.foranalysis.neurocognition(i,:) = NC.single_scores(log_temp, :);
            else
                input.data_complete.foranalysis.neurocognition(i,:) = array2table(nan(1, size(NC.single_scores,2)), 'VariableNames', NC.single_scores.Properties.VariableNames);
            end
        end
        input.data_complete.foranalysis.neurocognition.Properties.RowNames = input.data_complete.foranalysis.basic.Properties.RowNames;
        correlation_names = [correlation_names, input.data_complete.foranalysis.neurocognition.Properties.VariableNames];
    end
    
    if any(contains(IN.SD_selection, 'CTQ'))
        ctq_names = sort(input.data_complete.foranalysis.CTQ.Properties.VariableNames);
        ctq_data_old = input.data_complete.foranalysis.CTQ.Variables;
        CTQ.emotional_abuse = [3,8,14,18,25];
        CTQ.physical_abuse = [9,11,12,15,17];
        CTQ.sexual_abuse = [20,21,23,24,27];
        CTQ.emotional_neglect = [5,7,13,19,28];
        CTQ.physical_neglect = [1,2,4,6,26];
        CTQ.denial = [10,16,22];
        ctq_data = [sum(ctq_data_old(:, CTQ.emotional_abuse),2), sum(ctq_data_old(:, CTQ.physical_abuse),2), sum(ctq_data_old(:, CTQ.sexual_abuse),2), sum(ctq_data_old(:, CTQ.emotional_neglect),2), sum(ctq_data_old(:, CTQ.physical_neglect),2), sum(ctq_data_old(:, CTQ.denial),2)];
        input.data_complete.foranalysis.CTQ{:, fieldnames(CTQ)'} = ctq_data;
        correlation_names = [correlation_names, fieldnames(CTQ)'];
    end
    
    if any(contains(IN.SD_selection, 'LOS'))
        correlation_names = [correlation_names, 'Length_of_Storage_Days'];
    end
    
    collect_table = dp_find_tables(input.data_complete.foranalysis, correlation_names, input.final_ID);
    
    [RHO, P] = corr(collect_table.Variables, output.latent_scores_table{input.final_ID,:}, 'rows', 'pairwise', 'type', 'Spearman');
    output.post_hoc_correlations.RHO_table = array2table(RHO, 'RowNames', collect_table.Properties.VariableNames, 'VariableNames', output.latent_scores_table.Properties.VariableNames);
    output.post_hoc_correlations.P_table = array2table(P, 'RowNames', collect_table.Properties.VariableNames, 'VariableNames', output.latent_scores_table.Properties.VariableNames);
    p_adj = dp_FDR_adj(output.post_hoc_correlations.P_table.Variables);
    output.post_hoc_correlations.P_table_adj = output.post_hoc_correlations.P_table;
    output.post_hoc_correlations.P_table_adj.Variables = p_adj;
                    
    writetable(output.post_hoc_correlations.RHO_table, [p_folder, '/post_hoc_correlations.xlsx'], 'WriteRowNames', true, 'Sheet', 'RHO');
    writetable(output.post_hoc_correlations.P_table_adj, [p_folder, '/post_hoc_correlations.xlsx'], 'WriteRowNames', true, 'Sheet', 'P');
end


if val_log
    output.final_parameters = output.validation_results;
end

input.Y_names = strrep(input.Y_names, '_T0', '');
load(input.NM_structure);

%% prepare questionnaire data

% define column names for matrices so that you can access them later by
% indexing
% output.parameters_names = {'w', 'cu', 'cv', 'u', 'v', 'U_opt', 'S_opt' ,'V_opt', 'success', 'RHO', 'p'}; % names of parameters for optimal cu/cv combination of the w-th loop
opt_parameters_names = {'w', 'cu', 'cv', 'u', 'v', 'success', 'RHO', 'p'};
opt_parameters_names_long = {'w', 'cu', 'cv', 'u', 'v', 'U_opt', 'S_opt' ,'V_opt', 'success', 'RHO', 'p'}; % names of parameters for optimal cu/cv combination of the w-th loop

% indices for later use
opt_u = strcmp(output.parameters_names,'u');
opt_v = strcmp(output.parameters_names,'v');
opt_p = strcmp(output.parameters_names, 'p');
opt_RHO = strcmp(output.parameters_names, 'RHO');
index_epsilon = strcmp(output.parameters_names, 'epsilon');
index_omega = strcmp(output.parameters_names, 'omega');
index_epsilon_all = strcmp(output.parameters_names, 'epsilon_all');
index_omega_all = strcmp(output.parameters_names, 'omega_all');

if any(strcmp(IN.specific, 'images') | all_jobs)
    i_folder = [collection_folder, '/images'];
    mkdir(i_folder);
    cd(i_folder);
    % write brain vector to nifti file
    for i=1:size(output.final_parameters,1)
        nk_WriteVol(output.final_parameters{i,4}, ['brain_LV_final_' num2str(i)], 2, NM.brainmask{1,1}, cell2mat(NM.badcoords), 0, 'gt');
    end
end

if any(strcmp(IN.specific, 'atlas') | all_jobs)
    
    a_folder = [collection_folder, '/atlases'];
    mkdir(a_folder);
    writetable(output.min_max_table, [a_folder, '/min_max.xlsx'], 'WriteRowNames',true);
    
    IN.sample_size = size(input.Y_final,1);
    if contains(IN.overall_analysis, 'immune')
        IN.var_names = {'immune'};
    elseif contains(IN.overall_analysis, 'MDD')
        IN.var_names = {'MDD_Trauma'};
    else
        IN.var_names = input.analysis_identifier;
    end
    IN.atlases = dp_atlas_find(IN);
    
    IN.voxel_cutoff = 0;
    
    IN.analysis_name = input.name;
    
    for i=1:size(output.final_parameters,1)
        IN.vectors(i,:) = output.final_parameters{i,4};
    end
    
    IN.a_folder = a_folder;
    [output.atlas_readouts, output.atlas_readouts_filled] = dp_atlas_readout_new(IN);
    
    IN.atlas_readouts = output.atlas_readouts;
    output.atlas_table_readouts = dp_atlas_table_readout(IN);
    
    spider_folder = [collection_folder, '/spider_plots'];
    mkdir(spider_folder);
    save(IN.results_path, 'input', 'setup', 'output');
    [output.Table_Regions_Lateralized, output.Table_Regions_Vermis] = dp_atlases_vis_spiderplot(IN.results_path, spider_folder);
    writetable(output.Table_Regions_Lateralized, [spider_folder, '/spider_tables.xlsx'], 'Sheet', 'Table_Regions_Lateralized', 'WriteRowNames',true);
    writetable(output.Table_Regions_Vermis, [spider_folder, '/spider_tables.xlsx'], 'Sheet', 'Table_Regions_Vermis', 'WriteRowNames',true);
    save(IN.results_path, 'input', 'output', 'setup');
end

if contains(IN.overall_analysis, 'immune')
    try input = rmfield(input, 'selected_features');
    end
end

if any(strcmp(IN.specific, 'behavior') | all_jobs)
    
    b_folder = [collection_folder, '/behavior'];
    mkdir(b_folder);
    
    switch IN.analysis_origin
        case 1 % classical PRONIA analysis
            
            load('/volume/HCStress/Doc/all_questionnaires.mat', 'BS_questionnaire', 'CISS_questionnaire', 'CTQ_questionnaire', 'RSA_questionnaire', 'WSS_questionnaire');
            all_questionnaires = [CISS_questionnaire', RSA_questionnaire', CTQ_questionnaire', BS_questionnaire', WSS_questionnaire'];
            
            %% visualize behavior vector in barplot
            % CTQ.emotional_abuse)), (CTQ.physical_abuse)), (CTQ.sexual_abuse)), (CTQ.emotional_neglect)), (CTQ.physical_neglect)), (CTQ.denial))];
            % load([setup.analysis_folder, '/' setup.date, '_', name, '_data_collection.mat']);
            %     if sum(ismember(input.Y_names, 'Age')) + sum(ismember(input.Y_names, 'Sex'))==2
            %         selected_variables = [input.selected_features(1,:), 'Age', 'Sex'; input.selected_features(2,:), 1, 1];
            %     elseif sum(ismember(input.Y_names, 'Age'))
            %         selected_variables = [input.selected_features(1,:), 'Age'; input.selected_features(2,:), 1];
            %     elseif sum(ismember(input.Y_names, 'Sex'))
            %         selected_variables = [input.selected_features(1,:), 'Sex'; input.selected_features(2,:), 1];
            %     else
            %         selected_variables = input.selected_features;
            %     end
            
            log_f=[];
            if isfield(input, 'selected_features')
                for i=1:size(input.selected_features,2)
                    log_f(i,:) = contains(input.Y_names(1,:), input.selected_features{1,i});
                end
                
                if size(log_f,1)>1
                    log_c = sum(sum(log_f==0)==size(input.selected_features,2));
                else
                    log_c = sum((log_f==0)==size(input.selected_features,2));
                end
                
                selected_variables = [[input.selected_features(1,:), input.Y_names((end-(log_c-1)):end)]; [input.selected_features(2,:), num2cell(ones(1,log_c))]];
                
                count=0;
                for i=1:size(input.selected_features(2,:),2)
                    temp=size(fieldnames(input.selected_features{2,i}),1);
                    count=count+temp;
                end
                colorpattern_LV = colorcube((count + log_c));
                
                % collect items from vectors
                for i=1:(size(output.final_parameters,1))
                    x = output.final_parameters{i,opt_v};
                    output.questions_collection.(['LV_', num2str(i)]).items = input.Y_names(x~=0)';
                    output.questions_collection.(['LV_', num2str(i)]).subscales = input.subscales(x~=0)';
                    for q=1:size(output.questions_collection.(['LV_', num2str(i)]).items,1)
                        output.questions_collection.(['LV_', num2str(i)]).questions{q} = all_questionnaires(1,strcmp(all_questionnaires(2,:), output.questions_collection.(['LV_', num2str(i)]).items{q}));
                    end
                    
                    %     errors = output.CI_v{i};
                    f=figure();
                    nn=0; cc=0;
                    hold on
                    temp_legend=[]; temp_all_x = 0; temp_all_bars = 0; temp_handle=[]; %sel_temp_names=[];
                    for ii=1:(size(selected_variables,2))
                        switch class(selected_variables{2,ii})
                            case 'struct'
                                fields = fieldnames(input.selected_features{2,(strcmp(input.selected_features(1,:),input.selected_features{1,ii}))});
                                for iii=1:size(fields,1)
                                    temp_current=size(selected_variables{2,ii}.(fields{iii}),2);
                                    if IN.plot_all
                                        nn=nn+1; cc=cc+1;
                                        temp_handle(nn) = bar((temp_all_bars+1):(temp_all_bars+temp_current),x((temp_all_x+1):(temp_all_x+temp_current)), 'FaceColor', colorpattern_LV(cc,:));
                                        temp_all_bars = temp_all_bars+temp_current;
                                        temp_legend{nn} = [selected_variables{1,ii}, ' ', strrep(fields{iii}, '_', ' ')];
                                        hold on
                                    else
                                        log_nonzero = x((temp_all_x+1):(temp_all_x+temp_current))~=0;
                                        if any(log_nonzero)
                                            nn=nn+1; cc=cc+1;
                                            x_temp = x((temp_all_x+1):(temp_all_x+temp_current));
                                            temp_handle(nn) = bar((temp_all_bars+1):(temp_all_bars+sum(log_nonzero)), x_temp(log_nonzero), 'FaceColor', colorpattern_LV(cc,:));
                                            temp_all_bars = temp_all_bars+sum(log_nonzero);
                                            temp_legend{nn} = [selected_variables{1,ii}, ' ', strrep(fields{iii}, '_', ' ')];
                                            hold on
                                        else
                                            cc=cc+1;
                                        end
                                    end
                                    temp_all_x = temp_all_x+temp_current;
                                end
                            case 'double'
                                temp_current=size(selected_variables{2,ii},2);
                                if IN.plot_all
                                    nn=nn+1; cc=cc+1;
                                    temp_handle(nn) = bar((temp_all_bars+1):(temp_all_bars+temp_current),x((temp_all_x+1):(temp_all_x+temp_current)), 'FaceColor', colorpattern_LV(cc,:));
                                    temp_all_bars = temp_all_bars+temp_current;
                                    temp_legend{nn} = strrep(selected_variables{1,ii}, '_', ' ');
                                    hold on
                                else
                                    log_nonzero = x((temp_all_x+1):(temp_all_x+temp_current))~=0;
                                    if any(log_nonzero)
                                        nn=nn+1; cc=cc+1;
                                        x_temp = x((temp_all_x+1):(temp_all_x+temp_current));
                                        temp_handle(nn) = bar((temp_all_bars+1):(temp_all_bars+sum(log_nonzero)),x_temp(log_nonzero), 'FaceColor', colorpattern_LV(cc,:));
                                        temp_all_bars = temp_all_bars+sum(log_nonzero);
                                        temp_legend{nn} = strrep(selected_variables{1,ii}, '_', ' ');
                                        hold on
                                    else
                                        cc=cc+1;
                                    end
                                end
                                temp_all_x = temp_all_x+temp_current;
                        end
                    end
                    
                    if any(i==input.grid_dynamic.onset)
                        grid_x = input.grid_dynamic.(['LV_', num2str(i)]).x;
                        grid_y = input.grid_dynamic.(['LV_', num2str(i)]).y;
                    end
                    
                    first_line = strrep([input.name, ' grid_density_x=' num2str(grid_x.density), ', grid_density_y=' num2str(grid_y.density), ', LV ',num2str(i)], '_', ' ');
                    second_line = ['p-value (FDR-corrected) = ' num2str(output.final_parameters{i,opt_p}), ', Spearman''s RHO = ', num2str(output.final_parameters{i,opt_RHO})];
                    
                    if ~val_log
                        p_sig_threshold = output.pvalue_FDR(i);
                        fourth_line = 'main sample';
                    else
                        p_sig_threshold = 0.05;
                        fourth_line = 'validation sample';
                    end
                    
                    if output.final_parameters{i,opt_p}<=p_sig_threshold
                        third_line = 'significant';
                    else
                        third_line = 'not significant';
                    end
                    title({first_line; second_line; third_line; fourth_line}); % add third line
                    xlabel(LV_x, 'FontSize', font_size);
                    ylabel(LV_y, 'FontSize', font_size);
                    
                    legend(temp_handle, temp_legend, 'Location', 'bestoutside', 'FontSize', font_size);
                    hold all
                    
                    set(gcf,'Position', get(0,'Screensize'));
                    set(gcf,'PaperPositionMode','auto')
                    print(f, [b_folder, '/behavior_LV_' num2str(i)], '-dpng', '-r0');
                    saveas(f, [b_folder, '/behavior_LV' num2str(i), '.fig']);
                    %         saveas(f, [b_folder, '/behavior_LV' num2str(i), '.eps']);
                    saveas(f,[b_folder, '/behavior_LV' num2str(i)],'epsc');
                    
                    close(f);
                    output.questions_collection.(['LV_', num2str(i)]).final={};
                    for qq=1:size(output.questions_collection.(['LV_', num2str(i)]).items,1)
                        try output.questions_collection.(['LV_', num2str(i)]).final{qq} = [output.questions_collection.(['LV_', num2str(i)]).items{qq}, ': ', output.questions_collection.(['LV_', num2str(i)]).questions{1,qq}{1}];
                        catch
                            output.questions_collection.(['LV_', num2str(i)]).final{qq} = [output.questions_collection.(['LV_', num2str(i)]).items{qq}];
                        end
                    end
                    
                    output.questions_collection.(['LV_', num2str(i)]).combined = [output.questions_collection.(['LV_', num2str(i)]).final; output.questions_collection.(['LV_', num2str(i)]).subscales'];
                    
                    save(IN.results_path, 'input', 'output', 'setup');
                    
                end
                
                % print out questions
                fields = fieldnames(output.questions_collection);
                name_txt = [b_folder, '/collected_questions.txt'];
                writecell(cell(1,1), name_txt);
                
                for i=1:size(fields,1)
                    writecell(['LV_', num2str(i); output.questions_collection.(['LV_', num2str(i)]).combined(1,:)'; cell(2,1)],[b_folder, '/collected_questions.txt'],'WriteMode','append')
                end
                
                % print out post hoc correlations
                fields = fieldnames(output.post_hoc_correlations.all.correlations.test);
                fields = fields(contains(fields, 'table'));
                for i=1:size(fields,1)
                    writetable(output.post_hoc_correlations.all.correlations.test.(fields{i}), [b_folder, '/post_hoc_corr_test.xlsx'], 'WriteRowNames', true, 'Sheet', fields{i});
                end
                
                fields = fieldnames(output.post_hoc_correlations.all.correlations.validation);
                fields = fields(contains(fields, 'table'));
                for i=1:size(fields,1)
                    writetable(output.post_hoc_correlations.all.correlations.validation.(fields{i}), [b_folder, '/post_hoc_corr_val.xlsx'], 'WriteRowNames', true, 'Sheet', fields{i});
                end
                
                %% plot latent scores epsilon and omega, color code diagnostic groups (HC,
                % ROD, ROP, CHR)
                
                for i=1:size(input.selected_studygroups,2)
                    input.selected_studygroups{2,i} = strcmp(input.data_collection.Labels, input.selected_studygroups{1,i});
                end
                
                for i=1:size(output.CV.cv_outer_indices.TestInd,2)
                    %         temp_1 = zeros(size(input.selected_studygroups{2,1},1),size(input.selected_studygroups{2,1},2))>0;
                    %         temp_1(output.CV.cv_outer_indices.TestInd{1,i}) = true;
                    output.CV.cv_outer_indices.TestInd{2,i} = input.data_collection.Labels(output.CV.cv_outer_indices.TestInd{1,i});
                end
                
                % plot all latent scores according to LVs
                colorpattern_LS = hsv(size(input.selected_studygroups,2));
                for i=1:size(output.final_parameters,1)
                    f=figure();
                    for ii=1:size(input.selected_studygroups,2)
                        w = output.final_parameters{i,1};
                        x = rescale(output.final_parameters{i,index_epsilon}(strcmp(output.CV.cv_outer_indices.TestInd{2,w}, input.selected_studygroups{1,ii})), -1, 1);
                        y = rescale(output.final_parameters{i,index_omega}(strcmp(output.CV.cv_outer_indices.TestInd{2,w}, input.selected_studygroups{1,ii})), -1, 1);
                        plot(x,y,'.', 'MarkerSize', marker_size, 'color',colorpattern_LS(ii,:));
                        hold on
                    end
                    first_line = strrep([input.name, ' grid_density_x=' num2str(grid_x.density), ', grid_density_y=' num2str(grid_y.density), ', LV ',num2str(i)], '_', ' ');
                    second_line = ['p-value (FDR-corrected) = ' num2str(output.final_parameters{i,opt_p}), ', Spearman''s RHO = ', num2str(output.final_parameters{i,opt_RHO})];
                    
                    if ~val_log
                        p_sig_threshold = output.pvalue_FDR(i);
                        fourth_line = 'main sample';
                    else
                        p_sig_threshold = 0.05;
                        fourth_line = 'validation sample';
                    end
                    
                    if output.final_parameters{i,opt_p}<=p_sig_threshold
                        third_line = 'significant';
                    else
                        third_line = 'not significant';
                    end
                    
                    title({first_line; second_line; third_line; fourth_line}); % add third line
                    xlabel(LS_epsilon);
                    ylabel(LS_omega);
                    lsline
                    [~, lgd_data] = legend({input.selected_studygroups{1,:}}, 'FontSize', font_size);
                    lgd_line = findobj(lgd_data, 'type', 'line'); %// objects of legend of type line
                    set(lgd_line, 'Markersize', marker_size); %// set marker size as desired
                    
                    set(gcf, 'Position', get(0, 'Screensize'));
                    set(gcf,'PaperPositionMode','auto')
                    %     print(f, [results_folder, '/latent_scores_LV' num2str(i)], '-dpng', '-r0');
                    print(f, [b_folder, '/latent_scores_LV' num2str(i)], '-dpng', '-r0');
                    %     saveas(f, [results_folder, '/latent_scores_LV' num2str(i), '.fig']);
                    saveas(f, [b_folder, '/latent_scores_LV' num2str(i), '.fig']);
                    %         saveas(f, [b_folder, '/latent_scores_LV' num2str(i), '.eps']);
                    saveas(f,[b_folder, '/latent_scores_LV' num2str(i)],'epsc');
                    %     saveas(f, [results_folder, '/latent_scores_LV' num2str(i), '.fig']);
                    %     saveas(f, [results_folder, '/latent_scores_LV' num2str(i), '.png']);
                    close(f)
                end
                
                % plot latent scores all combined across significant LVs, colorcoded by LVs
                colorpattern_LS = hsv(size(output.final_parameters,1));
                f=figure(); temp_legend = [];
                for i=1:size(output.final_parameters,1)
                    x = rescale(output.final_parameters{i,index_epsilon}, -1, 1);
                    y = rescale(output.final_parameters{i,index_omega}, -1, 1);
                    plot(x,y,'.', 'MarkerSize', marker_size, 'color',colorpattern_LS(i,:));
                    
                    first_line = strrep([input.name, ', grid_density_x=' num2str(grid_x.density), ' grid_density_y=' num2str(grid_y.density)], '_', ' ');
                    second_line = 'latent scores from all LVs combined';
                    if output.final_parameters{i,opt_p}<=output.pvalue_FDR(i)
                        sig_line = 'significant';
                    else
                        sig_line = 'not significant';
                    end
                    title({first_line;second_line}); % add third line
                    hold on
                    temp_legend{i} = ['LV ', num2str(i), ': ', sig_line];
                end
                %     temp_legend{nn} = selected_variables{1,ii};
                xlabel(LS_epsilon);
                ylabel(LS_omega);
                lsline
                [~, lgd_data] = legend(temp_legend, 'Location', 'bestoutside', 'FontSize', font_size);
                lgd_line = findobj(lgd_data, 'type', 'line'); %// objects of legend of type line
                set(lgd_line, 'Markersize', marker_size); %// set marker size as desired
                
                set(gcf, 'Position', get(0, 'Screensize'));
                set(gcf,'PaperPositionMode','auto')
                print(f, [b_folder, '/latent_scores_combined_LV_color'], '-dpng', '-r0');
                saveas(f, [b_folder, '/latent_scores_combined_LV_color.fig']);
                %     saveas(f, [b_folder, '/latent_scores_combined_LV_color.eps']);
                saveas(f,[b_folder, '/latent_scores_combined_LV_color'],'epsc');
                close(f)
                
                %     % standardize latent scores and plot all of them in one graph, colorcoded
                %     % by diagnoses
                %     % first transpose latent scores so that they fit function, then standardize
                %     % feature-wise (per LV)
                %     output.epsilon_stand = (dp_standardize(output.epsilon'))';
                %     output.omega_stand = (dp_standardize(output.omega'))';
                
                % plot all latent scores according to LVs, colorcoded by diagnoses
                colorpattern_LS = hsv(size(input.selected_studygroups,2));
                f=figure();
                for i=1:size(output.final_parameters,1)
                    if output.final_parameters{i,opt_p}<=output.pvalue_FDR(i)
                        for ii=1:size(input.selected_studygroups,2)
                            w = output.final_parameters{i,1};
                            x = rescale(output.final_parameters{i,index_epsilon}(strcmp(output.CV.cv_outer_indices.TestInd{2,w}, input.selected_studygroups{1,ii})), -1, 1);
                            y = rescale(output.final_parameters{i,index_omega}(strcmp(output.CV.cv_outer_indices.TestInd{2,w}, input.selected_studygroups{1,ii})), -1 , 1);
                            plot(x,y,'.', 'MarkerSize', marker_size, 'color',colorpattern_LS(ii,:));
                            hold on
                        end
                    end
                end
                first_line = strrep([input.name, ', grid_density_x=' num2str(grid_x.density), ' grid_density_y=' num2str(grid_y.density)], '_', ' ');
                second_line = 'latent scores from all LVs combined';
                %             if output.final_parameters{i,opt_p}<=output.pvalue_FDR(i)
                %                 third_line = ['significant'];
                %             else
                %                 third_line = ['not significant'];
                %             end
                title({first_line; second_line}); % add third line
                xlabel(LS_epsilon);
                ylabel(LS_omega);
                [~, lgd_data] = legend({input.selected_studygroups{1,:}}, 'FontSize', font_size);
                lgd_line = findobj(lgd_data, 'type', 'line'); %// objects of legend of type line
                set(lgd_line, 'Markersize', marker_size); %// set marker size as desired
                
                set(gcf, 'Position', get(0, 'Screensize'));
                set(gcf,'PaperPositionMode','auto')
                %     print(f, [results_folder, '/latent_scores_LV' num2str(i)], '-dpng', '-r0');
                print(f, [b_folder, '/latent_scores_combined_diagnosis_color'], '-dpng', '-r0');
                %     saveas(f, [results_folder, '/latent_scores_LV' num2str(i), '.fig']);
                saveas(f, [b_folder, '/latent_scores_combined_diagnosis_color', '.fig']);
                %     saveas(f, [b_folder, '/latent_scores_combined_diagnosis_color', '.eps']);
                saveas(f,[b_folder, '/latent_scores_combined_diagnosis_color'],'epsc');
                %     saveas(f, [results_folder, '/latent_scores_LV' num2str(i), '.fig']);
                %     saveas(f, [results_folder, '/latent_scores_LV' num2str(i), '.png']);
                close(f);
            else
                dp_vis_general(b_folder, IN.results_path, {'v'}, IN.plot_all, false, {'abstinence_aktuell', 'abstinence_always'}, input.validation_set)
            end
        case 2
            dp_vis_general(b_folder, IN.results_path, {'v'}, IN.plot_all, false, {'abstinence_aktuell', 'abstinence_always'}, input.validation_set)
    end
end

if any(strcmp(IN.specific, 'detailed') | all_jobs)
    d_folder = [collection_folder, '/detailed'];
    mkdir(d_folder);
    
    %     detailed_results_folder = [IN.results_path(1:(strfind(IN.results_path, 'final')-1)), 'detailed_results'];
    cd(d_folder);
    
    for i=1:size(output.final_parameters,1)
        load([detailed_results_folder, '/opt_parameters_' num2str(i), '.mat']);
        if exist('opt_parameters')
            if size(opt_parameters,2) == 8
                opt_param_p = strcmp(opt_parameters_names, 'p');
                opt_param_v = strcmp(opt_parameters_names, 'v');
                opt_param_u = strcmp(opt_parameters_names, 'u');
                opt_param_RHO = strcmp(opt_parameters_names, 'RHO');
            elseif size(opt_parameters,2) == 11
                opt_param_p = strcmp(opt_parameters_names_long, 'p');
                opt_param_v = strcmp(opt_parameters_names_long, 'v');
                opt_param_u = strcmp(opt_parameters_names_long, 'u');
                opt_param_RHO = strcmp(opt_parameters_names_long, 'RHO');
            else
                disp('Something''s wrong!');
            end
            temp_opt_param = opt_parameters;
        elseif exist('opt_parameters_temp')
            if size(opt_parameters_temp,2) == 8 || size(opt_parameters_temp,2) == 10
                opt_param_p = strcmp(opt_parameters_names, 'p');
                opt_param_v = strcmp(opt_parameters_names, 'v');
                opt_param_u = strcmp(opt_parameters_names, 'u');
                opt_param_RHO = strcmp(opt_parameters_names, 'RHO');
            elseif size(opt_parameters_temp,2) == 11
                opt_param_p = strcmp(opt_parameters_names_long, 'p');
                opt_param_v = strcmp(opt_parameters_names_long, 'v');
                opt_param_u = strcmp(opt_parameters_names_long, 'u');
                opt_param_RHO = strcmp(opt_parameters_names_long, 'RHO');
            else
                disp('Something''s wrong!');
            end
            temp_opt_param = opt_parameters_temp;
        else
            disp('Something''s wrong!');
        end
        %     opt_dir = [detailed_results_folder, '/opt_parameters_' num2str(i)];
        %     mkdir(opt_dir);
        %     cd(opt_dir);
        
        d_folder_opt = [d_folder, '/opt_parameters_' num2str(i)];
        mkdir(d_folder_opt);
        cd(d_folder_opt);
        
        %     for ii=1:size(temp_opt_param,1)
        %         %% visualize results
        %         % write brain vector to nifti file
        %         nk_WriteVol(temp_opt_param{ii,opt_param_u}, ['brain_opt_parameters_' num2str(i), '_' num2str(temp_opt_param{ii,1})], 2, NM.brainmask{1,1}, cell2mat(NM.badcoords), 0, 'gt');
        % %         dp_resample_image([d_folder_opt, '/brain_opt_parameters_' num2str(i), '_' num2str(temp_opt_param{ii,1})], [1 1 1]);
        %     end
        
        %% visualize behavior vector in barplot
        
        for ii=1:(size(temp_opt_param,1))
            x = temp_opt_param{ii,opt_param_v};
            f=subplot(round(size(temp_opt_param,1)/2),2,ii);
            nn=0;
            hold on
            temp_legend=[]; temp_all = 0;
            for iii=1:size(selected_variables,2)
                switch class(selected_variables{2,iii})
                    case 'struct'
                        fields = fieldnames(input.selected_features{2,(strcmp(input.selected_features(1,:),input.selected_features{1,iii}))});
                        for iiii=1:size(fields,1)
                            nn=nn+1;
                            temp_current=size(selected_variables{2,iii}.(fields{iiii}),2);
                            bar((temp_all+1):(temp_all+temp_current),x((temp_all+1):(temp_all+temp_current)), 'FaceColor', colorpattern_LV(nn,:));
                            temp_all = temp_all+temp_current;
                            temp_legend{nn} = [selected_variables{1,iii}, ' ', strrep(fields{iiii}, '_', ' ')];
                            hold on
                        end
                    case 'double'
                        nn=nn+1;
                        temp_current=size(selected_variables{2,iii},2);
                        bar((temp_all+1):(temp_all+size(temp_current,2)),x((temp_all+1):(temp_all+size(temp_current,2))), 'FaceColor', colorpattern_LV(nn,:));
                        temp_all = temp_all+size(temp_current,2);
                        temp_legend{nn} = selected_variables{1,iii};
                        hold on
                end
            end
            axis([0 (size(temp_opt_param{ii,opt_param_v},1)+1) -1 1]);
            xlabel({'\color{black}clinical features'}, 'FontSize', font_size);
            ylabel('weights', 'FontSize', font_size);
            if temp_opt_param{ii,opt_param_p}<=output.pvalue_FDR(i)
                significance_opt = 'significant';
            else
                significance_opt = 'not significant';
            end
            subplot_title = {['Iteration ', num2str(ii), ', p-value (FDR-corrected) = ' num2str(temp_opt_param{ii,opt_param_p}), ', Spearman''s RHO = ', num2str(temp_opt_param{ii,opt_param_RHO}), significance_opt]};
            title(subplot_title, 'FontSize', font_size, 'FontWeight', 'normal');
            
        end
        %     set(gcf, 'Position', get(0, 'Screensize'));
        set(gcf, 'Position', get(0, 'Screensize'));
        set(gcf,'PaperPositionMode','auto')
        first_line = strrep([input.name, ' grid_density_x=' num2str(grid_x.density), ' grid_density_y=' num2str(grid_y.density), ', LV ',num2str(i)], '_', ' ');
        second_line = ['p-value (FDR-corrected) = ' num2str(output.final_parameters{i,opt_p}), ', Spearman''s RHO = ', num2str(output.final_parameters{i,opt_RHO})];
        if output.final_parameters{i,opt_p}<=output.pvalue_FDR(i)
            third_line = ['significant'];
        else
            third_line = ['not significant'];
        end
        suptitle({first_line; [second_line, ', ' third_line]}); % add third line
        
        %     print([detailed_results_folder, '/behavior_opt_parameters_',num2str(i)], '-dpng', '-r0');
        %     saveas(f, [detailed_results_folder, '/behavior_opt_parameters_',num2str(i)], 'png');
        
        print([d_folder_opt, '/behavior_opt_parameters_',num2str(i)], '-dpng', '-r0');
        saveas(f, [d_folder_opt, '/behavior_opt_parameters_',num2str(i)], 'png');
        %         saveas(f, [d_folder_opt, '/behavior_opt_parameters_',num2str(i)], 'eps');
        saveas(f,[d_folder_opt, '/behavior_opt_parameters_',num2str(i)],'epsc');
        %     saveas(f, [detailed_results_folder, '/behavior_opt_parameters_',num2str(i)], 'fig');
        %     saveas(f, [detailed_results_folder, '/behavior_opt_parameters_',num2str(i)], 'png');
        close all;
    end
    
end

if any(strcmp(IN.specific, 'correlation') | all_jobs)
    
    for i=1:size(output.CV.cv_outer_indices.TestInd,2)
        output.CV.cv_outer_indices.TestInd{2,i} = input.data_collection.Labels(output.CV.cv_outer_indices.TestInd{1,i});
    end
    
    % plot all latent scores according to LVs
    s_folder = [collection_folder, '/correlation'];
    mkdir(s_folder);
    latent_selection = {index_epsilon, index_omega; LS_epsilon, LS_omega; 'epsilon', 'omega'};
    colorpattern_LS = hsv(size(input.selected_studygroups,2));
    fields = fieldnames(output.post_hoc_correlations.correlations);
    fields(contains(fields, 'RHO')|contains(fields, 'p'))=[];
    for ff=1:size(fields,1)
        p_sig_collection = output.post_hoc_correlations.correlations.(fields{ff}).table_p.Variables < 0.05;
        for l=1:size(latent_selection,2)
            for i=1:size(output.post_hoc_correlations.correlations.(fields{ff}).table_p.Variables,1)
                for ii=1:size(output.post_hoc_correlations.correlations.(fields{ff}).table_p.Variables,2)
                    if p_sig_collection(i,ii)
                        f=figure();
                        for iii=1:size(input.selected_studygroups,2)
                            w = output.final_parameters{round(ii/2),1};
                            x = output.final_parameters{round(ii/2),latent_selection{1,l}}(strcmp(output.CV.cv_outer_indices.TestInd{2,w}, input.selected_studygroups{1,iii}));
                            temp=output.post_hoc_correlations.data_collection(:,i);
                            y = temp(strcmp(output.CV.cv_outer_indices.TestInd{2,w}, input.selected_studygroups{1,iii}));
                            plot(x,y,'.', 'MarkerSize', marker_size, 'color',colorpattern_LS(iii,:));
                            hold on
                        end
                        
                        xlabel(latent_selection{2,l}, 'FontSize', font_size);
                        y_temp = output.post_hoc_correlations.data_table.Properties.VariableNames(ii);
                        y_temp = strrep(y_temp, '_', ' ');
                        y_temp = strrep(y_temp, 'T0','');
                        y_temp = strrep(y_temp, 'Screening','');
                        ylabel(y_temp, 'FontSize', font_size);
                        title(['p-value = ' num2str(output.post_hoc_correlations.correlations.(fields{ff}).table_p{i,ii}), ', Spearman''s RHO = ', num2str(output.post_hoc_correlations.correlations.(fields{ff}).table_RHO{i,ii})], 'FontSize', font_size); % add third line
                        [~, lgd_data] = legend({input.selected_studygroups{1,:}}, 'FontSize', font_size);
                        lgd_line = findobj(lgd_data, 'type', 'line'); %// objects of legend of type line
                        set(lgd_line, 'Markersize', marker_size); %// set marker size as desired
                        hold on
                        
                        set(gcf, 'Position', get(0, 'Screensize'));
                        set(gcf,'PaperPositionMode','auto')
                        print(f, [s_folder, '/latent_scores_LV', num2str(round(ii/2)), '_', latent_selection{3,l}, '_', output.post_hoc_correlations.correlations.(fields{ff}).table_RHO.Properties.VariableNames{ii}], '-dpng', '-r0');
                        saveas(f, [s_folder, '/latent_scores_LV', num2str(round(ii/2)), '_', latent_selection{3,l}, '_', output.post_hoc_correlations.correlations.(fields{ff}).table_RHO.Properties.VariableNames{ii}, '.fig']);
                        %                         saveas(f, [s_folder, '/latent_scores_LV', num2str(round(ii/2)), '_', latent_selection{3,l}, '_', output.post_hoc_correlations.correlations.(fields{ff}).table_RHO.Properties.VariableNames{ii}, '.eps']);
                        saveas(f,[s_folder, '/latent_scores_LV', num2str(round(ii/2)), '_', latent_selection{3,l}, '_', output.post_hoc_correlations.correlations.(fields{ff}).table_RHO.Properties.VariableNames{ii}],'epsc');
                        close all;
                        
                    end
                end
            end
        end
    end
    save(IN.results_path, 'input', 'output', 'setup');
    
end


end
