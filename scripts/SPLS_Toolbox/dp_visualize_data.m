%% function to visualize SPLS output

function dp_visualize_data(IN)

switch IN.overall_analysis
    case 'Stress'
        overall_folder = '/volume/HCStress/Analysis/Stress';
    case 'Resilience'
        overall_folder = '/volume/HCStress/Analysis/Resilience';
end

if ~isfield(IN, 'specific')
    all_jobs = true;
    IN.specific = 'empty';
else
    all_jobs = false;
end

%% visualize results
load(IN.results_path);
collection_folder = [overall_folder, '/', input.name];
mkdir(collection_folder);

input.behavior_names = strrep(input.behavior_names, '_T0', '');
load(input.NM_structure);
load('/volume/HCStress/Doc/Stress_Resilience_questionnaires.mat');

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

% results_folder = IN.results_path(1:(strfind(IN.results_path, '/result.mat')-1));

% cd(results_folder);
%
% % write brain vector to nifti file
% for i=1:size(output.final_parameters,1)
%     nk_WriteVol(output.final_parameters{i,4}, ['brain_LV' num2str(i)], 2, NM.brainmask{1,1}, cell2mat(NM.badcoords), 0, 'gt');
% end

if any(strcmp(IN.specific, 'images') | all_jobs)
    cd(collection_folder);
    % write brain vector to nifti file
    for i=1:size(output.final_parameters,1)
        nk_WriteVol(output.final_parameters{i,4}, ['brain_LV' num2str(i)], 2, NM.brainmask{1,1}, cell2mat(NM.badcoords), 0, 'gt');
        %     dp_resample_image([collection_folder, '/brain_LV' num2str(i)], [1 1 1]);
    end
end


if any(strcmp(IN.specific, 'atlas') | all_jobs)
    
    % get clusters for brain regions using hammers and aal atlas
    % filepath hammers nifti: /opt/SPM/spm12_v6685_cat12_r1207/atlas/hammers.nii
    % filepath hammers description: /opt/SPM/spm12_v6685_cat12_r1207/atlas/labels_dartel_hammers.xml
    
    switch size(input.behavior,1)
        case 626
            hammers_path_full = '/volume/HCStress/Data/MRI/Hammers_Atlas/Hammers_mith-n30r95_spm12_full__626_NM_X.mat';
        case 627
            hammers_path_full = '/volume/HCStress/Data/MRI/Hammers_Atlas/Hammers_mith-n30r95_spm12_full__627_NM_X.mat';
        case 631
            hammers_path_full = '/volume/HCStress/Data/MRI/Hammers_Atlas/Hammers_mith-n30r95_spm12_full__631_NM_X.mat';
        case 634
            hammers_path_full = '/volume/HCStress/Data/MRI/Hammers_Atlas/Hammers_mith-n30r95_spm12_full__634_NM_X.mat';
    end
    
    load(hammers_path_full);
    [C_hammers, ~, ic_hammers] = unique(atlas_hammers_full_for_analysis);
    counts_hammers = accumarray(ic_hammers, 1);
    output.regions_hammers_count = [C_hammers', counts_hammers];
    load('/volume/HCStress/Data/MRI/Hammers_Atlas/Hammers_mith-n30-ancillary-data.tar/Hammers_mith-n30-ancillary-data/Hammers_mith-n30-ancillary-data/labels_hammerssmith_n30_ancillary_data.mat');
    output.regions.count_names = {'region_number', 'n_voxels', 'voxel_percentage', 'region_name', 'median_weights', 'mean_weights'};
    output.regions.count = struct('raw', [], 'voxels', [], 'weights', []);
    output.regions.count_names_sorting = {'positive', 'negative', 'overall'};
    voxel_cutoff = 10;
    fields = fieldnames(output.regions.count);
    for i=1:size(fields,1)
        for fff = 1:size(output.regions.count_names_sorting,2)
            FID = fopen([collection_folder, '/brain_regions_hammers_', output.regions.count_names_sorting{fff}, '_', fields{i}, '.txt'], 'w');
            fprintf(FID, [strrep(input.name,'_',' '), '\n', fields{i} ' sorted']);
            fclose(FID);
        end
    end
    
    for i=1:size(output.final_parameters,1)
        output.regions.log{i,1} = output.final_parameters{i,4}>0;
        output.regions.log{i,2} = output.final_parameters{i,4}<0;
        output.regions.log{i,3} = output.final_parameters{i,4}~=0;
        for ii=1:size(output.regions.log,2)
            output.regions.sum{i,ii} = atlas_hammers_full_for_analysis((atlas_hammers_full_for_analysis~=0)' & output.regions.log{i,ii});
            [C, ~, ic] = unique(output.regions.sum{i,ii});
            a_counts = accumarray(ic, 1);
            voxel_percentage = a_counts./(counts_hammers(ismember(C_hammers, C)));
            output.regions.count.raw{i,ii} = [num2cell(C'), num2cell(a_counts), num2cell(voxel_percentage), labels_regions_hammers((C'),2)];
            fix_end = size(output.regions.count.raw{i,ii},2);
            log_cutoff = cell2mat(output.regions.count.raw{i,ii}(:, find(strcmp(output.regions.count_names, 'n_voxels')))) <voxel_cutoff;
            output.regions.count.raw{i,ii}(log_cutoff,:) = [];
            if isempty(output.regions.count.raw{i,ii})
                output.regions.count.raw{i,ii} = [];
                output.regions.count.voxels{i,ii} = [];
                output.regions.count.weights{i,ii} = [];
            else
                for iii=1:size(output.regions.count.raw{i,ii},1)
                    if strcmp(output.regions.count_names_sorting{ii}, 'overall')
                        output.regions.count.raw{i,ii}{iii,fix_end+1} = median(abs(output.final_parameters{i,4}(output.regions.log{i,ii} & (atlas_hammers_full_for_analysis==output.regions.count.raw{i,ii}{iii,1})')));
                        output.regions.count.raw{i,ii}{iii,fix_end+2} = mean(abs(output.final_parameters{i,4}(output.regions.log{i,ii} & (atlas_hammers_full_for_analysis==output.regions.count.raw{i,ii}{iii,1})')));
%                         output.regions.count.raw{i,ii}{iii,6} = 
                    else
                        output.regions.count.raw{i,ii}{iii,fix_end+1} = median(output.final_parameters{i,4}(output.regions.log{i,ii} & (atlas_hammers_full_for_analysis==output.regions.count.raw{i,ii}{iii,1})'));
                        output.regions.count.raw{i,ii}{iii,fix_end+2} = mean(output.final_parameters{i,4}(output.regions.log{i,ii} & (atlas_hammers_full_for_analysis==output.regions.count.raw{i,ii}{iii,1})'));
                    end
                end
                
                if ~strcmp(output.regions.count_names_sorting{ii}, 'negative')
                    output.regions.count.voxels{i,ii} = flipud(sortrows((output.regions.count.raw{i,ii}), find(strcmp(output.regions.count_names, 'n_voxels'))));
                    output.regions.count.weights{i,ii} = flipud(sortrows((output.regions.count.raw{i,ii}), find(strcmp(output.regions.count_names, 'median_weights'))));
                else
                    output.regions.count.voxels{i,ii} = flipud(sortrows((output.regions.count.raw{i,ii}), find(strcmp(output.regions.count_names, 'n_voxels'))));
                    output.regions.count.weights{i,ii} = sortrows((output.regions.count.raw{i,ii}), find(strcmp(output.regions.count_names, 'median_weights')));
                end
            
            end
        end
        
        temp_add = [];
        temp_print = output.regions.count_names([1,2,4,5,3]);
        for f=1:size(temp_print,2)
            %         try output.regions.count_names{f} = num2str(temp_print{f});
            %         catch ME
            %         end
            temp_add = [temp_add, sprintf('\t'), temp_print{f}];
        end
        
        
        fields = fieldnames(output.regions.count);
        for ff=1:size(fields,1)
            for fff = 1:size(output.regions.count.(fields{ff}),2)
                FID = fopen([collection_folder, '/brain_regions_hammers_', output.regions.count_names_sorting{fff}, '_', fields{ff}, '.txt'], 'a');
                fprintf(FID, ['\n \n \n', sprintf('\t'), 'Latent Variable ' num2str(i), '\n \n', temp_add, '\n']);
                fclose(FID);
            end
        end
        
        fields = fieldnames(output.regions.count);
        for ff=1:size(fields,1)
            for fff = 1:size(output.regions.count.(fields{ff}),2)
                if ~isempty(output.regions.count.(fields{ff}){i,fff})
                    for r=1:size(output.regions.count.(fields{ff}){i,fff},1)
                        temp_print = output.regions.count.(fields{ff}){i,fff}(r,[1,2,4,5,3]);
                        temp_add = [];
                        for ii=1:size(temp_print,2)
                            try temp_print{ii} = num2str(temp_print{ii});
                            catch ME
                            end
                            temp_add = [temp_add, sprintf('\t'), temp_print{ii}];
                        end
                        
                        FID = fopen([collection_folder, '/brain_regions_hammers_', output.regions.count_names_sorting{fff}, '_', fields{ff}, '.txt'], 'a');
                        fprintf(FID, ['\n' temp_add]);
                        fclose(FID);
                    end
                end
            end
        end
        
        % %     dp_resample_image([collection_folder, '/brain_LV' num2str(i)], [1 1 1]);
    end
    
    save(IN.results_path, 'input', 'output', 'setup');
        
end

if any(strcmp(IN.specific, 'behavior') | all_jobs)
    
    %% visualize behavior vector in barplot
    % CTQ.emotional_abuse)), (CTQ.physical_abuse)), (CTQ.sexual_abuse)), (CTQ.emotional_neglect)), (CTQ.physical_neglect)), (CTQ.denial))];
    % load([setup.analysis_folder, '/' setup.date, '_', name, '_data_collection.mat']);
    if sum(ismember(input.behavior_names, 'Age')) + sum(ismember(input.behavior_names, 'Sex'))==2
        selected_variables = [input.selected_features(1,:), 'Age', 'Sex'; input.selected_features(2,:), 1, 1];
    elseif sum(ismember(input.behavior_names, 'Age'))
        selected_variables = [input.selected_features(1,:), 'Age'; input.selected_features(2,:), 1];
    elseif sum(ismember(input.behavior_names, 'Sex'))
        selected_variables = [input.selected_features(1,:), 'Sex'; input.selected_features(2,:), 1];
    else
        selected_variables = input.selected_features;
    end
    
    % compute measures for effects
    % Cohen's d
%     for i=1:size(output.final_parameters,1)
%         x1 = output.final_parameters{i, index_epsilon};
%         x2 = output.final_parameters{i, index_omega};
%         output.Cohen(i) = computeCohen_d(x1, x2, 'independent');
%         output.Spearman(i) = corr(x1, x2, 'Type', 'Spearman');
%         output.Kendall(i) = corr(x1, x2, 'Type', 'Kendall');
%         output.MI_kernel(i) = kernelmi(x1',x2');
%         output.MI_peng(i) = mutualinfo(x1,x2);
%     end
    
    count=0;
    for i=1:size(input.selected_features(2,:),2)
        temp=size(fieldnames(input.selected_features{2,i}),1);
        count=count+temp;
    end
    colorpattern_LV = hsv((count + 2));
    
    for i=1:(size(output.final_parameters,1))
        questions_collection=[];
        x = output.final_parameters{i,opt_v};
        questions_collection.items = input.behavior_names(x~=0)';
        Resilience_Stress_questionnaires = [CISS_questionnaire', RSA_questionnaire', CTQ_questionnaire', BS_questionnaire'];
        for q=1:size(questions_collection.items,1)
            questions_collection.questions{q} = Resilience_Stress_questionnaires(1,strcmp(Resilience_Stress_questionnaires(2,:), questions_collection.items{q}));
        end
        %     errors = output.CI_v{i};
        f=figure();
        nn=0;
        hold on
        temp_legend=[]; temp_all = 0; %sel_temp_names=[];
        for ii=1:(size(selected_variables,2))
            switch class(selected_variables{2,ii})
                case 'struct'
                    fields = fieldnames(input.selected_features{2,(strcmp(input.selected_features(1,:),input.selected_features{1,ii}))});
                    for iii=1:size(fields,1)
                        nn=nn+1;
                        temp_current=size(selected_variables{2,ii}.(fields{iii}),2);
                        try temp_handle(nn) = dp_barwitherr(errors([(temp_all+1):(temp_all+temp_current)],:),(temp_all+1):(temp_all+temp_current),x((temp_all+1):(temp_all+temp_current)), 'FaceColor', colorpattern_LV(nn,:));
                        catch 
                            temp_handle(nn) = bar((temp_all+1):(temp_all+temp_current),x((temp_all+1):(temp_all+temp_current)), 'FaceColor', colorpattern_LV(nn,:));
                        end
                        %                     temp_names = input.behavior_names((temp_all+1):(temp_all+temp_current));
                        %                     sel_temp_names = [sel_temp_names; temp_names(x((temp_all+1):(temp_all+temp_current))~=0)];
                        temp_all = temp_all+temp_current;
                        temp_legend{nn} = [selected_variables{1,ii}, ' ', strrep(fields{iii}, '_', ' ')];
                        hold on
                    end
                case 'double'
                    nn=nn+1;
                    temp_current=size(selected_variables{2,ii},2);
                    try temp_handle(nn) = dp_barwitherr(errors([(temp_all+1):(temp_all+temp_current)],:),(temp_all+1):(temp_all+temp_current),x((temp_all+1):(temp_all+temp_current)), 'FaceColor', colorpattern_LV(nn,:));
                    catch 
                        temp_handle(nn) = bar((temp_all+1):(temp_all+temp_current),x((temp_all+1):(temp_all+temp_current)), 'FaceColor', colorpattern_LV(nn,:));
                    end
                    temp_all = temp_all+temp_current;
                    temp_legend{nn} = selected_variables{1,ii};
                    hold on
            end
        end
        first_line = strrep([input.name, ' grid_density=' num2str(input.grid_density), ', LV ',num2str(i)], '_', ' ');
        second_line = ['p-value (FDR-corrected) = ' num2str(output.final_parameters{i,opt_p}), ', Spearman''s RHO = ', num2str(output.final_parameters{i,opt_RHO})];
        if output.final_parameters{i,opt_p}<=output.pvalue_FDR(i)
            third_line = ['significant'];
        else
            third_line = ['not significant'];
        end
        title({first_line; second_line; third_line}); % add third line
        xlabel('weight vector v', 'FontSize', 12);
        ylabel('score', 'FontSize', 12);
        axis([0 (size(output.final_parameters{i,opt_v},1)+1) -1 1]);
        legend(temp_handle, temp_legend, 'Location', 'bestoutside', 'FontSize', 12);
        hold all
        %     legend('test', 'Location', 'bestoutside');
        
        x_pos = 1:size(x,1);
        y_value = x;
        ygap = 0.05;  % Specify vertical gap between the bar and label
        ylimits = get(gca,'YLim');
        %     set(gca,'YLim',[ylim(1),ylim(2)+0.2*max(y)]); % Increase y limit for labels
        
        
        for xi=1:size(x_pos,2) % Loop over each bar
            %         xpos = x_pos(i);        % Set x position for the text label
            if xi==1 %|| i==size(x_pos,2)
                if y_value(xi)<0
                    ypos(xi) = ygap;
                elseif y_value(xi)>0
                    ypos(xi) = -ygap;% Set y position, including gap
                else
                    ypos(xi)=0;
                end
            else
                if y_value(xi)~=0
                    if y_value(xi)>0 && y_value(xi-1)>0 %abs(y_value(i)-y_value(i-1))<=ygap
                        ypos(xi) = ypos(xi-1) - ygap;
                    elseif y_value(xi)<0 && y_value(xi-1)<0
                        ypos(xi) = ypos(xi-1) + ygap;
                    elseif y_value(xi) > 0
                        ypos(xi) = -ygap;
                    elseif y_value(xi) < 0
                        ypos(xi) = ygap;
                    end
                else
                    ypos(xi)=0;
                end
            end
            if y_value(xi)~=0
                htext = text(x_pos(xi),ypos(xi),strrep(input.behavior_names{xi},'_',' '));  % Add text label
                set(htext,'VerticalAlignment','bottom','HorizontalAlignment','center', 'FontSize', 8); % Adjust properties
            end
        end
        
        for qq=1:size(questions_collection.items,1)
            try questions_collection.final{qq} = [questions_collection.items{qq}, ': ', questions_collection.questions{1,qq}{1}];
            catch
                questions_collection.final{qq} = [questions_collection.items{qq}];
            end
        end
        
        try annotation(f, 'textbox', [0.79, 0.2, 0.16, 0.4], 'string', strrep(questions_collection.final, '_', ' '), 'FontSize', 8, 'FitHeightToText', 'on');
        catch
        end
        
        set(gcf,'Position', get(0,'Screensize'));
        set(gcf,'PaperPositionMode','auto')
        print(f, [collection_folder, '/behavior_LV_' num2str(i)], '-dpng', '-r0');
        saveas(f, [collection_folder, '/behavior_LV' num2str(i), '.fig']);
        
        close(f);
    end
    
    %% plot latent scores epsilon and omega, color code diagnostic groups (HC,
    % ROD, ROP, CHR)
    
    for i=1:size(input.selected_studygroups,2)
        input.selected_studygroups{2,i} = strcmp(input.data_collection.Labels, input.selected_studygroups{1,i});
    end
    
    for i=1:size(output.cv_outer.TestInd,2)
%         temp_1 = zeros(size(input.selected_studygroups{2,1},1),size(input.selected_studygroups{2,1},2))>0; 
%         temp_1(output.cv_outer.TestInd{1,i}) = true;
        output.cv_outer.TestInd{2,i} = input.data_collection.Labels(output.cv_outer.TestInd{1,i});
    end
    
    % plot all latent scores according to LVs
    colorpattern_LS = hsv(size(input.selected_studygroups,2));
    for i=1:size(output.final_parameters,1)
        f=figure();
        for ii=1:size(input.selected_studygroups,2)
            w = output.final_parameters{i,1};
            x = output.final_parameters{i,index_epsilon}(strcmp(output.cv_outer.TestInd{2,w}, input.selected_studygroups{1,ii}));
            y = output.final_parameters{i,index_omega}(strcmp(output.cv_outer.TestInd{2,w}, input.selected_studygroups{1,ii}));
            plot(x,y,'.', 'MarkerSize', 20, 'color',colorpattern_LS(ii,:));
            hold on
        end
        first_line = strrep([input.name, ' grid_density=' num2str(input.grid_density), ', LV ',num2str(i)], '_', ' ');
        second_line = ['p-value (FDR-corrected) = ' num2str(output.final_parameters{i,opt_p}), ', Spearman''s RHO = ', num2str(output.final_parameters{i,opt_RHO})];
        if output.final_parameters{i,opt_p}<=output.pvalue_FDR(i)
            third_line = 'significant';
        else
            third_line = 'not significant';
        end
        title({first_line; second_line; third_line}); % add third line
        xlabel('epsilon(brain score)');
        ylabel('omega(behavior score)');
        [~, lgd_data] = legend({input.selected_studygroups{1,:}}, 'FontSize', 18);
        lgd_line = findobj(lgd_data, 'type', 'line'); %// objects of legend of type line
        set(lgd_line, 'Markersize', 20); %// set marker size as desired
        
        set(gcf, 'Position', get(0, 'Screensize'));
        set(gcf,'PaperPositionMode','auto')
        %     print(f, [results_folder, '/latent_scores_LV' num2str(i)], '-dpng', '-r0');
        print(f, [collection_folder, '/latent_scores_LV' num2str(i)], '-dpng', '-r0');
        %     saveas(f, [results_folder, '/latent_scores_LV' num2str(i), '.fig']);
        saveas(f, [collection_folder, '/latent_scores_LV' num2str(i), '.fig']);
        
        %     saveas(f, [results_folder, '/latent_scores_LV' num2str(i), '.fig']);
        %     saveas(f, [results_folder, '/latent_scores_LV' num2str(i), '.png']);
        close(f)
    end
    
    % plot latent scores all combined across significant LVs, colorcoded by LVs
    colorpattern_LS = hsv(size(output.final_parameters,1));
    f=figure(); temp_legend = [];
    for i=1:size(output.final_parameters,1)
        x = output.final_parameters{i,index_epsilon};
        y = output.final_parameters{i,index_omega};
        plot(x,y,'.', 'MarkerSize', 20, 'color',colorpattern_LS(i,:));
        first_line = strrep([input.name, ' grid_density=' num2str(input.grid_density)], '_', ' ');
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
    xlabel('epsilon(brain score)');
    ylabel('omega(behavior score)');
    [~, lgd_data] = legend(temp_legend, 'Location', 'bestoutside', 'FontSize', 18);
    lgd_line = findobj(lgd_data, 'type', 'line'); %// objects of legend of type line
    set(lgd_line, 'Markersize', 20); %// set marker size as desired
    
    set(gcf, 'Position', get(0, 'Screensize'));
    set(gcf,'PaperPositionMode','auto')
    print(f, [collection_folder, '/latent_scores_combined_LV_color'], '-dpng', '-r0');
    saveas(f, [collection_folder, '/latent_scores_combined_LV_color.fig']);
    
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
                x = output.final_parameters{i,index_epsilon}(strcmp(output.cv_outer.TestInd{2,w}, input.selected_studygroups{1,ii}));
                y = output.final_parameters{i,index_omega}(strcmp(output.cv_outer.TestInd{2,w}, input.selected_studygroups{1,ii}));
                plot(x,y,'.', 'MarkerSize', 20, 'color',colorpattern_LS(ii,:));
                hold on
            end
        end
    end
    first_line = strrep([input.name, ', grid_density=' num2str(input.grid_density)], '_', ' ');
    second_line = 'latent scores from all LVs combined';
    %             if output.final_parameters{i,opt_p}<=output.pvalue_FDR(i)
    %                 third_line = ['significant'];
    %             else
    %                 third_line = ['not significant'];
    %             end
    title({first_line; second_line}); % add third line
    xlabel('epsilon(brain score)');
    ylabel('omega(behavior score)');
    [~, lgd_data] = legend({input.selected_studygroups{1,:}}, 'FontSize', 18);
    lgd_line = findobj(lgd_data, 'type', 'line'); %// objects of legend of type line
    set(lgd_line, 'Markersize', 20); %// set marker size as desired
    
    set(gcf, 'Position', get(0, 'Screensize'));
    set(gcf,'PaperPositionMode','auto')
    %     print(f, [results_folder, '/latent_scores_LV' num2str(i)], '-dpng', '-r0');
    print(f, [collection_folder, '/latent_scores_combined_diagnosis_color'], '-dpng', '-r0');
    %     saveas(f, [results_folder, '/latent_scores_LV' num2str(i), '.fig']);
    saveas(f, [collection_folder, '/latent_scores_combined_diagnosis_color', '.fig']);
    
    %     saveas(f, [results_folder, '/latent_scores_LV' num2str(i), '.fig']);
    %     saveas(f, [results_folder, '/latent_scores_LV' num2str(i), '.png']);
    close(f);
    
end


if any(strcmp(IN.specific, 'detailed') | all_jobs)
    
    detailed_results_folder = [IN.results_path(1:(strfind(IN.results_path, 'final')-1)), 'detailed_results'];
    cd(detailed_results_folder);
    
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
        
        collection_folder_opt = [collection_folder, '/opt_parameters_' num2str(i)];
        mkdir(collection_folder_opt);
        cd(collection_folder_opt);
        
        %     for ii=1:size(temp_opt_param,1)
        %         %% visualize results
        %         % write brain vector to nifti file
        %         nk_WriteVol(temp_opt_param{ii,opt_param_u}, ['brain_opt_parameters_' num2str(i), '_' num2str(temp_opt_param{ii,1})], 2, NM.brainmask{1,1}, cell2mat(NM.badcoords), 0, 'gt');
        % %         dp_resample_image([collection_folder_opt, '/brain_opt_parameters_' num2str(i), '_' num2str(temp_opt_param{ii,1})], [1 1 1]);
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
            xlabel({'\color{black}weight vector v'}, 'FontSize', 6);
            ylabel('score', 'FontSize', 6);
            if temp_opt_param{ii,opt_param_p}<=output.pvalue_FDR(i)
                significance_opt = 'significant';
            else
                significance_opt = 'not significant';
            end
            subplot_title = {['Iteration ', num2str(ii), ', p-value (FDR-corrected) = ' num2str(temp_opt_param{ii,opt_param_p}), ', Spearman''s RHO = ', num2str(temp_opt_param{ii,opt_param_RHO}), significance_opt]};
            title(subplot_title, 'FontSize', 6, 'FontWeight', 'normal');
            
        end
        %     set(gcf, 'Position', get(0, 'Screensize'));
        set(gcf, 'Position', get(0, 'Screensize'));
        set(gcf,'PaperPositionMode','auto')
        first_line = strrep([input.name, ', grid_density=' num2str(input.grid_density), ', LV ',num2str(i)], '_', ' ');
        second_line = ['p-value (FDR-corrected) = ' num2str(output.final_parameters{i,opt_p}), ', Spearman''s RHO = ', num2str(output.final_parameters{i,opt_RHO})];
        if output.final_parameters{i,opt_p}<=output.pvalue_FDR(i)
            third_line = ['significant'];
        else
            third_line = ['not significant'];
        end
        suptitle({first_line; [second_line, ', ' third_line]}); % add third line
        
        %     print([detailed_results_folder, '/behavior_opt_parameters_',num2str(i)], '-dpng', '-r0');
        %     saveas(f, [detailed_results_folder, '/behavior_opt_parameters_',num2str(i)], 'png');
        
        print([collection_folder_opt, '/behavior_opt_parameters_',num2str(i)], '-dpng', '-r0');
        saveas(f, [collection_folder_opt, '/behavior_opt_parameters_',num2str(i)], 'png');
        
        %     saveas(f, [detailed_results_folder, '/behavior_opt_parameters_',num2str(i)], 'fig');
        %     saveas(f, [detailed_results_folder, '/behavior_opt_parameters_',num2str(i)], 'png');
        close all;
    end
    
end

end