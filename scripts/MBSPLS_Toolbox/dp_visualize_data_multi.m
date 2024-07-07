%% function to visualize SPLS output

function dp_visualize_data_multi(IN)

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

marker_size = 14;
font_size = 14;
LV_x = 'clinical features';
LV_y = 'weights';
LS_epsilon = 'brain score';
LS_omega = 'behavior score';

%% visualize results
load(IN.results_path);
if strfind(IN.results_path, '2018')
    analysis_date = '2018';
elseif strfind(IN.results_path, '2019')
    analysis_date = '2019';
elseif strfind(IN.results_path, '2020')
    analysis_date = '2020';
end
folder_name = [setup.date, '-', IN.results_path(5+strfind(IN.results_path, analysis_date):(strfind(IN.results_path, '/final_results')-1))];
collection_folder = [overall_folder, '/', folder_name];
mkdir(collection_folder);

%% compute sociodemographic and clinical outcome correlations
if any(strcmp(IN.specific, 'sociodemographic') | all_jobs)
    [input, output, setup] = dp_sociodemographic(IN);
    
    subsets = {'all', 'hold_out'};
    
    corr_folder = [collection_folder, '/correlations'];
    mkdir(corr_folder);
    
    for s=1:size(subsets,2)
        fields=fieldnames(output.(['tables_', subsets{s}, '_Rho_p']));
        for i=1:size(fields,1)
            latent_scores = fieldnames(output.(['tables_', subsets{s}, '_Rho_p']).(fields{i}));
            for ii=1:size(latent_scores,1)
                %             dp_txt_write(s_folder, ['/corr_RHO_', fields{i}, '_', latent_scores{ii}], output.tables_ho_RHO.(fields{i}).(latent_scores{ii})', '%.3f \t %.3f \t %.3f \t %.3f \t %.3f \n');
                %             dp_txt_write(s_folder, ['/corr_p_', fields{i}, '_', latent_scores{ii}], output.tables_ho_p.(fields{i}).(latent_scores{ii})', '%.3f \t %.3f \t %.3f \t %.3f \t %.3f \n');
                dp_txt_write(corr_folder, ['/corr_all_', fields{i}, '_', latent_scores{ii}], output.(['tables_', subsets{s}, '_Rho_p']).(fields{i}).(latent_scores{ii})', '%.3f \t %.3f \t %.3f \t %.3f \t %.3f \n');
            end
        end
    end
    
end

input.behavior_names = strrep(input.behavior_names, '_T0', '');
load(input.NM_structure);

%% prepare questionnaire data
load('/volume/HCStress/Doc/Stress_Resilience_questionnaires.mat');
Resilience_Stress_questionnaires = [CISS_questionnaire', RSA_questionnaire', CTQ_questionnaire', BS_questionnaire'];

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
    i_folder = [collection_folder, '/images'];
    mkdir(i_folder);
    cd(i_folder);
    % write brain vector to nifti file
    for i=1:size(output.final_parameters,1)
        nk_WriteVol(output.final_parameters{i,4}, ['brain_LV' num2str(i)], 2, NM.brainmask{1,1}, cell2mat(NM.badcoords), 0, 'gt');
        %     dp_resample_image([collection_folder, '/brain_LV' num2str(i)], [1 1 1]);
    end
end


if any(strcmp(IN.specific, 'atlas') | all_jobs)
    a_folder = [collection_folder, '/atlases'];
    mkdir(a_folder);
    atlases = {'brainnetome', 'cerebellum'};
    % get clusters for brain regions using hammers and aal atlas
    % filepath hammers nifti: /opt/SPM/spm12_v6685_cat12_r1207/atlas/hammers.nii
    % filepath hammers description: /opt/SPM/spm12_v6685_cat12_r1207/atlas/labels_dartel_hammers.xml
    
    for aa=1:size(atlases,2)
        
        switch atlases{aa}
            case 'hammers'
                switch size(input.behavior,1)
                    case 626
                        atlas_path_full = '/volume/HCStress/Data/MRI/Hammers_Atlas/Hammers_mith-n30r95_spm12_full__626_NM_X.mat';
                    case 627
                        atlas_path_full = '/volume/HCStress/Data/MRI/Hammers_Atlas/Hammers_mith-n30r95_spm12_full__627_NM_X.mat';
                    case 630
                        atlas_path_full = '/volume/HCStress/Data/MRI/Hammers_Atlas/Hammers_mith-n30r95_spm12_full__630_CISS_RSA_NM_X.mat';
                        a_number = 1;
                    case 631
                        atlas_path_full = '/volume/HCStress/Data/MRI/Hammers_Atlas/Hammers_mith-n30r95_spm12_full__631_NM_X.mat';
                    case 634
                        atlas_path_full = '/volume/HCStress/Data/MRI/Hammers_Atlas/Hammers_mith-n30r95_spm12_full__634_NM_X.mat';
                    case 649
                        atlas_path_full = '/volume/HCStress/Data/MRI/hammers_3mm_649_X.mat';
                        a_number = 1;
                end
                temp = load('/volume/HCStress/Data/MRI/Hammers_Atlas/Hammers_mith-n30-ancillary-data.tar/Hammers_mith-n30-ancillary-data/Hammers_mith-n30-ancillary-data/labels_hammerssmith_n30_ancillary_data.mat');
                fields = fieldnames(temp);
                labels_regions = temp.(fields{1});
                
                for i=1:size(labels_regions,1)
                    strf = strfind(labels_regions{i,2}, ' ');
                    labels_regions{i,3} = labels_regions{i,2}(1:(strf-1));
                end
                
            case 'juelich'
                switch size(input.behavior,1)
                    case 630
                        atlas_path_full = '/volume/HCStress/Data/MRI/Julich_Atlas/juelich_maxprob_0_r3mm_630_CISS_RSA_NM_X.mat';
                        a_number = 1;
                    case 627
                        atlas_path_full = '/volume/HCStress/Data/MRI/Julich_Atlas/juelich_maxprob_0_r3mm_627_NM_X.mat';
                    case 634
                        atlas_path_full = '/volume/HCStress/Data/MRI/Julich_Atlas/juelich_maxprob_0_r3mm_634_NM_X.mat';
                    case 649
                        atlas_path_full = '/volume/HCStress/Data/MRI/BNA_hammers_juelich_3mm_649_X.mat';
                        a_number = 3;
                end
                temp = load('/volume/HCStress/Data/MRI/Julich_Atlas/juelich_indices_labels.mat');
                fields = fieldnames(temp);
                labels_regions = temp.(fields{1});
                
                for i=1:size(labels_regions,1)
                    strf = strfind(labels_regions{i,2}, ' ');
                    labels_regions{i,3} = labels_regions{i,2}(1:(strf-1));
                end
                
            case 'brainnetome'
                switch size(input.behavior,1)
                    case 621
                        atlas_path_full = '/volume/HCStress/Data/MRI/BNA_CISSRSAPAS_3mm_621_X.mat';
                    case 627
                        atlas_path_full = '/volume/HCStress/Data/MRI/Brainnetome_Atlas/brainnetome_3mm_627_CTQ_BS_NM_X.mat';
                    case 630
                        atlas_path_full = '/volume/HCStress/Data/MRI/Brainnetome_Atlas/brainnetome_3mm_630_CISS_RSA_NM_X.mat';
                        a_number = 1;
                    case 634
                        atlas_path_full = '/volume/HCStress/Data/MRI/Brainnetome_Atlas/brainnetome_3mm_634_CISS_RSA_NM_X.mat';
                    case 649
                        atlas_path_full = '/volume/HCStress/Data/MRI/BNA_hammers_juelich_3mm_649_X.mat';
                        a_number = 1;
                end
                temp = load('/volume/HCStress/Data/MRI/Brainnetome_Atlas/brainnetome_indices.mat');
                fields = fieldnames(temp);
                labels_regions = temp.(fields{1});
                
                %                 for i=1:size(labels_regions,1)
                %                     strf = strfind(labels_regions{i,2}, '_');
                %                     labels_regions{i,3} = labels_regions{i,2}(1:(strf-1));
                %                 end
            case 'cerebellum'
                switch size(input.behavior,1)
                    case 621
                        atlas_path_full = '/volume/HCStress/Data/MRI/Cerebellum-MNIflirt-MRICroN/Cerebellum-MNIflirt_621_CISSRSAPAS_X.mat';
                    case 649
                        atlas_path_full = '/volume/HCStress/Data/MRI/Cerebellum-MNIflirt-MRICroN/Cerebellum-MNIflirt_649_CTQ_X.mat';
                end
                temp = load('/volume/HCStress/Data/MRI/Cerebellum-MNIflirt-MRICroN/CerebellumMNIflirt_indices.mat');
                fields = fieldnames(temp);
                labels_regions = temp.(fields{1});
                
                
        end
        
        temp=load(atlas_path_full);
        fields = fieldnames(temp);
        if exist('a_number', 'var')
            atlas_for_analysis = round(temp.(fields{1})(a_number,:));
        else
            atlas_for_analysis = round(temp.(fields{1}));
        end
        [C_atlas, ~, ic_atlas] = unique(atlas_for_analysis);
        counts_atlas = accumarray(ic_atlas, 1);
        output.regions.(atlases{aa}).regions_count = [C_atlas', counts_atlas];
        output.regions.(atlases{aa}).count_names = {'region_number', 'n_voxels', 'voxel_percentage', 'region_name', 'median_weights', 'mean_weights'};
        output.regions.(atlases{aa}).count = struct('raw', [], 'voxels', [], 'weights', []);
        output.regions.(atlases{aa}).count_names_sorting = {'positive', 'negative', 'overall'};
        voxel_cutoff = 0;
        fields = fieldnames(output.regions.(atlases{aa}).count);
        for i=1:size(fields,1)
            for fff = 1:size(output.regions.(atlases{aa}).count_names_sorting,2)
                FID = fopen([a_folder, '/brain_regions_', atlases{aa}, '_', output.regions.(atlases{aa}).count_names_sorting{fff}, '_', fields{i}, '.txt'], 'w');
                fprintf(FID, [strrep(input.name,'_',' '), '\n', fields{i} ' sorted']);
                fclose(FID);
            end
        end
        
        for i=1:size(output.final_parameters,1)
            output.regions.(atlases{aa}).log{i,1} = output.final_parameters{i,4}>0;
            output.regions.(atlases{aa}).log{i,2} = output.final_parameters{i,4}<0;
            output.regions.(atlases{aa}).log{i,3} = output.final_parameters{i,4}~=0;
            for ii=1:size(output.regions.(atlases{aa}).log,2)
                output.regions.(atlases{aa}).sum{i,ii} = atlas_for_analysis((atlas_for_analysis~=0)' & output.regions.(atlases{aa}).log{i,ii});
                [C, ~, ic] = unique(output.regions.(atlases{aa}).sum{i,ii});
                a_counts = accumarray(ic, 1);
                voxel_percentage = a_counts./(counts_atlas(ismember(C_atlas, C)));
                output.regions.(atlases{aa}).count.raw{i,ii} = [num2cell(C'), num2cell(a_counts), num2cell(voxel_percentage), labels_regions((C'),2)];
                fix_end = size(output.regions.(atlases{aa}).count.raw{i,ii},2);
                log_cutoff = cell2mat(output.regions.(atlases{aa}).count.raw{i,ii}(:, find(strcmp(output.regions.(atlases{aa}).count_names, 'n_voxels')))) <voxel_cutoff;
                output.regions.(atlases{aa}).count.raw{i,ii}(log_cutoff,:) = [];
                if isempty(output.regions.(atlases{aa}).count.raw{i,ii})
                    output.regions.(atlases{aa}).count.raw{i,ii} = [];
                    output.regions.(atlases{aa}).count.voxels{i,ii} = [];
                    output.regions.(atlases{aa}).count.weights{i,ii} = [];
                else
                    for iii=1:size(output.regions.(atlases{aa}).count.raw{i,ii},1)
                        if strcmp(output.regions.(atlases{aa}).count_names_sorting{ii}, 'overall')
                            output.regions.(atlases{aa}).count.raw{i,ii}{iii,fix_end+1} = median(abs(output.final_parameters{i,4}(output.regions.(atlases{aa}).log{i,ii} & (atlas_for_analysis==output.regions.(atlases{aa}).count.raw{i,ii}{iii,1})')));
                            output.regions.(atlases{aa}).count.raw{i,ii}{iii,fix_end+2} = mean(abs(output.final_parameters{i,4}(output.regions.(atlases{aa}).log{i,ii} & (atlas_for_analysis==output.regions.(atlases{aa}).count.raw{i,ii}{iii,1})')));
                            %                         output.regions.(atlases{aa}).count.raw{i,ii}{iii,6} =
                        else
                            output.regions.(atlases{aa}).count.raw{i,ii}{iii,fix_end+1} = median(output.final_parameters{i,4}(output.regions.(atlases{aa}).log{i,ii} & (atlas_for_analysis==output.regions.(atlases{aa}).count.raw{i,ii}{iii,1})'));
                            output.regions.(atlases{aa}).count.raw{i,ii}{iii,fix_end+2} = mean(output.final_parameters{i,4}(output.regions.(atlases{aa}).log{i,ii} & (atlas_for_analysis==output.regions.(atlases{aa}).count.raw{i,ii}{iii,1})'));
                        end
                    end
                    
                    if ~strcmp(output.regions.(atlases{aa}).count_names_sorting{ii}, 'negative')
                        output.regions.(atlases{aa}).count.voxels{i,ii} = flipud(sortrows((output.regions.(atlases{aa}).count.raw{i,ii}), find(strcmp(output.regions.(atlases{aa}).count_names, 'n_voxels'))));
                        output.regions.(atlases{aa}).count.weights{i,ii} = flipud(sortrows((output.regions.(atlases{aa}).count.raw{i,ii}), find(strcmp(output.regions.(atlases{aa}).count_names, 'median_weights'))));
                    else
                        output.regions.(atlases{aa}).count.voxels{i,ii} = flipud(sortrows((output.regions.(atlases{aa}).count.raw{i,ii}), find(strcmp(output.regions.(atlases{aa}).count_names, 'n_voxels'))));
                        output.regions.(atlases{aa}).count.weights{i,ii} = sortrows((output.regions.(atlases{aa}).count.raw{i,ii}), find(strcmp(output.regions.(atlases{aa}).count_names, 'median_weights')));
                    end
                    
                end
            end
            
            temp_add = [];
            temp_print = output.regions.(atlases{aa}).count_names([1,2,4,5,3]);
            for f=1:size(temp_print,2)
                %         try output.regions.(atlases{aa}).count_names{f} = num2str(temp_print{f});
                %         catch ME
                %         end
                temp_add = [temp_add, sprintf('\t'), temp_print{f}];
            end
            
            
            fields = fieldnames(output.regions.(atlases{aa}).count);
            for ff=1:size(fields,1)
                for fff = 1:size(output.regions.(atlases{aa}).count.(fields{ff}),2)
                    FID = fopen([a_folder, '/brain_regions_', atlases{aa}, '_', output.regions.(atlases{aa}).count_names_sorting{fff}, '_', fields{ff}, '.txt'], 'a');
                    fprintf(FID, ['\n \n \n', sprintf('\t'), 'Latent Variable ' num2str(i), '\n \n', temp_add, '\n']);
                    fclose(FID);
                end
            end
            
            
            fields = fieldnames(output.regions.(atlases{aa}).count);
            for ff=1:size(fields,1)
                for fff = 1:size(output.regions.(atlases{aa}).count.(fields{ff}),2)
                    if ~isempty(output.regions.(atlases{aa}).count.(fields{ff}){i,fff})
                        for r=1:size(output.regions.(atlases{aa}).count.(fields{ff}){i,fff},1)
                            temp_print = output.regions.(atlases{aa}).count.(fields{ff}){i,fff}(r,[1,2,4,5,3]);
                            temp_add = [];
                            for ii=1:size(temp_print,2)
                                try temp_print{ii} = num2str(temp_print{ii});
                                catch ME
                                end
                                temp_add = [temp_add, sprintf('\t'), temp_print{ii}];
                            end
                            
                            FID = fopen([a_folder, '/brain_regions_', atlases{aa}, '_', output.regions.(atlases{aa}).count_names_sorting{fff}, '_', fields{ff}, '.txt'], 'a');
                            fprintf(FID, ['\n' temp_add]);
                            fclose(FID);
                        end
                    end
                end
            end
            
            % %     dp_resample_image([collection_folder, '/brain_LV' num2str(i)], [1 1 1]);
            
            % just deactivated for paper preparation, needs to be activated
            % later on
            %             if i==size(output.final_parameters,1)
            %                 temp =  output.regions.(atlases{aa}).count.voxels;
            %                 for b=1:size(temp,1)
            %                     for bb=1:size(temp,2)
            %                         mat = temp{b,bb};
            %                         temp_mat=[];
            %                         while size(mat,1)>0
            %                             strf = strfind(mat{1,4}, ' ');
            %                             temp_prefix = mat{1,4}(1:(strf-1));
            %                             strf_log = ~cellfun(@isempty, strfind(mat(:,4), temp_prefix));
            %                             temp_mat = [temp_mat; mat(strf_log,:)];
            %                             output.regions.(atlases{aa}).count.alphabet{b,bb} = temp_mat;
            %                             mat(strf_log,:)=[];
            %                         end
            %                     end
            %                 end
            %
            %
            %
            %
            %                 for fi=1:size(output.regions.(atlases{aa}).count.alphabet,1)
            %                     temp_add = [];
            %                     temp_print = output.regions.(atlases{aa}).count_names([1,2,4,5,3]);
            %                     for f=1:size(temp_print,2)
            %                         %         try output.regions.(atlases{aa}).count_names{f} = num2str(temp_print{f});
            %                         %         catch ME
            %                         %         end
            %                         temp_add = [temp_add, sprintf('\t'), temp_print{f}];
            %                     end
            %                     %                     fields = fieldnames(output.regions.(atlases{aa}).count);
            %                     %                     for ff=1:size(fields,1)
            %                     for fff = 1:size(output.regions.(atlases{aa}).count.alphabet,2)
            %                         FID = fopen([a_folder, '/brain_regions_', atlases{aa}, '_', output.regions.(atlases{aa}).count_names_sorting{fff}, '_alphabet.txt'], 'a');
            %                         fprintf(FID, ['\n \n \n', sprintf('\t'), 'Latent Variable ' num2str(fi), '\n \n', temp_add, '\n']);
            %                         fclose(FID);
            %                     end
            %                     %                     end
            %
            %                     %                     fields = fieldnames(output.regions.(atlases{aa}).count);
            %                     %                     for ff=1:size(fields,1)
            %                     for fff = 1:size(output.regions.(atlases{aa}).count.alphabet,2)
            %                         if ~isempty(output.regions.(atlases{aa}).count.alphabet{fi,fff})
            %                             for r=1:size(output.regions.(atlases{aa}).count.alphabet{fi,fff},1)
            %                                 temp_print = output.regions.(atlases{aa}).count.alphabet{fi,fff}(r,[1,2,4,5,3]);
            %                                 temp_add = [];
            %                                 for ii=1:size(temp_print,2)
            %                                     try temp_print{ii} = num2str(temp_print{ii});
            %                                     catch ME
            %                                     end
            %                                     temp_add = [temp_add, sprintf('\t'), temp_print{ii}];
            %                                 end
            %
            %                                 FID = fopen([a_folder, '/brain_regions_', atlases{aa}, '_', output.regions.(atlases{aa}).count_names_sorting{fff}, '_alphabet.txt'], 'a');
            %                                 fprintf(FID, ['\n' temp_add]);
            %                                 fclose(FID);
            %                             end
            %                         end
            %                     end
            %                     %                     end
            %                 end
            %
            %             end
            
        end
        % just deactivated for paper preparation, needs to be reactivated
        % later on
        %         load(input.MRI)
        %         MRI_volumes = MRI_for_analysis;
        %         for i=1:size(output.final_parameters,1)
        %             u=output.final_parameters{i,4};
        %             for ii=1:size(output.regions.(atlases{aa}).log,2)
        %                 log_weights = output.regions.(atlases{aa}).log{i,ii}';
        %                 for iii=1:size(output.regions.(atlases{aa}).count.alphabet{i,ii},1)
        %                     log_region = atlas_for_analysis == output.regions.(atlases{aa}).count.alphabet{i,ii}{iii,1};
        %                     output.volumes.(atlases{aa}).names{i,ii} = [output.regions.(atlases{aa}).count.alphabet{i,ii}(:,1)';output.regions.(atlases{aa}).count.alphabet{i,ii}(:,4)'];
        %                     for s=1:size(MRI_volumes,1)
        %                         output.volumes.(atlases{aa}).raw{i,ii}(s,iii) = sum(MRI_volumes(s,log_region&log_weights));
        %                         %                     output.volumes.(['LV_', num2str(i)]).weighted.(atlases{aa}){s,iii} = sum(output.final_parameters{i,4}(log_region&log_weights)' .* MRI_volumes(s,log_region&log_weights));
        %                     end
        %                 end
        %             end
        %             MRI_volumes = MRI_volumes - (MRI_volumes*u)*u';
        %         end
        
    end
    
    save(IN.results_path, 'input', 'output', 'setup');
    %             end
    
end


if any(strcmp(IN.specific, 'behavior') | all_jobs)
    
    b_folder = [collection_folder, '/behavior'];
    mkdir(b_folder);
    
    %% visualize behavior vector in barplot
    % CTQ.emotional_abuse)), (CTQ.physical_abuse)), (CTQ.sexual_abuse)), (CTQ.emotional_neglect)), (CTQ.physical_neglect)), (CTQ.denial))];
    % load([setup.analysis_folder, '/' setup.date, '_', name, '_data_collection.mat']);
    %     if sum(ismember(input.behavior_names, 'Age')) + sum(ismember(input.behavior_names, 'Sex'))==2
    %         selected_variables = [input.selected_features(1,:), 'Age', 'Sex'; input.selected_features(2,:), 1, 1];
    %     elseif sum(ismember(input.behavior_names, 'Age'))
    %         selected_variables = [input.selected_features(1,:), 'Age'; input.selected_features(2,:), 1];
    %     elseif sum(ismember(input.behavior_names, 'Sex'))
    %         selected_variables = [input.selected_features(1,:), 'Sex'; input.selected_features(2,:), 1];
    %     else
    %         selected_variables = input.selected_features;
    %     end
    
    log_f=[];
    for i=1:size(input.selected_features,2)
        log_f(i,:) = ~cellfun(@isempty,(strfind(input.behavior_names(1,:), input.selected_features{1,i})));
    end
    
    if size(log_f,1)>1
        log_c = sum(sum(log_f==0)==size(input.selected_features,2));
    else
        log_c = sum((log_f==0)==size(input.selected_features,2));
    end
    
    selected_variables = [[input.selected_features(1,:), input.behavior_names((end-(log_c-1)):end)]; [input.selected_features(2,:), num2cell(ones(1,log_c))]];
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
    colorpattern_LV = colorcube((count + log_c));
    
    for i=1:(size(output.final_parameters,1))
        x = output.final_parameters{i,opt_v};
        output.questions_collection.(['LV_', num2str(i)]).items = input.behavior_names(x~=0)';
        output.questions_collection.(['LV_', num2str(i)]).subscales = input.subscales(x~=0)';
        for q=1:size(output.questions_collection.(['LV_', num2str(i)]).items,1)
            output.questions_collection.(['LV_', num2str(i)]).questions{q} = Resilience_Stress_questionnaires(1,strcmp(Resilience_Stress_questionnaires(2,:), output.questions_collection.(['LV_', num2str(i)]).items{q}));
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
                    temp_legend{nn} = strrep(selected_variables{1,ii}, '_', ' ');
                    hold on
            end
        end
        
        try grid_x = input.grid_dynamic.(['LV_', num2str(i)]).x;
            grid_x = input.grid_dynamic.(['LV_', num2str(i)]).x;
            grid_y = input.grid_dynamic.(['LV_', num2str(i)]).y;
        catch
            try grid_x = input.grid_dynamic.(['LV_', num2str(i-1)]).x;
                grid_y = input.grid_dynamic.(['LV_', num2str(i-1)]).y;
            catch
                input.grid_dynamic.(['LV_', num2str(i-1)]).x = input.grid_dynamic.(['LV_', num2str(i-2)]).x;
                input.grid_dynamic.(['LV_', num2str(i-1)]).y = input.grid_dynamic.(['LV_', num2str(i-2)]).y;
                grid_x = input.grid_dynamic.(['LV_', num2str(i-1)]).x;
                grid_y = input.grid_dynamic.(['LV_', num2str(i-1)]).y;
            end
        end
            
            first_line = strrep([input.name, ' grid_density_x=' num2str(grid_x.density), ', grid_density_y=' num2str(grid_y.density), ', LV ',num2str(i)], '_', ' ');
            second_line = ['p-value (FDR-corrected) = ' num2str(output.final_parameters{i,opt_p}), ', Spearman''s RHO = ', num2str(output.final_parameters{i,opt_RHO})];
            if output.final_parameters{i,opt_p}<=output.pvalue_FDR(i)
                third_line = ['significant'];
            else
                third_line = ['not significant'];
            end
            title({first_line; second_line; third_line}); % add third line
            xlabel(LV_x, 'FontSize', font_size);
            ylabel(LV_y, 'FontSize', font_size);
            %         if min(x)>=0
            %             axis([0 (size(output.final_parameters{i,opt_v},1)+1) 0 1]);
            %         elseif max(x)<=0
            %             axis([0 (size(output.final_parameters{i,opt_v},1)+1) 0 -1]);
            %         end
            legend(temp_handle, temp_legend, 'Location', 'bestoutside', 'FontSize', font_size);
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
                    %                 htext = text(x_pos(xi),ypos(xi),strrep(input.behavior_names{xi},'_',' '));  % Add text label
                    %                 set(htext,'VerticalAlignment','bottom','HorizontalAlignment','center', 'FontSize', 8); % Adjust properties
                end
            end
            
            
            %         try annotation(f, 'textbox', [0.79, 0.2, 0.16, 0.4], 'string', strrep(output.questions_collection.(['LV_', num2str(i)]).final, '_', ' '), 'FontSize', 8, 'FitHeightToText', 'on');
            %         catch
            %         end
            
            set(gcf,'Position', get(0,'Screensize'));
            set(gcf,'PaperPositionMode','auto')
            print(f, [b_folder, '/behavior_LV_' num2str(i)], '-dpng', '-r0');
            saveas(f, [b_folder, '/behavior_LV' num2str(i), '.fig']);
            
            close(f);
            
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
        FID = fopen([b_folder, '/collected_questions.txt'], 'w');
        fclose(FID);
        
        for i=1:size(fields,1)
            FID = fopen([b_folder, '/collected_questions.txt'], 'a');
            fprintf(FID, ['\n \n \nLV_', num2str(i)]);
            fclose(FID);
            temp1 = output.questions_collection.(['LV_', num2str(i)]).combined(1,:);
            temp2 = output.questions_collection.(['LV_', num2str(i)]).combined(2,:);
            temp_unique = unique(output.questions_collection.LV_1.combined(2,:));
            for ii=1:size(temp_unique,2)
                FID = fopen([b_folder, '/collected_questions.txt'], 'a');
                fprintf(FID, ['\n \n', temp_unique{ii}]);
                fclose(FID);
                temp_unique_log = strcmp(temp2, temp_unique{ii});
                temp_unique_coll = temp1(temp_unique_log);
                for iii=1:size(temp_unique_coll,2)
                    FID = fopen([b_folder, '/collected_questions.txt'], 'a');
                    fprintf(FID, ['\n', temp_unique_coll{iii}]);
                    fclose(FID);
                end
            end
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
                x = output.final_parameters{i,index_epsilon}(strcmp(output.CV.cv_outer_indices.TestInd{2,w}, input.selected_studygroups{1,ii}));
                y = output.final_parameters{i,index_omega}(strcmp(output.CV.cv_outer_indices.TestInd{2,w}, input.selected_studygroups{1,ii}));
                plot(x,y,'.', 'MarkerSize', marker_size, 'color',colorpattern_LS(ii,:));
                hold on
            end
            first_line = strrep([input.name, ' grid_density_x=' num2str(grid_x.density), ', grid_density_y=' num2str(grid_y.density), ', LV ',num2str(i)], '_', ' ');
            second_line = ['p-value (FDR-corrected) = ' num2str(output.final_parameters{i,opt_p}), ', Spearman''s RHO = ', num2str(output.final_parameters{i,opt_RHO})];
            if output.final_parameters{i,opt_p}<=output.pvalue_FDR(i)
                third_line = 'significant';
            else
                third_line = 'not significant';
            end
            title({first_line; second_line; third_line}); % add third line
            xlabel(LS_epsilon);
            ylabel(LS_omega);
            [~, lgd_data] = legend({input.selected_studygroups{1,:}}, 'FontSize', font_size);
            lgd_line = findobj(lgd_data, 'type', 'line'); %// objects of legend of type line
            set(lgd_line, 'Markersize', marker_size); %// set marker size as desired
            
            set(gcf, 'Position', get(0, 'Screensize'));
            set(gcf,'PaperPositionMode','auto')
            %     print(f, [results_folder, '/latent_scores_LV' num2str(i)], '-dpng', '-r0');
            print(f, [b_folder, '/latent_scores_LV' num2str(i)], '-dpng', '-r0');
            %     saveas(f, [results_folder, '/latent_scores_LV' num2str(i), '.fig']);
            saveas(f, [b_folder, '/latent_scores_LV' num2str(i), '.fig']);
            
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
        [~, lgd_data] = legend(temp_legend, 'Location', 'bestoutside', 'FontSize', font_size);
        lgd_line = findobj(lgd_data, 'type', 'line'); %// objects of legend of type line
        set(lgd_line, 'Markersize', marker_size); %// set marker size as desired
        
        set(gcf, 'Position', get(0, 'Screensize'));
        set(gcf,'PaperPositionMode','auto')
        print(f, [b_folder, '/latent_scores_combined_LV_color'], '-dpng', '-r0');
        saveas(f, [b_folder, '/latent_scores_combined_LV_color.fig']);
        
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
                    x = output.final_parameters{i,index_epsilon}(strcmp(output.CV.cv_outer_indices.TestInd{2,w}, input.selected_studygroups{1,ii}));
                    y = output.final_parameters{i,index_omega}(strcmp(output.CV.cv_outer_indices.TestInd{2,w}, input.selected_studygroups{1,ii}));
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
        
        %     saveas(f, [results_folder, '/latent_scores_LV' num2str(i), '.fig']);
        %     saveas(f, [results_folder, '/latent_scores_LV' num2str(i), '.png']);
        close(f);
        
        
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
            
            %     saveas(f, [detailed_results_folder, '/behavior_opt_parameters_',num2str(i)], 'fig');
            %     saveas(f, [detailed_results_folder, '/behavior_opt_parameters_',num2str(i)], 'png');
            close all;
        end
        
    end
    
    if any(strcmp(IN.specific, 'correlation') | all_jobs)
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
        s_folder = [collection_folder, '/correlation/' subsets{s}];
        mkdir(s_folder);
        
        
        for i=1:size(output.final_parameters,1)
            fields = fieldnames(output.post_hoc_correlations.data_collection]).(['LV_', num2str(i)]));
            for ii=1:size(input.selected_studygroups,2)
                p_value=output.([subsets{s}, '_correlations']).(['LV_', num2str(i)]).epsilon.(input.selected_studygroups{1,ii}).(fields{iii})(2);
                
                if p_value<=FDR_value
                    %                 h1 = lsline;
                    %                 h1.Color = 'k';
                    %                 h1.LineWidth = 2;
                    FDR_value = output.([subsets{s}, '_FDR_values']).(['LV_', num2str(i)]).epsilon.(input.selected_studygroups{1,ii});
                    RHO_value = output.([subsets{s}, '_correlations']).(['LV_', num2str(i)]).epsilon.(input.selected_studygroups{1,ii}).(fields{iii})(1);
                    first_line = strrep([input.name, ' grid_density_x=' num2str(grid_x.density), ' grid_density_x=' num2str(grid_x.density), ' grid_density_y=' num2str(grid_y.density), ', LV ',num2str(i)], '_', ' ');
                    second_line = strrep([input.selected_studygroups{1,ii}, ', epsilon x ', fields{iii}, ', p-value = ' num2str(p_value), ', Spearman''s RHO = ', num2str(RHO_value)], '_', ' ');
                    if output.final_parameters{i,opt_p}<=output.pvalue_FDR(i)
                        third_line1 = 'LV significant';
                    else
                        third_line1 = 'LV not significant';
                    end
                    f=figure();
                    w = output.final_parameters{i,1};
                    x = output.final_parameters{i,index_epsilon}(strcmp(output.CV.cv_outer_indices.TestInd{2,w}, input.selected_studygroups{1,ii}));
                    y = output.([subsets{s}, '_corr_data']).(['LV_', num2str(i)]).(fields{iii})(strcmp(output.CV.cv_outer_indices.TestInd{2,w}, input.selected_studygroups{1,ii}));
                    third_line2 = 'correlation significant';
                    plot(x,y,'.', 'MarkerSize', marker_size, 'color','k');
                    title(['p-value = ' num2str(p_value), ', Spearman''s RHO = ', num2str(RHO_value)], 'FontSize', font_size); % add third line
                    xlabel(LS_epsilon, 'FontSize', font_size);
                    y_temp = strrep(fields{iii}, '_', ' ');
                    y_temp = strrep(y_temp, 'T0','');
                    y_temp = strrep(y_temp, 'Screening','');
                    ylabel(y_temp, 'FontSize', font_size);
                    [~, lgd_data] = legend({input.selected_studygroups{1,ii}}, 'FontSize', font_size);
                    lgd_line = findobj(lgd_data, 'type', 'line'); %// objects of legend of type line
                    set(lgd_line, 'Markersize', marker_size); %// set marker size as desired
                    
                    set(gcf, 'Position', get(0, 'Screensize'));
                    set(gcf,'PaperPositionMode','auto')
                    %     print(f, [results_folder, '/latent_scores_LV' num2str(i)], '-dpng', '-r0');
                    print(f, [s_folder, '/latent_scores_LV', num2str(i), '_', subsets{s},  '_', input.selected_studygroups{1,ii}, '_epsilon_', fields{iii}], '-dpng', '-r0');
                    %     saveas(f, [results_folder, '/latent_scores_LV' num2str(i), '.fig']);
                    %                         saveas(f, [s_folder, '/latent_scores_LV', num2str(i), '_', subsets{s},  '_', input.selected_studygroups{1,ii}, '_epsilon_', fields{iii}, '.fig']);
                    close(f);
                    
                else
                    third_line2 = 'correlation not significant';
                end
                
            end
            
            %             g=figure();
            %             a = output.final_parameters{i, index_epsilon};
            %             b = output.([subsets{s}, '_corr_data']).(['LV_', num2str(i)]).(fields{iii});
            %             plot(a,b,'.', 'MarkerSize', marker_size, 'color','k');
            %             hold on
            %             h1 = lsline;
            %             h1.Color = 'k';
            %             h1.LineWidth = 2;
            %             k1=h1;
            
            
            p_value=output.([subsets{s}, '_correlations']).(['LV_', num2str(i)]).epsilon.all.(fields{iii})(2);
            FDR_value = output.([subsets{s}, '_FDR_values']).(['LV_', num2str(i)]).epsilon.all;
            RHO_value = output.([subsets{s}, '_correlations']).(['LV_', num2str(i)]).epsilon.all.(fields{iii})(1);
            first_line = strrep([input.name, ' grid_density_x=' num2str(grid_x.density), ' grid_density_x=' num2str(grid_x.density), ' grid_density_y=' num2str(grid_y.density), ', LV ',num2str(i)], '_', ' ');
            second_line = strrep(['all groups, epsilon x ', fields{iii}, ', p-value = ' num2str(p_value), ', Spearman''s RHO = ', num2str(RHO_value)], '_', ' ');
            if output.final_parameters{i,opt_p}<=output.pvalue_FDR(i)
                third_line1 = 'LV significant';
            else
                third_line1 = 'LV not significant';
            end
            
            if p_value<=FDR_value
                third_line2 = 'correlation significant';
                f=figure();
                for ii=1:size(input.selected_studygroups,2)
                    w = output.final_parameters{i,1};
                    x = output.final_parameters{i,index_epsilon}(strcmp(output.CV.cv_outer_indices.TestInd{2,w}, input.selected_studygroups{1,ii}));
                    y = output.([subsets{s}, '_corr_data']).(['LV_', num2str(i)]).(fields{iii})(strcmp(output.CV.cv_outer_indices.TestInd{2,w}, input.selected_studygroups{1,ii}));
                    plot(x,y,'.', 'MarkerSize', marker_size, 'color',colorpattern_LS(ii,:));
                    %                 ax1.Visible = 'off';
                    hold on
                    %                     title({first_line; second_line; [third_line1, '   ', third_line2]}); % add third line
                    
                end
                xlabel(LS_epsilon, 'FontSize', font_size);
                y_temp = strrep(fields{iii}, '_', ' ');
                y_temp = strrep(y_temp, 'T0','');
                y_temp = strrep(y_temp, 'Screening','');
                ylabel(y_temp, 'FontSize', font_size);
                title(['p-value = ' num2str(p_value), ', Spearman''s RHO = ', num2str(RHO_value)], 'FontSize', font_size); % add third line
                [~, lgd_data] = legend({input.selected_studygroups{1,:}}, 'FontSize', font_size);
                lgd_line = findobj(lgd_data, 'type', 'line'); %// objects of legend of type line
                set(lgd_line, 'Markersize', marker_size); %// set marker size as desired
                hold on
                %             line(k1.XData, k1.YData, 'Color','k','LineWidth',2);
                
                set(gcf, 'Position', get(0, 'Screensize'));
                set(gcf,'PaperPositionMode','auto')
                %     print(f, [results_folder, '/latent_scores_LV' num2str(i)], '-dpng', '-r0');
                print(f, [s_folder, '/latent_scores_LV', num2str(i), '_', subsets{s},  '_', 'all_epsilon_', fields{iii}], '-dpng', '-r0');
                %     saveas(f, [results_folder, '/latent_scores_LV' num2str(i), '.fig']);
                %                     saveas(f, [s_folder, '/latent_scores_LV', num2str(i), '_', subsets{s},  '_', 'all_epsilon_', fields{iii}, '.fig']);
                
                %     saveas(f, [results_folder, '/latent_scores_LV' num2str(i), '.fig']);
                %     saveas(f, [results_folder, '/latent_scores_LV' num2str(i), '.png']);
                close all;
                
                
            else
                third_line2 = 'correlation not significant';
            end
            
            
            for iii=1:size(fields,1)
                for ii=1:size(input.selected_studygroups,2)
                    %                 h1 = lsline;
                    %                 h1.Color = 'k';
                    %                 h1.LineWidth = 2;
                    p_value = output.([subsets{s}, '_correlations']).(['LV_', num2str(i)]).omega.(input.selected_studygroups{1,ii}).(fields{iii})(2);
                    FDR_value = output.([subsets{s}, '_FDR_values']).(['LV_', num2str(i)]).omega.(input.selected_studygroups{1,ii});
                    RHO_value = output.([subsets{s}, '_correlations']).(['LV_', num2str(i)]).omega.(input.selected_studygroups{1,ii}).(fields{iii})(1);
                    first_line = strrep([input.name, ' grid_density_x=' num2str(grid_x.density), ', grid_density_y=' num2str(grid_y.density), ', LV ',num2str(i)], '_', ' ');
                    second_line = strrep([input.selected_studygroups{1,ii}, ', omega x ', fields{iii}, ', p-value = ' num2str(p_value), ', Spearman''s RHO = ', num2str(RHO_value)], '_', ' ');
                    if output.final_parameters{i,opt_p}<=output.pvalue_FDR(i)
                        third_line1 = 'LV significant';
                    else
                        third_line1 = 'LV not significant';
                    end
                    
                    if p_value<=FDR_value
                        third_line2 = 'correlation significant';
                        f=figure();
                        w = output.final_parameters{i,1};
                        x = output.final_parameters{i,index_omega}(strcmp(output.CV.cv_outer_indices.TestInd{2,w}, input.selected_studygroups{1,ii}));
                        y = output.([subsets{s}, '_corr_data']).(['LV_', num2str(i)]).(fields{iii})(strcmp(output.CV.cv_outer_indices.TestInd{2,w}, input.selected_studygroups{1,ii}));
                        plot(x,y,'.', 'MarkerSize', marker_size, 'color','k');
                        title(['p-value = ' num2str(p_value), ', Spearman''s RHO = ', num2str(RHO_value)], 'FontSize', font_size); % add third line
                        xlabel(LS_omega, 'FontSize', font_size);
                        y_temp = strrep(fields{iii}, '_', ' ');
                        y_temp = strrep(y_temp, 'T0','');
                        y_temp = strrep(y_temp, 'Screening','');
                        ylabel(y_temp, 'FontSize', font_size);                    [~, lgd_data] = legend({input.selected_studygroups{1,ii}}, 'FontSize', font_size);
                        lgd_line = findobj(lgd_data, 'type', 'line'); %// objects of legend of type line
                        set(lgd_line, 'Markersize', marker_size); %// set marker size as desired
                        
                        set(gcf, 'Position', get(0, 'Screensize'));
                        set(gcf,'PaperPositionMode','auto')
                        %     print(f, [results_folder, '/latent_scores_LV' num2str(i)], '-dpng', '-r0');
                        print(f, [s_folder, '/latent_scores_LV', num2str(i), '_', subsets{s},  '_', input.selected_studygroups{1,ii}, '_omega_', fields{iii}], '-dpng', '-r0');
                        %     saveas(f, [results_folder, '/latent_scores_LV' num2str(i), '.fig']);
                        %                         saveas(f, [s_folder, '/latent_scores_LV', num2str(i), '_', subsets{s},  '_', input.selected_studygroups{1,ii}, '_omega_', fields{iii}, '.fig']);
                        close(f);
                        
                    else
                        third_line2 = 'correlation not significant';
                    end
                    
                    
                end
                
                %             g=figure();
                %             a = output.final_parameters{i, index_omega};
                %             b = output.([subsets{s}, '_corr_data']).(['LV_', num2str(i)]).(fields{iii});
                %             plot(a,b,'.', 'MarkerSize', 20, 'color','k');
                %             hold on
                %             h1 = lsline;
                %             h1.Color = 'k';
                %             h1.LineWidth = 2;
                %             k1=h1;
                
                
                p_value=output.([subsets{s}, '_correlations']).(['LV_', num2str(i)]).omega.all.(fields{iii})(2);
                FDR_value = output.([subsets{s}, '_FDR_values']).(['LV_', num2str(i)]).omega.all;
                RHO_value = output.([subsets{s}, '_correlations']).(['LV_', num2str(i)]).omega.all.(fields{iii})(1);
                first_line = strrep([input.name, ' grid_density_x=' num2str(grid_x.density), ', grid_density_y=' num2str(grid_y.density), ', LV ',num2str(i)], '_', ' ');
                second_line = strrep(['all groups, omega x ', fields{iii}, ', p-value = ' num2str(p_value), ', Spearman''s RHO = ', num2str(RHO_value)], '_', ' ');
                if output.final_parameters{i,opt_p}<=output.pvalue_FDR(i)
                    third_line1 = 'LV significant';
                else
                    third_line1 = 'LV not significant';
                end
                
                if p_value<=FDR_value
                    third_line2 = 'correlation significant';
                    f=figure();
                    for ii=1:size(input.selected_studygroups,2)
                        w = output.final_parameters{i,1};
                        x = output.final_parameters{i,index_omega}(strcmp(output.CV.cv_outer_indices.TestInd{2,w}, input.selected_studygroups{1,ii}));
                        y = output.([subsets{s}, '_corr_data']).(['LV_', num2str(i)]).(fields{iii})(strcmp(output.CV.cv_outer_indices.TestInd{2,w}, input.selected_studygroups{1,ii}));
                        plot(x,y,'.', 'MarkerSize', marker_size, 'color',colorpattern_LS(ii,:));
                        hold on
                    end
                    
                    title(['p-value = ' num2str(p_value), ', Spearman''s RHO = ', num2str(RHO_value)], 'FontSize', font_size); % add third line
                    xlabel(LS_omega, 'FontSize', font_size);
                    y_temp = strrep(fields{iii}, '_', ' ');
                    y_temp = strrep(y_temp, 'T0','');
                    y_temp = strrep(y_temp, 'Screening','');
                    ylabel(y_temp, 'FontSize', font_size);
                    [~, lgd_data] = legend({input.selected_studygroups{1,:}}, 'FontSize', font_size);
                    lgd_line = findobj(lgd_data, 'type', 'line'); %// objects of legend of type line
                    set(lgd_line, 'Markersize', marker_size); %// set marker size as desired
                    hold on
                    %             line(k1.XData, k1.YData, 'Color','k','LineWidth',2);
                    
                    set(gcf, 'Position', get(0, 'Screensize'));
                    set(gcf,'PaperPositionMode','auto')
                    %     print(f, [results_folder, '/latent_scores_LV' num2str(i)], '-dpng', '-r0');
                    print(f, [s_folder, '/latent_scores_LV', num2str(i), '_', subsets{s},  '_', 'all_omega_', fields{iii}], '-dpng', '-r0');
                    %     saveas(f, [results_folder, '/latent_scores_LV' num2str(i), '.fig']);
                    %                     saveas(f, [s_folder, '/latent_scores_LV', num2str(i), '_', subsets{s},  '_', 'all_omega_', fields{iii}, '.fig']);
                    
                    %     saveas(f, [results_folder, '/latent_scores_LV' num2str(i), '.fig']);
                    %     saveas(f, [results_folder, '/latent_scores_LV' num2str(i), '.png']);
                    close all;
                    
                else
                    third_line2 = 'correlation not significant';
                end
            end
            
            % plot correlations according to main dimension, all groups
            % combined, colorcoded for dimensions
            % GAF
            %             fields = fieldnames(output.hold_out_corr_data.LV_1);
            %             dimensions_selection = {'GAF', 'GF', 'BDI', 'WHO'};
            %             for d=1:size(dimensions_selection,2)
            %                 leg_selection = ~cellfun(@isempty, strfind(fields, dimensions_selection{d}));
            %                 colorpattern_DD = hsv(sum(leg_selection));
            %                 for dd=1:size(fields,1)
            %                     if leg_selection(dd)
            %                         x1 = output.final_parameters{i, index_omega};
            %                         y = output.hold_out_corr_data.(['LV_', num2str(i)]).(fields{dd});
            %                         plot(x1,y,'.', 'MarkerSize', marker_size, 'color',colorpattern_DD(dd,:));
            %                         hold on
            %                     end
            %
            %                 end
            %             end
            
            %             p_value=output.([subsets{s}, '_correlations']).(['LV_', num2str(i)]).omega.all.(fields{iii})(2);
            %             FDR_value = output.([subsets{s}, '_FDR_values']).(['LV_', num2str(i)]).omega.all;
            %             RHO_value = output.([subsets{s}, '_correlations']).(['LV_', num2str(i)]).omega.all.(fields{iii})(1);
            %             first_line = strrep([input.name, ' grid_density_x=' num2str(grid_x.density), ', grid_density_y=' num2str(grid_y.density), ', LV ',num2str(i)], '_', ' ');
            %             second_line = strrep(['all groups, omega x ', fields{iii}, ', p-value = ' num2str(p_value), ', Spearman''s RHO = ', num2str(RHO_value)], '_', ' ');
            %             if output.final_parameters{i,opt_p}<=output.pvalue_FDR(i)
            %                 third_line1 = 'LV significant';
            %             else
            %                 third_line1 = 'LV not significant';
            %             end
            %
            %             if p_value<=FDR_value
            %                 third_line2 = 'correlation significant';
            %                 f=figure();
            %                 for ii=1:size(input.selected_studygroups,2)
            %                     w = output.final_parameters{i,1};
            %                     x = output.final_parameters{i,index_omega};
            %                     y = output.([subsets{s}, '_corr_data']).(['LV_', num2str(i)]).(fields{iii})(strcmp(output.CV.cv_outer_indices.TestInd{2,w}, input.selected_studygroups{1,ii}));
            %                     plot(x,y,'.', 'MarkerSize', marker_size, 'color',colorpattern_LS(ii,:));
            %                     hold on
            %                 end
            %
            %                 title(['p-value = ' num2str(p_value), ', Spearman''s RHO = ', num2str(RHO_value)], 'FontSize', font_size); % add third line
            %                 xlabel(LS_omega, 'FontSize', font_size);
            %                 y_temp = strrep(fields{iii}, '_', ' ');
            %                 y_temp = strrep(y_temp, 'T0','');
            %                 y_temp = strrep(y_temp, 'Screening','');
            %                 ylabel(y_temp, 'FontSize', font_size);
            %                 [~, lgd_data] = legend({input.selected_studygroups{1,:}}, 'FontSize', font_size);
            %                 lgd_line = findobj(lgd_data, 'type', 'line'); %// objects of legend of type line
            %                 set(lgd_line, 'Markersize', marker_size); %// set marker size as desired
            %                 hold on
            %                 %             line(k1.XData, k1.YData, 'Color','k','LineWidth',2);
            %
            %                 set(gcf, 'Position', get(0, 'Screensize'));
            %                 set(gcf,'PaperPositionMode','auto')
            %                 %     print(f, [results_folder, '/latent_scores_LV' num2str(i)], '-dpng', '-r0');
            %                 print(f, [s_folder, '/latent_scores_LV', num2str(i), '_', subsets{s},  '_', 'all_omega_', fields{iii}], '-dpng', '-r0');
            %                 %     saveas(f, [results_folder, '/latent_scores_LV' num2str(i), '.fig']);
            %                 %                     saveas(f, [s_folder, '/latent_scores_LV', num2str(i), '_', subsets{s},  '_', 'all_omega_', fields{iii}, '.fig']);
            %
            %                 %     saveas(f, [results_folder, '/latent_scores_LV' num2str(i), '.fig']);
            %                 %     saveas(f, [results_folder, '/latent_scores_LV' num2str(i), '.png']);
            %                 close all;
            
            
        end
        
        
        save(IN.results_path, 'input', 'output', 'setup');
        
    end
    
end