%% function to visualize SPLS output

function dp_visualize_data_Dev(IN)

load(IN.results_path);

switch IN.overall_analysis
    case 'Stress'
        overall_folder = '/volume/HCStress/Analysis/Stress';
    case 'Resilience'
        overall_folder = '/volume/HCStress/Analysis/Resilience';
    case 'Schizotypy'
        overall_folder = '/volume/MU_Pronia_PLS_Schizotypy/Analysis/SPLS/Schizotypy';
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

%% compute sociodemographic and clinical outcome correlations
if any(strcmp(IN.specific, 'sociodemographic') | all_jobs)
    [input, output, setup] = dp_sociodemographic_2020_full(IN);
end

if val_log
    output.final_parameters = output.validation_results;
end

input.Y_names = strrep(input.Y_names, '_T0', '');
load(input.NM_structure);

%% prepare questionnaire data
load('/volume/HCStress/Doc/all_questionnaires.mat', 'BS_questionnaire', 'CISS_questionnaire', 'CTQ_questionnaire', 'RSA_questionnaire', 'WSS_questionnaire');
all_questionnaires = [CISS_questionnaire', RSA_questionnaire', CTQ_questionnaire', BS_questionnaire', WSS_questionnaire'];

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
    atlases = {'brainnetome', 'cerebellum', 'glasser', 'yeo', 'buckner'};
    % get clusters for brain regions using hammers and aal atlas
    % filepath hammers nifti: /opt/SPM/spm12_v6685_cat12_r1207/atlas/hammers.nii
    % filepath hammers description: /opt/SPM/spm12_v6685_cat12_r1207/atlas/labels_dartel_hammers.xml
    
    for aa=1:size(atlases,2)
        
        switch atlases{aa}
            case 'glasser'
                switch size(input.behavior,1)
                    case 636
                        atlas_path_full = '/volume/HCStress/Data/MRI/Glasser_Atlas/glasser_atlas_CISS_636_NM_X.mat';
                    case 652
                        atlas_path_full = '/volume/HCStress/Data/MRI/Glasser_Atlas/glasser_atlas_WSS_652_NM_X.mat';
                end
                temp = load('/volume/HCStress/Data/MRI/Glasser_Atlas/glasser_indices.mat');
                fields = fieldnames(temp);
                labels_regions = temp.(fields{1});
            case 'yeo'
                switch size(input.behavior,1)
                    case 636
                        atlas_path_full = '/volume/HCStress/Data/MRI/Yeo_Atlas/Yeo_17strict_CISS_636_NM_X.mat';
                    case 652
                        atlas_path_full = '/volume/HCStress/Data/MRI/Yeo_Atlas/Yeo_17strict_WSS_652_NM_X.mat';
                end
                temp = load('/volume/HCStress/Data/MRI/Yeo_Atlas/yeo_buckner_indices.mat');
                fields = fieldnames(temp);
                labels_regions = temp.(fields{1});
                
            case 'buckner'
                switch size(input.behavior,1)
                    case 636
                        atlas_path_full = '/volume/HCStress/Data/MRI/Buckner_Atlas/Buckner_17tight_CISS_636_NM_X.mat';
                    case 652
                        atlas_path_full = '/volume/HCStress/Data/MRI/Buckner_Atlas/Buckner_17tight_WSS_652_NM_X.mat';
                end
                temp = load('/volume/HCStress/Data/MRI/Yeo_Atlas/yeo_buckner_indices.mat');
                fields = fieldnames(temp);
                labels_regions = temp.(fields{1});

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
                    case 636
                        atlas_path_full = '/volume/HCStress/Data/MRI/Brainnetome_Atlas/brainnetome_3mm_636_CISS_NM_X.mat';
                    case 649
                        atlas_path_full = '/volume/HCStress/Data/MRI/BNA_hammers_juelich_3mm_649_X.mat';
                        a_number = 1;
                    case 652
                        atlas_path_full = '/volume/HCStress/Data/MRI/Brainnetome_Atlas/brainnetome_3mm_652_WSS_NM_X.mat';
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
                    case 636
                        atlas_path_full = '/volume/HCStress/Data/MRI/Cerebellum-MNIflirt-MRICroN/Cerebellum-MNIflirt_636_CISS_X.mat';
                    case 649
                        atlas_path_full = '/volume/HCStress/Data/MRI/Cerebellum-MNIflirt-MRICroN/Cerebellum-MNIflirt_649_CTQ_X.mat';
                    case 634
                        atlas_path_full = '/volume/HCStress/Data/MRI/Cerebellum-MNIflirt-MRICroN/Cerebellum-MNIflirt_634_CISSRSA_X.mat';
                    case 652
                        atlas_path_full = '/volume/HCStress/Data/MRI/Cerebellum-MNIflirt-MRICroN/Cerebellum-MNIflirt_652_WSS_X.mat';
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
        
        if any(contains('buckner', atlases{aa}))
            log_no_regions = ismember(atlas_for_analysis, [1, 5, 14]);
            atlas_for_analysis(log_no_regions) = 0;
            [C_atlas, ~, ic_atlas] = unique(atlas_for_analysis);
            counts_atlas = accumarray(ic_atlas, 1);
            temp=[];
            for bb=1:17
                if any(bb==C_atlas)
                    temp(bb,:) = [C_atlas(bb==C_atlas), counts_atlas(bb==C_atlas)];
                else
                    temp(bb,:) = [bb, 0];
                end
            end
            temp = [C_atlas(1), counts_atlas(1); temp];
            C_atlas = temp(:,1)';
            counts_atlas = temp(:,2);
        end
        
        output.regions_max.(atlases{aa}) = [C_atlas(2:end)', counts_atlas(2:end)];

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
            
            save(IN.results_path, 'input', 'output', 'setup');
            %     dp_resample_image([collection_folder, '/brain_LV' num2str(i)], [1 1 1]);
            
            %             just deactivated for paper preparation, needs to be activated
            %             later on
            
            if i==size(output.final_parameters,1)
                if contains(atlases{aa}, 'brainnetome')
                    %                 if contains(atlases{aa}, 'cerebellum')
                    temp =  output.regions.(atlases{aa}).count.voxels;
                    for b=1:size(temp,1)
                        for bb=1:size(temp,2)
                            mat = temp{b,bb};
                            temp_mat=[];
                            while size(mat,1)>0
                                strf = strfind(mat{1,4}, ' ');
                                temp_prefix = mat{1,4}(1:(strf-1));
                                strf_log = ~cellfun(@isempty, strfind(mat(:,4), temp_prefix));
                                temp_mat = [temp_mat; mat(strf_log,:)];
                                output.regions.(atlases{aa}).count.alphabet{b,bb} = temp_mat;
                                mat(strf_log,:)=[];
                            end
                        end
                    end
                    
                else
                    temp =  output.regions.(atlases{aa}).count.voxels;
                    for b=1:size(temp,1)
                        for bb=1:size(temp,2)
                            mat = temp{b,bb};
                            if ~isempty(mat)
                                [~,idu] = sort([mat(:,4)]');
                                output.regions.(atlases{aa}).count.alphabet{b,bb} = mat(idu,:);
                            else
                                output.regions.(atlases{aa}).count.alphabet{b,bb} = mat;
                            end
                        end
                    end
                    
                end
                
                for fi=1:size(output.regions.(atlases{aa}).count.alphabet,1)
                    temp_add = [];
                    temp_print = output.regions.(atlases{aa}).count_names([1,2,4,5,3]);
                    for f=1:size(temp_print,2)
                        %         try output.regions.(atlases{aa}).count_names{f} = num2str(temp_print{f});
                        %         catch ME
                        %         end
                        temp_add = [temp_add, sprintf('\t'), temp_print{f}];
                    end
                    %                     fields = fieldnames(output.regions.(atlases{aa}).count);
                    %                     for ff=1:size(fields,1)
                    for fff = 1:size(output.regions.(atlases{aa}).count.alphabet,2)
                        FID = fopen([a_folder, '/brain_regions_', atlases{aa}, '_', output.regions.(atlases{aa}).count_names_sorting{fff}, '_alphabet.txt'], 'a');
                        fprintf(FID, ['\n \n \n', sprintf('\t'), 'Latent Variable ' num2str(fi), '\n \n', temp_add, '\n']);
                        fclose(FID);
                    end
                    %                     end
                    
                    %                     fields = fieldnames(output.regions.(atlases{aa}).count);
                    %                     for ff=1:size(fields,1)
                    for fff = 1:size(output.regions.(atlases{aa}).count.alphabet,2)
                        if ~isempty(output.regions.(atlases{aa}).count.alphabet{fi,fff})
                            for r=1:size(output.regions.(atlases{aa}).count.alphabet{fi,fff},1)
                                temp_print = output.regions.(atlases{aa}).count.alphabet{fi,fff}(r,[1,2,4,5,3]);
                                temp_add = [];
                                for ii=1:size(temp_print,2)
                                    try temp_print{ii} = num2str(temp_print{ii});
                                    catch ME
                                    end
                                    temp_add = [temp_add, sprintf('\t'), temp_print{ii}];
                                end
                                
                                FID = fopen([a_folder, '/brain_regions_', atlases{aa}, '_', output.regions.(atlases{aa}).count_names_sorting{fff}, '_alphabet.txt'], 'a');
                                fprintf(FID, ['\n' temp_add]);
                                fclose(FID);
                            end
                        end
                    end
                    %                     end
                end
                
            end
            
        end
        %         just deactivated for paper preparation, needs to be reactivated
        %         later on
        try temp = load(input.X);
        catch temp = load(input.MRI);
        end
        field = fieldnames(temp);
        MRI_volumes = temp.(field{1});
        for i=1:size(output.final_parameters,1)
            u=output.final_parameters{i,4};
            %             log_u{1,1} = u>0;
            %             log_u{1,2} = u<0;
            for ii=1:(size(output.regions.(atlases{aa}).log,2)-1)
                log_weights = output.regions.(atlases{aa}).log{i,ii}';
                for iii=1:size(output.regions.(atlases{aa}).count.alphabet{i,ii},1)
                    log_region = atlas_for_analysis == output.regions.(atlases{aa}).count.alphabet{i,ii}{iii,1};
                    output.volumes.(atlases{aa}).names{i,ii} = [output.regions.(atlases{aa}).count.alphabet{i,ii}(:,1)';output.regions.(atlases{aa}).count.alphabet{i,ii}(:,4)'];
                    %                     for s=1:size(MRI_volumes,1)
                    output.volumes.(atlases{aa}).raw{i,ii}(:,iii) = sum(MRI_volumes(:,log_region&log_weights),2);
                    %                     output.volumes.(['LV_', num2str(i)]).weighted.(atlases{aa}){s,iii} = sum(output.final_parameters{i,4}(log_region&log_weights)' .* MRI_volumes(s,log_region&log_weights));
                    %                     end
                end
            end
            %             MRI_volumes = MRI_volumes - (MRI_volumes*u)*u';
        end
        
    end
    
    temp=[]; names_atlases = {'yeo', 'buckner'};
    for an=1:size(names_atlases,2)
        for rr=1:size(output.regions.(names_atlases{an}).count.raw,1)
            for cc=1:size(output.regions.(names_atlases{an}).count.raw,2)
                temp = output.regions.(names_atlases{an}).count.raw{rr,cc};
                temp_new={};
                for nr=1:size(labels_regions,1)
                    try log_find = ismember([temp{:,1}], labels_regions{nr,1});
                        if any(log_find)
                            temp_new(nr,:) = [temp(log_find,:), temp{log_find,2}/output.regions_max.(names_atlases{an})(nr, 2)];
                        else
                            temp_new(nr,:) = [nr, 0, 0, labels_regions((nr),2), 0, 0, 0];
                            
                        end
                    catch
                        temp_new(nr,:) = [nr, 0, 0, labels_regions((nr),2), 0, 0, 0];
                    end
                end
                output.regions.yeo_buckner_collection.(names_atlases{an}){rr,cc} = temp_new;
            end
        end
    end
    
    save(IN.results_path, 'input', 'output', 'setup');

    % get max voxels for cerebrum and cerebellum to compose yeo and buckner
    % into one matrix

    temp = output.regions.yeo_buckner_collection.yeo; new_collection={};
    for cc=1:size(temp,2)
        type_voxels = output.regions.(atlases{aa}).count_names_sorting{cc};
        name_txt = [a_folder, '/yeo_buckner_combined_', type_voxels, '.txt'];
        for rr=1:size(temp,1)
            writecell({['LV_', num2str(rr)]}, name_txt, 'WriteMode', 'append');
            temp_sum = [temp{rr,cc}{:,2}] + [output.regions.yeo_buckner_collection.buckner{rr,cc}{:,2}];
            temp_ratio = temp_sum'./(output.regions_max.yeo(:,2)+output.regions_max.buckner(:,2));
            new_temp_sum = sum(temp_sum(9:10));
            new_temp_ratio = new_temp_sum/sum((sum(output.regions_max.yeo(9:10,2), sum(output.regions_max.buckner(9:10,2)))));
            temp_sum(9) = new_temp_sum;
            temp_ratio(9) = new_temp_ratio;
            temp_sum(10) = [];
            temp_ratio(10) = [];
            temp_array = [[temp{rr,cc}{[1:9, 11:17],1}]', temp_sum', temp_ratio];
            output.regions.yeo_buckner_collection.combined{rr,cc} = array2table(temp_array, 'VariableNames', {'region_number', 'voxels_sum', 'voxels_ratio'}, 'RowNames', labels_regions([1:9, 11:17],2));
            writetable(output.regions.yeo_buckner_collection.combined{rr,cc}, name_txt, 'WriteMode', 'append', 'WriteRowNames', true);
        end
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
    for i=1:size(input.selected_features,2)
        log_f(i,:) = contains(input.Y_names(1,:), input.selected_features{1,i});
    end
    
    if size(log_f,1)>1
        log_c = sum(sum(log_f==0)==size(input.selected_features,2));
    else
        log_c = sum((log_f==0)==size(input.selected_features,2));
    end
    
    selected_variables = [[input.selected_features(1,:), input.Y_names((end-(log_c-1)):end)]; [input.selected_features(2,:), num2cell(ones(1,log_c))]];
    % compute measures for effects
    
    log_coll = {index_epsilon, index_epsilon_all; index_omega, index_omega_all; 'test', 'all'};
    for l=1:size(log_coll,2)
        for i=1:size(output.final_parameters,1)
            x1 = output.final_parameters{i, log_coll{1, l}};
            x2 = output.final_parameters{i, log_coll{2, l}};
            output.effect_sizes.(log_coll{3,l}).Cohen(i) = computeCohen_d(x1, x2, 'independent');
            output.effect_sizes.(log_coll{3,l}).Spearman(i) = corr(x1, x2, 'Type', 'Spearman');
            output.effect_sizes.(log_coll{3,l}).Kendall(i) = corr(x1, x2, 'Type', 'Kendall');
            output.effect_sizes.(log_coll{3,l}).MI_kernel(i) = kernelmi(x1',x2');
            output.effect_sizes.(log_coll{3,l}).MI_peng(i) = mutualinfo(x1,x2);
            mdl = fitlm(array2table([x1, x2], 'VariableNames', {'epsilon', 'omega'}));
            output.effect_sizes.(log_coll{3,l}).R2(i) = mdl.Rsquared.Adjusted;
        end
    end
    
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
                        %                     temp_names = input.Y_names((temp_all+1):(temp_all+temp_current));
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
