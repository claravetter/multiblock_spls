%% DP script for retrieving sociodemographic data
function [input, output, setup] = dp_sociodemographic(IN)

load(IN.results_path)

if IN.type ==1
    if any(strfind(IN.results_path, 'CTQ'))
        overall_folder = '/volume/HCStress/Analysis/Stress';
        load('/volume/HCStress/Data/Stress_SPLS_DP/Stress_SPLS_DP/DATA/17-Jul-2018/Stress_SPLS_DP_data_table_NM_17-Jul-2018.mat')
    elseif any(strfind(IN.results_path, 'CISS'))
        overall_folder = '/volume/HCStress/Analysis/Resilience';
        load('/volume/HCStress/Data/BHAM/All_Birmingham_IncludingBrain_data_table_NM_01-Mar-2018.mat')
        load([overall_folder, '/education_years.mat']);
        for i=1:size(ID_name,1)
            try
                education_collection(i,1) = education_file{strcmp(education_file(:,1), ID_name{i}),2};
            catch
                education_collection(i,1) = NaN;
            end
        end
        data_table_names_NM = [data_table_names_NM, 'DEMOG_T0T1T2_31AA_EducationYears_T0'];
        data_table_NM = [data_table_NM, education_collection];
    else
        disp('Something''s wrong!');
    end
    
    if strfind(IN.results_path, '2018')
        analysis_date = '2018';
    elseif strfind(IN.results_path, '2019')
        analysis_date = '2019';
    elseif strfind(IN.results_path, '2020')
        analysis_date = '2020';
    end
    
    folder_name = IN.results_path(5+strfind(IN.results_path, analysis_date):(strfind(IN.results_path, '/final_results')-1));
    collection_folder = [overall_folder, '/', folder_name];
    mkdir(collection_folder);
    
    s_folder = [collection_folder, '/sociodemographic'];
    mkdir(s_folder);
    
    % IN.results_path = '/volume/HCStress/Analysis/15-Feb-2019/DP_CTQ_allgroups_649_GM_80PI_1GO_1X0_20GD_1Y0_20GD_correct10_10_XY_SC_3_Diag_1/final_results/result.mat';
    load(IN.results_path);
    load('/volume/HCStress/Data/EHI_data_for_req_PSN.mat');
    load('/volume/DP_FEF/WHOQOL_data.mat');
    
    % load('/volume/data/PRONIA/DataDump/29-Jan-2019/table_export/table_export_pruned/Self_Rating_Instruments/complete_information_T0PAT_EHI_SV_pruned.mat');
    temp_EHI_01=[];
    
    for i=1:size(input.data_collection.PSN_BOG,1)
        try
            temp_EHI_01(i,1)=cell2mat(data_table.EHI_01(strcmp(data_table.PATIENT_ID, input.data_collection.PSN_BOG{i,1})));
        catch
            disp([input.data_collection.PSN_BOG{i}, ' is missing.']);
            temp_EHI_01(i,1)=NaN;
        end
    end
    
    % temp_EHI_01{cellfun(@isempty, temp_EHI_01)}=NaN;
    EHI_01 = temp_EHI_01;
    EHI_01(EHI_01>0)=1;
    EHI_01(EHI_01<0)=0;
    data_table_names_NM = [data_table_names_NM, 'EHI_01'];
    addpath(genpath('/opt/NM/NeuroMiner_Release/'));
    addpath('/opt/SPM/spm12_v6685_cat12_r1207/');
    addpath(genpath('/volume/DP_FEF/ScrFun/ScriptsRepository'));
    
    load(input.MRI);
    
    % % remove subjects that have less than 80% complete GAF, impute the rest
    % input.data_collection.PSN_BOG(GAF_T1_toremove&pat_study,:)=[];
    
    %% get data
    % for i=1:size(input.data_collection.PSN_BOG,1)
    data_table_selected=[];
    for i=1:size(input.data_collection.PSN_BOG,1)
        try
            data_table_selected(i,:)=data_table_NM(strcmp(ID_name, input.data_collection.PSN_BOG{i}),:);
        catch
            disp([input.data_collection.PSN_BOG{i}, ' is missing.']);
            data_table_selected(i,:)=NaN;
        end
    end
    
    % data_table_selected = data_table_NM(ismember(ID_name, input.data_collection.PSN_BOG(i,1)), :);
    % end
    
    data_table_selected = [data_table_selected, EHI_01];
    
    SD_variables1 = {'AGE_T0', 'SEX_T0', 'DEMOG_T0T1T2_31AA_EducationYears_T0', 'GF_S_1_Current_T0',...
        'GF_R_1_Current_T0', 'GAF_S_PastMonth_Screening', 'GAF_DI_PastMonth_Screening', 'EHI_01'};
    
    % SD_variables_screen ={'GAF', 'GF'};
    % temp=[];
    % for i=1:size(SD_variables_screen,2)
    %     temp1 = data_table_names_NM(~cellfun(@isempty,(strfind(data_table_names_NM, SD_variables_screen{i}))));
    %     temp=[temp,temp1];
    % end
    
    % data_table_names_NM_selected = cellfun(data_table_names_NM, @(x) ~isempty(strfind(x, SD_variables{1})));
    
    temp=[];
    for i=1:size(SD_variables1,2)
        temp1 = data_table_names_NM(strcmp(data_table_names_NM, SD_variables1{i}));
        temp=[temp,temp1];
    end
    
    names_selected1 = temp;
    
    names_selected1 = temp;
    for i=1:size(names_selected1,2)
        data_table_selected_1(:,i) = data_table_selected(:,strcmp(data_table_names_NM, names_selected1(i)));
    end
    
    SD_variables2 = {'PANSS', 'BDI'};
    temp=[];
    
    for i=1:size(SD_variables2,2)
        temp1 = data_table_names_NM(~cellfun(@isempty,(strfind(data_table_names_NM, SD_variables2{i}))));
        temp=[temp,temp1];
    end
    
    names_selected2 = temp(~cellfun(@isempty,(strfind(temp, 'T0'))));
    
    names_selected_all = [names_selected1, names_selected2];
    
    % compute PANSS scores
    panss_total_names = {'PANSS_P1_T0','PANSS_P2_T0','PANSS_P3_T0','PANSS_P4_T0','PANSS_P5_T0',...
        'PANSS_P6_T0','PANSS_P7_T0','PANSS_N1_T0','PANSS_N2_T0','PANSS_N3_T0',...
        'PANSS_N4_T0','PANSS_N5_T0','PANSS_N6_T0','PANSS_N7_T0','PANSS_G01_T0',...
        'PANSS_G02_T0','PANSS_G03_T0','PANSS_G04_T0','PANSS_G05_T0','PANSS_G06_T0',...
        'PANSS_G07_T0','PANSS_G08_T0','PANSS_G09_T0','PANSS_G10_T0','PANSS_G11_T0',...
        'PANSS_G12_T0','PANSS_G13_T0','PANSS_G14_T0','PANSS_G15_T0','PANSS_G16_T0'};
    panss_pos_names = {'PANSS_P1_T0','PANSS_P2_T0','PANSS_P3_T0','PANSS_P4_T0','PANSS_P5_T0',...
        'PANSS_P6_T0','PANSS_P7_T0'};
    panss_neg_names = {'PANSS_N1_T0','PANSS_N2_T0','PANSS_N3_T0',...
        'PANSS_N4_T0','PANSS_N5_T0','PANSS_N6_T0','PANSS_N7_T0'};
    panss_gen_names = {'PANSS_G01_T0',...
        'PANSS_G02_T0','PANSS_G03_T0','PANSS_G04_T0','PANSS_G05_T0','PANSS_G06_T0',...
        'PANSS_G07_T0','PANSS_G08_T0','PANSS_G09_T0','PANSS_G10_T0','PANSS_G11_T0',...
        'PANSS_G12_T0','PANSS_G13_T0','PANSS_G14_T0','PANSS_G15_T0','PANSS_G16_T0'};
    
    % panss positive scores
    panss_total = sum(data_table_selected(:,ismember(data_table_names_NM, panss_total_names)),2);
    panss_pos = sum(data_table_selected(:,ismember(data_table_names_NM, panss_pos_names)),2);
    panss_neg = sum(data_table_selected(:,ismember(data_table_names_NM, panss_neg_names)),2);
    panss_gen = sum(data_table_selected(:,ismember(data_table_names_NM, panss_gen_names)),2);
    
    % compute BDI scores
    BDI_names = {'BDI2_01_T0','BDI2_02_T0','BDI2_03_T0','BDI2_04_T0','BDI2_05_T0',...
        'BDI2_06_T0','BDI2_07_T0','BDI2_08_T0','BDI2_09_T0','BDI2_10_T0','BDI2_11_T0',...
        'BDI2_12_T0','BDI2_13_T0','BDI2_14_T0','BDI2_15_T0','BDI2_16_T0','BDI2_17_T0',...
        'BDI2_18_T0','BDI2_19_T0','BDI2_20_T0','BDI2_21_T0'};
    
    BDI_scores = sum(data_table_selected(:,ismember(data_table_names_NM, BDI_names)),2);
    
    data_table_selected_2 = [panss_total, panss_pos, panss_neg, panss_gen, BDI_scores];
    data_table_selected_2_names = {'PANSS_total', 'PANSS_pos', 'PANSS_neg', 'PANSS_gen', 'BDI_scores'};
    


output.data_table_study.all = [data_table_selected_1, data_table_selected_2];
output.data_table_study_names = [names_selected1, data_table_selected_2_names];

end


% if any(strfind(IN.results_path, 'CTQ'))
%     overall_folder = '/volume/HCStress/Analysis/Stress';
% elseif any(strfind(IN.results_path, 'CISS'))
%     overall_folder = '/volume/HCStress/Analysis/Resilience';
% end
% 
% if strfind(IN.results_path, '2018')
%     analysis_date = '2018';
% elseif strfind(IN.results_path, '2019')
%     analysis_date = '2019';
% else
%     disp('Something is wrong with the file');
% end
% 
% folder_name = IN.results_path(5+strfind(IN.results_path, analysis_date):(strfind(IN.results_path, '/final_results')-1));
% collection_folder = [overall_folder, '/', folder_name];
% mkdir(collection_folder);
% 
% s_folder = [collection_folder, '/sociodemographic'];
% mkdir(s_folder);

% sort according to diagnostic groups
studygroups = {'HC', 'ROD', 'CHR', 'ROP'};
for i=1:size(studygroups,2)
    output.data_table_study.(studygroups{i}) = output.data_table_study.all(strcmp(input.data_collection.Labels, studygroups{i}),:);
end

fields = fieldnames (output.data_table_study);
output.data_table_final=struct;
for i=1:size(fields,1)
    output.data_table_final.means(:,i) = (nanmean(output.data_table_study.(fields{i}),1))';
    output.data_table_final.std(:,i) = (nanstd(output.data_table_study.(fields{i}),1))';
end

output.data_table_final.names = output.data_table_study_names';
output.data_table_final.labels = ['all', studygroups];

% get site numbers for diagnoses
studygroups = {'HC', 'ROD', 'CHR', 'ROP'};
for i=1:size(input.sites,1)
    sites_all(i,1) = find(input.sites(i,:)==1);
end

sex_all = output.data_table_study.all(:,strcmp(output.data_table_study_names, 'SEX_T0'));
hand_all = output.data_table_study.all(:,strcmp(output.data_table_study_names, 'EHI_01'));

for i=1:size(studygroups,2)
    output.sites_numbers.(studygroups{i}).raw = sites_all(strcmp(input.data_collection.Labels, studygroups{i}),:);
    output.sites_numbers.(studygroups{i}).raw = output.sites_numbers.(studygroups{i}).raw(~isnan(output.sites_numbers.(studygroups{i}).raw));
    output.sex_numbers.(studygroups{i}).raw = sex_all(strcmp(input.data_collection.Labels, studygroups{i}),:);
    output.sex_numbers.(studygroups{i}).raw = output.sex_numbers.(studygroups{i}).raw(~isnan(output.sex_numbers.(studygroups{i}).raw));
    output.hand_numbers.(studygroups{i}).raw = hand_all(strcmp(input.data_collection.Labels, studygroups{i}),:);
    output.hand_numbers.(studygroups{i}).raw = output.hand_numbers.(studygroups{i}).raw(~isnan(output.hand_numbers.(studygroups{i}).raw));
end

for i=1:size(studygroups,2)
    [C, ~, ic] = unique(output.sites_numbers.(studygroups{i}).raw);
    counts = accumarray(ic, 1);
    output.sites_numbers.(studygroups{i}).count = [C, counts];
%     dp_txt_write(s_folder, ['sites_', studygroups{i}], output.sites_numbers.(studygroups{i}).count', '%d \t %d \n');
    [C, ~, ic] = unique(output.sex_numbers.(studygroups{i}).raw);
    counts = accumarray(ic, 1);
    output.sex_numbers.(studygroups{i}).count = [C, counts];
    [C, ~, ic] = unique(output.hand_numbers.(studygroups{i}).raw);
    counts = accumarray(ic, 1);
    output.hand_numbers.(studygroups{i}).count = [C, counts];
end

% test for imbalances
output.collected_sites = [output.sites_numbers.HC.count(:,2), output.sites_numbers.ROD.count(:,2),output.sites_numbers.CHR.count(:,2),output.sites_numbers.ROP.count(:,2)];
output.collected_sex = [output.sex_numbers.HC.count(:,2), output.sex_numbers.ROD.count(:,2),output.sex_numbers.CHR.count(:,2),output.sex_numbers.ROP.count(:,2)];
output.collected_hand = [output.hand_numbers.HC.count(:,2), output.hand_numbers.ROD.count(:,2),output.hand_numbers.CHR.count(:,2),output.hand_numbers.ROP.count(:,2)];

% [output.sites.collection, output.sites.collection_names, output.sites.h, output.sites.p, output.sites.stats]=dp_chi2(output.collected_sites, 'absolute');
% [output.sex.collection, output.sex.collection_names, output.sex.h, output.sex.p, output.sex.stats]=dp_chi2(output.collected_sex, 'absolute');
% [output.hand.collection, output.hand.collection_names, output.hand.h, output.hand.p, output.hand.stats]=dp_chi2(sum(output.collected_hand,2), 'absolute');

nn=1;
for i=1:size(output.data_table_final.means,1)
    output.data_table_final.complete(nn,:) = output.data_table_final.means(i,:);
    nn=nn+1;
    output.data_table_final.complete(nn,:) = output.data_table_final.std(i,:);
    nn=nn+1;
end

% dp_txt_write('/volume/HCStress/Analysis/Stress/', 'CTQ_BS_627_means', output.data_table_final.means', '%.1f \t %.1f \t %.1f \t %.1f \t %.1f \t  \n');
% dp_txt_write('/volume/HCStress/Analysis/Stress/', 'CTQ_BS_627_std', output.data_table_final.std', '%.1f \t %.1f \t %.1f \t %.1f \t %.1f \t  \n');
% % dp_txt_write('/volume/HCStress/Analysis/Stress/', 'CTQ_BS_627_names', output.data_table_final.names, '%s \n');
% % dp_txt_write('/volume/HCStress/Analysis/Stress/', 'CTQ_BS_627_labels', output.data_table_final.labels, '%s');
dp_txt_write(s_folder, 'complete_table', output.data_table_final.complete', '%.2f \t %.2f \t %.2f \t %.2f \t %.2f \t  \n');

%% test for significant differences between groups
% first test with ANOVA for overall differences
groups_to_choose = {'ROD', 'CHR', 'ROP'};
log_groups_to_choose =  ismember(input.data_collection.Labels, groups_to_choose);

for i=1:size(output.data_table_study_names,2)
    [p,tbl,stats] = kruskalwallis(output.data_table_study.all(log_groups_to_choose,i), input.data_collection.Labels(log_groups_to_choose));
    output.KW_results(i,:) = [p,tbl{2,5}];
    [output.Dunn_results{i,1}, ~, ~, output.Dunn_results_labels] = multcompare(stats, 'CType', 'dunn-sidak', 'Estimate', 'kruskalwallis');
    close all
    dp_txt_write(s_folder, ['Dunn_results_', num2str(i)], output.Dunn_results{i,1}', '%.3f \t %.3f \t %.3f \t %.3f \t %.3f \t %.3f \t \n \n');
end

dp_txt_write(s_folder, 'KW_results', output.KW_results', '%.3f \t %.3f \n \n');
dp_txt_write(s_folder, 'Dunn_labels', '', '%s');
for i=1:size(output.Dunn_results_labels,1)
    FID = fopen([s_folder, '/Dunn_labels.txt'], 'a');
    fprintf(FID, '%s\n', output.Dunn_results_labels{i,1});
    fclose(FID);
    dp_txt_write(s_folder, 'Dunn_labels', [output.Dunn_results_labels{:}], '%s \n');
end

% compute differences regarding latent scores between study groups
groups_to_choose = {'HC', 'ROD', 'CHR', 'ROP'};
log_groups_to_choose =  ismember(input.data_collection.Labels, groups_to_choose);
output.latent_scores = struct;

for i=1:size(output.final_parameters,1)
%     output.latent_scores.IN.(['LV_', num2str(i)]).epsilon=struct;
%     output.latent_scores.data.(['LV_', num2str(i)])(:,1)=dp_standardize(output.final_parameters{i,strcmp(output.parameters_names, 'epsilon')}, output.latent_scores.IN.(['LV_', num2str(i)]).epsilon);
    output.latent_scores.data.(['LV_', num2str(i)])(:,1)=output.final_parameters{i,strcmp(output.parameters_names, 'epsilon')};
    output.latent_scores.median_iqr.groups = {'HC', 'ROD', 'CHR', 'ROP'};
    for ii=1:size(output.latent_scores.median_iqr.groups,2)
        output.latent_scores.median_iqr.(['LV_', num2str(i)]).epsilon(1,ii) = median(output.latent_scores.data.(['LV_', num2str(i)])(strcmp(input.data_collection.Labels(output.CV.cv_outer_indices.TestInd{1,output.final_parameters{i,1}},1), output.latent_scores.median_iqr.groups{ii}),1));
        output.latent_scores.median_iqr.(['LV_', num2str(i)]).epsilon(2,ii) = iqr(output.latent_scores.data.(['LV_', num2str(i)])(strcmp(input.data_collection.Labels(output.CV.cv_outer_indices.TestInd{1,output.final_parameters{i,1}},1), output.latent_scores.median_iqr.groups{ii}),1));
    end
%     output.latent_scores.IN.(['LV_', num2str(i)]).omega=struct;
%     output.latent_scores.data.(['LV_', num2str(i)])(:,2)=dp_standardize(output.final_parameters{i,strcmp(output.parameters_names, 'omega')}, output.latent_scores.IN.(['LV_', num2str(i)]).omega);
    output.latent_scores.data.(['LV_', num2str(i)])(:,2)=output.final_parameters{i,strcmp(output.parameters_names, 'omega')};
    for ii=1:size(output.latent_scores.median_iqr.groups,2)
        output.latent_scores.median_iqr.(['LV_', num2str(i)]).omega(1,ii) = median(output.latent_scores.data.(['LV_', num2str(i)])(strcmp(input.data_collection.Labels(output.CV.cv_outer_indices.TestInd{1,output.final_parameters{i,1}},1), output.latent_scores.median_iqr.groups{ii}),2));
        output.latent_scores.median_iqr.(['LV_', num2str(i)]).omega(2,ii) = iqr(output.latent_scores.data.(['LV_', num2str(i)])(strcmp(input.data_collection.Labels(output.CV.cv_outer_indices.TestInd{1,output.final_parameters{i,1}},1), output.latent_scores.median_iqr.groups{ii}),2));
    end
end

% output.latent_scores = struct;
for i=1:size(output.final_parameters,1)
%     output.latent_scores.data.(['LV_', num2str(i)])(:,1)=output.final_parameters{i,strcmp(output.parameters_names, 'epsilon')};
%     output.latent_scores.data.(['LV_', num2str(i)])(:,2)=output.final_parameters{i,strcmp(output.parameters_names, 'omega')};
    for ii=1:size(output.latent_scores.data.(['LV_', num2str(i)]),2)
        [p,tbl,stats] = kruskalwallis(output.latent_scores.data.(['LV_', num2str(i)])(:,ii), input.data_collection.Labels(output.CV.cv_outer_indices.TestInd{1,output.final_parameters{i,1}},1));
        output.latent_scores.KW_results.(['LV_', num2str(i)])(:,ii) = [p,tbl{2,5}]';
        [output.latent_scores.Dunn_results.(['LV_', num2str(i)]){1,ii}, ~, ~, output.latent_scores.Dunn_results_labels] = multcompare(stats, 'CType', 'dunn-sidak', 'Estimate', 'kruskalwallis');
        close all
    end
%     dp_txt_write(s_folder, ['Dunn_results_', num2str(i)], output.latent_scores_differences.Dunn_results{i,1}', '%.3f \t %.3f \t %.3f \t %.3f \t %.3f \t %.3f \t \n \n');
end

% if ANOVA shows significant differences, then do binary tests for all
% groups, while correcting for multiple testing

%% look for specific correlations
corr_data=[];
corr_variables=[];

if any(strcmp(IN.SD_selection, 'GAF'))
    
    % get GAF data
    
    GAF_names = {'GAF_S_LifeTime_Screening','GAF_S_PastYearT0_Screening',...
        'GAF_S_PastMonth_Screening','GAF_DI_LifeTime_Screening','GAF_DI_PastYear_Screening',...
        'GAF_DI_PastMonth_Screening', 'GF_S_1_Current_T0','GF_S_2_LowPastYearT0_T0',...
        'GF_S_3_HighPastYearT0_T0','GF_S_4_HighLifetimeT0_T0', 'GF_R_1_Current_T0',...
        'GF_R_2_LowPastYearT0_T0','GF_R_3_HighPastYearT0_T0','GF_R_4_HighLifetimeT0_T0'};
    
    
    % GAF_variables = data_table_names_NM(~cellfun(@isempty,(strfind(data_table_names_NM, 'GAF'))));
    GAF_data=[];
    for i=1:size(GAF_names,2)
        temp = data_table_selected(:, strcmp(data_table_names_NM, GAF_names{i}));
        GAF_data=[GAF_data,temp];
    end
    corr_variables = [corr_variables, GAF_names];
    corr_data = [corr_data, GAF_data];
end

if any(strcmp(IN.SD_selection, 'BDI'))
    
    % get BDI data
    corr_variables = [corr_variables, 'BDI'];
    corr_data = [corr_data, BDI_scores];
    
end

if any(strcmp(IN.SD_selection, 'NEO'))
    
    % get NEO FFI data
    neo_log = ~cellfun(@isempty,(strfind(data_table_names_NM, 'NEO')));
    neo_names = data_table_names_NM(neo_log);
    neo_data_temp = data_table_selected(:, neo_log);
    
    NEO.neuroticism.negative_affect = [1, 11, 16, 31, 46];
    NEO.neuroticism.self_reproach = [6, 21, 26, 36, 41, 51, 56];
    NEO.extraversion.positive_affect = [7, 12, 37, 42];
    NEO.extraversion.sociability = [2, 17, 27, 57];
    NEO.extraversion.activity = [22, 32, 47, 52];
    NEO.openness.aesthetic_interests = [13, 23, 43];
    NEO.openness.intellectual_interests = [33, 48, 53, 58];
    NEO.openness.unconventionality = [3,8,18, 28, 38];
    NEO.agreeableness.nonantagonistic_orientation = [9,14,19,24,29,44,54,59];
    NEO.agreeableness.prosocial_orientation = [4, 34, 39, 49];
    NEO.conscientiousness.orderliness = [5,10,15,30,55];
    NEO.conscientiousness.goal_striving = [25, 35, 60];
    NEO.conscientiousness.dependability = [20, 40, 45, 50];
    NEO_inverse_questions = [1,16,31,46,12,42,27,57,23,33,48,3,8,18,38,9,14,24,29,44,54,59,39,15,30,55,45];
    NEO_inverse_algorithm = [1:1:5;5:-1:1];
    NEO_poorly_functioning = [6,12,27,42,3,8,28,38,9,19,24,29,34,15];
    
    %exclude poorly functioning questions (optional)
    fields=fieldnames(NEO);
    for i=1:size(fields,1)
        fields1=fieldnames(NEO.(fields{i}));
        for ii=1:size(fields1,1)
            temp = NEO.(fields{i}).(fields1{ii});
            log = ismember(temp, NEO_poorly_functioning);
            temp(log)=[];
            NEO.(fields{i}).(fields1{ii})=temp;
        end
    end
    
    neo_data_inv=[];
    for i=1:size(neo_data_temp,1)
        temp_row = neo_data_temp(i,:);
        for ii=1:size(temp_row,2)
            if sum(ii==NEO_inverse_questions)>0
                try
                    temp_row(ii)=NEO_inverse_algorithm(2,temp_row(ii)==NEO_inverse_algorithm(1,:));
                catch
                    temp_row(ii)=NaN;
                end
            end
        end
        neo_data_inv(i,:)=temp_row;
    end
    
    % impute data with at least 80% complete questions
    log_impute = sum(isnan(neo_data_inv),2)<=(0.2*size(neo_data_inv,2));
    temp_neo_data_inv = dp_impute(neo_data_inv(log_impute,:), 'euclidean');
    neo_data_inv_imp = neo_data_inv;
    neo_data_inv_imp(log_impute,:) = temp_neo_data_inv;
    
    % compute sum scores
    fields = fieldnames(NEO);
    for i=1:size(fields,1)
        fields1 = fieldnames(NEO.(fields{i}));
        temp_collect.(fields{i})=[];
        for ii=1:size(fields1,1)
            temp_collect.(fields{i}) = [temp_collect.(fields{i}), NEO.(fields{i}).(fields1{ii})];
        end
    end
    
    fields=fieldnames(temp_collect);
    neo_names=[];
    for i=1:size(fields,1)
        neo_names=[neo_names, {['NEO_', fields{i}]}];
    end
    
    neo_data_foranalysis = neo_data_inv_imp;
    for i=1:size(neo_data_foranalysis,1)
        temp_sum=[];
        for ii=1:size(fields,1)
            temp_sum = [temp_sum, sum(neo_data_foranalysis(i,temp_collect.(fields{ii})))];
        end
        neo_data_sum(i,:)=temp_sum;
    end
    
    % add NEO FFI data to corr data
    
    corr_variables = [corr_variables,neo_names];
    corr_data = [corr_data,neo_data_sum];
end


% get QUOL data
if any(strcmp(IN.SD_selection, 'QOL'))
    WHO_table = [HC;PAT];
    
    for i=1:26
        if i<10
            try
                WHO_raw(:,i) = WHO_table.(['WHOQOL_0', num2str(i)]);
            catch
                temp = WHO_table.(['WHOQOL_0', num2str(i)]);
                temp(cellfun(@isempty, temp))={'NaN'};
                WHO_raw(:,i)=temp;
            end
        else
            try
                WHO_raw(:,i) = WHO_table.(['WHOQOL_', num2str(i)]);
            catch
                temp = WHO_table.(['WHOQOL_', num2str(i)]);
                temp(cellfun(@isempty, temp))={'NaN'};
                WHO_raw(:,i)=temp;
            end
        end
    end
    
    for i=1:size(WHO_raw,1)
        temp=WHO_raw(i,:);
        log=cellfun(@isempty, temp);
        temp(log)={num2str(NaN)};
        WHO_raw(i,:)=temp;
        for ii=1:size(temp,2)
            try
                temp_new(i,ii)=str2num(temp{ii});
            catch
                temp_new(i,ii)=NaN;
            end
        end
    end
    
    WHO_new=temp_new;
    
    WHO_ID=[HC.PATIENT_ID; PAT.PATIENT_ID];
    
    WHO_selected=[];
    for i=1:size(input.data_collection.PSN_BOG,1)
        try
            WHO_selected(i,:)=WHO_new(strcmp(WHO_ID, input.data_collection.PSN_BOG{i}),:);
        catch
            disp([input.data_collection.PSN_BOG{i}, ' is missing.']);
            WHO_selected(i,:)=NaN;
        end
    end
    
    % impute data with at least 80% complete questions
    log_impute = sum(isnan(WHO_selected),2)<=(0.2*size(WHO_selected,2));
    temp_WHO_selected = dp_impute(WHO_selected(log_impute,:), 'euclidean');
    WHO_selected_imp = WHO_selected;
    WHO_selected_imp(log_impute,:) = temp_WHO_selected;
    
    % compute sum scores
    WHO.physical = [3, 4, 10, 15, 16, 17, 18];
    WHO.psychosocial = [5, 6, 7, 11, 19, 26];
    WHO.social_relationship = [20, 21, 22];
    WHO.environment = [8, 9, 12, 13, 14, 23, 24, 25];
    
    WHO_sum = [sum(WHO_selected_imp(:,WHO.physical),2), sum(WHO_selected_imp(:,WHO.psychosocial),2),...
        sum(WHO_selected_imp(:,WHO.social_relationship),2), sum(WHO_selected_imp(:,WHO.environment),2)];
    
    WHO_names = {'WHO_physical', 'WHO_psychosocial', 'WHO_social_relationship', 'WHO_environment'};
    
    corr_variables = [corr_variables,WHO_names];
    corr_data = [corr_data,WHO_sum];
end



%% other stuff
% corr_tokeep = sum(isnan(corr_data),2)==0;
%
% MRI_for_analysis_pruned = MRI_for_analysis(corr_tokeep,:);
% fields=fieldnames(input.data_collection);
% for i=1:size(fields,1)
%     temp = input.data_collection.(fields{i});
%     try
%         input.data_collection.(fields{i}) = temp(corr_tokeep,:);
%     end
% end
% behavior_pruned = input.behavior(corr_tokeep,:);
% corr_data_pruned = corr_data(corr_tokeep,:);
%
% % correct MRI data, behavioral data and correlation data for site effects
% % just like in main analysis
% IN = struct;
% IN.TrCovars = input.sites(corr_tokeep,:);
% [MRI_for_analysis_c, ~] = nk_PartialCorrelationsObj(MRI_for_analysis_pruned, IN);
% IN = struct;
% IN.TrCovars = input.sites(corr_tokeep,:);
% [behavior_c, ~] = nk_PartialCorrelationsObj(behavior_pruned, IN);
% IN = struct;
% % IN.TrCovars = input.sites(corr_tokeep,:);
% % [corr_data, ~] = nk_PartialCorrelationsObj(corr_data, IN);
%
% studygroups = {'all' 'HC', 'ROD', 'CHR', 'ROP'};
% for i=1:size(studygroups,2)
%     switch studygroups{i}
%         case 'all'
%             corr_table_study.(studygroups{i}) = corr_data_pruned;
%             X_table_study.(studygroups{i}) = MRI_for_analysis_c;
%             Y_table_study.(studygroups{i}) = behavior_c;
%         otherwise
%             corr_table_study.(studygroups{i}) = corr_data_pruned(strcmp(input.data_collection.Labels, studygroups{i}),:);
%             X_table_study.(studygroups{i}) = MRI_for_analysis_c(strcmp(input.data_collection.Labels, studygroups{i}),:);
%             Y_table_study.(studygroups{i}) = behavior_c(strcmp(input.data_collection.Labels, studygroups{i}),:);
%     end
% %     corr_table_study.(studygroups{i}) = dp_standardize GAF_data;
% %     X_table_study.(studygroups{i}) = MRI_for_analysis;
% %     Y_table_study.(studygroups{i}) = input.behavior;
% end

% % assign X and Y matrices for all subjects
% corr_type = 'Spearman';
% for i=1:size(output.final_parameters,1)
% %     IN=struct;
% %     IN.method = 'mean-centering';
% %     X = dp_standardize(X, IN);
% %     IN=struct;
% %     IN.method = 'mean-centering';
% %     Y = dp_standardize(Y, IN);
%     u = output.final_parameters{i,4};
%     v = output.final_parameters{i,5};
% %     epsilon = X*u;
% %     omega = Y*v;
%     for ii=1:size(studygroups,2)
%         corr_scores = corr_table_study.(studygroups{ii});
%         for iii=1:size(corr_variables,2)
%             log_corr = strcmp(corr_variables, corr_variables{iii});
%             epsilon = X_table_study.(studygroups{ii})*u;
%             omega = Y_table_study.(studygroups{ii})*v;
%             corr_data = corr_scores(:,log_corr);
%             [RHO, p] = corr(corr_data, epsilon, 'Type', corr_type, 'Rows', 'complete');
%             corr_data_correlation.(['LV_', num2str(i)]).epsilon.(studygroups{ii}).(corr_variables{iii}) = [RHO,p];
%             [RHO, p] = corr(corr_data, omega, 'Type', corr_type, 'Rows', 'complete');
%             corr_data_correlation.(['LV_', num2str(i)]).omega.(studygroups{ii}).(corr_variables{iii}) = [RHO, p];
%         end
%         [X_table_study.(studygroups{ii}),Y_table_study.(studygroups{ii})] = proj_def(X_table_study.(studygroups{ii}), Y_table_study.(studygroups{ii}), u, v);
%     end
% end

% test correlations for entire sample and for held out subjects: GAF, GF,
% BDI


%% proper correlations
log_epsilon_opt = strcmp(output.parameters_names, 'epsilon');
log_omega_opt = strcmp(output.parameters_names, 'omega');
log_epsilon_all = strcmp(output.parameters_names, 'epsilon_all');
log_omega_all = strcmp(output.parameters_names, 'omega_all');
output.hold_out_corr_data = [];
output.hold_out_correlations = [];
output.all_corr_data = [];
output.all_correlations=[];
output.hold_out_sig=[];
output.hold_out_FDR_values=[];
output.all_sig=[];
output.all_FDR_values=[];

for i=1:size(output.final_parameters,1)
    
    hold_out_corr_data = corr_data(output.CV.cv_outer_indices.TestInd{1,output.final_parameters{i,1}},:);
    hold_out_epsilon_opt = output.final_parameters{i,log_epsilon_opt};
    hold_out_omega_opt = output.final_parameters{i,log_omega_opt};
    hold_out_labels = input.data_collection.Labels(output.CV.cv_outer_indices.TestInd{1,output.final_parameters{i,1}},1);
    
    all_corr_data = corr_data;
    all_epsilon = output.final_parameters{i,log_epsilon_all};
    all_omega = output.final_parameters{i,log_omega_all};
    labels_all = input.data_collection.Labels;
    
    for ii=1:size(corr_variables,2)
        
        % hold out data
        hold_out_corr_temp = hold_out_corr_data(:,ii);
        output.hold_out_corr_data.(['LV_', num2str(i)]).(corr_variables{ii}) = hold_out_corr_temp;
        [RHO,p] = corr(hold_out_epsilon_opt, hold_out_corr_temp, 'Rows', 'complete', 'Type', 'Spearman');
        output.hold_out_correlations.(['LV_', num2str(i)]).epsilon.all.(corr_variables{ii}) = [RHO, p];
        [RHO,p] = corr(hold_out_omega_opt, hold_out_corr_temp, 'Rows', 'complete', 'Type', 'Spearman');
        output.hold_out_correlations.(['LV_', num2str(i)]).omega.all.(corr_variables{ii}) = [RHO, p];
        
        % all data
        all_corr_temp = all_corr_data(:,ii);
        output.all_corr_data.(['LV_', num2str(i)]).(corr_variables{ii}) = all_corr_temp;
        [RHO,p] = corr(all_epsilon, all_corr_temp, 'Rows', 'complete', 'Type', 'Spearman');
        output.all_correlations.(['LV_', num2str(i)]).epsilon.all.(corr_variables{ii}) = [RHO, p];
        [RHO,p] = corr(all_omega, all_corr_temp, 'Rows', 'complete', 'Type', 'Spearman');
        output.all_correlations.(['LV_', num2str(i)]).omega.all.(corr_variables{ii}) = [RHO, p];
        
        for iii=1:size(input.selected_studygroups,2)
            %hold out data
            log_group = strcmp(hold_out_labels, input.selected_studygroups{1,iii});
            hold_out_corr_data_group = hold_out_corr_temp(log_group);
            hold_out_epsilon_group = hold_out_epsilon_opt(log_group);
            hold_out_omega_group = hold_out_omega_opt(log_group);
            [RHO,p] = corr(hold_out_epsilon_group, hold_out_corr_data_group, 'Rows', 'complete', 'Type', 'Spearman');
            output.hold_out_correlations.(['LV_', num2str(i)]).epsilon.(input.selected_studygroups{1,iii}).(corr_variables{ii}) = [RHO, p];
            [RHO,p] = corr(hold_out_omega_group, hold_out_corr_data_group, 'Rows', 'complete', 'Type', 'Spearman');
            output.hold_out_correlations.(['LV_', num2str(i)]).omega.(input.selected_studygroups{1,iii}).(corr_variables{ii}) = [RHO, p];
            
            % all data
            log_group = strcmp(labels_all, input.selected_studygroups{1,iii});
            all_corr_data_group = all_corr_temp(log_group);
            all_epsilon_group = all_epsilon(log_group);
            all_omega_group = all_omega(log_group);
            [RHO,p] = corr(all_epsilon_group, all_corr_data_group, 'Rows', 'complete', 'Type', 'Spearman');
            output.all_correlations.(['LV_', num2str(i)]).epsilon.(input.selected_studygroups{1,iii}).(corr_variables{ii}) = [RHO, p];
            [RHO,p] = corr(all_omega_group, all_corr_data_group, 'Rows', 'complete', 'Type', 'Spearman');
            output.all_correlations.(['LV_', num2str(i)]).omega.(input.selected_studygroups{1,iii}).(corr_variables{ii}) = [RHO, p];
            
        end
    end
end

% test for significance
fields1=fieldnames(output.hold_out_correlations);
for i=1:size(fields1,1)
    fields2=fieldnames(output.hold_out_correlations.(fields1{i}));
    for ii=1:size(fields2,1)
        fields3 = fieldnames(output.hold_out_correlations.(fields1{i}).(fields2{ii}));
        for iii=1:size(fields3,1)
            fields4 = fieldnames(output.hold_out_correlations.(fields1{i}).(fields2{ii}).(fields3{iii}));
            temp_p=[];
            for iiii=1:size(fields4,1)
                pp=output.hold_out_correlations.(fields1{i}).(fields2{ii}).(fields3{iii}).(fields4{iiii})(1,2);
                temp_p=[temp_p;pp];
            end
            output.hold_out_FDR_values.(fields1{i}).(fields2{ii}).(fields3{iii}) = dp_FDR(temp_p, 0.05);
            if output.hold_out_FDR_values.(fields1{i}).(fields2{ii}).(fields3{iii})>0
                fields4 = fieldnames(output.hold_out_correlations.(fields1{i}).(fields2{ii}).(fields3{iii}));
                temp_p=[];
                for iiii=1:size(fields4,1)
                    pp=output.hold_out_correlations.(fields1{i}).(fields2{ii}).(fields3{iii}).(fields4{iiii})(1,2);
                    if pp <= output.hold_out_FDR_values.(fields1{i}).(fields2{ii}).(fields3{iii})
                        output.hold_out_sig.(fields1{i}).(fields2{ii}).(fields3{iii}).(fields4{iiii}) =  output.hold_out_correlations.(fields1{i}).(fields2{ii}).(fields3{iii}).(fields4{iiii});
                    end
                end
            end
        end
        
    end
end

fields1=fieldnames(output.all_correlations);
for i=1:size(fields1,1)
    fields2=fieldnames(output.all_correlations.(fields1{i}));
    for ii=1:size(fields2,1)
        fields3 = fieldnames(output.all_correlations.(fields1{i}).(fields2{ii}));
        for iii=1:size(fields3,1)
            fields4 = fieldnames(output.all_correlations.(fields1{i}).(fields2{ii}).(fields3{iii}));
            temp_p=[];
            for iiii=1:size(fields4,1)
                pp=output.all_correlations.(fields1{i}).(fields2{ii}).(fields3{iii}).(fields4{iiii})(1,2);
                temp_p=[temp_p;pp];
            end
            output.all_FDR_values.(fields1{i}).(fields2{ii}).(fields3{iii}) = dp_FDR(temp_p, 0.05);
            if output.all_FDR_values.(fields1{i}).(fields2{ii}).(fields3{iii})>0
                fields4 = fieldnames(output.all_correlations.(fields1{i}).(fields2{ii}).(fields3{iii}));
                temp_p=[];
                for iiii=1:size(fields4,1)
                    pp=output.all_correlations.(fields1{i}).(fields2{ii}).(fields3{iii}).(fields4{iiii})(1,2);
                    if pp <= output.all_FDR_values.(fields1{i}).(fields2{ii}).(fields3{iii})
                        output.all_sig.(fields1{i}).(fields2{ii}).(fields3{iii}).(fields4{iiii}) =  output.all_correlations.(fields1{i}).(fields2{ii}).(fields3{iii}).(fields4{iiii});
                    end
                end
            end
        end
        
    end
end

save(IN.results_path, 'input', 'output', 'setup');

%% collect correlations

subsets = {'all', 'hold_out'};

for s=1:size(subsets,2)
    LVs = fieldnames(output.([subsets{s}, '_correlations']));
    for i=1:size(LVs,1)
        latent_scores = fieldnames(output.([subsets{s}, '_correlations']).(LVs{i}));
        for ii=1:size(latent_scores,1)
            groups = fieldnames(output.([subsets{s}, '_correlations']).(LVs{i}).(latent_scores{ii}));
            output.(['tables_', subsets{s}, '_Rho_p']).(LVs{i}).(latent_scores{ii})=[];
            for iii=1:size(groups,1)
                subfields = fieldnames(output.([subsets{s}, '_correlations']).(LVs{i}).(latent_scores{ii}).(groups{iii}));
                nn=1;
                for iiii=1:size(subfields,1)
                    output.(['tables_', subsets{s}, '_Rho_p']).(LVs{i}).(latent_scores{ii})(nn,iii) = output.([subsets{s}, '_correlations']).(LVs{i}).(latent_scores{ii}).(groups{iii}).(subfields{iiii})(1);
                    nn=nn+1;
                    output.(['tables_', subsets{s}, '_Rho_p']).(LVs{i}).(latent_scores{ii})(nn,iii) = output.([subsets{s}, '_correlations']).(LVs{i}).(latent_scores{ii}).(groups{iii}).(subfields{iiii})(2);
                    nn=nn+1;
                end
            end
        end
    end
    output.(['tables_', subsets{s}, '_Rho_p_labels']) = subfields;
    output.(['tables_', subsets{s}, '_Rho_p_names']) = groups';
end


save(IN.results_path, 'input', 'output', 'setup');

end




