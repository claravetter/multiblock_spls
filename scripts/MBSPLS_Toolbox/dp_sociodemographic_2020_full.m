%% DP script for retrieving sociodemographic data
function [input, output, setup]= dp_sociodemographic_2020_full(IN)

load(IN.results_path);

addpath(genpath('/opt/NM/NeuroMiner_Release/'));
addpath('/opt/SPM/spm12_v6685_cat12_r1207/');
addpath(genpath('/volume/DP_FEF/ScrFun/ScriptsRepository'));

if contains(IN.results_path, 'CTQ')
    load('/volume/HCStress/Data/Stress_SPLS_DP/Stress_SPLS_DP/DATA/17-Jul-2018/Stress_SPLS_DP_data_table_NM_17-Jul-2018.mat')
elseif contains(IN.results_path, 'CISS')
    load('/volume/data/PRONIA/DataDump/03-Apr-2020/QueryData/Munich/All_Birmingham_IncludingBrain_Disc/All_Birmingham_IncludingBrain_Disc/DATA/03-Aug-2020/All_Birmingham_IncludingBrain_Disc_data_table_NM_03-Aug-2020.mat')
elseif contains(IN.results_path, 'WSS')
    load('/volume/data/PRONIA/DataDump/03-Apr-2020/QueryData/Munich/PRONIAQueryTemplate_v3_1_MU_DP_WSS_PLS_Disc/PRONIAQueryTemplate_v3_1_MU_DP_WSS_PLS_Disc/DATA/19-May-2020/PRONIAQueryTemplate_v3_1_MU_DP_WSS_PLS_Disc_data_table_NM_19-May-2020.mat')
elseif contains(IN.results_path, 'immune', 'IgnoreCase', true)
    
else
    disp('Something''s wrong!');
end

load('/volume/HCStress/Data/EHI_data_for_req_PSN.mat');
load('/volume/DP_FEF/WHOQOL_data.mat');
EHI_01=[];

%% solve issues with w iteration for final parameters
if contains(IN.results_path, 'final_vis')
    for i=1:size(output.final_parameters,1)
        log_find = [output.opt_parameters.(['LV_', num2str(i)]){:,7}] == output.final_parameters{i,7};
        if sum(log_find)>1
            temp = output.opt_parameters.(['LV_', num2str(i)])(log_find,:);
            log_find = [temp{:,2}] == output.final_parameters{i,2};
            output.final_parameters{i,1} = temp{log_find,1};
        else
            output.final_parameters{i,1} = output.opt_parameters.(['LV_', num2str(i)]){log_find,1};
        end
    end
end

% compute epsilon and omega all by projecting u and v onto X and Y
try temp = load(input.X);
    temp_names = fieldnames(temp);
    X = temp.(temp_names{1});
catch
    X = input.X;
end

Y = input.Y;

RHO=[]; epsilon_all={}; omega_all={};

for i=1:size(output.final_parameters,1)
    
    if ~input.corrected_log(i)
        Covars = nan(size(input.Y,1),1);
        correction_target = 3;
    else
        Covars = input.covariates;
    end
    
    [OUT_x, OUT_y] = dp_master_correctscale(X, Y, Covars, input.scaling_method, input.correction_target);
    
    u = output.final_parameters{i, strcmp(output.parameters_names, 'u')};
    v = output.final_parameters{i, strcmp(output.parameters_names, 'v')};
    
    [RHO(i,1), epsilon_all{i,1}, omega_all{i,1}, ~, ~] = dp_projection(OUT_x, OUT_y, u, v, input.correlation_method);
    
    mdl = fitlm(epsilon_all{i,1}, omega_all{i,1});
    output.measures.linear_models{i} = mdl;
    output.measures.Rsquared_Ordinary(i) = mdl.Rsquared;
    
    [X, Y] = proj_def(X, Y, u, v);
    
end

output.final_parameters = [output.final_parameters, epsilon_all, omega_all];
output.parameters_names = [output.parameters_names, 'epsilon_all', 'omega_all'];

%% load and add additional handedness data

for i=1:size(input.data_collection.PSN_BOG,1)
    try
        EHI_01(i,1)=cell2mat(data_table.EHI_01(str2num(cell2mat(data_table{:, 'PATIENT_ID'})) == str2num(input.data_collection.PSN_BOG{i,1})));
    catch
        disp([input.data_collection.PSN_BOG{i}, ' is missing.']);
        EHI_01(i,1)=NaN; % EHI am Ende hinzufügen
    end
end

EHI_01 = EHI_01==2;

try temp = load(input.MRI);
catch
    temp = load(input.X);
end
field = fieldnames(temp);
MRI_for_analysis = temp.(field{1});

%% get data
% for i=1:size(input.data_collection.PSN_BOG,1)
data_table_selected=[];
for i=1:size(input.data_collection.PSN_BOG,1)
    try data_table_selected(i,:)=data_table_NM(contains(ID_name, input.data_collection.PSN_BOG{i,1}),:);
    catch
        disp([input.data_collection.PSN_BOG{i}, ' is missing.']);
        data_table_selected(i,:)=NaN;
    end
end

sociodem_variables = {'DEMOG_T0T1T2_31AA_EducationYears_T0', 'GF_S_1_Current_T0',...
    'GF_R_1_Current_T0', 'GAF_S_PastMonth_Screening', 'GAF_DI_PastMonth_Screening'};

sociodem_collection=[]; sociodem_collection_names={};
for i=1:size(sociodem_variables,2)
    temp = data_table_selected(:,contains(data_table_names_NM, sociodem_variables{i}));
    sociodem_collection=[sociodem_collection,temp];
    temp_names = data_table_names_NM(contains(data_table_names_NM, sociodem_variables{i}));
    sociodem_collection_names=[sociodem_collection_names,temp_names];
end

try sociodem_collection = [input.data_collection.age, input.data_collection.sex(:,2), input.data_collection.IQR, EHI_01, sociodem_collection];
catch
    input.data_collection.age = input.data_complete.foranalysis.basic.age;
    input.data_collection.sex = input.data_complete.foranalysis.basic{:, {'male_sex', 'female_sex'}};
    input.data_collection.IQR = input.data_complete.foranalysis.mri.Cat12_IQR_sMRI_T0;
    input.data_collection.Diagfull_names = input.selected_studygroups;
    input.data_collection.Labels = input.data_complete.foranalysis.basic{input.final_PSN, 'Labels'};
    sociodem_collection = [input.data_collection.age, input.data_collection.sex(:,2), input.data_collection.IQR, EHI_01, sociodem_collection];
end
sociodem_collection_names = ['age', 'female_sex', 'IQR', 'handedness', sociodem_collection_names];

% clinical_collection = clinical_collection(:,contains(clinical_collection_names, 'T0'));
% clinical_collection_names = clinical_collection_names(:,contains(clinical_collection_names, 'T0'));

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

clinical_collection = [panss_total, panss_pos, panss_neg, panss_gen, BDI_scores];
clinical_collection_names = {'PANSS_total', 'PANSS_pos', 'PANSS_neg', 'PANSS_gen', 'BDI_scores'};

output.socio_clin_data.collection = [sociodem_collection, clinical_collection];
output.socio_clin_data.collection_names = [sociodem_collection_names, clinical_collection_names];

output.socio_clin_data.means_std=[];

for i=1:size(unique(input.data_collection.Diag),1)
    output.socio_clin_data.means_std(1:2:size(output.socio_clin_data.collection,2)*2,i) = nanmean(output.socio_clin_data.collection(input.data_collection.Diag==i,:),1)';
    output.socio_clin_data.means_std(2:2:size(output.socio_clin_data.collection,2)*2,i) = nanstd(output.socio_clin_data.collection(input.data_collection.Diag==i,:),1)';
end

temp_all=[];
temp_all(1:2:size(output.socio_clin_data.collection,2)*2,1) = nanmean(output.socio_clin_data.collection,1)';
temp_all(2:2:size(output.socio_clin_data.collection,2)*2,1) = nanstd(output.socio_clin_data.collection,1)';

output.socio_clin_data.means_std = [temp_all, output.socio_clin_data.means_std];
output.socio_clin_data.names = output.socio_clin_data.collection_names';
output.socio_clin_data.labels = ['all', input.data_collection.Diagfull_names];

% get site numbers for diagnoses
output.sites_data.counts=[];
for i=1:size(input.sites,2)
    [~, ~, ic_sites] = unique(input.data_collection.Diag(input.sites(:,i)>0));
    output.sites_data.counts(i,:) = accumarray(ic_sites, 1)';
end

[~, ~, ic_sites] = unique(input.data_collection.Diag);
temp_all = accumarray(ic_sites, 1)';

output.sites_data.counts = [temp_all; output.sites_data.counts];
output.sites_data.names = input.data_collection.Diagfull_names;
output.sites_data.labels = input.sites_names;

% test for site imbalances
bins = 0:(size(output.sites_data.labels,2)-1);
x = sum(output.sites_data.counts(2:end,:),2);
obsCounts = x;
expCounts = ones(7,1)*round(mean(obsCounts));
try [h,p,st] = chi2gof(bins,'Ctrs',bins, 'Frequency',obsCounts, 'Expected',expCounts,'NParams',0);
catch
    difference = sum(expCounts)-sum(obsCounts);
    obsCounts(end) = obsCounts(end)+difference;
    [h,p,st] = chi2gof(bins,'Ctrs',bins, 'Frequency',obsCounts, 'Expected',expCounts,'NParams',0);
end

output.sites_data.main_stats = [h,p,st.chi2stat];
labels = output.sites_data.labels; nn=1;
for ii=1:size(labels,2)
    for iii=(ii+1):size(labels,2)
        if ii~=iii
            bins = 0:1;
            x_1 = x(contains(labels, labels(ii)));
            x_2 = x(contains(labels, labels(iii)));
            obsCounts = [sum(x_1), sum(x_2)];
            expCounts = [round(mean(obsCounts)), sum([x_1; x_2])-round(mean(obsCounts))];
            [h,p,st] = chi2gof(bins,'Ctrs',bins, 'Frequency',obsCounts, 'Expected',expCounts,'NParams',0);
            temp_sites(nn,:) = [st.chi2stat, p];
            row_names{nn,1} = [labels{ii}, '_', labels{iii}];
        else
            temp_sites(nn,:) = [NaN, NaN];
            row_names{nn,1} = ['NaN_', num2str(nn)];
        end
        nn=nn+1;
    end
end

output.sites_data.multcomp_stats = array2table(temp_sites, 'RowNames', strrep(row_names, 'INSTITUTE_SHORTNAME_sMRI_T0_', ''), 'VariableNames', {'Chi2Stat', 'p'});
output.sites_data.multcomp_stats.p = dp_FDR_adj(output.sites_data.multcomp_stats.p);

%% get medication and IQ
IN.groups = {'ROD', 'CHR', 'ROP'};
output.meds_and_iq = dp_pronia_meds_and_IQ(IN);

%% test for significant differences between groups
% first test with ANOVA for overall differences
output.socio_clin_data.sociodem_results=[]; output.socio_clin_data.sociodem_results_ext={};
output.socio_clin_data.chi2square_mult=[]; output.socio_clin_data.chi2square_mult_ext={};
for i=1:size(output.socio_clin_data.names,1)
    if size(unique(output.socio_clin_data.collection(:,i)),1)<4
        bins = 0:1;
        x = output.socio_clin_data.collection(:,i);
        obsCounts = [sum(x), size(x,1)-sum(x)];
        expCounts = [round(mean(obsCounts)), size(x,1)-round(mean(obsCounts))];
        [h,p,st] = chi2gof(bins,'Ctrs',bins, 'Frequency',obsCounts, 'Expected',expCounts,'NParams',0);
        temp_sociodem_results(i,:) = [st.chi2stat, p];
        output.socio_clin_data.sociodem_results_ext(i,:) = {p,h,st};
        labels = unique(input.data_collection.Labels);
        nn=1; temp_mult=[]; row_names={};
        for ii=1:size(labels,1)
            for iii=(ii+1):size(labels,1)
                if ii~=iii
                    bins = 0:1;
                    x_1 = output.socio_clin_data.collection(contains(input.data_collection.Labels, labels(ii)),i);
                    x_2 = output.socio_clin_data.collection(contains(input.data_collection.Labels, labels(iii)),i);
                    obsCounts = [sum(x_1), sum(x_2)];
                    expCounts = [round(mean(obsCounts)), sum([x_1; x_2])-round(mean(obsCounts))];
                    [h,p,st] = chi2gof(bins,'Ctrs',bins, 'Frequency',obsCounts, 'Expected',expCounts,'NParams',0);
                    temp_mult(nn,:) = [st.chi2stat, p];
                    output.socio_clin_data.chi2square_mult_ext.(output.socio_clin_data.names{i})(nn,:)= {p,h,st};
                    row_names{nn,1} = [labels{ii}, '_', labels{iii}];
                else
                    temp_mult(nn,:) = [NaN, NaN];
                    output.socio_clin_data.chi2square_mult_ext.(output.socio_clin_data.names{i})(nn,:)= {NaN,NaN,NaN};
                    row_names{nn,1} = ['NaN_', num2str(nn)];
                end
                nn=nn+1;
            end
        end
        output.socio_clin_data.chi2square_mult.(output.socio_clin_data.names{i}) = array2table(temp_mult, 'RowNames', row_names, 'VariableNames', {'Chi2Stat', 'p'});
        output.socio_clin_data.chi2square_mult.(output.socio_clin_data.names{i}).p = dp_FDR_adj(output.socio_clin_data.chi2square_mult.(output.socio_clin_data.names{i}).p);
    else
        [p,tbl,stats] = kruskalwallis(output.socio_clin_data.collection(:,i), input.data_collection.Labels);
        close all
        temp_sociodem_results(i,:) = [tbl{2,5}, p];
        output.socio_clin_data.sociodem_results_ext(i,:) = {p,tbl,stats};
        [output.socio_clin_data.Dunn_results{i,1}, ~, ~, output.socio_clin_data.Dunn_results_labels] = multcompare(stats, 'CType', 'dunn-sidak', 'Estimate', 'kruskalwallis');
        close all
    end
end

output.socio_clin_data.sociodem_results = array2table(temp_sociodem_results, 'RowNames', output.socio_clin_data.names, 'VariableNames', {'Chi2Stat', 'p'});
output.socio_clin_data.sociodem_results.p = dp_FDR_adj(output.socio_clin_data.sociodem_results.p);

%% look for specific correlations

output.post_hoc_correlations.data_collection=[];
output.post_hoc_correlations.names=[];

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
        temp = data_table_selected(:, contains(data_table_names_NM, GAF_names{i}));
        GAF_data=[GAF_data,temp];
    end
    output.post_hoc_correlations.names = [output.post_hoc_correlations.names, GAF_names];
    output.post_hoc_correlations.data_collection = [output.post_hoc_correlations.data_collection, GAF_data];
end

if any(strcmp(IN.SD_selection, 'BDI'))
    
    % get BDI data
    output.post_hoc_correlations.names = [output.post_hoc_correlations.names, 'BDI'];
    output.post_hoc_correlations.data_collection = [output.post_hoc_correlations.data_collection, BDI_scores];
    
end

if any(strcmp(IN.SD_selection, 'neurocognition'))
    
    % get extraction data
    load('/volume/RU_DP_immune/Data/Immune_megapaper_request_Disc/Immune_megapaper_request_Disc/DATA/18-Dec-2020/Immune_megapaper_request_Disc_Data_all_18-Dec-2020.mat');
    temp_data_table = data_table_all;
    
    load('/volume/RU_DP_immune/Data/Immune_megapaper_request_Repl/Immune_megapaper_request_Repl/DATA/18-Dec-2020/Immune_megapaper_request_Repl_Data_all_18-Dec-2020.mat');
    temp_names = temp_data_table.Properties.VariableNames(matches(temp_data_table.Properties.VariableNames, data_table_all.Properties.VariableNames));
    data_table_all = [temp_data_table(:, temp_names); data_table_all(:, temp_names)];
    
    NC = rs_neurocognition(data_table_all);
    output.neurocognition_data=table;
    for i=1:size(input.data_collection.PSN_BOG,1)
        log_temp = str2num(cell2mat(data_table_all.PSN)) == str2num(input.data_collection.PSN_BOG{i,1});
        if sum(log_temp)~=0
            output.neurocognition_data(i,:) = NC.single_scores(log_temp, :);
        else
            output.neurocognition_data(i,:) = array2table(nan(1, size(NC.single_scores,2)), 'VariableNames', NC.single_scores.Properties.VariableNames);
        end
    end
    output.neurocognition_data.Properties.RowNames = input.data_collection.PSN_BOG;
    
    output.post_hoc_correlations.names = [output.post_hoc_correlations.names, output.neurocognition_data.Properties.VariableNames];
    output.post_hoc_correlations.data_collection = [output.post_hoc_correlations.data_collection, output.neurocognition_data.Variables];
    
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
    
    output.post_hoc_correlations.names = [output.post_hoc_correlations.names,neo_names];
    output.post_hoc_correlations.data_collection = [output.post_hoc_correlations.data_collection,neo_data_sum];
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
    
    output.post_hoc_correlations.names = [output.post_hoc_correlations.names,WHO_names];
    output.post_hoc_correlations.data_collection = [output.post_hoc_correlations.data_collection,WHO_sum];
end

output.post_hoc_correlations.data_table = array2table(output.post_hoc_correlations.data_collection, 'RowNames', input.data_collection.PSN_BOG(:,1));
output.post_hoc_correlations.data_table.Properties.VariableNames = output.post_hoc_correlations.names;

%% proper correlations
log_epsilon = strcmp(output.parameters_names, 'epsilon');
log_omega = strcmp(output.parameters_names, 'omega');
log_epsilon_all = strcmp(output.parameters_names, 'epsilon_all');
log_omega_all = strcmp(output.parameters_names, 'omega_all');
log_coll = {log_epsilon, log_epsilon_all; log_omega, log_omega_all; 'test', 'all'};

for l=1:size(log_coll,2)
    output.post_hoc_correlations.(log_coll{3,l}).correlations.test.RHO=[];
    output.post_hoc_correlations.(log_coll{3,l}).correlations.test.p=[];
    output.post_hoc_correlations.(log_coll{3,l}).correlations.validation.RHO=[];
    output.post_hoc_correlations.(log_coll{3,l}).correlations.validation.p=[];
    output.post_hoc_mdl.(log_coll{3,l}).test.R2=[];
    output.post_hoc_mdl.(log_coll{3,l}).test.p=[];
    output.post_hoc_mdl.(log_coll{3,l}).validation.R2=[];
    output.post_hoc_mdl.(log_coll{3,l}).validation.p=[];
    
    for i=1:(size(output.final_parameters,1)-1)
        
        epsilon = output.final_parameters{i,log_coll{1,l}};
        omega = output.final_parameters{i,log_coll{2,l}};
        X = [epsilon, omega];
        
        if input.selection_train == 1
            if strcmp(log_coll{3,l}, 'test')
                y = output.post_hoc_correlations.data_collection(output.CV.cv_outer_indices.TestInd{1,output.final_parameters{i,1}},:);
            elseif strcmp(log_coll{3,l}, 'all')
                idx = [];
                for ii=1:size(output.CV.cv_outer_indices.TestInd,2)
                    idx = [idx; output.CV.cv_outer_indices.TestInd{1,ii}];
                end
                X = X(idx,:);
                y = output.post_hoc_correlations.data_collection(idx,:);
            end
        elseif input.selection_train == 2
            idx = [];
            for ii=1:size(output.CV.cv_outer_indices.TestInd,2)
                idx = [idx; output.CV.cv_outer_indices.TestInd{1,ii}];
            end
            y = output.post_hoc_correlations.data_collection(idx,:);
        end
        
        [RHO,p] = corr(X, y, 'Type', 'Spearman', 'rows', 'complete');
        %     [p(1,:), ~] = dp_FDR_adj(p(1,:));
        %     [p(2,:), ~] = dp_FDR_adj(p(2,:));
        
        output.post_hoc_correlations.(log_coll{3,l}).correlations.test.RHO = [output.post_hoc_correlations.(log_coll{3,l}).correlations.test.RHO, RHO'];
        output.post_hoc_correlations.(log_coll{3,l}).correlations.test.p = [output.post_hoc_correlations.(log_coll{3,l}).correlations.test.p, p'];
        
        p=[]; R2=[]; mdl_coll={};
        
        for v=1:size(y,2)
            mdl = fitlm(X, y(:,v));
            p(v,1) = mdl.coefTest;
            R2(v,1) = mdl.Rsquared.Adjusted;
            %         mdl_coll = [mdl_coll; {mdl}];
        end
        
        output.post_hoc_mdl.(log_coll{3,l}).test.R2 = [output.post_hoc_mdl.(log_coll{3,l}).test.R2, R2];
        output.post_hoc_mdl.(log_coll{3,l}).test.p = [output.post_hoc_mdl.(log_coll{3,l}).test.p, p];
        %     output.post_hoc_mdl.(log_coll{3,l}).test.mdl = [output.post_hoc_mdl.(log_coll{3,l}).test.mdl, mdl_coll];
        
        %     [output.post_hoc_mdl.(log_coll{3,l}).test.p, ~] = dp_FDR_adj(output.post_hoc_mdl.(log_coll{3,l}).test.p);
        
        if ~islogical(input.validation_set)
            epsilon = output.validation_results{i,log_epsilon};
            omega = output.validation_results{i,log_omega};
            X = [epsilon, omega];
            y = output.post_hoc_correlations.data_collection(output.validation_indices.TestInd{1,1},:);
            [RHO,p] = corr(X, y, 'Type', 'Spearman', 'rows', 'complete');
            %         [p(1,:), ~] = dp_FDR_adj(p(1,:));
            %         [p(2,:), ~] = dp_FDR_adj(p(2,:));
            
            output.post_hoc_correlations.correlations.validation.RHO = [output.post_hoc_correlations.correlations.validation.RHO, RHO'];
            output.post_hoc_correlations.correlations.validation.p = [output.post_hoc_correlations.correlations.validation.p, p'];
            
            p=[]; R2=[];mdl_coll={};
            
            for v=1:size(y,2)
                mdl = fitlm(X, y(:,v));
                p(v,1) = mdl.coefTest;
                R2(v,1) = mdl.Rsquared.Adjusted;
                %             mdl_coll=[mdl_coll; {mdl}];
            end
            
            output.post_hoc_mdl.validation.R2 = [output.post_hoc_mdl.validation.R2, R2];
            output.post_hoc_mdl.validation.p = [output.post_hoc_mdl.validation.p, p];
            %         output.post_hoc_mdl.(log_coll{3,l}).validation.mdl = [output.post_hoc_mdl.(log_coll{3,l}).validation.mdl, mdl_coll];
            
            %         [output.post_hoc_mdl.(log_coll{3,l}).validation.p, ~] = dp_FDR_adj(output.post_hoc_mdl.(log_coll{3,l}).validation.p);
            output.post_hoc_mdl.validation.p = dp_FDR_adj(output.post_hoc_mdl.validation.p);
        end
        
    end
    
    output.post_hoc_correlations.(log_coll{3,l}).correlations.test.p = dp_FDR_adj(output.post_hoc_correlations.(log_coll{3,l}).correlations.test.p);
    output.post_hoc_mdl.(log_coll{3,l}).test.p = dp_FDR_adj(output.post_hoc_mdl.(log_coll{3,l}).test.p);
    
    output.post_hoc_correlations.(log_coll{3,l}).correlations.test.table_RHO = array2table(output.post_hoc_correlations.(log_coll{3,l}).correlations.test.RHO, 'RowNames', output.post_hoc_correlations.names);
    output.post_hoc_correlations.(log_coll{3,l}).correlations.test.table_p = array2table(output.post_hoc_correlations.(log_coll{3,l}).correlations.test.p, 'RowNames', output.post_hoc_correlations.names);
    output.post_hoc_mdl.(log_coll{3,l}).test.table_p = array2table(output.post_hoc_mdl.(log_coll{3,l}).test.p, 'RowNames', output.post_hoc_correlations.names);
    output.post_hoc_mdl.(log_coll{3,l}).test.table_R2 = array2table(output.post_hoc_mdl.(log_coll{3,l}).test.R2, 'RowNames', output.post_hoc_correlations.names);
    
    temp=[];temp_names={};nn=1;
    for rr=1:size(output.post_hoc_correlations.(log_coll{3,l}).correlations.test.RHO,1)
        temp(nn,:) =  output.post_hoc_correlations.(log_coll{3,l}).correlations.test.RHO(rr,:);
        temp_names{nn,1} = [output.post_hoc_correlations.names{rr}, '_RHO'];
        nn=nn+1;
        temp(nn,:) =  output.post_hoc_correlations.(log_coll{3,l}).correlations.test.p(rr,:);
        temp_names{nn,1} = [output.post_hoc_correlations.names{rr}, '_p'];
        nn=nn+1;
    end
    
    output.post_hoc_correlations.(log_coll{3,l}).correlations.test.table_RHO_p = array2table(temp, 'RowNames', temp_names);
    
    if ~islogical(input.validation_set)
        output.post_hoc_correlations.(log_coll{3,l}).correlations.validation.p = dp_FDR_adj(output.post_hoc_correlations.(log_coll{3,l}).correlations.validation.p);
        output.post_hoc_correlations.(log_coll{3,l}).correlations.validation.table_RHO = array2table(output.post_hoc_correlations.correlations.validation.RHO, 'RowNames', output.post_hoc_correlations.names);
        output.post_hoc_correlations.(log_coll{3,l}).correlations.validation.table_p = array2table(output.post_hoc_correlations.correlations.validation.p, 'RowNames', output.post_hoc_correlations.names);
        output.post_hoc_mdl.validation.table_p = array2table(output.post_hoc_mdl.validation.p, 'RowNames', output.post_hoc_correlations.names);
        output.post_hoc_mdl.validation.table_R2 = array2table(output.post_hoc_mdl.validation.R2, 'RowNames', output.post_hoc_correlations.names);
        
        temp=[];temp_names={};nn=1;
        for rr=1:size(output.post_hoc_correlations.correlations.validation.RHO,1)
            temp(nn,:) =  output.post_hoc_correlations.correlations.validation.RHO(rr,:);
            temp_names{nn,1} = [output.post_hoc_correlations.names{rr}, '_RHO'];
            nn=nn+1;
            temp(nn,:) =  output.post_hoc_correlations.correlations.validation.p(rr,:);
            temp_names{nn,1} = [output.post_hoc_correlations.names{rr}, '_p'];
            nn=nn+1;
        end
        
        output.post_hoc_correlations.correlations.validation.table_RHO_p = array2table(temp, 'RowNames', temp_names);
    end
    
    nn=1;temp_vars={};
    for i=1:size(output.post_hoc_correlations.(log_coll{3,l}).correlations.test.table_RHO,2)/2
        temp_vars{1,nn} = ['epsilon_LV', num2str(i)];
        nn=nn+1;
        temp_vars{1,nn} = ['omega_LV', num2str(i)];
        nn=nn+1;
    end
    
    temp_vars_names={};
    for i=1:size(output.post_hoc_correlations.(log_coll{3,l}).correlations.test.table_RHO,2)/2
        temp_vars_names{1,i} = ['LV', num2str(i), '_latent_scores'];
    end
    
    output.post_hoc_correlations.(log_coll{3,l}).correlations.test.table_RHO.Properties.VariableNames = temp_vars;
    output.post_hoc_correlations.(log_coll{3,l}).correlations.test.table_p.Properties.VariableNames = temp_vars;
    output.post_hoc_correlations.(log_coll{3,l}).correlations.test.table_RHO_p.Properties.VariableNames = temp_vars;
    output.post_hoc_mdl.(log_coll{3,l}).test.table_p.Properties.VariableNames = temp_vars_names;
    output.post_hoc_mdl.(log_coll{3,l}).test.table_R2.Properties.VariableNames = temp_vars_names;
    
    if ~islogical(input.validation_set)
        
        output.post_hoc_correlations.correlations.validation.table_RHO.Properties.VariableNames = temp_vars;
        output.post_hoc_correlations.correlations.validation.table_p.Properties.VariableNames = temp_vars;
        output.post_hoc_correlations.correlations.validation.table_RHO_p.Properties.VariableNames = temp_vars;
        output.post_hoc_mdl.validation.table_p.Properties.VariableNames = temp_vars_names;
        output.post_hoc_mdl.validation.table_R2.Properties.VariableNames = temp_vars_names;
        
    end
    
    % save(IN.results_path, 'input', 'output', 'setup');
    
end

end


