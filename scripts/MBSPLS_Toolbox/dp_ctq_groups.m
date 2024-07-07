%% new script to compute CTQ group differences

m_discovery = matfile('/volume/HCStress/Analysis/27-Mar-2019/DP_CTQ_allgroups_649_GM_80PI_12GO_110X060_40GD_110Y010_40GD_correct10_10_XY_SC_3_Diag/final_results/result_new.mat');
m_replication = matfile('/volume/HCStress/Analysis/25-Feb-2020_CTQ_627_replication_sample_CTQ_fixed_datafile.mat');

input_discovery = m_discovery.input;
input_replication = m_replication.input;

% get CTQ data and compute total and subscale scores
CTQ.total = 1:25;
CTQ.emotional_abuse = 1:5; %[3,8,14,18,25];
CTQ.physical_abuse = 6:10; %[9,11,12,15,17];
CTQ.sexual_abuse = 11:15; %[20,21,23,24,27];
CTQ.emotional_neglect = 16:20; %[5,7,13,19,28];
CTQ.physical_neglect = 21:25; %[1,2,4,6,26];
CTQ.denial = 26:28; %[10,16,22];
CTQ_discovery=[];CTQ_replication=[];
fields = fieldnames(CTQ);
for i=1:size(fields,1)
    CTQ_discovery(:,i) = sum(input_discovery.data_collection.data(:, CTQ.(fields{i})),2);
    CTQ_replication(:,i) = sum(input_replication.data_collection.data(:, CTQ.(fields{i})),2);
end

labels={'HC', 'ROD', 'CHR', 'ROP'};
labels_disc = zeros(size(input_discovery.data_collection.Labels,1),1);
labels_rep = ones(size(input_replication.data_collection.Labels,1),1);
CTQ_results=[];temp_results=[];

for i=1:size(CTQ_discovery,2)
    temp_results=[];
    temp_data = [CTQ_discovery(:,i); CTQ_replication(:,i)];
    [p,tbl,stats] = kruskalwallis(temp_data, [labels_disc; labels_rep]);
    temp_results = [temp_results, p];
    close all
    for ii=1:size(labels,2)
        log_temp_d = contains(input_discovery.data_collection.Labels, labels{ii});
        log_temp_r = contains(input_replication.data_collection.Labels, labels{ii});
        temp_data = [CTQ_discovery(log_temp_d,i); CTQ_replication(log_temp_r,i)];
        [p,tbl,stats] = kruskalwallis(temp_data, [labels_disc(log_temp_d); labels_rep(log_temp_r)]);
        temp_results = [temp_results, p];
    end
    CTQ_results = [CTQ_results; temp_results];
    close all
end



temp = input.data_collection.data;
for i=1:size(temp,1)
    CTQ_data_subscales(i,:) = [sum(temp(i,CTQ.emotional_abuse)), sum(temp(i,CTQ.physical_abuse)), sum(temp(i,CTQ.sexual_abuse)), sum(temp(i,CTQ.emotional_neglect)), sum(temp(i,CTQ.physical_neglect)), sum(temp(i,CTQ.denial))]; 
end

CTQ_data = [CTQ_data_subscales, CTQ_data_total];
CTQ_data_names = {'emotional_abuse', 'physical_abuse', 'sexual_abuse', 'emotional_neglect', 'physical_neglect', 'denial', 'total'};
groups_to_choose = {'HC', 'ROD', 'CHR', 'ROP'};
log_groups_to_choose =  ismember(input.data_collection.Labels, groups_to_choose);

s_folder = '/volume/HCStress/Analysis/27-Mar-2019/DP_CTQ_allgroups_649_GM_80PI_12GO_110X060_40GD_110Y010_40GD_correct10_10_XY_SC_3_Diag/final_results';

for i=1:size(CTQ_data,2)
    [p,tbl,stats] = kruskalwallis(CTQ_data(log_groups_to_choose,i), input.data_collection.Labels(log_groups_to_choose));
    output.CTQ.KW_results(i,:) = [tbl{2,5}, p];
    output.CTQ.KW_results_ext(i,:) = {p,tbl,stats};
    [output.CTQ.Dunn_results{i,1}, ~, ~, output.CTQ.Dunn_results_labels] = multcompare(stats, 'CType', 'dunn-sidak', 'Estimate', 'kruskalwallis');
    close all
%     dp_txt_write(s_folder, ['Dunn_results_CTQ_', num2str(i)], output.CTQ.Dunn_results{i,1}', '%.3f \t %.3f \t %.3f \t %.3f \t %.3f \t %.3f \t \n \n');
end

% dp_txt_write(s_folder, 'KW_CTQ_results', output.CTQ.KW_results', '%.3f \t %.3f \n \n');
% dp_txt_write(s_folder, 'Dunn_CTQ_labels', '', '%s');
% for i=1:size(output.CTQ.Dunn_results_labels,1)
%     FID = fopen([s_folder, '/Dunn_CTQ_labels.txt'], 'a');
%     fprintf(FID, '%s\n', output.CTQ.Dunn_results_labels{i,1});
%     fclose(FID);
%     dp_txt_write(s_folder, 'Dunn_CTQ_labels', [output.CTQ.Dunn_results_labels{:}], '%s \n');
% end

% adjust p values for multiple testing
[output.CTQ.KW_results_FDR_threshold,~,output.CTQ.KW_results(:,2)] = fdr(output.CTQ.KW_results(:,2));

% get mean and SD values for CTQ data
descriptive_groups = {'HC', 'ROD', 'CHR', 'ROP'};
log_groups_to_choose =  ismember(input.data_collection.Labels, descriptive_groups);
temp_descriptive_stats=[];
temp_all_descriptive_stats=[];
for i=1:size(CTQ_data,2)
    temp_mean_all = mean(CTQ_data(:,i));
    temp_std_all = std(CTQ_data(:,i));
    temp_all_descriptive_stats=[temp_all_descriptive_stats;temp_mean_all; temp_std_all];
    for ii=1:size(descriptive_groups,2)
       temp_means(1,ii) = mean(CTQ_data(ismember(input.data_collection.Labels, descriptive_groups{ii}),i));
       temp_std(1,ii) = std(CTQ_data(ismember(input.data_collection.Labels, descriptive_groups{ii}),i));
    end
    temp_descriptive_stats=[temp_descriptive_stats;temp_means; temp_std];
end

output.CTQ.descriptive_stats = [temp_all_descriptive_stats, temp_descriptive_stats];

save(results_path, 'input', 'setup', 'output');

% [pthr,pcor,padj] = fdr(pvals)
% 
% pvals=[0.127357313939555;0.223286174661224;0.987758230437541; 0.203652541916409;0.255893623925517;0.998512378851066;...
%     0.994306327317029;0.608632053506984;0.767610213762934;0.997076468205630;0.558523555856795;0.691558985229927;...
%     0.0583756086858885;0.609101209753551;0.511036107759435;0.0412122143925899;0.0374061120934872;0.999999972063131;...
%     0.170652060156671;0.0674167894612940;0.977288374963093]
% padj'

log_groups_to_choose =  input.data_collection.sex>0;

s_folder = '/volume/HCStress/Analysis/27-Mar-2019/DP_CTQ_allgroups_649_GM_80PI_12GO_110X060_40GD_110Y010_40GD_correct10_10_XY_SC_3_Diag/final_results';

for i=1:size(CTQ_data,2)
    [p,h,stats] = ranksum(CTQ_data(log_groups_to_choose(:,1),i), CTQ_data(log_groups_to_choose(:,2),i));
    output.CTQ.male_female_MW_results(i,:) = [h, p];
    output.CTQ.male_female_MW_results_ext(i,:) = {p,h,stats};
%     [output.CTQ.Dunn_results{i,1}, ~, ~, output.CTQ.Dunn_results_labels] = multcompare(stats, 'CType', 'dunn-sidak', 'Estimate', 'kruskalwallis');
    close all
%     dp_txt_write(s_folder, ['Dunn_results_CTQ_', num2str(i)], output.CTQ.Dunn_results{i,1}', '%.3f \t %.3f \t %.3f \t %.3f \t %.3f \t %.3f \t \n \n');
end

dp_txt_write(s_folder, 'male_female_MW_results', output.CTQ.male_female_MW_results', '%.3f \t %.3f \n \n');
% dp_txt_write(s_folder, 'Dunn_CTQ_labels', '', '%s');
% for i=1:size(output.CTQ.Dunn_results_labels,1)
%     FID = fopen([s_folder, '/Dunn_CTQ_labels.txt'], 'a');
%     fprintf(FID, '%s\n', output.CTQ.Dunn_results_labels{i,1});
%     fclose(FID);
%     dp_txt_write(s_folder, 'Dunn_CTQ_labels', [output.CTQ.Dunn_results_labels{:}], '%s \n');
% end

% adjust p values for multiple testing
[output.CTQ.male_female_MW_results_FDR_threshold,~,output.CTQ.male_female_MW_results(:,2)] = fdr(output.CTQ.male_female_MW_results(:,2));

% get mean and SD values for CTQ data
descriptive_groups = {'male', 'female'};
log_groups_to_choose =  input.data_collection.sex>0;
temp_descriptive_stats=[];
temp_all_descriptive_stats=[];
for i=1:size(CTQ_data,2)
    temp_mean_all = mean(CTQ_data(:,i));
    temp_std_all = std(CTQ_data(:,i));
    temp_all_descriptive_stats=[temp_all_descriptive_stats;temp_mean_all; temp_std_all];
    for ii=1:size(descriptive_groups,2)
       temp_means(1,ii) = mean(CTQ_data(log_groups_to_choose(:,ii),i));
       temp_std(1,ii) = std(CTQ_data(log_groups_to_choose(:,ii),i));
    end
    temp_descriptive_stats=[temp_descriptive_stats;temp_means; temp_std];
end



save(results_path, 'input', 'setup', 'output');

