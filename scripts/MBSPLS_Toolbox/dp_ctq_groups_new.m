%% new script to compute CTQ group differences

% m_discovery = matfile('/volume/HCStress/Analysis/27-Mar-2019/DP_CTQ_allgroups_649_GM_80PI_12GO_110X060_40GD_110Y010_40GD_correct10_10_XY_SC_3_Diag/final_results/result_new.mat', 'Writable', true);
% m_replication = matfile('/volume/HCStress/Analysis/25-Feb-2020_CTQ_627_replication_sample_CTQ_fixed_datafile.mat', 'Writable', true);
% 
% input_discovery = m_discovery.input;
% input_replication = m_replication.input;
% 
% output_discovery = m_discovery.output;

load('/volume/HCStress/Analysis/replication_results_SPLS_Stress/validation_replication_detailed_results_comb_app_1.mat')

% update neglect subscales
to_add1 = {'CTQ_02', 'CTQ_05', 'CTQ_07', 'CTQ_13', 'CTQ_19', 'CTQ_26', 'CTQ_28'};
CTQ_collection_discovery = input_discovery.data_collection.data;
CTQ_collection_discovery(:, contains(input_discovery.data_collection.names, to_add1)) = CTQ_collection_discovery(:, contains(input_discovery.data_collection.names, to_add1))+1;

CTQ_collection_replication = input_replication.data_collection.data;
CTQ_collection_replication(:, contains(input_replication.data_collection.names, to_add1)) = CTQ_collection_replication(:, contains(input_replication.data_collection.names, to_add1))+1;

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
    CTQ_discovery(:,i) = sum(CTQ_collection_discovery(:, CTQ.(fields{i})),2);
    CTQ_replication(:,i) = sum(CTQ_collection_replication(:, CTQ.(fields{i})),2);
end

input_replication.CTQ_collection = CTQ_replication;
input_replication.CTQ_collection_names = fieldnames(CTQ)';
output_discovery.CTQ_collection = CTQ_discovery;
output_discovery.CTQ_collection_names = fieldnames(CTQ)';

m_discovery.input = input_discovery;
m_replication.input = input_replication;
m_discovery.output = output_discovery;


labels={'HC', 'ROD', 'CHR', 'ROP'};
labels_disc = zeros(size(input_discovery.data_collection.Labels,1),1);
labels_rep = ones(size(input_replication.data_collection.Labels,1),1);
CTQ_results=[];temp_results=[];mean_std_results=[];temp_mean_std=[];

for i=1:size(CTQ_discovery,2)
    temp_results=[];temp_mean_std=[];
    temp_data = [CTQ_discovery(:,i); CTQ_replication(:,i)];    
    [p,tbl,stats] = kruskalwallis(temp_data, [labels_disc; labels_rep]);
    temp_results = [temp_results, p];
    temp_mean_std = [nanmean(CTQ_discovery(:,i)), nanstd(CTQ_discovery(:,i)); nanmean(CTQ_replication(:,i)), nanstd(CTQ_replication(:,i))];
    close all
    for ii=1:size(labels,2)
        log_temp_d = contains(input_discovery.data_collection.Labels, labels{ii});
        log_temp_r = contains(input_replication.data_collection.Labels, labels{ii});
        temp_data = [CTQ_discovery(log_temp_d,i); CTQ_replication(log_temp_r,i)];
        temp_mean_std = [temp_mean_std, [nanmean(CTQ_discovery(log_temp_d,i)), nanstd(CTQ_discovery(log_temp_d,i)); nanmean(CTQ_replication(log_temp_r,i)), nanstd(CTQ_replication(log_temp_r,i))]];
        [p,tbl,stats] = kruskalwallis(temp_data, [labels_disc(log_temp_d); labels_rep(log_temp_r)]);
        temp_results = [temp_results, p];
    end
    CTQ_results = [CTQ_results; temp_results];
    mean_std_results = [mean_std_results; temp_mean_std];
    close all
end

[CTQ_results, ~] = dp_FDR_adj(CTQ_results);

groups_to_choose = {'ROD', 'CHR', 'ROP'};
log_groups_disc = contains(input_discovery.data_collection.Labels, groups_to_choose);
log_groups_rep = contains(input_replication.data_collection.Labels, groups_to_choose);
CTQ_discovery_results_KW=[];CTQ_discovery_results_Dunn=[];
CTQ_replication_results_KW=[];CTQ_replication_results_Dunn=[];

for i=1:size(CTQ_discovery,2)
    [p,tbl,stats] = kruskalwallis(CTQ_discovery(log_groups_disc,i), input_discovery.data_collection.Labels(log_groups_disc));
    CTQ_discovery_results_KW(i,:) = [p, tbl{2,5}, tbl{2,3}];
    CTQ_discovery_results_Dunn{i,1} = multcompare(stats, 'Estimate', 'kruskalwallis', 'CType', 'dunn-sidak');
    close all
    
    [p,tbl,stats] = kruskalwallis(CTQ_replication(log_groups_rep,i), input_replication.data_collection.Labels(log_groups_rep));
    CTQ_replication_results_KW(i,:) = [p, tbl{2,5}, tbl{2,3}];
    CTQ_replication_results_Dunn{i,1} = multcompare(stats, 'Estimate', 'kruskalwallis', 'CType', 'dunn-sidak');
    
    close all
    
    %    CTQ_discovery_results = [CTQ_discovery_results,
end

[CTQ_discovery_results_KW(:,1), ~] = dp_FDR_adj(CTQ_discovery_results_KW(:,1));
[CTQ_replication_results_KW(:,1), ~] = dp_FDR_adj(CTQ_replication_results_KW(:,1));
    
    
    
    
    

