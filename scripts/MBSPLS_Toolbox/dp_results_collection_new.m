%% DP function to analyze results in varying grid densities
%% prepare indices

output.parameters_names = {'w', 'cu', 'cv', 'u', 'v', 'U_opt', 'S_opt' ,'V_opt', 'success', 'RHO', 'p'}; % names of parameters for optimal cu/cv combination of the w-th loop
opt_u = strcmp(output.parameters_names,'u');
opt_v = strcmp(output.parameters_names,'v');
opt_p = strcmp(output.parameters_names,'p');
opt_RHO = strcmp(output.parameters_names,'RHO');
opt_cu = strcmp(output.parameters_names,'cu');
opt_cv = strcmp(output.parameters_names,'cv');

opt_parameters_names = {'w', 'cu', 'cv', 'u', 'v', 'success', 'RHO', 'p'};
opt_param_p = strcmp(opt_parameters_names, 'p');
opt_param_v = strcmp(opt_parameters_names, 'v');
opt_param_u = strcmp(opt_parameters_names, 'u');
opt_param_RHO = strcmp(opt_parameters_names, 'RHO');
opt_param_cu = strcmp(opt_parameters_names, 'cu');
opt_param_cv = strcmp(opt_parameters_names, 'cv');

load('/volume/HCStress/Doc/Stress_Resilience_questionnaires.mat');

%% % variables to assess: significant_LVs, range_RHO, range_p
results_paths_stress = {'/volume/HCStress/Analysis/20-Sep-2018/DP_CTQ_BS_HC_ROD_CHR_ROP_626_GM_80PI_4GD_correct_1/final_results/result.mat',...
    '/volume/HCStress/Analysis/22-Oct-2018/DP_CTQ_BS_HC_ROD_CHR_ROP_627_GM_80PI_3GD_correct_1/final_results/result.mat',...
    '/volume/HCStress/Analysis/20-Sep-2018/DP_CTQ_BS_HC_ROD_CHR_ROP_626_GM_80PI_6GD_correct/final_results/result.mat',...
    '/volume/HCStress/Analysis/22-Oct-2018/DP_CTQ_BS_HC_ROD_CHR_ROP_627_GM_80PI_5GD_correct/final_results/result.mat',...
    '/volume/HCStress/Analysis/22-Sep-2018/DP_CTQ_BS_HC_ROD_CHR_ROP_626_GM_80PI_8GD_correct_2/final_results/result.mat',...
    '/volume/HCStress/Analysis/23-Oct-2018/DP_CTQ_BS_HC_ROD_CHR_ROP_627_GM_80PI_7GD_correct/final_results/result.mat',...
    '/volume/HCStress/Analysis/24-Sep-2018/DP_CTQ_BS_HC_ROD_CHR_ROP_626_GM_80PI_8GD_correct/final_results/result.mat',...
    '/volume/HCStress/Analysis/23-Oct-2018/DP_CTQ_BS_HC_ROD_CHR_ROP_627_GM_80PI_9GD_correct/final_results/result.mat',...
    '/volume/HCStress/Analysis/25-Sep-2018/DP_CTQ_BS_HC_ROD_CHR_ROP_626_GM_80PI_10GD_correct/final_results/result.mat',...
    '/volume/HCStress/Analysis/27-Sep-2018/DP_CTQ_BS_HC_ROD_CHR_ROP_626_GM_80PI_12GD_correct/final_results/result.mat',...
    '/volume/HCStress/Analysis/01-Oct-2018/DP_CTQ_BS_HC_ROD_CHR_ROP_626_GM_80PI_14GD_correct/final_results/result.mat'};

results_paths_resilience = {'/volume/HCStress/Analysis/21-Sep-2018/DP_CISS_RSA_HC_ROD_CHR_ROP_631_GM_80PI_4GD_correct/final_results/result.mat',...
    '/volume/HCStress/Analysis/23-Oct-2018/DP_CISS_RSA_HC_ROD_CHR_ROP_634_GM_80PI_3GD_correct/final_results/result.mat',...
    '/volume/HCStress/Analysis/22-Sep-2018/DP_CISS_RSA_HC_ROD_CHR_ROP_631_GM_80PI_6GD_correct/final_results/result.mat',...
    '/volume/HCStress/Analysis/23-Oct-2018/DP_CISS_RSA_HC_ROD_CHR_ROP_634_GM_80PI_5GD_correct/final_results/result.mat',...
    '/volume/HCStress/Analysis/22-Sep-2018/DP_CISS_RSA_HC_ROD_CHR_ROP_631_GM_80PI_8GD_correct_3/final_results/result.mat',...
    '/volume/HCStress/Analysis/26-Oct-2018/DP_CISS_RSA_HC_ROD_CHR_ROP_634_GM_80PI_7GD_correct/final_results/result.mat',...
    '/volume/HCStress/Analysis/24-Sep-2018/DP_CISS_RSA_HC_ROD_CHR_ROP_631_GM_80PI_8GD_correct/final_results/result.mat',...
    '/volume/HCStress/Analysis/26-Oct-2018/DP_CISS_RSA_HC_ROD_CHR_ROP_634_GM_80PI_9GD_correct/final_results/result.mat',...
    '/volume/HCStress/Analysis/25-Sep-2018/DP_CISS_RSA_HC_ROD_CHR_ROP_631_GM_80PI_10GD_correct/final_results/result.mat',...
    '/volume/HCStress/Analysis/27-Sep-2018/DP_CISS_RSA_HC_ROD_CHR_ROP_631_GM_80PI_12GD_correct/final_results/result.mat',...
    '/volume/HCStress/Analysis/01-Oct-2018/DP_CISS_RSA_HC_ROD_CHR_ROP_631_GM_80PI_14GD_correct/final_results/result.mat'};

results_paths = {'/volume/HCStress/Analysis/12-Nov-2018/DP_CTQ_BS_HC_ROD_CHR_ROP_627_GM_80PI_10GD_correct/final_results/result.mat',...
    '/volume/HCStress/Analysis/12-Nov-2018/DP_CTQ_BS_HC_ROD_CHR_ROP_627_GM_80PI_10GD_correct/final_results/result_pruned_0_5onlysig_top10.mat',...
    '/volume/HCStress/Analysis/12-Nov-2018/DP_CTQ_BS_HC_ROD_CHR_ROP_627_GM_80PI_10GD_correct/final_results/result_top10.mat',...
    '/volume/HCStress/Analysis/09-Nov-2018/DP_CTQ_BS_HC_ROD_CHR_ROP_627_GM_80PI_20GD_correct/final_results/result.mat',...
    '/volume/HCStress/Analysis/09-Nov-2018/DP_CTQ_BS_HC_ROD_CHR_ROP_627_GM_80PI_20GD_correct/final_results/result_pruned_0_5onlysig.mat',...
    '/volume/HCStress/Analysis/09-Nov-2018/DP_CTQ_BS_HC_ROD_CHR_ROP_627_GM_80PI_20GD_correct/final_results/result_top10.mat'};

load(results_paths{i});
input.name=[input.name, '_new'];
save(results_paths{1}, 'input', 'output', 'setup');

% CTQ and BS combined
results_paths = {'/volume/HCStress/Analysis/14-Feb-2019/DP_CTQ_BS_allgroups_627_GM_80PI_1GO_10X10_20GD_10Y10_20GD_correct10_10_XY_SC_3_Diag/final_results/result.mat',...
    '/volume/HCStress/Analysis/25-Feb-2019/DP_CTQ_BS_allgroups_627_GM_80PI_1GO_10X15_20GD_10Y10_20GD_correct10_10_XY_SC_3_Diag/final_results/result.mat',...
    '/volume/HCStress/Analysis/26-Feb-2019/DP_CTQ_BS_allgroups_627_GM_80PI_12GO_520X_03020GD_110Y_01020GD_correct10_10_XY_SC_3_Diag/final_results/result.mat',...
    '/volume/HCStress/Analysis/26-Feb-2019/DP_CTQ_BS_allgroups_627_GM_80PI_1GO_10X20_20GD_10Y10_20GD_correct10_10_XY_SC_3_Diag/final_results/result.mat',...
    '/volume/HCStress/Analysis/27-Feb-2019/DP_CTQ_BS_allgroups_627_GM_80PI_12GO_5010X_06010GD_110Y_01010GD_correct10_10_XY_SC_3_Diag_1/final_results/result.mat',...
    '/volume/HCStress/Analysis/28-Feb-2019/DP_CTQ_BS_allgroups_627_GM_80PI_12GO_110X_06020GD_110Y_01020GD_correct10_5_XY_SC_3_Diag/final_results/result.mat',...
    '/volume/HCStress/Analysis/28-Feb-2019/DP_CTQ_BS_allgroups_627_GM_80PI_1GO_1X_020GD_1Y_020GD_correct10_5_XY_SC_3_Diag/final_results/result.mat',...
    '/volume/HCStress/Analysis/28-Feb-2019/DP_CTQ_BS_allgroups_627_GM_80PI_1GO_1X_020GD_1Y_020GD_correct5_5_XY_SC_3_Diag/final_results/result.mat',...
    '/volume/HCStress/Analysis/28-Feb-2019/DP_CTQ_BS_allgroups_627_GM_80PI_12GO_15X_06020GD_15Y_01020GD_correct5_10_XY_SC_3_Diag/final_results/result.mat',...
    '/volume/HCStress/Analysis/01-Mar-2019/DP_CTQ_BS_allgroups_627_GM_80PI_12GO_110X_07020GD_110Y_01020GD_correct10_5_XY_SC_3_Diag/final_results/result.mat',...
    '/volume/HCStress/Analysis/04-Mar-2019/DP_CTQ_BS_allgroups_627_GM_80PI_12GO_11X_07020GD_11Y_01020GD_correct10_5_XY_SC_3_Diag/final_results/result.mat'};
% only CTQ
results_paths = {'/volume/HCStress/Analysis/15-Feb-2019/DP_CTQ_allgroups_649_GM_80PI_1GO_1X0_20GD_1Y0_20GD_correct10_10_XY_SC_3_Diag_1/final_results/result.mat',...
    '/volume/HCStress/Analysis/17-Feb-2019/DP_CTQ_allgroups_649_GM_80PI_1GO_5X5_20GD_5Y5_20GD_correct10_10_XY_SC_3_Diag/final_results/result.mat',...
    '/volume/HCStress/Analysis/17-Feb-2019/DP_CTQ_allgroups_649_GM_80PI_2GO_5X5_20GD_5Y5_20GD_correct10_10_XY_SC_3_Diag/final_results/result.mat',...
    '/volume/HCStress/Analysis/17-Feb-2019/DP_CTQ_allgroups_649_GM_80PI_1GO_5X5_20GD_5Y5_20GD_correct10_10_XY_SC_3_Diag/final_results/result.mat'};

% single
results_paths = {'/volume/HCStress/Analysis/17-Feb-2020/DP_CTQ_allgroups_649_GM_80PI_1GO_1X0_20GD_1Y0_20GD_correct5_5_XY_SC_3_Diag/final_results/result.mat'};

% try rmdir('/home/dpopovic/.mcrCache8.5', 's')
% catch ME
% end

IN.specific = {'correlation', 'atlas', 'behavior', 'images', 'sociodemographic'}; % 'correlation', 'atlas', 'behavior', 'images', 'sociodemographic', 'outcome' 
IN.SD_selection = {'BDI', 'GAF', 'NEO', 'QOL'}; %'BDI', 'GAF', 'NEO', 'QOL'
IN.type = 1; % 1:new, 2:old

for i=1:numel(results_paths)
    IN.results_path = results_paths{i};
    if any(strfind(IN.results_path, 'CTQ'))
        IN.overall_analysis = 'Stress';
    elseif any(strfind(IN.results_path, 'CISS'))
        IN.overall_analysis = 'Resilience';
    else
        disp('Something''s wrong!');
    end
    dp_visualize_data_multi_2020(IN);
%     dp_brain_regions(results_paths{i}, 'Stress');
end

% compute adjusted p values (FDR)
output = dp_fdr_posthoc_adjust(results_paths{1});

%% single groups 4GD
detailed_paths=[];
detailed_paths.GD_4 = {'/volume/HCStress/Analysis/05-Oct-2018/DP_CTQ_BS_HC_251_GM_80PI_4GD_correct_2/final_results/result.mat',...
    '/volume/HCStress/Analysis/05-Oct-2018/DP_CTQ_BS_ROD_123_GM_80PI_4GD_correct/final_results/result.mat',...
    '/volume/HCStress/Analysis/05-Oct-2018/DP_CTQ_BS_CHR_115_GM_80PI_4GD_correct/final_results/result.mat',...
    '/volume/HCStress/Analysis/05-Oct-2018/DP_CTQ_BS_ROP_124_GM_80PI_4GD_correct/final_results/result.mat',...
    '/volume/HCStress/Analysis/06-Oct-2018/DP_CTQ_BS_CHR_ROD_ROP_373_GM_80PI_4GD_correct/final_results/result.mat'};

%% single groups 2GD
detailed_paths.GD_2 = {'/volume/HCStress/Analysis/07-Oct-2018/DP_CTQ_BS_HC_251_GM_80PI_2GD_correct/final_results/result.mat',...
    '/volume/HCStress/Analysis/07-Oct-2018/DP_CTQ_BS_ROD_123_GM_80PI_2GD_correct/final_results/result.mat',...
    '/volume/HCStress/Analysis/07-Oct-2018/DP_CTQ_BS_CHR_115_GM_80PI_2GD_correct/final_results/result.mat',...
    '/volume/HCStress/Analysis/07-Oct-2018/DP_CTQ_BS_ROP_124_GM_80PI_2GD_correct/final_results/result.mat'};

%% single groups 6GD
detailed_paths.GD_6 = {'/volume/HCStress/Analysis/07-Oct-2018/DP_CTQ_BS_HC_251_GM_80PI_6GD_correct/final_results/result.mat',...
    '/volume/HCStress/Analysis/07-Oct-2018/DP_CTQ_BS_ROD_123_GM_80PI_6GD_correct/final_results/result.mat',...
    '/volume/HCStress/Analysis/07-Oct-2018/DP_CTQ_BS_CHR_115_GM_80PI_6GD_correct/final_results/result.mat',...
    '/volume/HCStress/Analysis/08-Oct-2018/DP_CTQ_BS_ROP_124_GM_80PI_6GD_correct/final_results/result.mat'};

%% single groups 8GD
detailed_paths.GD_8={'/volume/HCStress/Analysis/09-Oct-2018/DP_CTQ_BS_HC_251_GM_80PI_8GD_correct/final_results/result.mat',...
    '/volume/HCStress/Analysis/09-Oct-2018/DP_CTQ_BS_ROD_123_GM_80PI_8GD_correct/final_results/result.mat',...
    '/volume/HCStress/Analysis/12-Oct-2018/DP_CTQ_BS_CHR_115_GM_80PI_8GD_correct/final_results/result.mat',...
    '/volume/HCStress/Analysis/12-Oct-2018/DP_CTQ_BS_ROP_124_GM_80PI_8GD_correct/final_results/result.mat'};

%% single groups 10GD
detailed_paths.GD_10={'/volume/HCStress/Analysis/12-Oct-2018/DP_CTQ_BS_HC_251_GM_80PI_10GD_correct/final_results/result.mat',...
    '/volume/HCStress/Analysis/12-Oct-2018/DP_CTQ_BS_ROD_123_GM_80PI_10GD_correct/final_results/result.mat',...
    '/volume/HCStress/Analysis/13-Oct-2018/DP_CTQ_BS_CHR_115_GM_80PI_10GD_correct/final_results/result.mat',...
    '/volume/HCStress/Analysis/13-Oct-2018/DP_CTQ_BS_ROP_124_GM_80PI_10GD_correct/final_results/result.mat'};

fields = fieldnames(detailed_paths);

for i=1:numel(fields)
    for ii=1:numel(detailed_paths.(fields{i}))
        dp_visualize_data(detailed_paths.(fields{i}){ii}, 'Stress');
    end
end

for i=1:size(results_paths,2)
    load(results_paths{i});
    results_collection.grid_density(i) = input.grid_density(1);
    detailed_results_folder = [results_paths{i}(1:(strfind(results_paths{i}, 'final')-1)), 'detailed_results'];
    for o=1:size(output.final_parameters,1)
        load([detailed_results_folder, '/opt_parameters_' num2str(o), '.mat']);
        if exist('opt_parameters', 'var')
            if size(opt_parameters,2) == 11
                index_RHO = opt_RHO;
                index_p = opt_p;
            elseif size(opt_parameters,2) == 8
                index_RHO = opt_param_RHO;
                index_p = opt_param_p;
            else
                disp('Something''s wrong');
            end
            results_collection.range_RHO.(['GD_',num2str(input.grid_density(1,1))])(:,o) = [opt_parameters{:,index_RHO}];
            results_collection.range_p.(['GD_',num2str(input.grid_density(1,1))])(:,o) = [opt_parameters{:,index_p}];
        elseif exist('opt_parameters_temp', 'var')
            if size(opt_parameters_temp,2) == 11
                index_RHO = opt_RHO;
                index_p = opt_p;
            elseif size(opt_parameters_temp,2) == 8
                index_RHO = opt_param_RHO;
                index_p = opt_param_p;
            else
                disp('Something''s wrong');
            end
            results_collection.range_RHO.(['GD_',num2str(input.grid_density(1,1))])(:,o) = [opt_parameters_temp{:,index_RHO}];
            results_collection.range_p.(['GD_',num2str(input.grid_density(1,1))])(:,o) = [opt_parameters_temp{:,index_p}];
        end
        clear opt_parameters opt_parameters_temp;
    end
    
    if size(output.final_parameters,2) == 11
        index_RHO = opt_RHO;
        index_p = opt_p;
        index_cu = opt_cu;
        index_cv = opt_cv;
    elseif size(output.final_parameters,2) == 8
        index_RHO = opt_param_RHO;
        index_p = opt_param_p;
        index_cu = opt_param_cu;
        index_cv = opt_param_cv;
    else
        disp('Something''s wrong');
    end
    
    if size(output.pvalue_FDR,2)>1
        output.pvalue_FDR = output.pvalue_FDR';
    end
        
    try log_sig = [output.final_parameters{:,index_p}]<=output.pvalue_FDR;
    catch
        log_sig = [output.final_parameters{:,index_p}]'<=output.pvalue_FDR;
    end
    
    results_collection.sig_LV(i) = sum(log_sig);
    
    temp_cu = cell2mat(output.final_parameters(:, index_cu));
    temp_cv = cell2mat(output.final_parameters(:, index_cv));
        
    results_collection.cu_mean(i) = mean(temp_cu(log_sig));
    results_collection.cv_mean(i) = mean(temp_cv(log_sig));
    
end

fields = fieldnames(results_collection);

% results_collection_temp.cu_mean = results_collection.cu_mean/sqrt(size(output.final_parameters{1,index_u},1));
% results_collection_temp.cv_mean = results_collection.cv_mean/sqrt(size(output.final_parameters{1,index_v},1));

for r=1:size(fields,1)
    temp = results_collection.(fields{r});
    if ~ strcmp(fields{r},'grid_density')
        try results_collection_temp.(fields{r}) = temp/(max(temp));
        catch
            results_collection_temp.(fields{r}) = temp;
        end
    else
        results_collection_temp.(fields{r}) = temp;
    end
end

%% plot results
collection_folder = '/volume/HCStress/Analysis/Resilience';
if ~exist(collection_folder)
    mkdir(collection_folder);
end
f = figure();
plot(results_collection.grid_density, results_collection.sig_LV, '--o');
name = strrep(input.name(1:(strfind(input.name, 'GD')-3)), '_', ' ');
title([name, ' grid density plot']); % add third line
axis([0 max(results_collection.grid_density)+2 0 max(results_collection.sig_LV)+2]);
set(gca, 'XTick', [0 : 20])%1 : results_collection.grid_density(end)]);
set(gca, 'YTick', [0 : 20]);% : results_collection.sig_LV(end)]);
xlabel('grid density');
ylabel('number of significant LV');
set(gcf, 'Position', get(0, 'Screensize'));
set(gcf,'PaperPositionMode','auto')
print(f, [collection_folder, '/sig_LV_GD'], '-dpng', '-r0');
saveas(f, [collection_folder, '/sig_LV_GD.fig']);

% range p
f = figure();
fields = fieldnames(results_collection.range_p);
for i=1:size(fields,1)
    subplot(round(size(fields,1)/2), 2, i);
    boxplot(results_collection.range_p.(fields{i}));
    xlabel('LV No.');
    ylabel('p value');
    title(strrep(fields{i},'_',' '));
end
set(gcf, 'Position', get(0, 'Screensize'));
set(gcf,'PaperPositionMode','auto')
print([collection_folder, '/range_p_GD'], '-dpng', '-r0');
saveas(f, [collection_folder, '/range_p_GD.fig']);

% range RHO
f = figure();
fields = fieldnames(results_collection.range_RHO);
for i=1:size(fields,1)
    subplot(round(size(fields,1)/2), 2, i);
    boxplot(results_collection.range_RHO.(fields{i}));
    xlabel('LV No.');
    ylabel('RHO value');
    title(strrep(fields{i},'_',' '));
end
set(gcf, 'Position', get(0, 'Screensize'));
set(gcf,'PaperPositionMode','auto')
print([collection_folder, '/range_RHO_GD'], '-dpng', '-r0');
saveas(f, [collection_folder, '/range_RHO_GD.fig']);

% plotting GD, sig LV, cu and cv means as a decision function
f = figure();
temp = results_collection_temp;
x = temp.grid_density;
y1 = temp.sig_LV;
y2 = temp.cu_mean;
y3 = temp.cv_mean;
plot(x,y1, 'b--o', x,y2, 'g--o', x,y3, 'r--o');
legend('significant LVs','cu mean','cv mean');
set(gcf, 'Position', get(0, 'Screensize'));
set(gcf,'PaperPositionMode','auto')
print(f, [collection_folder, '/sig_LV_cu_cv'], '-dpng', '-r0');
saveas(f, [collection_folder, '/sig_LV_cu_cv.fig']);

close all;
