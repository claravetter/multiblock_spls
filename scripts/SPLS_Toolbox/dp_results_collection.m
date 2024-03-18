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
results_paths = {'/volume/HCStress/Analysis/20-Sep-2018/DP_CTQ_BS_HC_ROD_CHR_ROP_626_GM_80PI_4GD_correct_1/final_results/result.mat',...
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

for i=1:numel(results_paths)
    dp_visualize_data(results_paths{i}, 'Resilience');
%     dp_brain_regions(results_paths{i}, 'Stress');
end

results_paths = {'/volume/HCStress/Analysis/21-Sep-2018/DP_CISS_RSA_HC_ROD_CHR_ROP_631_GM_80PI_4GD_correct/final_results/result.mat',...
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

% temp = input.behavior_names;
% 
% for i=1:size(results_paths,2)
%     load(results_paths{i});
%     input.behavior_names = temp;
% %     input.grid_density = input.grid_density-2;
% %     new_name = [input.name(1:(strfind(input.name, 'GD')-3)), num2str(input.grid_density), input.name(strfind(input.name, 'GD'):end)];
% %     input.name = new_name;
%     save(results_paths{i}, 'input', 'output', 'setup');
% end

% for i=1:size(results_paths,2)
%     load(results_paths{i});
%     detailed_results_folder = [results_paths{i}(1:(strfind(results_paths{i}, 'final')-1)), 'detailed_results'];
%     if ~isfield(output, 'pvalue_FDR')
%         output.pvalue_FDR =[];
%         FDRvalue=0.05;
%         for o=1:size(output.final_parameters,1)
%             load([detailed_results_folder, '/opt_parameters_' num2str(o), '.mat']);
%             if exist('opt_parameters', 'var')
%                 if size(opt_parameters,2) == 11
%                     index_RHO = opt_RHO;
%                     index_p = opt_p;
%                 elseif size(opt_parameters,2) == 8
%                     index_RHO = opt_param_RHO;
%                     index_p = opt_param_p;
%                 else
%                     disp('Something''s wrong');
%                 end
%                 output.pvalue_FDR(o,1) = dp_FDR([opt_parameters{:,index_p}], FDRvalue);
%             elseif exist('opt_parameters_temp', 'var')
%                 if size(opt_parameters_temp,2) == 11
%                     index_RHO = opt_RHO;
%                     index_p = opt_p;
%                 elseif size(opt_parameters_temp,2) == 8
%                     index_RHO = opt_param_RHO;
%                     index_p = opt_param_p;
%                 else
%                     disp('Something''s wrong');
%                 end
%                 output.pvalue_FDR(o,1) = dp_FDR([opt_parameters_temp{:,index_p}], FDRvalue);
%             end
%         end
%         save(results_paths{i}, 'input', 'output', 'setup');
%     end
% end

for i=1:size(results_paths,2)
    load(results_paths{i});
    results_collection.grid_density(i) = input.grid_density(1);
    detailed_results_folder = [results_paths{i}(1:(strfind(results_paths{i}, 'final')-1)), 'detailed_results'];
    for o=1:size(output.final_parameters,1)
%         temp_pos = output.final_parameters{o,4}(output.final_parameters{o,4}>0);
%         temp_neg = output.final_parameters{o,4}(output.final_parameters{o,4}<0);
%         prctile(temp_neg, 25);
%         prctile(temp_neg, 75);
%         results_collection.limits_MRI.(['GD_',num2str(input.grid_density(1,1))])(o,:) = [max(temp_pos), min(temp_pos), max(temp_neg), min(temp_neg)];

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
%             results_collection.range_RHO{i,o} = [min([opt_parameters{:,index_RHO}]), max([opt_parameters{:,index_RHO}])];
%             results_collection.range_p{i,o} = [min([opt_parameters{:,index_p}]), max([opt_parameters{:,index_p}])];
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

% results_collection.coeff=[];
for ff=1:size(results_paths,2)
%     coeff_temp=[];
    load(results_paths{ff})
    if exist('output.parameters_names', 'var')
        index_u = strcmp(output.parameters_names,'u');
        index_v = strcmp(output.parameters_names,'v');
        index_p = strcmp(output.parameters_names,'p');
        index_RHO = strcmp(output.parameters_names,'RHO');
        index_cu = strcmp(output.parameters_names,'cu');
        index_cv = strcmp(output.parameters_names,'cv');
    else
        output.parameters_names = {'w', 'cu', 'cv', 'u', 'v', 'U_opt', 'S_opt' ,'V_opt', 'success', 'RHO', 'p'}; % names of parameters for optimal cu/cv combination of the w-th loop
        index_p = strcmp(output.parameters_names, 'p');
        index_v = strcmp(output.parameters_names, 'v');
        index_u = strcmp(output.parameters_names, 'u');
        index_RHO = strcmp(output.parameters_names, 'RHO');  
        index_cu = strcmp(output.parameters_names, 'cu');
        index_cv = strcmp(output.parameters_names, 'cv');
    end
    
    if isfield(input, 'behavior_names')
        index_age = strcmp(input.behavior_names,'Age');
        index_sex = strcmp(input.behavior_names,'Sex'); 
    else
        disp('Something''s wrong!');
    end
% 
%     log_as=true;
%     for i=1:size(output.final_parameters,1)
%         v_temp = output.final_parameters{i, index_v};%(~index_age&~index_sex);
%         age_sex_coeff = [abs(v_temp(index_age)), abs(v_temp(index_sex))];
%         coeff_temp(i,1) = output.final_parameters{i, index_RHO}./(sum(v_temp~=0)/numel(v_temp));
%         if max(age_sex_coeff)>0.8
%             log_as(i,1) = false;
%         else
%             log_as(i,1) = true;
%         end
%     end
    
    try log_sig = [output.final_parameters{:,index_p}]<=output.pvalue_FDR;
    catch
        log_sig = [output.final_parameters{:,index_p}]'<=output.pvalue_FDR;
    end
    
    temp_cu = cell2mat(output.final_parameters(:, index_cu));
    temp_cv = cell2mat(output.final_parameters(:, index_cv));
    
%     try
        results_collection.cu_mean(ff) = mean(temp_cu(log_sig));
        results_collection.cv_mean(ff) = mean(temp_cv(log_sig));
%     catch 
%         results_collection.cu_mean(ff) = median(temp_cu(log_sig));
%         results_collection.cv_mean(ff) = mean(temp_cv(log_sig));
%     end
    
%     if isnan(results_collection.coeff_median(ff))
%         results_collection.coeff_median(ff) = 0;
%     end
    
%     if isnan(results_collection.coeff_mean(ff))
%         results_collection.coeff_mean(ff) = 0;
%     end
    
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

% coefficient
f = figure();
[hAx, hline1, hline2] = plotyy(results_collection.grid_density, results_collection.sig_LV, results_collection.grid_density, results_collection.coeff_median);
name = strrep(input.name(1:(strfind(input.name, 'GD')-3)), '_', ' ');
title([name, ' grid density plot']); % add third line
xlabel('grid density');
ylabel(hAx(1), 'number of significant LV');
ylabel(hAx(2), 'LV coefficient median');
set(hAx(1), 'ylim', [0 max(results_collection.sig_LV)+2]);
set(hAx(2), 'ylim', [0 round(max(results_collection.coeff_median))+2]);
set(hAx, 'xlim', [0 max(results_collection.grid_density)+2]);
set(hAx(1), 'YTick', [0 : 20]);%1 : results_collection.grid_density(end)]);
set(hAx(2), 'YTick', [0 : 20]);%1 : results_collection.grid_density(end)]);
set(hline1, 'LineStyle', '--', 'Marker', 'o', 'Color', 'b');
set(hline2, 'LineStyle', '--', 'Marker', 'o', 'Color', 'r');

set(gcf, 'Position', get(0, 'Screensize'));
set(gcf,'PaperPositionMode','auto')
print(f, [collection_folder, '/sig_LV_coeff_median'], '-dpng', '-r0');
saveas(f, [collection_folder, '/sig_LV_coeff_median.fig']);

close all;

% coefficient
f = figure();
[hAx, hline1, hline2] = plotyy(results_collection.grid_density, results_collection.sig_LV, results_collection.grid_density, results_collection.coeff_mean);
name = strrep(input.name(1:(strfind(input.name, 'GD')-3)), '_', ' ');
title([name, ' grid density plot']); % add third line
xlabel('grid density');
ylabel(hAx(1), 'number of significant LV');
ylabel(hAx(2), 'LV coefficient mean');
set(hAx(1), 'ylim', [0 max(results_collection.sig_LV)+2]);
set(hAx(2), 'ylim', [0 round(max(results_collection.coeff_mean))+2]);
set(hAx, 'xlim', [0 max(results_collection.grid_density)+2]);
set(hAx(1), 'YTick', [0 : 20]);%1 : results_collection.grid_density(end)]);
set(hAx(2), 'YTick', [0 : 20]);%1 : results_collection.grid_density(end)]);
set(hline1, 'LineStyle', '--', 'Marker', 'o', 'Color', 'b');
set(hline2, 'LineStyle', '--', 'Marker', 'o', 'Color', 'r');

set(gcf, 'Position', get(0, 'Screensize'));
set(gcf,'PaperPositionMode','auto')
print(f, [collection_folder, '/sig_LV_coeff_mean'], '-dpng', '-r0');
saveas(f, [collection_folder, '/sig_LV_coeff_mean.fig']);

close all;

%% plot all LV of GD

% for i=1:size(results_paths,2)
%     results_folder = results_paths{i}(1:(strfind(results_paths{i}, '/result.mat')-1));
%     h1 = openfig([results_folder, '/behavior_LV1.fig'],'reuse'); % open figure
%     ax1 = gca; % get handle to axes of figure
%     h2 = openfig('test2.fig','reuse');
%     ax2 = gca;
%     h3 = figure; %create new figure
%     s1 = subplot(2,1,1); %create and get handle to the subplot axes
%     s2 = subplot(2,1,2);
%     fig1 = get(ax1,'children'); %get handle to all the children in the figure
%     fig2 = get(ax2,'children');
%     copyobj(fig1,s1); %copy children to new parent axes i.e. the subplot axes
%     copyobj(fig2,s2);
% end
