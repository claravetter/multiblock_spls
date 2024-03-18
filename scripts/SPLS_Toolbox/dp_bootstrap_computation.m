%% DP function for bootstrapping mean and SD for detailed results

final_results_path = '/volume/HCStress/Analysis/09-Nov-2018/DP_CTQ_BS_HC_ROD_CHR_ROP_627_GM_80PI_20GD_correct/final_results/';
file = [final_results_path, 'result.mat'];
load(file);
IN.specific = {'atlas'};

for i=1:numel(results_paths)
    IN.results_path = file;
    if any(strfind(IN.results_path, 'CTQ'))
        IN.overall_analysis = 'Stress';
    elseif any(strfind(IN.results_path, 'CISS'))
        IN.overall_analysis = 'Resilience';
    else
        disp('Something''s wrong!');
    end
    dp_visualize_data(IN);
%     dp_brain_regions(results_paths{i}, 'Stress');
end

bootstrap_folder = '/volume/HCStress/Analysis/09-Nov-2018/DP_CTQ_BS_HC_ROD_CHR_ROP_627_GM_80PI_20GD_correct/bootstrap';
% collect all brain and behavior data across LVs
LV_dir = '/volume/HCStress/Analysis/12-Nov-2018/DP_CTQ_BS_HC_ROD_CHR_ROP_627_GM_80PI_10GD_correct/detailed_results/opt_parameters_';
LV_selected = 1:4;

for i=1:size(LV_selected,2)
    inv=struct;
    load([LV_dir, num2str(i), '.mat']);
    inv.to_extract = 'v';
    [~, inv] = dp_inverse_array(opt_parameters_temp, inv, opt_parameters_names);
    % take only significant iterations
    log_sig = cell2mat(opt_parameters_temp(:,strcmp(opt_parameters_names, 'p'))) <= output.pvalue_FDR(i);
    temp = opt_parameters_temp(log_sig, :);
    for ii=1:size(temp,1)
        if inv.log(ii)
            temp{ii,4} = temp{ii,4}*(-1);
            temp{ii,5} = temp{ii,5}*(-1);
        end
    end
    
    opt_parameters_temp1 = temp;
    
    for ii=1:size(opt_parameters_temp1,1)
        LV_collection.brain.(['LV_', num2str(i)])(ii,:) = opt_parameters_temp1{ii,4};
        LV_collection.behavior.(['LV_', num2str(i)])(ii,:) = opt_parameters_temp1{ii,5};
    end
    
end

% find most reliable brain and behavior features across opt_parameters
fields1 = fieldnames(LV_collection); % paramci
threshold = 0.5;
for i=1:size(fields1,1)
    fields2 = fieldnames(LV_collection.(fields1{i}));
    for ii=1:size(fields2,1)
        temp=LV_collection.(fields1{i}).(fields2{ii});
        for iii=1:size(temp,2)
            DC_collection.(fields1{i}).(fields2{ii}).log(iii) = (sum(temp(:,iii)~=0))>=(threshold*size(LV_collection.(fields1{i}).(fields2{ii}),1));
        end
%         DC_collection.(fields1{i}).(fields2{ii}).log = DC_collection.(fields1{i}).(fields2{ii}).log>0;
        switch fields1{i}
            case 'behavior'
                n=5;
            case 'brain'
                n=4;
        end
        output.final_parameters{ii,n}(~DC_collection.(fields1{i}).(fields2{ii}).log) = 0;
    end
end

input.name = [input.name, '_pruned_0_5onlysig'];
name= [final_results_path, 'result_pruned_0_5onlysig.mat'];
save(name, 'input', 'output', 'setup');

% apply logicals to output files and write out atlases
results_paths = {name};
IN.specific = {'atlas'};

for i=1:numel(results_paths)
    IN.results_path = results_paths{i};
    if any(strfind(IN.results_path, 'CTQ'))
        IN.overall_analysis = 'Stress';
    elseif any(strfind(IN.results_path, 'CISS'))
        IN.overall_analysis = 'Resilience';
    else
        disp('Something''s wrong!');
    end
    dp_visualize_data(IN);
%     dp_brain_regions(results_paths{i}, 'Stress');
end

% take top 10 regions if possible and apply them to mri images
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

for i=1:4
    
    try output.regions.top12.weights{i,1} = cell2mat(output.regions.count.weights{i,1}(1:12,1));
    catch
        try output.regions.top12.weights{i,1} = cell2mat(output.regions.count.weights{i,1}(1:end,1));
        catch
            output.regions.top12.weights{i,1} = output.regions.count.weights{i,1};
        end
    end
    
    try output.regions.top12.weights{i,2} = cell2mat(output.regions.count.weights{i,2}(1:12,1));
    catch
        try output.regions.top12.weights{i,2} = cell2mat(output.regions.count.weights{i,2}(1:end,1));
        catch
            output.regions.top12.weights{i,2} = output.regions.count.weights{i,2};
        end
    end
    
end

for i=1:4
    log = ismember(atlas_hammers_full_for_analysis, [output.regions.top12.weights{i,1};output.regions.top12.weights{i,2}]);
    output.final_parameters{i,4}(~log)=0;
end

input.name = [input.name, '_top12.mat'];
name= [final_results_path, 'result_top12.mat'];
save(name, 'input', 'output', 'setup');

% apply logicals to output files and write out atlases
results_paths = {name};
IN.specific = {'atlas', 'images', 'behavior'};

for i=1:numel(results_paths)
    IN.results_path = results_paths{i};
    if any(strfind(IN.results_path, 'CTQ'))
        IN.overall_analysis = 'Stress';
    elseif any(strfind(IN.results_path, 'CISS'))
        IN.overall_analysis = 'Resilience';
    else
        disp('Something''s wrong!');
    end
    dp_visualize_data(IN);
%     dp_brain_regions(results_paths{i}, 'Stress');
end


for i=1:size(fields1,1)
    fields2 = fieldnames(LV_collection.(fields1{i}));
    switch fields1{i}
        case 'brain'
            for ii=1:size(fields2,1)
                %         for iii=1:size(LV_collection.(fields1{i}).(fields2{ii}),2)
                temp = LV_collection.(fields1{i}).(fields2{ii});
                mem_total           = 80;   % max: 40
                max_sim_jobs        = 40;   % max: 60
                CI_collection.(fields1{i}).(fields2{ii}).log_unreliable = dp_bootstrap_opt_para(setup.spls_standalone_path, bootstrap_folder, 'bootstrap', mem_total, max_sim_jobs, setup.queue_name_slave, temp, input.bootstrap);
                
            end
        case 'behavior'
            for ii=1:size(fields2,1)
                temp = LV_collection.(fields1{i}).(fields2{ii});
                CI_collection.(fields1{i}).(fields2{ii}).log_unreliable = dp_bootci(temp);
            end
    end
end

save('/volume/HCStress/Analysis/09-Nov-2018/DP_CTQ_BS_HC_ROD_CHR_ROP_627_GM_80PI_20GD_correct/final_results/CI_collection_improved.mat', 'CI_collection');

output.final_parameters_prebts = output.final_parameters;
output.bts_final_parameters = output.final_parameters;

for i=1:size(output.final_parameters,1)
    temp_u = output.final_parameters{i,4};
    temp_u(CI_collection.brain.(['LV_', num2str(i)]).log_unreliable>0) = 0;
    output.bts_final_parameters{i,4} = temp_u;
    temp_v = output.final_parameters{i,5};
    temp_v(CI_collection.behavior.(['LV_', num2str(i)]).log_unreliable>0) = 0;
    output.bts_final_parameters{i,5} = temp_v;
end

output.final_parameters = output.bts_final_parameters;
input.name = [input.name, '_BTS'];

output.bts_final_parameters = [];

save('/volume/HCStress/Analysis/09-Nov-2018/DP_CTQ_BS_HC_ROD_CHR_ROP_627_GM_80PI_20GD_correct/final_results/result_bts.mat', 'input', 'output', 'setup');







