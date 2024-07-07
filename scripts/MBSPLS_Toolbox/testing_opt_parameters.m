%% testing of opt_parameters
addpath(genpath('/opt/NM/NeuroMiner_Release/'));
addpath('/opt/SPM/spm12_v6685_cat12_r1207/');
addpath(genpath('/volume/DP_FEF/ScrFun/ScriptsRepository'));

%% visualize results
results_path = '/volume/HCStress/Analysis/04-Sep-2018/DP_CTQBS_HCRODCHRROP_625_GM_mask_80P/final_results/result.mat';
NM_structure = '/volume/HCStress/Analysis/01-Sep-2018/allgroups_625_NM.mat';
load(NM_structure);
load(results_path);
detailed_results_folder = [results_path(1:(strfind(results_path, 'final')-1)), 'detailed_results'];
name = data_collection_path(strfind(data_collection_path, 'DP'):(strfind(data_collection_path, 'data_collection')-2));
cd(detailed_results_folder);
FDRvalue = 0.05;

for i=1:size(final_parameters,1)
    load([detailed_results_folder, '/opt_parameters_' num2str(i), '.mat']);
    temp = [opt_parameters{:,8}]';
    pvalue_FDR = dp_FDR(temp, FDRvalue);
    temp_sig = opt_parameters;
    opt_dir = [detailed_results_folder, '/opt_parameters_' num2str(i)];
    mkdir(opt_dir);
    cd(opt_dir);
    for ii=1:size(temp_sig,1)
        %% visualize results
        % write brain vector to nifti file
        nk_WriteVol(temp_sig{ii,4}, ['brain_opt_parameters_' num2str(i), '_' num2str(temp_sig{ii,1})], 2, NM.brainmask{1,1}, cell2mat(NM.badcoords), 0, 'gt');
    end
    
    %% visualize behavior vector in barplot
    % CTQ.emotional_abuse)), (CTQ.physical_abuse)), (CTQ.sexual_abuse)), (CTQ.emotional_neglect)), (CTQ.physical_neglect)), (CTQ.denial))];
    selected_variables = [AllVarNames, 'Age', 'Sex'];
    % colorpattern = {'by', 'm', 'c', 'r', 'g', 'b', 'k', 'y', 'm', 'c', 'r', 'g', 'b', 'k'};
    colorpattern = hsv(sum([size(fieldnames(CTQ),1), size(fieldnames(BS),1), 2]));
    
    for ii=1:(size(temp_sig,1))
        x = temp_sig{ii,5};
        f=figure();
        nn=0;
        hold on
        temp_legend=[];
        for iii=1:size(selected_variables,2)
            switch selected_variables{iii}
                case 'CTQ'
                    fields = fieldnames(CTQ);
                    temp_all = 0;
                    for iiii=1:size(fields,1)
                        nn=nn+1;
                        temp_current=size(CTQ.(fields{iiii}),2);
                        bar((temp_all+1):(temp_all+temp_current),x((temp_all+1):(temp_all+temp_current)), 'FaceColor', colorpattern(nn,:));
                        temp_all = temp_all+temp_current;
                        temp_legend{nn} = [selected_variables{iii}, ' ', strrep(fields{iiii}, '_', ' ')];
                        hold on
                    end
                case 'BS'
                    fields = fieldnames(BS);
                    for iiii=1:size(fields,1)
                        nn=nn+1;
                        temp_current=BS.(fields{iiii});
                        bar((temp_all+1):(temp_all+size(temp_current,2)),x((temp_all+1):(temp_all+size(temp_current,2))), 'FaceColor', colorpattern(nn,:));
                        temp_all = temp_all+size(temp_current,2);
                        temp_legend{nn} = [selected_variables{iii}, ' ', strrep(fields{iiii}, '_', ' ')];
                    end
                case 'Age'
                    nn=nn+1;
                    temp_current=data_collection.foranalysis.Age_Sex(:,1);
                    bar((temp_all+1):(temp_all+size(temp_current,2)),x((temp_all+1):(temp_all+size(temp_current,2))), 'FaceColor', colorpattern(nn,:));
                    temp_all = temp_all+size(temp_current,2);
                    temp_legend{nn} = selected_variables{iii};
                    hold on
                case 'Sex'
                    nn=nn+1;
                    temp_current=data_collection.foranalysis.Age_Sex(:,2);
                    bar((temp_all+1):(temp_all+size(temp_current,2)),x((temp_all+1):(temp_all+size(temp_current,2))), 'FaceColor', colorpattern(nn,:));
                    temp_all = temp_all+size(temp_current,2);
                    temp_legend{nn} = selected_variables{iii};
                    hold on
            end
        end
        first_line = strrep([name, ' grid_density=' num2str(RP), ' opt_parameters-',num2str(i), '-' num2str(temp_sig{ii,1})], '_', ' ');
        second_line = ['p-value = ' num2str(temp_sig{ii,8})];
        if temp_sig{ii,8}<=pvalue_FDR
            third_line = ['significant'];
        else
            third_line = ['not significant'];
        end
        title({first_line; second_line; third_line});
        xlabel('weight vector v');
        ylabel('score');
        axis([0 (size(temp_sig{ii,5},1)+1) -1 1]);
        legend(temp_legend);
        saveas(f, [opt_dir, '/behavior_opt_parameters_',num2str(i), '_' num2str(temp_sig{ii,1}), '.fig']);
        saveas(f, [opt_dir, '/behavior_opt_parameters_',num2str(i), '_' num2str(temp_sig{ii,1}), '.png']);
        close(f)
        %     set(get(subplot(size(temp_sig,1)-1,1,i), 'XLabel'), 'String', 'weight vector v');
        %     set(get(subplot(size(temp_sig,1)-1,1,i), 'YLabel'), 'String', 'score');
    end
    
%     %% plot latent scores epsilon and omega, color code diagnostic groups (HC,
%     % ROD, ROP, CHR)
%     
%     for ii=1:size(Studygroups_selected,2)
%         Studygroups_selected{2,ii} = strcmp(data_collection.foranalysis.Labels, Studygroups_selected{1,ii});
%     end
%     
%     colorpattern = hsv(size(Studygroups_selected,2));
%     for ii=1:size(temp_sig,1)
%         f=figure();
%         for iii=1:size(Studygroups_selected,2)
%             %     subplot(2,1,i)
%             %     for ii=1:size(sites_names,2)
%             %         temp = epsilon(i,:);
%             %         x = temp(logicals_sites{ii}');
%             x = epsilon(ii,Studygroups_selected{2,iii});
%             %         temp = omega(i,:);
%             y = omega(ii,Studygroups_selected{2,iii});
%             plot(x,y,'.', 'color',colorpattern(ii,:));
%             %     end
%             title(strrep([name, ' grid_density=' num2str(RP), ' LV ',num2str(ii)], '_', ' '));
%             xlabel('epsilon(brain score)');
%             ylabel('omega(behavior score)');
%             %         legend(Studygroups_selected{ii});
%             hold on
%             %         set(get(subplot(2,1,i), 'XLabel'), 'String', 'epsilon(brain score)');
%             %         set(get(subplot(2,1,i), 'YLabel'), 'String', 'omega(behavior score)');
%         end
%         legend(Studygroups_selected{1,:});
%         saveas(f, [results_folder, '/latent_scores_LV' num2str(ii), '.fig']);
%         saveas(f, [results_folder, '/latent_scores_LV' num2str(ii), '.png']);
%         close(f)
%     end
    
end

