%% rest of plotting script

test=[5 4 7];
errors=[1,2,3]
bar(test, 'Facecolor', 'g');
hold on
errorbar([1:3], test, errors, 'b.');

load(results_path);
load(input.NM_structure);
results_folder = results_path(1:(strfind(results_path, '/result.mat')-1));
cd(results_folder);

% write brain vector to nifti file
for i=1:size(output.final_parameters,1)
    nk_WriteVol(output.final_parameters{i,4}, ['brain_LV' num2str(i)], 2, NM.brainmask{1,1}, cell2mat(NM.badcoords), 0, 'gt');
end

%% visualize behavior vector in barplot
% CTQ.emotional_abuse)), (CTQ.physical_abuse)), (CTQ.sexual_abuse)), (CTQ.emotional_neglect)), (CTQ.physical_neglect)), (CTQ.denial))];
% load([setup.analysis_folder, '/' setup.date, '_', name, '_data_collection.mat']);
selected_variables = [input.selected_features(1,:), 'Age', 'Sex'; input.selected_features(2,:), 1, 1];

% colorpattern = {'by', 'm', 'c', 'r', 'g', 'b', 'k', 'y', 'm', 'c', 'r', 'g', 'b', 'k'};
count=0;
for i=1:size(input.selected_features(2,:),2)
    temp=size(fieldnames(input.selected_features{2,i}),1);
    count=count+temp;
end

colorpattern_LV = hsv((count + 2));
% fieldnames(input.selected_features{2,:})
for i=1:(size(output.final_parameters,1))
    x = output.final_parameters{i,5};
    f=figure();
    nn=0;
    hold on
    temp_legend=[]; temp_all = 0;
    for ii=1:size(selected_variables,2)
        switch class(selected_variables{2,ii})
            case 'struct'
                fields = fieldnames(input.selected_features{2,(strcmp(input.selected_features(1,:),input.selected_features{1,ii}))});
                for iii=1:size(fields,1)
                    nn=nn+1;
                    temp_current=size(selected_variables{2,ii}.(fields{iii}),2);
                    bar((temp_all+1):(temp_all+temp_current),x((temp_all+1):(temp_all+temp_current)), 'FaceColor', colorpattern_LV(nn,:));
                    temp_all = temp_all+temp_current;
                    temp_legend{nn} = [selected_variables{1,ii}, ' ', strrep(fields{iii}, '_', ' ')];
                    hold on
                end
                
            case 'double'
                nn=nn+1;
                temp_current=size(selected_variables{2,ii},2);
                bar((temp_all+1):(temp_all+size(temp_current,2)),x((temp_all+1):(temp_all+size(temp_current,2))), 'FaceColor', colorpattern_LV(nn,:));
                temp_all = temp_all+size(temp_current,2);
                temp_legend{nn} = selected_variables{1,ii};
                hold on
        end
    end
    title(strrep([input.name, ' grid_density=' num2str(input.grid_density), ' LV ',num2str(i)], '_', ' '));    
    xlabel('weight vector v');
    ylabel('score');
    axis([0 (size(output.final_parameters{i,5},1)+1) -1 1]);
    legend(temp_legend, 'Location', 'bestoutside');
    saveas(f, [results_folder, '/behavior_LV' num2str(i), '.fig']);
    saveas(f, [results_folder, '/behavior_LV' num2str(i), '.png']);
    close(f)
%     set(get(subplot(size(final_parameters,1)-1,1,i), 'XLabel'), 'String', 'weight vector v');
%     set(get(subplot(size(final_parameters,1)-1,1,i), 'YLabel'), 'String', 'score');
end

%% plot latent scores epsilon and omega, color code diagnostic groups (HC,
% ROD, ROP, CHR)

for i=1:size(input.selected_studygroups,2)
    input.selected_studygroups{2,i} = strcmp(input.data_collection.Labels, input.selected_studygroups{1,i});
end

colorpattern_LS = hsv(size(input.selected_studygroups,2));
for i=1:size(output.final_parameters,1)
    f=figure();
    for ii=1:size(input.selected_studygroups,2)
        %     subplot(2,1,i)
        %     for ii=1:size(sites_names,2)
%         temp = epsilon(i,:);
%         x = temp(logicals_sites{ii}');
        x = output.epsilon(i,input.selected_studygroups{2,ii});
%         temp = omega(i,:);
        y = output.omega(i,input.selected_studygroups{2,ii});
        plot(x,y,'.', 'color',colorpattern_LS(ii,:));
        %     end
        title(strrep([input.name, ' grid_density=' num2str(input.grid_density), ' LV ',num2str(i)], '_', ' '));    
        xlabel('epsilon(brain score)');
        ylabel('omega(behavior score)');
%         legend(Studygroups_selected{ii});
        hold on
%         set(get(subplot(2,1,i), 'XLabel'), 'String', 'epsilon(brain score)');
%         set(get(subplot(2,1,i), 'YLabel'), 'String', 'omega(behavior score)');
    end
    legend(input.selected_studygroups{1,:});
    saveas(f, [results_folder, '/latent_scores_LV' num2str(i), '.fig']);
    saveas(f, [results_folder, '/latent_scores_LV' num2str(i), '.png']);
    close(f)
end

detailed_results_folder = [results_path(1:(strfind(results_path, 'final')-1)), 'detailed_results'];
% name = data_collection_path(strfind(data_collection_path, 'DP'):(strfind(data_collection_path, 'data_collection')-2));
cd(detailed_results_folder);
FDRvalue = 0.05;

for i=1:size(output.final_parameters,1)
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
%     selected_variables = [AllVarNames, 'Age', 'Sex'];
    % colorpattern = {'by', 'm', 'c', 'r', 'g', 'b', 'k', 'y', 'm', 'c', 'r', 'g', 'b', 'k'};
%     colorpattern = hsv(sum([size(fieldnames(CTQ),1), size(fieldnames(BS),1), 2]));
    
    for ii=1:(size(temp_sig,1))
        x = temp_sig{ii,5};
        f=figure();
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
                        num2str(iiii);
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
        first_line = strrep([input.name, ' grid_density=' num2str(input.grid_density), ' opt_parameters-',num2str(i), '-' num2str(temp_sig{ii,1})], '_', ' ');
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
        legend(temp_legend, 'Location', 'bestoutside');
        saveas(f, [opt_dir, '/behavior_opt_parameters_',num2str(i), '_' num2str(temp_sig{ii,1}), '.fig']);
        saveas(f, [opt_dir, '/behavior_opt_parameters_',num2str(i), '_' num2str(temp_sig{ii,1}), '.png']);
        close(f)
        %     set(get(subplot(size(temp_sig,1)-1,1,i), 'XLabel'), 'String', 'weight vector v');
        %     set(get(subplot(size(temp_sig,1)-1,1,i), 'YLabel'), 'String', 'score');
    end
    
end



