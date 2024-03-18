%% additional FDR testing

function output = dp_fdr_posthoc_adjust(results_path)

load(results_path);

fields1 = fieldnames(output.hold_out_correlations);
output.adj_tables_holdout_RHO_p = output.tables_hold_out_Rho_p;
output.adj_FDR_values = output.hold_out_FDR_values;
for i=1:size(fields1,1)
    fields2 = fieldnames(output.hold_out_correlations.(['LV_', num2str(i)]));
    for ii=1:size(fields2,1)
        fields3 = fieldnames(output.hold_out_correlations.(['LV_', num2str(i)]).(fields2{ii}));
        for iii=1:size(fields3,1)
            fields4 = fieldnames(output.hold_out_correlations.(['LV_', num2str(i)]).(fields2{ii}).(fields3{iii}));
            temp_p_list.(fields1{i}).(fields2{ii}).(fields3{iii}) = [];
            for iiii=1:size(fields4,1)
                temp_p = output.hold_out_correlations.(['LV_', num2str(i)]).(fields2{ii}).(fields3{iii}).(fields4{iiii})(1,2);
                temp_p_list.(fields1{i}).(fields2{ii}).(fields3{iii}) = [temp_p_list.(fields1{i}).(fields2{ii}).(fields3{iii}); temp_p];
            end
            [output.adj_FDR_values.(fields1{i}).(fields2{ii}).(fields3{iii}),~,output.adj_tables_holdout_RHO_p.(fields1{i}).(fields2{ii})(2:2:end,iii)] = fdr(temp_p_list.(fields1{i}).(fields2{ii}).(fields3{iii}));
        end
        
    end
    
end

save(results_path, 'setup', 'input', 'output');

if any(strfind(results_path, 'CTQ'))
    overall_folder = '/volume/HCStress/Analysis/Stress';
elseif any(strfind(results_path, 'CISS'))
    overall_folder = '/volume/HCStress/Analysis/Resilience';
end

folder_name = results_path(5+strfind(results_path, '2019'):(strfind(results_path, '/final_results')-1));
collection_folder = [overall_folder, '/', folder_name];
corr_folder = [collection_folder, '/correlations'];
mkdir(corr_folder);


fields_adj=fieldnames(output.adj_tables_holdout_RHO_p);
for i=1:size(fields_adj,1)
    latent_scores = fieldnames(output.adj_tables_holdout_RHO_p.(fields_adj{i}));
    for ii=1:size(latent_scores,1)
        %             dp_txt_write(s_folder, ['/corr_RHO_', fields{i}, '_', latent_scores{ii}], output.tables_ho_RHO.(fields{i}).(latent_scores{ii})', '%.3f \t %.3f \t %.3f \t %.3f \t %.3f \n');
        %             dp_txt_write(s_folder, ['/corr_p_', fields{i}, '_', latent_scores{ii}], output.tables_ho_p.(fields{i}).(latent_scores{ii})', '%.3f \t %.3f \t %.3f \t %.3f \t %.3f \n');
        dp_txt_write(corr_folder, ['adj_corr_all_', fields_adj{i}, '_', latent_scores{ii}], output.adj_tables_holdout_RHO_p.(fields_adj{i}).(latent_scores{ii})', '%.3f \t %.3f \t %.3f \t %.3f \t %.3f \n');
    end
end


end