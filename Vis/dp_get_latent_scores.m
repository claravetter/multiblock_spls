%% DP function to get latent scores from SPLS results
function latent_scores_table = dp_get_latent_scores(results_file, group_select, correct)
load(results_file);
% group_select = 'all'; % all OR fold
% correct = 1; % 1=yes, 2=no
latent_scores_all =[];
latent_scores_names = [];
RHO_all=[];

switch group_select
    case 'all'
        try temp = load(input.X);
            name_temp = fieldnames(temp);
            X = temp.(name_temp{1});
        catch
            X = input.X;
        end    
        
        Y = input.Y;
        
        for i=1:(size(output.final_parameters,1)-1)
            
            IN_x.train = X;
            IN_x.test = X;
            IN_y.train = Y;
            IN_y.test = Y;
            
            if correct

                % with covariate correction
                if i==1
                    COV.train = input.covariates;
                    COV.test = input.covariates;
                else
                    COV.train = nan(size(input.covariates,1),1);
                    COV.test = nan(size(input.covariates,1),1);
                    input.correction_target = 3;
                end
            else % without covariate correction
                COV.train = nan(size(input.X,1),1);
                COV.test = nan(size(input.X,1),1);
                input.correction_target = 3;
            end
            
            if ~isempty(input.cs_method.correction_subgroup)
                try labels_temp = input.data_complete.foranalysis.basic{input.Y_final.Properties.RowNames, 'Labels'};
                catch
                    labels_temp = input.DiagNames;
                end
                cs_method.correction_subgroup = input.cs_method.correction_subgroup;
                cs_method.method = input.cs_method.method;
                cs_method.subgroup_train = contains(labels_temp, cs_method.correction_subgroup);
                cs_method.subgroup_test = contains(labels_temp, cs_method.correction_subgroup);
            else
                cs_method.correction_subgroup = [];
                cs_method.method = 'mean-centering';
                cs_method.subgroup_train = [];
                cs_method.subgroup_test = [];
            end
            
            [OUT_x, OUT_y] = dp_master_correctscale(IN_x, IN_y, COV, cs_method, input.correction_target);
                       
            log_u = matches(output.parameters_names, 'u');
            log_v = matches(output.parameters_names, 'v');
            u = output.final_parameters{i,log_u};
            v = output.final_parameters{i,log_v};
            
            [RHO, epsilon, omega, u, v] = dp_projection(OUT_x.train, OUT_y.train, u, v, input.correlation_method);
            
            [X,Y] = proj_def(X, Y, u, v);
            
            RHO_all = [RHO_all, RHO];
            latent_scores_all = [latent_scores_all, omega, epsilon];
            latent_scores_names = [latent_scores_names, {['omega_LV', num2str(i)]}, {['epsilon_LV', num2str(i)]}];
            
        end
        
        try latent_scores_table = array2table(latent_scores_all, 'RowNames', input.final_PSN, 'VariableNames', latent_scores_names);
        catch
            latent_scores_table = array2table(latent_scores_all, 'RowNames', input.final_ID, 'VariableNames', latent_scores_names);
        end
    case 'fold'
        
        latent_scores_table = array2table([omega_all{3}, epsilon_all{3}, omega_all{4}, epsilon_all{4}], 'RowNames', input.final_PSN, 'VariableNames', {'omega_LV3', 'epsilon_LV3', 'omega_LV4', 'epsilon_LV4'});
        
end