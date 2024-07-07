%% DP function for bootstrap sampling

function OUT = dp_bootstrap_sampling(results_file)

load(results_file);

boot_size = inp;
[~,bootsam] = bootstrp(boot_size,[],input.final_PSN);

try temp = load(input.X);
    fieldnames_temp = fields(temp);
    X = temp.(fieldnames_temp{1});
catch
    X = input.X;
end

for i=1:size(bootsam,2)
    X_boot = X(bootsam(:,i),:);
    Y_boot = input.Y(bootsam(:,i),:);
    labels_boot = input.DiagNames(bootsam(:,i),:);
    
    for ii=1:size(output.final_parameters,1)
        if ii == 1
            covars_boot = input.covariates(bootsam(:,i),:);
        else
            covars_boot = nan(size(input.covariates,1),1);
        end
        input.cs_method.subgroup_train = matches(labels_boot, input.cs_method.correction_subgroup);
        input.cs_method.subgroup_test = matches(labels_boot, input.cs_method.correction_subgroup);
        [OUT_x, OUT_y] = dp_master_correctscale(X_boot, Y_boot, covars_boot, input.cs_method, input.correction_target);
        cu = output.final_parameters{ii, matches(output.opt_parameters_names, 'cu')};
        cv = output.final_parameters{ii, matches(output.opt_parameters_names, 'cv')};
        
        [RHO_boot.(['LV_', num2str(ii)])(1,i), u_boot.(['LV_', num2str(ii)])(:,i), v_boot.(['LV_', num2str(ii)])(:,i), ~, ~, ~] = dp_spls_full(OUT_x,OUT_y,OUT_x, OUT_y, cu, cv, input.correlation_method);
        %         [u_boot.(['LV_', num2str(ii)])(:,i), v_boot.(['LV_', num2str(ii)])(:,i), ~, ~, ~, ~] = dp_spls(OUT_x, OUT_y, cu, cv);
        
    end
    i
end

output.old_final_parameters = output.final_parameters;

log_application = 'BS';

for i=1:size(output.final_parameters,1)
    lv_name = ['LV_', num2str(i)];
%     RHO_boot.(lv_name) = output.bootsampling_results.(lv_name).RHO_sample;
    RHO = output.final_parameters{i, matches(output.opt_parameters_names, 'RHO')};
    RHO_sample = RHO_boot.(lv_name);
    RHO_mean = mean(RHO_sample);
    RHO_SE = std(RHO_sample)/(sqrt(boot_size));
    ci_RHO = [RHO_mean - 1.96 * RHO_SE, RHO_mean + 1.96 * RHO_SE];
    bs_ratio_RHO = RHO/RHO_SE;
    output.bootsampling_results.(lv_name).ci_RHO = ci_RHO;
    output.bootsampling_results.(lv_name).bs_ratio_RHO = bs_ratio_RHO;
    output.bootsampling_results.(lv_name).RHO_sample = RHO_sample;
    
%     u_boot.(lv_name) = output.bootsampling_results.(lv_name).u_sample;
    u_sample = u_boot.(lv_name);
    u_analysis = output.final_parameters{i, matches(output.opt_parameters_names, 'u')};
    u_mean = mean(u_sample,2);
    u_SE = std(u_sample,0,2)/(sqrt(boot_size));
    bs_ratio_u = u_analysis./u_SE;
    bs_ratio_u(isnan(bs_ratio_u)) = 0;
    bs_ratio_u(bs_ratio_u == Inf) = 0;
    log_bs_u = abs(bs_ratio_u)<=2;
    ci_u = [u_mean - 1.96 * u_SE, u_mean + 1.96 * u_SE];
    log_ci_u = ((sum(ci_u>0, 2) == 2) + (sum(ci_u<0, 2) == 2)) == 0;
    u = output.final_parameters{i, matches(output.opt_parameters_names, 'u')};
    switch log_application
        case 'CI'
                u(log_ci_u) = 0;
        case 'BS'
                u(log_bs_u) = 0;
    end
    output.final_parameters{i, matches(output.opt_parameters_names, 'u')} = u;
    output.bootsampling_results.(lv_name).ci_u = ci_u;
    output.bootsampling_results.(lv_name).bs_ratio_u = bs_ratio_u;
    output.bootsampling_results.(lv_name).u_sample = u_sample;
    output.bootsampling_results.(lv_name).log_bs_u = log_bs_u;
    output.bootsampling_results.(lv_name).sum_bs_u = sum(log_bs_u);
    output.bootsampling_results.(lv_name).log_ci_u = log_ci_u;
    output.bootsampling_results.(lv_name).sum_ci_u = sum(log_ci_u);
    
%     v_boot.(lv_name) = output.bootsampling_results.(lv_name).v_sample;
    v_sample = v_boot.(lv_name);
    v_analysis = output.final_parameters{i, matches(output.opt_parameters_names, 'v')};
    v_mean = mean(v_sample,2);
    v_SE = std(v_sample,0,2)/(sqrt(boot_size));
    bs_ratio_v = v_analysis./v_SE;
    bs_ratio_v(isnan(bs_ratio_v)) = 0;
    bs_ratio_v(bs_ratio_v == Inf) = 0;
    log_bs_v = abs(bs_ratio_v)<=2;
    ci_v = [v_mean - 1.96 * v_SE, v_mean + 1.96 * v_SE];
    log_ci_v = ((sum(ci_v>0, 2) == 2) + (sum(ci_v<0, 2) == 2)) == 0;
    v = output.final_parameters{i, matches(output.opt_parameters_names, 'v')};
    switch log_application
        case 'CI'
            v(log_ci_v) = 0;
        case 'BS'
            v(log_bs_v) = 0;
    end
    output.final_parameters{i, matches(output.opt_parameters_names, 'v')} = v;
    output.bootsampling_results.(lv_name).ci_v = ci_v;
    output.bootsampling_results.(lv_name).bs_ratio_v = bs_ratio_v;
    output.bootsampling_results.(lv_name).v_sample = v_sample;
    output.bootsampling_results.(lv_name).log_bs_v = log_bs_v;
    output.bootsampling_results.(lv_name).sum_bs_v = sum(log_bs_v);
    output.bootsampling_results.(lv_name).log_ci_v = log_ci_v;
    output.bootsampling_results.(lv_name).sum_ci_v = sum(log_ci_v);
    
    % testing
    output.testing.(lv_name).u_comparison = sum(log_bs_u == log_ci_u)/size(log_bs_u,1);
    output.testing.(lv_name).v_comparison = sum(log_bs_v == log_ci_v)/size(log_bs_v,1);
end

input.name = [input.name, 'boot_', log_application];

save(['result_bootstrapping_', log_application, '.mat'], 'input', 'output', 'setup');


end
