%% DP function for bootstrap sampling

function boot_results_file = dp_bootstrap_pruning(results_file, log_application)

load(results_file);

output.old_final_parameters = output.final_parameters;

for i=1:size(output.final_parameters,1)
    lv_name = ['LV_', num2str(i)];
    u = output.final_parameters{i, matches(output.opt_parameters_names, 'u')};
    switch log_application
        case 'CI'
            log_ci_u = output.bootstrap_results.(lv_name).log_ci_u;
            u(log_ci_u) = 0;
        case 'BS'
            log_bs_u = output.bootstrap_results.(lv_name).log_bs_u;
            u(log_bs_u) = 0;
    end
    output.final_parameters{i, matches(output.opt_parameters_names, 'u')} = u;
    
    v = output.final_parameters{i, matches(output.opt_parameters_names, 'v')};
    switch log_application
        case 'CI'
            log_ci_v = output.bootstrap_results.(lv_name).log_ci_v;
            v(log_ci_v) = 0;
        case 'BS'
            log_bs_v = output.bootstrap_results.(lv_name).log_bs_v;
            v(log_bs_v) = 0;
    end
    output.final_parameters{i, matches(output.opt_parameters_names, 'v')} = v;
    
end

input.name = [input.name, 'boot_', log_application];
boot_results_file = strrep(results_file, '.mat', ['_', log_application, '.mat']);
save(boot_results_file, 'input', 'output', 'setup');

end
