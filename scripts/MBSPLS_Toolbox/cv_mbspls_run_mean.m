function cv_mbspls_run_mean(datafile)
any_sig = true; 



while any_sig 

   [lv_filepath, dir] = cv_mbspls_main_1LV(datafile);

    mean_lv_filepath = cv_mbspls_outer_mean_average(lv_filepath, dir);

    deflated_datafile = cv_mbspls_outer_mean_deflate(mean_lv_filepath, dir);

    load(deflated_datafile);
    temp = output.opt_parameters.(['LV_' num2str(input.n_LV)]);
    log_sig = [temp{:,matches(output.opt_parameters_names, 'p')}] <= output.pvalue_FDR(input.n_LV,1);
    if ~any(log_sig)
        any_sig = false; 
        disp(['checkpoint: LV' input.n_LV ' not significant. Finished.'])
    else
    datafile = deflated_datafile; 
    end
    disp(['checkpoint: LV' input.n_LV ' done'])
end
end
