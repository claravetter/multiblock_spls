    %% Projection of the original data onto the SPLS latent space
    
    % output.epsilon and output.omega are the projections of the original data onto the
    % SPLS latent space, spanned by the weight vector pair u and v of this
    % LV. Therefore each LV spans a specific latent space onto wihich the
    % samples can be projected. These latent spaces can then be used to
    % stratify patients.
%     output.epsilon(ff,:) = X_original * output.final_parameters{ff,opt_u};
%     output.omega(ff,:) = Y_original * output.final_parameters{ff,opt_v};
    
    
%     output.epsilon_stand(ff,:) = zscore(X_original) * output.final_parameters{ff,opt_u};
%     output.omega_stand(ff,:) = zscore(Y_original) * output.final_parameters{ff,opt_v};
    
    
%     %% compute CI for item weights using bootstrapping of 10 opt_parameter iterations
%     
%     % first find inverted weight vectors in turn them around in order to
%     % avoid problems with mean computation during bootstrapping
%     inv=struct; inv.to_extract = 'v';
%     [opt_datamatrix_temp, inv] = dp_inverse_array(opt_parameters, inv, output.parameters_names);
%     inv.to_extract = 'u';
%     [opt_datamatrix, ~] = dp_inverse_array(opt_datamatrix_temp, inv, output.parameters_names);
%     
%     disp('checkpoint inversion');
%     
%     % set features to zero, which are also zero in the final LV
%     temp = opt_datamatrix(:,(opt_u + opt_v >0));
%     log_zero = {output.final_parameters{ff,opt_u}==0, output.final_parameters{ff,opt_v}==0};
%     for i=1:size(temp,1)
%         temp{i,1}(log_zero{1}) = 0;
%         temp{i,2}(log_zero{2}) = 0;
%     end
%     
%     opt_datamatrix(:,(opt_u + opt_v >0))=temp;
%     
%     disp('checkpoint zero');
%     
%     % vector u
%     temp=opt_datamatrix(:,opt_u); temp_u=[]; ci_collection_temp=[]; bootstat_collection_temp=[];
%     for i=1:size(temp,1)
%         for ii=1:size(opt_datamatrix{1,opt_u},1)
%             temp_u(i,ii) = opt_datamatrix{i,opt_u}(ii);
%         end
%     end
%     
%     mem_total           = 80;   % max: 40
%     max_sim_jobs        = 40;   % max: 60
%     ci_collection_temp = dp_bootstrap_opt_para(setup.spls_standalone_path, bootstrap_folder, 'bootstrap', mem_total, max_sim_jobs, setup.queue_name_slave, temp_u, input.bootstrap);
%     output.CI_u{ff} =[abs(ci_collection_temp(:,1)-output.final_parameters{ff,opt_u}),abs(ci_collection_temp(:,2)-output.final_parameters{ff,opt_u})];
%     
%     disp('checkpoint bootstrap u');
%     
%     % vector v
%     temp=opt_datamatrix(:,opt_v); temp_v=[]; ci_collection_temp=[]; bootstat_collection_temp=[];
%     for i=1:size(temp,1)
%         for ii=1:size(opt_datamatrix{1,opt_v},1)
%             temp_v(i,ii) = opt_datamatrix{i,opt_v}(ii);
%         end
%     end
%     
%     mem_total           = 80;   % max: 40
%     max_sim_jobs        = 40;   % max: 60
%     ci_collection_temp = dp_bootstrap_opt_para(setup.spls_standalone_path, bootstrap_folder, 'bootstrap', mem_total, max_sim_jobs, setup.queue_name_slave, temp_v, input.bootstrap);
%     output.CI_v{ff} =[abs(ci_collection_temp(:,1)-output.final_parameters{ff,opt_v}),abs(ci_collection_temp(:,2)-output.final_parameters{ff,opt_v})];
%     
%     disp('checkpoint bootstrap v');