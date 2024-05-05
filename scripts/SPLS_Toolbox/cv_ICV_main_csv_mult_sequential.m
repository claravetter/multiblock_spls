%% function for RHO parallel computing

function [RHO_ICV, weights_ICV] = cv_ICV_main_csv_mult_sequential(analysis_folder, type_analysis, total_jobs, cache_path, num_matrices, xy_spls_flag)

%RHO_bash = dp_ICV_bash_job_mult(spls_standalone_path, queue_name, analysis_folder, type_analysis, total_jobs, parallel_jobs, mem_request, matlab_version, compilation_subpath, cache_path);

cd(analysis_folder);
%system(['qsub ' RHO_bash]);

switch type_analysis
    case 'hyperopt'
        if xy_spls_flag
            dp_ICV_hyperopt_csv(num2str(1), analysis_folder)
        else
            cv_ICV_hyperopt_csv(1, analysis_folder)
        end
    case 'permutation'
        if xy_spls_flag
            dp_ICV_permutation_csv(num2str(1), analysis_folder)
        else
            cv_ICV_permutation_csv(1, analysis_folder)
        end
    case 'bootstrap'
        if xy_spls_flag
            dp_ICV_bootstrap_csv(num2str(1), analysis_folder)
        else
            cv_ICV_bootstrap_csv(1, analysis_folder)
        end
end


% see hyperparameter optimization
%target = [analysis_folder '/RHO_results_*.csv'];
%mydir = size(dir(target),1);

%while mydir<total_jobs
%    mydir = size(dir(target),1);
%end

%pause(30);

switch type_analysis
    case 'bootstrap'

        target = [analysis_folder '/weights_results_*.csv'];
        mydir = size(dir(target),1);
        RHO_collection = [];

        if xy_spls_flag
            u_collection = [];
            v_collection = [];
        else
            weights_collection{num_matrices} = [];
        end

        for i=1:total_jobs
            if xy_spls_flag
                path_file = [analysis_folder '/u_results_' num2str(i),'.csv'];
                u_collection_temp = readmatrix(path_file);
                u_collection = [u_collection, u_collection_temp];

                u_ICV = u_collection;

                path_file = [analysis_folder '/v_results_' num2str(i),'.csv'];
                v_collection_temp = readmatrix(path_file);
                v_collection = [v_collection, v_collection_temp];

                v_ICV = v_collection;
                delete(path_file);

                weights_ICV = {u_ICV, v_ICV};
            else

                for j=1:num_matrices
                    if isnumeric(i)
                        path_file = [analysis_folder '/weights_' num2str(j) '_results_' num2str(i),'.csv'];
                    else
                        path_file = [analysis_folder '/weights_' num2str(j) '_results_' i,'.csv'];
                    end
                    weights_collection_temp = readmatrix(path_file);
                    weights_collection{j} = [weights_collection{j}, weights_collection_temp];

                    weights_ICV = weights_collection;
                    delete(path_file);
                end

            end

        end


        for i=1:total_jobs
            if isnumeric(i)
                path_file = [analysis_folder '/RHO_results_' num2str(i),'.csv'];
            else
                path_file = [analysis_folder '/RHO_results_' i,'.csv'];
            end
            RHO_collection_temp = load(path_file);
            RHO_collection = [RHO_collection, RHO_collection_temp];
            delete(path_file);
        end

        RHO_ICV = RHO_collection;
    otherwise
        RHO_collection = [];


        for i=1:total_jobs
            if isnumeric(i)
                path_file = [analysis_folder '/RHO_results_' num2str(i),'.csv'];
            else
                path_file = [analysis_folder '/RHO_results_' i,'.csv'];
            end
            RHO_collection_temp = readmatrix(path_file);
            RHO_collection = [RHO_collection; RHO_collection_temp];
            delete(path_file);
        end

        RHO_ICV = RHO_collection;

end

end