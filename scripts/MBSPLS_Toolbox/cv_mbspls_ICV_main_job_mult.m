%% function for RHO parallel computing

function [RHO_ICV, weights_ICV] = cv_mbspls_ICV_main_job_mult(analysis_folder, type_analysis, total_jobs, setup, num_matrices) %...



switch setup.cluster
    case {'SLURM', 'SGE'}
        RHO_bash = cv_mbspls_ICV_job_mult(analysis_folder, type_analysis, total_jobs, setup);
    case 'sequential'
        switch type_analysis
            case 'hyperopt'
                cv_mbspls_ICV_hyperopt_csv(1, analysis_folder)
            case 'permutation'
                cv_mbspls_ICV_permutation_csv(1, analysis_folder)
            case 'bootstrap'
                cv_mbspls_ICV_bootstrap_csv(1, analysis_folder)
        end
end

cd(analysis_folder);

switch setup.cluster
    case 'SLURM'
        system(['sbatch ' RHO_bash]);
    case 'SGE'
        system(['qsub ' RHO_bash]);
end


% see hyperparameter optimization
disp(analysis_folder);
target = [analysis_folder '/RHO_results_*.csv'];
disp(target);

if ~strcmp(setup.cluster, 'sequential') 
    mydir = size(dir(target),1);
    while mydir<total_jobs
        mydir = size(dir(target),1);
    end
    pause(30);
end


switch type_analysis
    case 'bootstrap'

        RHO_collection = [];

    
        target = [analysis_folder '/weights_results_*.csv'];
        mydir = size(dir(target),1);
        weights_collection{num_matrices} = [];


        for i=1:total_jobs


            for j=1:num_matrices
                path_file = [analysis_folder '/weights_' num2str(j) '_results_' num2str(i),'.csv'];
                weights_collection_temp = readmatrix(path_file);
                weights_collection{j} = [weights_collection{j}, weights_collection_temp];

                weights_ICV = weights_collection;
                delete(path_file);
            end



        end


        for i=1:total_jobs
            path_file = [analysis_folder '/RHO_results_' num2str(i),'.csv'];
            RHO_collection_temp = load(path_file);
            RHO_collection = [RHO_collection, RHO_collection_temp];
            delete(path_file);
        end

        RHO_ICV = RHO_collection;
    otherwise
        RHO_collection = [];


        for i=1:total_jobs
            path_file = [analysis_folder '/RHO_results_' num2str(i),'.csv'];
            RHO_collection_temp = readmatrix(path_file);
            RHO_collection = [RHO_collection; RHO_collection_temp];
            delete(path_file);
        end

        RHO_ICV = RHO_collection;

end

end