%% function for permutation testing

function cv_mbspls_mean_ICV_permutation_csv(i, analysis_folder)



% permute all matrices
for num_m=1:length(IN_matrices.train)
    for pp=1:str2double(i)
        permmat = nk_PermInd(m_setup.size_sets_permutation, m_data.train_Diag);
    end
    permmats{num_m} = permmat;
end
clear permmat

RHO_collection_ICV_mean = nan(m_setup.size_sets_permutation,1);

for ii=1:m_setup.size_sets_permutation
    RHO_collection_ICV = nan(W,1);
    for ob=1:OB
        for w=1:W

            for num_m=1:num_matrices
                test_covariates{num_m} = Covars_ana{num_m}(output.CV.cv_outer_indices.TestInd{ob,w},:);
                train_covariates{num_m} = Covars_ana{num_m}(output.CV.cv_outer_indices.TrainInd{ob,w},:);
            end


            wrapper_test = @(x) get_split_matrix(x, output.CV.cv_outer_indices.TestInd{ob,w}) ;
            wrapper_train = @(x) get_split_matrix(x, output.CV.cv_outer_indices.TrainInd{ob,w}) ;
            train_data_matrices = cellfun(wrapper_train, matrices_ana, UniformOutput=false);

            test_data_matrices = cellfun(wrapper_test, matrices_ana, UniformOutput=false);

            train_Diag = Diag_ana(output.CV.cv_outer_indices.TrainInd{ob,w},:);

            IN_matrices.train = train_data_matrices;
            IN_matrices.test = test_data_matrices;

            train_covariates = train_covariates;
            test_covariates = test_covariates;
            for num_m=1:size(IN_matrices.train,2)
                COV{num_m}.train = train_covariates{num_m};
                COV{num_m}.test = test_covariates{num_m};
            end

            cs_method_permutation = input.cs_method;

            for num_m=1:length(IN_matrices.train)-1
                IN_matrices.train{num_m+1} = train_data_matrices{num_m+1}(permmats{num_m}(ii,:),:);
            end


            [OUT_matrices] = cv_master_correctscale(IN_matrices, COV, cs_method_permutation, input.correction_target);

            if ~islogical(m_opt.Vs_opt) && m_setup.procrustes_flag
                RHO_collection_ICV(w,1) = cv_mbspls_wrapper(OUT_matrices.train, ...
                    OUT_matrices.test, ...
                    m_opt.c_weights_opt, ...
                    m_setup.mbspls_params, ...
                    m_setup.correlation_method, ...
                    m_setup.matrix_norm, ...
                    m_opt.Vs_opt); % slim
            else
                RHO_collection_ICV(w,1) = cv_mbspls_wrapper(OUT_matrices.train, ...
                    OUT_matrices.test, ...
                    m_opt.c_weights_opt, ...
                    m_setup.mbspls_params, ...
                    m_setup.correlation_method, ...
                    m_setup.matrix_norm, ...
                    []); % slim
            end

        end
    end
    RHO_collection_ICV_mean(ii,1) = mean(RHO_collection_ICV);
end

RHO_collection = [];


for i=1:total_jobs
    path_file = [analysis_folder '/RHO_results_' num2str(i),'.csv'];
    RHO_collection_temp = readmatrix(path_file);
    RHO_collection = [RHO_collection; RHO_collection_temp];
    delete(path_file);
end

RHO_ICV = RHO_collection;


end