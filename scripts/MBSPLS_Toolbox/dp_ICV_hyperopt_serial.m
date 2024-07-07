%% new function for permutation testing

function [RHO_collection_ICV, u_collection_ICV, v_collection_ICV] = dp_ICV_hyperopt_serial(X, Y, Covariates, B, K, TrainInd, TestInd, cu_cv_combination, scaling_method, correlation_method)

u_collection_ICV = cell(size(cu_cv_combination,1),1);
v_collection_ICV = cell(size(cu_cv_combination,1),1);
RHO_collection_ICV = nan(size(cu_cv_combination,1), B*K);

for i=1:size(cu_cv_combination,1)    
    cu = cu_cv_combination(i,1);
    cv = cu_cv_combination(i,2);
    nn=1;
    for b=1:B
        for k=1:K
            
            IN_x.train = X(TrainInd{b,k},:);
            IN_x.test = X(TestInd{b,k},:);
            
            IN_y.train = Y(TrainInd{b,k},:);
            IN_y.test = Y(TestInd{b,k},:);

            COV.train = Covariates(TrainInd{b,k},:);
            COV.test = Covariates(TestInd{b,k},:);
            
            OUT_x = dp_correctscale(IN_x,COV,scaling_method);
            OUT_y = dp_correctscale(IN_y,COV,scaling_method);
            
            [RHO_collection_ICV(i,nn), u_collection_ICV{i,1}(:,nn), v_collection_ICV{i,1}(:,nn), ~, ~, ~] = dp_spls_full(OUT_x.train,OUT_y.train,OUT_x.test,OUT_y.test, cu, cv, correlation_method);
            nn=nn+1;
        end
    end
end

end