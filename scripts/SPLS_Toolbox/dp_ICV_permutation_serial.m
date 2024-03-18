%% function for permutation testing

function [RHO_collection_ICV, u_collection_ICV, v_collection_ICV] = dp_ICV_permutation_serial(X, Y, COV, PARAM, DIAG)

switch PARAM.train
    case 1
        IN_x.train = X.train;
        IN_x.test = X.test;
        IN_y.train = Y.train;
        IN_y.test = Y.test;
        COV.train = COV.train;
        COV.test = COV.test;
        permmat = nk_PermInd2(PARAM.B, DIAG.train);
end

% 1) retrain on single test splits within
% folds, then merge, 2) retrain on all inner folds
% separately, then merge with mean or median, 3)
% retrain on entirety of inner folds, 4) use already
% existing u and v from inner folds without retraining

RHO_collection_ICV = nan(PARAM.B,1);
u_collection_ICV = nan(PARAM.B, size(IN_x.train,2));
v_collection_ICV = nan(PARAM.B, size(IN_y.train,2));

for ii=1:PARAM.B
    
    switch PARAM.train
        case 2
            find_log = ii==PARAM.perm_coding(3,:);
            if any(find_log)
                w = PARAM.perm_coding(2, find_log);
                ob=1;
                COV.test = COV(PARAM.cv.TestInd{ob,w},:);
                COV.train = COV(PARAM.cv.TrainInd{ob,w},:);
                
                IN_x.train = X(PARAM.cv.TrainInd{ob,w},:);
                IN_x.test = X(PARAM.cv.TestInd{ob,w},:);
                
                IN_y.train = Y(PARAM.cv.TrainInd{ob,w},:);
                IN_y.test = Y(PARAM.cv.TestInd{ob,w},:);
                
                DIAG.train = DIAG(PARAM.cv.TrainInd{ob,w},:);
                
                permmat = nk_PermInd2(PARAM.B, DIAG.train);
            end
    end
    
    % perform procrustean transformation to minimize rotation effects of
    % permutated y matrix, if V_opt available
    IN_y.train = IN_y.train(permmat(ii,:),:);
    
    OUT_x = dp_correctscale(IN_x,COV,PARAM.scale);
    OUT_y = dp_correctscale(IN_y,COV,PARAM.scale);
    
    if ~islogical(PARAM.V)
        [RHO_collection_ICV(ii,1), u_collection_ICV(ii,:), v_collection_ICV(ii,:), ~, ~, ~] = dp_spls_full(OUT_x.train,OUT_y.train,OUT_x.test, OUT_y.test, PARAM.cu, PARAM.cv, PARAM.correlate, PARAM.V);
    else
        [RHO_collection_ICV(ii,1), u_collection_ICV(ii,:), v_collection_ICV(ii,:), ~, ~, ~] = dp_spls_full(OUT_x.train,OUT_y.train,OUT_x.test, OUT_y.test, PARAM.cu, PARAM.cv, PARAM.correlate);
    end
    
end

end