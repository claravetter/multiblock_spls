%% function for RHO_avg BO

function [RHO_avg] = dp_RHO_avg_serial(K, tp, keep_in_data_x, keep_in_data_y, cu, cv) 

RHO_collection = nan(K,1);

for i=1:K
    [test_data_x, training_data_x, test_data_y, training_data_y] = dp_partition_holdout(tp, keep_in_data_x, keep_in_data_y);

    %perform SPLS on the training data using the current cu/cv combination
    [u,v,~] = spls_suppressed_display(training_data_x,training_data_y,cu,cv); 

    %compute the correlation between the projections of the
    %training and test matrices onto the SPLS latent space
    %spanned by the weight vectors
    RHO_collection(i,1) = abs(corr(test_data_x*u,test_data_y*v)); 
end

RHO_avg = mean(RHO_collection);

end