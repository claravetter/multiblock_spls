%% new function for cu/cv combination and 100 splits

function [RHO_avg] = dp_RHO_avg_k_AR(i)

load('/volume/DP_FEF/Analysis/28-Mar-2018/cvcu/cvcu.mat')
%% start of the function
cu_cv_combination   = cu_cv_combination(:,i);

RHO_collection = zeros(K,1);

for k=1:K
            
    % separate the keep_in data into training and test data according to
    % the chosen test percentage tp
    [test_data_x, training_data_x, test_data_y, training_data_y] = dp_partition_holdout(tp, keep_in_data_x, keep_in_data_y);
        
    %perform SPLS on the training data using the current cu/cv combination
    [u, v, ~] = spls(training_data_x,training_data_y,cu_cv_combination(1),cu_cv_combination(2)); 
                    
    %compute the correlation between the projections of the training and
    %test matrices onto the SPLS latent space spanned by the weight vectors
    RHO = abs(corr(test_data_x*u,test_data_y*v)); 
                
    %store the results  
    RHO_collection(k) = RHO;
                
end


% computing the mean value of all RHO calculated for the k iterations for
% one specific cu/cv combination
            
RHO_avg = mean(RHO_collection); 
% M  = matfile('/volume/DP_FEF/Analysis/28-Mar-2018/cvcu/RHO_avg.mat');
disp(num2str(RHO_avg))
save(['/volume/DP_FEF/Analysis/28-Mar-2018/cvcu/RHO_avg_',num2str(i),'.mat'],'RHO_avg');

end
