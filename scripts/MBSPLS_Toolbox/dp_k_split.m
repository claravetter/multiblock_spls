%% DP function for one k split

function [RHO] = dp_k_split(training_data_x,training_data_y,test_data_x, test_data_y, cu, cv, correlation_method)

% % initialize variables
% test_data_x = []; training_data_x=[]; test_data_y=[]; training_data_y=[];
% u=[]; v=[]; RHO=[];     
% 
% % separate the keep_in data into training and test data according to
% % the chosen test percentage tp
% [test_data_x, training_data_x, test_data_y, training_data_y] = dp_partition_holdout(tp, keep_in_data_x, keep_in_data_y);
        
%perform SPLS on the training data using the current cu/cv combination
[u, v, ~] = spls_suppressed_display(training_data_x,training_data_y,cu,cv); 
                    
%compute the correlation between the projections of the training and
%test matrices onto the SPLS latent space spanned by the weight vectors
RHO = abs(corr(test_data_x*u,test_data_y*v, 'Type', correlation_method)); 
 
end