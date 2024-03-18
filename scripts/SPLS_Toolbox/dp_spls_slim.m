%% DP function for one k split

function RHO = dp_spls_slim(training_data_x,training_data_y,test_data_x, test_data_y, cu, cv, correlation_method, V_original)
        
%perform SPLS on the training data using the current cu/cv combination
if exist('V_original', 'var')
    [u, v, ~, ~, ~, ~] = dp_spls(training_data_x, training_data_y, cu, cv, V_original);
else
    [u, v, ~, ~, ~, ~] = dp_spls(training_data_x, training_data_y, cu, cv);
end

%compute the correlation between the projections of the training and
%test matrices onto the SPLS latent space spanned by the weight vectors

RHO = dp_projection(test_data_x, test_data_y, u, v, correlation_method);
 
end