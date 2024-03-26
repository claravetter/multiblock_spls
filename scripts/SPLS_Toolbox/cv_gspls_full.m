%% DP function for one k split

function [RHO, weights, covariances, lVs] = cv_gspls_full(train_test_data, cs, correlation_method, printflag, V_original)

% CV: implement as user input
num_matrices = length(train_test_data);
gs = (1/(num_matrices-1)) * (ones(num_matrices) - eye(num_matrices));
e = 1e-5;
itr_lim = 1000;

%perform SPLS on the training data using the current cu/cv combination
if exist('V_original', 'var')
    train_data = cellfun(@(x) x.train, train_test_data);
    test_data = cellfun(@(x) (x.test), train_test_data);
    [weights, covariances, success, spls_itr, diff, dim_Ms] = cv_gspls(train_data, cs, gs, e, itr_lim, printflag, V_original);
            
else
    [weights, covariances, success, spls_itr, diff, dim_Ms] = cv_gspls(training_data_x, cs, gs, e, itr_lim, printflag);
end

%compute the correlation between the projections of the training and
%test matrices onto the SPLS latent space spanned by the weight vectors

[RHO, lVs, weights] = cv_projection(test_data, weights, correlation_method);
 
end