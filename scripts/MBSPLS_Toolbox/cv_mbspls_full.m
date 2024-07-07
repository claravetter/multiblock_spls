%% DP function for one k split

function [RHO, weights, Vs, lVs] = cv_mbspls_full(training_data, test_data, cs, mbspls_params, correlation_method, matrix_norm, Vs_original)

% CV: implement as user input
num_matrices = length(training_data);
gs = mbspls_params.gs; %(1/(num_matrices-1)) * (ones(num_matrices) - eye(num_matrices));
e = mbspls_params.e;% 1e-5;
itr_lim = mbspls_params.itr_lim;% 1000;
printflag = mbspls_params.printflag;%0; 
%perform SPLS on the training data using the current cu/cv combination
if exist('Vs_original', 'var')
    [weights, ~, Vs, ~, ~, ~, ~] = cv_mbspls(training_data, cs, gs, e, itr_lim, printflag, Vs_original);
            
else
    [weights, ~, Vs, ~, ~, ~, ~] = cv_mbspls(training_data, cs, gs, e, itr_lim, printflag);
end

%compute the correlation between the projections of the training and
%test matrices onto the SPLS latent space spanned by the weight vectors

[RHO, lVs, weights] = cv_mbspls_projection(test_data, weights, correlation_method, matrix_norm);
 
end