
function [RHO] = cv_gspls_slim(training_data, test_data, cs, correlation_method, Vs_original, matrix_norm)

% CV: implement as user input
num_matrices = length(training_data);
gs = (1/(num_matrices-1)) * (ones(num_matrices) - eye(num_matrices));
e = 1e-5;
itr_lim = 1000;
printflag = 0; 
%perform SPLS on the training data using the current cu/cv combination
if exist('Vs_original', 'var')
    [weights, ~, ~, ~, ~, ~, ~] = cv_gspls(training_data, cs, gs, e, itr_lim, printflag, Vs_original);
            
else
    [weights, ~, ~, ~, ~, ~, ~] = cv_gspls(training_data, cs, gs, e, itr_lim, printflag);
end

%compute the correlation between the projections of the training and
%test matrices onto the SPLS latent space spanned by the weight vectors

RHO = cv_projection(test_data, weights, correlation_method, matrix_norm);
 
end