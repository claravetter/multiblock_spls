function [p] = test_validation_lv_significance (input, output, Xs_validation, weights, RHO)
%input.validation_perm_flag = true;
for num_m=1:length(Xs_validation)
    permmat = cv_simplePermInd(input.validation_perm, size(input.Diag_validation,1));
    permmats{num_m} = permmat;

end

RHO_perms = [];
Xs_validation_ii = Xs_validation;
for ii=1:input.validation_perm
    % perform procrustean transformation to minimize rotation effects of
    % permutated y matrix, if V_opt available
    for num_m=1:length(Xs_validation)-1
        Xs_validation_ii{num_m+1} = Xs_validation{num_m+1}(permmats{num_m}(ii,:),:);%CV: does it matter which matrix is being permuted? Only one?
        % CV 10.4.2024: alle anderen Matrizen außer einer müssen permutiert werden --> permmat{} pro Matrix
    end
    %input.validation_perm_flag = false
    [RHO_perm, ~, ~] = cv_mbspls_projection(Xs_validation_ii, weights, input.correlation_method, input.matrix_norm);


    RHO_perms = [RHO_perms, RHO_perm];
end



switch input.statistical_testing
    case 1 % counting
        RHO_count_b = sum(RHO_perms > RHO);
        p = (RHO_count_b+1)/(input.validation_perm+1);
    case 2 % AUC testing
        IN.dist = RHO_perms;
        IN.val = RHO;
        IN.testing_precision = input.permutation_testing_precision;
        [p, ~] = dp_auc_testing(IN);
end



end


function permmat = cv_simplePermInd(nPerms, N)
% SIMPLEPERMIND  Return nPerms full‐size permutations of 1:N
%
%   permmat = simplePermInd(nPerms, N) returns an nPerms-by-N matrix
%   whose i-th row is randperm(N).

    permmat = zeros(nPerms, N);
    for j = 1:nPerms
        permmat(j,:) = randperm(N);
    end
end
