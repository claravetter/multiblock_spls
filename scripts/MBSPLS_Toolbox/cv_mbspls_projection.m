%% DP function for projection of u and v onto test data
function [RHO, lVs, weights, corrlVs] = cv_mbspls_projection(data, weights, correlation_method, matrix_norm)

for i=1:length(data)
    lVs(:,i) = data{i}*weights{i};
end

% TO DO: user defined input what to use (also max pairwise correlation etc.

if size(lVs,2) == 2 && matrix_norm == 0
    RHO= corr(lVs(:, 1), lVs(:, 2), 'Type', correlation_method);
    f_invert = @(x)(-1*x);
    if RHO<0 % if correlation negative, invert one weight vector --> correlation positive
        weights{2} = f_invert(weights{2});
        lVs(:, 2) = data{2}*weights{2};
        RHO = corr(lVs(:, 1), lVs(:, 2), 'Type', correlation_method);
    end

    if sum(weights{2})<0 % easier interpretation; v = phenotypic vector
        weights{1} = f_invert(weights{1});
        weights{2} = f_invert(weights{2});
        lVs(:, 1) = data{1}*weights{1};
        lVs(:, 2) = data{2}*weights{2};
        RHO = corr(lVs(:, 1), lVs(:, 2), 'Type', correlation_method);
    end
else
    corr_lVs = corr(lVs,'Type', correlation_method );
    if exist('matrix_norm', 'var')
        mn = norm(corr_lVs, matrix_norm); % 'fro', 1, 2, Inf
    else
        mn = norm(corr_lVs, 2);
    end
    RHO = mn;
end


end