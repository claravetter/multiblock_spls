%% DP function for projection of u and v onto test data
function [RHO, lVs, weights] = cv_projection(data, weights, correlation_method)

for i=1:length(data)
    lVs(:,i) = data{i}*weights{i};
end

% TO DO: user defined input what to use (also max pairwise correlation etc.

if size(lVs,2) == 2
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
    mnfrob = norm(corr_lVs, 'fro');
    mn1 = norm(corr_lVs, 1);
    mn2 = norm(corr_lVs, 2);
    mnInf = norm(corr_lVs, Inf);
    RHO = mn2;
end


end