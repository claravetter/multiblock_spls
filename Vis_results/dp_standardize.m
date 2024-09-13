%% DP function for standardizing

function [Z, IN] = dp_standardize(X, IN)

Z = nan(size(X,1),size(X,2));

if isfield(IN, 'method')
    method = IN.method;
else
    method = 'mean-centering';
end

switch method
    case 'min_max'
        if isfield(IN, 'min')
            for i=1:size(X,2)
                Z(:,i) = (X(:,i)-IN.min(i))./(IN.max(i)-IN.min(i));
            end
        else
            for i=1:size(X,2)
                IN.max(i) = max(X(:,i));
                IN.min(i) = min(X(:,i));
            end
            for i=1:size(X,2)
                Z(:,i) = (X(:,i)-IN.min(i))./(IN.max(i)-IN.min(i));
            end
        end
    case 'mean-centering'
        if isfield(IN, 'means')
            for i=1:size(X,2)
                Z(:,i) = (X(:,i)-IN.means(i))/IN.stds(i);
                if any(isnan(Z(:,i)))
                    Z(:,i) = X(:,i);
                end
            end
        else
            IN.means = mean(X,1);
            IN.stds = std(X,1);
            for i=1:size(X,2)
                Z(:,i) = (X(:,i)-IN.means(i))/IN.stds(i);
                if any(isnan(Z(:,i)))
                    Z(:,i) = X(:,i);
                end
            end
        end
end

end