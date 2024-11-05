function explained_variance_per_component = compute_explained_variance_X(X, T, P)
% Computes the explained variance in block X explained by the latent variables
% X: data matrix for the block (n_samples x n_variables)
% T: latent variable matrix (n_samples x n_LVs)
% P: loading vector matrix (n_variables x n_LVs)

% Step 1: Compute the total variance in the block X
total_variance_X = trace(X' * X);

% Step 2: Initialize variable to store explained variance per component
n_LVs = size(T, 2);  % Number of latent variables
explained_variance_per_component = zeros(1, n_LVs);

% Step 3: Loop through each latent variable to compute explained variance
for h = 1:n_LVs
    % Reconstruct the data using the h-th latent variable
    X_hat = T(:, h) * P(:, h)';  % Reconstructed data for component h
    
    % Compute the variance explained by the h-th latent variable
    explained_variance_per_component(h) = trace(X_hat' * X_hat) / total_variance_X;
end

% Display results
disp('Explained variance per latent variable (in X):');
disp(explained_variance_per_component);

end
