% explained covariance 

function [explained_covariance_per_component,total_explained_covariance]  = cv_compute_explained_covariance_LV(matrices, T_matrices)
% matrices = {X1, X2, X3};
% T_matrices = {T1, T2, T3}; dimensions: n_samples * n_LVs
n_Xs = length(matrices);  % Total number of data blocks (matrices)
n_LVs = size(T_matrices{1}, 2);  % Assuming T1 has the latent components

% Initialize a cell array to store covariance matrices
covariances = cell(n_Xs, n_Xs);

% Loop through each pair of matrices to compute the covariance matrices
for i = 1:n_Xs
    for j = 1:n_Xs
        % Ensure the number of samples (rows) in the two blocks are the same
        if size(matrices{i}, 1) == size(matrices{j}, 1)
            % Compute covariance matrix between block i and block j
            covariances{i, j} = (matrices{i}' * matrices{j}) / size(matrices{i}, 1);  % Normalized by sample size
        else
            error('Number of rows (samples) in matrices %d and %d do not match', i, j);
        end
    end
end

% Now, calculate the explained covariance for each latent component across the blocks
explained_covariance = zeros(n_LVs, n_Xs, n_Xs);
% Assuming you have the following data:
% matrices: cell array where each entry is a data block (e.g., X1, X2, X3, etc.)
% T_matrices: cell array where each entry contains the latent components (e.g., T1, T2, T3, etc.)
% n_Xs: number of blocks (e.g., 3 for X1, X2, X3)
% n_LVs: number of latent components

n_Xs = length(matrices);  % Total number of data blocks (matrices)
n_LVs = size(T_matrices{1}, 2);  % Assuming T1 has the latent components

% Initialize a cell array to store covariance matrices
covariances = cell(n_Xs, n_Xs);

% Loop through each pair of matrices to compute the covariance matrices
for i = 1:n_Xs
    for j = 1:n_Xs
        % Ensure the number of samples (rows) in the two blocks are the same
        if size(matrices{i}, 1) == size(matrices{j}, 1)
            % Compute covariance matrix between block i and block j
            covariances{i, j} = (matrices{i}' * matrices{j}); 
        else
            error('Number of rows (samples) in matrices %d and %d do not match', i, j);
        end
    end
end

% Now, calculate the explained covariance for each latent component across the blocks
explained_covariance = zeros(n_LVs, n_Xs, n_Xs);

% Loop over latent components to compute explained covariance for each pair of blocks
for h = 1:n_LVs
    for i = 1:n_Xs
        for j = 1:n_Xs
            if i ~= j  % We are interested in the covariance between different blocks
                % Get the latent components for the current block pair
                T_i = T_matrices{i}(:, h);  % Latent component h for block i
                T_j = T_matrices{j}(:, h);  % Latent component h for block j
                
                % Compute the covariance between the latent components T_i and T_j
                cov_Ti_Tj = cov(T_i, T_j);  % This gives a 2x2 covariance matrix
                
                % Use SVD on the covariance matrix between the blocks
                [~, S, ~] = svd(covariances{i, j});  % SVD of the covariance matrix
                
                % Sum the singular values of the covariance matrix
                sum_singular_values = sum(diag(S));
                
                % Compute explained covariance using sum of singular values
                explained_covariance(h, i, j) = cov_Ti_Tj(1, 2) / sum_singular_values;
            end
        end
    end
end

% Calculate the average explained covariance across all block pairs for each latent component
explained_covariance_per_component = zeros(1, n_LVs);
for h = 1:n_LVs
    total_cov = 0;
    count = 0;
    for i = 1:n_Xs
        for j = 1:n_Xs
            if i ~= j  % Avoid self-covariance
                total_cov = total_cov + explained_covariance(h, i, j);
                count = count + 1;
            end
        end
    end
    explained_covariance_per_component(h) = total_cov / count;  % Average explained covariance
end

% Display results
disp('Explained covariance per latent component (SVD-based):');
disp(explained_covariance_per_component);


% Initialize variables to store total covariance
total_shared_covariance = 0;

% Calculate total shared covariance across all block pairs (using SVD)
for i = 1:n_Xs
    for j = 1:n_Xs
        if i ~= j
            % SVD on the covariance matrix between blocks i and j
            [~, S, ~] = svd(covariances{i, j});
            % Sum the singular values (total covariance between blocks)
            total_shared_covariance = total_shared_covariance + sum(diag(S));
        end
    end
end

% Now, normalize the explained covariance for each latent component
explained_covariance_per_component = zeros(1, n_LVs);
for h = 1:n_LVs
    total_cov = 0;
    count = 0;
    for i = 1:n_Xs
        for j = 1:n_Xs
            if i ~= j  % Avoid self-covariance
                total_cov = total_cov + explained_covariance(h, i, j);
                count = count + 1;
            end
        end
    end
    % Normalize by the total shared covariance
    explained_covariance_per_component(h) = total_cov / total_shared_covariance;
end

% Check the total explained covariance across all components (should be ≤ 1)
total_explained_covariance = sum(explained_covariance_per_component);

% Display results
disp('Normalized explained covariance per latent component (SVD-based):');
disp(explained_covariance_per_component);
disp(['Total explained covariance (as proportion of total shared covariance): ', num2str(total_explained_covariance)]);

end

