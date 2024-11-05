function explained_variance_all_blocks = compute_explained_variance_all_blocks(latentscores_file, resultpath)
% Computes explained variance for all blocks in Xs using latent variables in Ts
% Xs: Cell array of data blocks (e.g., {X1, X2, X3, ...})
% Ts: Cell array of latent variables for each block (e.g., {T1, T2, T3, ...})
% output: Structure containing final parameters (loadings) for each block
% explained_variance_all_blocks: Cell array storing explained variance for each block

load(resultpath)
Xs = input.Xs; 

% Center and scale Xs; replace with covariate and standardisation as in
% training
for i = 1:size(Xs,2)
    Xs{i} = (Xs{i} - mean(Xs{i}, 1)) ./ std(Xs{i}, 0, 1);
end
% Get all sheet names in the Excel file
lv_sheets = sheetnames(latentscores_file);

n_LVs = size(lv_sheets,1);
% Loop through each LV
Ts = {};
for i = 1:n_LVs
    lv_sheet = lv_sheets{i};
    LV = readmatrix(latentscores_file, 'Sheet', lv_sheet);
    for num_m = 1:size(Xs,2)
        Ts{num_m}(:,i) = LV(:,1+num_m);
        temp_weights = output.final_parameters(i, 3); 
        Ps{num_m}(:,i) = temp_weights{1,1}{1, num_m};  % Load the loading vectors for block i
    
    end
end

n_blocks = length(Xs);  % Number of data blocks
explained_variance_all_blocks = cell(1, n_blocks);  % Initialize cell array to store results

% Loop through each block
for i = 1:n_blocks
    X = Xs{i};  % Get the data block
    T = Ts{i};  % Get the corresponding latent variables (scores)
    P = Ps{i};
    % Load the corresponding loading matrix P from output.final_parameters
    % Assuming the first dimension in output.final_parameters corresponds to latent variables
    % and the third element in output.final_parameters contains loadings for the blocks
    
    % Call the function to compute explained variance for the current block
    explained_variance_all_blocks{i} = compute_explained_variance_X(X, T, P);
end

% Display the explained variance for all blocks
disp('Explained variance for all blocks:');
disp(explained_variance_all_blocks);

end
