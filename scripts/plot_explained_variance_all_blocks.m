function plot_explained_variance_all_blocks(explained_variance_all_blocks)
% Visualize explained variance for all blocks in a single grouped bar plot
% explained_variance_all_blocks: Cell array storing explained variance for each block

n_blocks = length(explained_variance_all_blocks);  % Number of blocks
n_LVs = length(explained_variance_all_blocks{1});  % Number of latent variables

% Combine explained variance into one matrix for grouped bar plot
explained_variance_matrix = zeros(n_blocks, n_LVs);
for i = 1:n_blocks
    explained_variance_matrix(i, :) = explained_variance_all_blocks{i};
end

% Create grouped bar plot
figure;
b = bar(explained_variance_matrix, 'grouped');  % Grouped bar plot

% Set colors for the different latent variables
colormap(parula(n_LVs));  % Use a colormap with as many colors as there are LVs

% Add labels and title
xlabel('Blocks');
ylabel('Explained Variance');
title('Explained Variance Grouped by Block and Latent Variable');

% Add legend for the latent variables
legend(arrayfun(@(x) ['LV' num2str(x)], 1:n_LVs, 'UniformOutput', false), 'Location', 'Best');

% Display grid for clarity
grid on;

end

