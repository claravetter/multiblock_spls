function cv_plot_aligned_matrices(original_matrices, aligned_matrices)
    % Plot original and aligned matrices for visualization.

    % Parameters:
    %     original_matrices: cell array of matrices
    %         Original matrices.
    %     aligned_matrices: cell array of matrices
    %         Aligned matrices.

    num_matrices = length(original_matrices);
    for i = 1:num_matrices
        subplot(2, num_matrices, i);
        imshow(original_matrices{i}, []);
        title(['Original Matrix ', num2str(i)]);
        
        subplot(2, num_matrices, num_matrices + i);
        imshow(aligned_matrices{i}, []);
        title(['Aligned Matrix ', num2str(i)]);
    end
    sgtitle('Original and Aligned Matrices');
end
