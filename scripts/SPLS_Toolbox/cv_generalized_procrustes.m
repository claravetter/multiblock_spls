function [aligned_matrices, converged] = cv_generalized_procrustes(varargin)
    % Generalized Procrustes analysis for aligning multiple matrices.

    % Parameters:
    %     varargin: cell array of matrices
    %         Input matrices to be aligned.
    %     max_iter: int, optional
    %         Maximum number of iterations for Procrustes alignment. Default is 100.
    %     tol: float, optional
    %         Convergence tolerance. Default is 1e-6.

    % Returns:
    %     aligned_matrices: cell array of matrices
    %         Aligned matrices.
    %     converged: logical
    %         True if the algorithm converged within tolerance, False otherwise.

    matrices = varargin;
    num_matrices = length(matrices);
    shape = size(matrices{1});
    
    % Check if all matrices have the same shape
    for i = 2:num_matrices
        if ~isequal(size(matrices{i}), shape)
            error('All matrices must have the same shape.');
        end
    end
    
    aligned_matrices = matrices;
    rotations = cell(1, num_matrices);
    for i = 1:num_matrices
        rotations{i} = eye(shape(2)); % Initialize rotation matrices randomly
    end
    
    max_iter = 100;
    tol = 1e-6;
    converged = false;
    for iter = 1:max_iter
        prev_rotations = rotations;
        
        % Align each matrix to the mean of the other matrices
        for i = 1:num_matrices
            mean_matrix = mean(cat(3, aligned_matrices{[1:i-1, i+1:end]}), 3);
            [U, ~, V] = svd(mean_matrix' * aligned_matrices{i});
            rotations{i} = U * V';
            aligned_matrices{i} = aligned_matrices{i} * rotations{i}';
        end
        
        % Check for convergence
        if all(cellfun(@(x, y) norm(x - y, 'fro'), rotations, prev_rotations) < tol)
            converged = true;
            break;
        end
    end
    
    aligned_matrices = aligned_matrices;
end
