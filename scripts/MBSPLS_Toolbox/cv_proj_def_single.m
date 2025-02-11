function deflated_matrix = cv_proj_def_single(matrix, weights)
    % NaN-aware matrix deflation: Deflates the matrix while ignoring NaNs.
    %
    % INPUTS:
    % matrix  - Input matrix (can contain NaN values)
    % weights - Column vector of weights (same length as number of columns in matrix)
    %
    % OUTPUT:
    % deflated_matrix - Deflated matrix with NaNs preserved

    [rows, cols] = size(matrix);
    
    % Ensure weights is a column vector
    weights = weights(:);

    % Create a logical mask for valid (non-NaN) values
    valid_mask = ~isnan(matrix);

    if any(~valid_mask(:))  % If there are NaNs in the matrix
        projection = NaN(rows, 1);  % Preallocate with NaNs

        % Compute projection term, ignoring NaNs
        for i = 1:rows
            valid_entries = matrix(i, valid_mask(i, :));  % Extract valid values
            valid_weights = weights(valid_mask(i, :));  % Extract corresponding weights
            
            if ~isempty(valid_entries)
                projection(i) = sum(valid_entries .* valid_weights') / sum(valid_weights); 
            end
        end

        % Compute the deflated matrix while preserving NaNs
        deflated_matrix = matrix;
        for i = 1:rows
            for j = 1:cols
                if valid_mask(i, j)  % Only modify non-NaN values
                    deflated_matrix(i, j) = matrix(i, j) - projection(i) * weights(j);
                end
            end
        end
    else
        % If there are no NaNs, apply direct matrix deflation
        deflated_matrix = matrix - (matrix * weights) * weights';
    end
end
