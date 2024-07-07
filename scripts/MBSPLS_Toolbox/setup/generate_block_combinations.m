function [combinations, combination_names, combination_feature_names] = generate_block_combinations(data_matrices, data_names, feature_names)
    % Initialize cell arrays for combinations
    all_Xs = {};
    all_Xs_names = {};
    
    % Number of modalities
    num_modalities = length(data_matrices);
    
    % Generate all possible combinations
    combinations = {};
    combination_names = {};
    combination_feature_names = {};
    n_comb = 0;
    for num_comb = 1:num_modalities
        comb_indices = nchoosek(1:num_modalities, num_comb);
       
        for i = 1:size(comb_indices, 1)
            combination = {};
            names = {};
            feat_names = {};
            if num_comb == 1
                combination = {data_matrices{comb_indices(i)}, data_matrices{comb_indices(i)}};
                names = {data_names{comb_indices(i)},data_names{comb_indices(i)}};
                feat_names = {feature_names{comb_indices(i)},feature_names{comb_indices(i)}};
            else
                for j = 1:num_comb
                    combination{j} = data_matrices{comb_indices(i,j)};
                    names{j} = data_names{comb_indices(i,j)};
                    feat_names{j} = feature_names{comb_indices(i,j)};
                end
            end
            combinations{n_comb+i} = combination;
            combination_names{n_comb+i} = names;
            combination_feature_names{n_comb+i} = feat_names;
        end
        n_comb = size(combinations,2);

    end
end