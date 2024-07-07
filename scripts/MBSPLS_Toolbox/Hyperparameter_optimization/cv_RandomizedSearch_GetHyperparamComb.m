function hyperparameterCombinations = cv_RandomizedSearch_GetHyperparamComb(hyperparameterDistributions, numCombinations, seed)
    % Set the random seed
    if exist('seed', 'var')
        rng(seed);
    end
    
    % Number of hyperparameters
    numHyperparameters = numel(hyperparameterDistributions);
    
    % Initialize array to store hyperparameter combinations
    hyperparameterCombinations = zeros(numCombinations, numHyperparameters);
    
    % Generate random samples from the distributions for each hyperparameter
    for i = 1:numHyperparameters
        distribution = hyperparameterDistributions{i};
        hyperparameterCombinations(:, i) = random(distribution, numCombinations, 1);
    end
end
