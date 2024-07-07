function modalities = simulate_data(sample_size, nr_variables, nr_active_variables, nr_modalities, epsilon)
    N = sample_size;
    p = nr_variables; % Number of variables in each modality

    % Shared latent variable
    t = randn(N, 1);

    % Initialize weight vectors and modalities cell array
    W = cell(1, nr_modalities);
    modalities = cell(1, nr_modalities);

    % Generate weight vectors and data for each modality except X
    for i = 1:nr_modalities-1
        % Generate weight vector for current modality
        W{i} = zeros(p(i), 1);
        active_indices = randperm(p(i), nr_active_variables(i));
        W{i}(active_indices) = 1;

        % Generate error matrix for current modality
        E = epsilon * randn(N, p(i));

        % Generate data for current modality
        modalities{i} = t * W{i}' + E;
    end

    % Generate data for X (last modality)
    W{nr_modalities} = zeros(p(nr_modalities), 1);
    active_indices = randperm(p(nr_modalities), nr_active_variables(nr_modalities));
    W{nr_modalities}(active_indices) = 1;
    E = epsilon * randn(N, p(nr_modalities));
    modalities{nr_modalities} = t * W{nr_modalities}' + E;
end
