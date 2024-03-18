%% function to find best opt_param combination

function best_row = dp_find_opt(IN)

matrix_names = IN.matrix_names;
matrix = IN.matrix;
type_analysis = IN.type_analysis;
p = IN.p;

opt_RHO = strcmp(matrix_names, 'RHO');
opt_p = strcmp(matrix_names, 'p');
opt_cu = strcmp(matrix_names, 'cu');
opt_cv = strcmp(matrix_names, 'cv');
opt_u = strcmp(matrix_names, 'u');
opt_v = strcmp(matrix_names, 'v');

% find all significant rows
log_sig = cell2mat(matrix(:,opt_p))<= p;
matrix(~log_sig,:)=[];

max_cu = sqrt(size(matrix{1, opt_u},1));
max_cv = sqrt(size(matrix{1, opt_v},1));

% scale hyperparameters
matrix_temp = matrix;
matrix_temp(:,opt_cu) = num2cell(cell2mat(matrix_temp(:, opt_cu))./max_cu);
matrix_temp(:,opt_cv) = num2cell(cell2mat(matrix_temp(:, opt_cv))./max_cv);


switch type_analysis
    case 1
        % choose the sparsest solution from all significant iterations
        log_target = (cell2mat(matrix_temp(:, opt_cu))+cell2mat(matrix_temp(:, opt_cv))) == min((cell2mat(matrix_temp(:, opt_cu))+cell2mat(matrix_temp(:, opt_cv))));
        if sum(log_target)>1
            new_matrix = matrix(log_target,:);
            log_target1 = abs(cell2mat(new_matrix(:, opt_RHO))) == max(abs(cell2mat(new_matrix(:,opt_RHO))));
            best_row = new_matrix(log_target1,:);
        else
            best_row = matrix(log_target,:);
        end
    case 2
        % choose the solution with the highest RHO from all significant
        % iterations
        log_target = abs(cell2mat(matrix_temp(:, opt_RHO))) == max(abs(cell2mat(matrix_temp(:,opt_RHO))));
        if sum(log_target)>1
            new_matrix = matrix_temp(log_target,:);
            log_target1 = (cell2mat(new_matrix(:, opt_cu))+cell2mat(new_matrix(:, opt_cv))) == min((cell2mat(new_matrix(:, opt_cu))+cell2mat(new_matrix(:, opt_cv))));
            new_matrix1 = matrix(log_target,:);
            best_row = new_matrix1(log_target1,:);
        else
            best_row = matrix(log_target,:);
        end
    case 3
        % choose the lowest p value and if needed the one with
        % addiationally the highest RHO value from all significant
        % iterations
        log_target = cell2mat(matrix_temp(:, opt_p)) == min(cell2mat(matrix_temp(:, opt_p)));
        if sum(log_target)>1
            new_matrix = matrix(log_target,:);
            log_target1 = abs(cell2mat(new_matrix(:, opt_RHO))) == max(abs(cell2mat(new_matrix(:,opt_RHO))));
            best_row = new_matrix(log_target1,:);
        else
            best_row = matrix(log_target,:);
        end
end

end