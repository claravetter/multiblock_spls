%% DP SPLS standalone function

function dp_spls_standalone_cleanup(spls_standalone_path, analysis_folder, queue_name, variables_for_analysis)

%% 1. prepare analysis folders
mkdir([analysis_folder '/hyperopt']);
hyperopt_folder = [analysis_folder '/hyperopt']; % folder for hyperparameter optimization
mkdir([analysis_folder '/permutation']);
permutation_folder = [analysis_folder '/permutation']; % folder for permutation testing
mkdir([analysis_folder '/utilities']);
utilities_folder = [analysis_folder '/utilities']; % folder containing general settings => Section 3 parameter settings

%% 2. set parameters for queue submission
mem_total           = 40;   % max: 40
max_sim_jobs        = 60;   % max: 60

%% 3. set parameters for SPLS analysis
hp  = 10;       % hold-out percentage, Default: 10
tp  = 20;       % test percentage, Default: 20
B   = 10000;    % number of permutations for statistical evaluation, Default: 10000
W   = 10;       % number of iterations for one LV, Default: 10
K   = 100;      % number of K-splits during hyperparameter optimization, Default: 100
RP  = 40;       % number of points for grid search of cu/cv combinations, Default: 40
FDRvalue = 0.05; % significance level for FDR testing

utilities_names = {'hyperopt_folder', 'permutation_folder', 'mem_total', 'max_sim_jobs', 'hp', 'tp', 'B', 'W', 'K',... 
    'RP', 'spls_standalone_path', 'analysis_folder', 'queue_name', 'variables_for_analysis'};
utilities = {hyperopt_folder, permutation_folder, mem_total, max_sim_jobs, hp, tp, B, W, K, RP, spls_standalone_path,... 
    analysis_folder, queue_name, variables_for_analysis};

for i=1:size(utilities,2)
    switch class(utilities{i})
        case 'char'
            dp_txt_write(utilities_folder, utilities_names{i}, utilities{i}, '%s');
        case 'double'
            dp_txt_write(utilities_folder, utilities_names{i}, utilities{i}, '%d');
    end
end

%% 4. prepare accessory variables, indices and initialize variables
load(variables_for_analysis);
X_original = zscore(X);
Y_original = zscore(Y);

% define the data matrices X and Y for input into the SPLS algorithm
X = X_original;
Y = Y_original;

%% 5. LV loop

% this loop is repeated for the next vector pair u and v as long as the
% previous vector pair was significant. If the previous vector was not
% significant, then the loop will stop

% pLV is a starting p value set to 1, so that the while loop starts. Later
% on it will be updated at the end of each LV iteration. Then it serves the
% purpose to stop the while loop as soon as an LV was not significant.
p_LV = 1; 

% ff counts the iteration through the LV loops, so that the results for each
% LV can be stored in separate rows
ff = 1;  

% define column names for matrices so that you can access them later by
% indexing
parameters_names = {'w', 'cu', 'cv', 'u', 'v', 'success', 'RHO', 'p'}; % names of parameters for optimal cu/cv combination of the w-th loop
cu_cv_combination_names = {'cu';'cv'};

% indices for later use
opt_u = strcmp(parameters_names,'u');
opt_v = strcmp(parameters_names,'v');
comb_cu = strcmp(cu_cv_combination_names,'cu');
comb_cv = strcmp(cu_cv_combination_names,'cv');
opt_RHO = strcmp(parameters_names,'RHO');
opt_p = strcmp(parameters_names, 'p');  

% cu_range and cv_range define the range for the grid search. The grid
% search is performed along 40 points from 1 to sqrt(number of variables)
% as proposed in Monteiro et al. 2016.
cu_range_points = linspace(1,sqrt(size(X,2)),RP);
cv_range_points = linspace(1,sqrt(size(Y,2)),RP);

% compile a matrix with separate row for all possible cu and cv
% combinations by taking cu and repeating every single element 40 times and
% then takin cv and repeating the entire vector 40 times
cu_range = repelem(cu_range_points,numel(cu_range_points));
cu_range(:,1:RP)=[]; %adapted for better performance
cv_range = repmat(cv_range_points,1,numel(cv_range_points));
cv_range(:,1:RP)=[]; %adapted for better performance
cu_cv_combination = [cu_range;cv_range]; % this matrix contains all the possible combinations between cu and cv

for i=1:size(cu_range,2)
    dp_txt_write(hyperopt_folder, ['cu_' num2str(i)], cu_range(i), '%.4f')
end

for i=1:size(cv_range,2)
    dp_txt_write(hyperopt_folder, ['cv_' num2str(i)], cv_range(i), '%.4f')
end

% preallocate placeholder matrices for higher speed
final_parameters = num2cell(nan(size(Y,2), numel(parameters_names))); % cell array to store final parameters for each LV iteration
epsilon = nan(size(final_parameters,1),size(X_original,1)); %projection of weight vector u on corresponding matrix X
omega = nan(size(final_parameters,1),size(Y_original,1)); % projection of weight vector v on corresponding matrix Y

% here starts the outer loop for each single LV, the next LV is only
% computed if the previous LV was significant (with FWE correction)
while p_LV > 0
 
%% 3. Wide loop
% repeats the entire process W times to obtain W different final p values,
% which go into the omnibus hypothesis

% matrix to store the optimal parameters of the w-th loop
opt_parameters = num2cell(nan(W,numel(parameters_names)));

    for w=1:W

        %% hyper-parameter optimisation
        % remove hp% of the data randomly and keep it in a hold-out
        % dataset, the rest is called the keep_in data for training/testing
        % the algorithm
        [hold_out_data_x, keep_in_data_x, hold_out_data_y, keep_in_data_y] = dp_partition_holdout(hp, X, Y);

        % save the hyperparameter information for optimization step
        save([hyperopt_folder '/keep_in_partition.mat'], 'keep_in_data_x','keep_in_data_y');
                
        % write a bash script for hyperparameter optimization, which can
        % later be called
        no_jobs = size(cu_cv_combination,2);
        hyperopt_bash = dp_bash_para_new(spls_standalone_path, hyperopt_folder, utilities_folder, 'hyperopt', mem_total, no_jobs, max_sim_jobs, queue_name);
        
        % access the temp_RHO folder and initialize the bash script for
        % parallel computing of the hyperparameter optimization
        cd(hyperopt_folder);
        system(['qsub ' hyperopt_bash]);

        % set a variable which counts the files in the RHO folder
        filecount = size(dir([hyperopt_folder '/RHO_avg_*.txt']),1);
        RHO_avg_collection = nan(size(cu_cv_combination,2),1);

        % use a while loop which recounts all the files in the RHO folder
        % until all computations from the parallelized hyperparameter
        % optimization step are completed and all files are saved in that
        % folder, the while loop will be exited as soon as all files are
        % saved in the folder
        while filecount < size(cu_cv_combination,2)
            filecount = size(dir([hyperopt_folder '/RHO_avg_*.txt']),1);
        end
        
        % load all RHO values which were saved in mat files into a
        % RHO_collection vector
        for i=1:size(cu_cv_combination,2)
            if exist([hyperopt_folder '/RHO_avg_' num2str(i) '.txt'],'file')
                RHO_avg_collection(i,1) = dp_txtscan([hyperopt_folder '/RHO_avg_' num2str(i) '.txt'], '%f');
                delete([hyperopt_folder '/RHO_avg_' num2str(i) '.txt']); 
            else
                RHO_avg_collection(i,1) = NaN;
            end
        end
        
%     select the hyperparameter combination with the highest average
%     correlation, ie the one with the highest RHO_avg
        [~, index_RHO] = max(RHO_avg_collection); 
    
%     find the corresponding cu and cv values of the max RHO element, this
%     gives you the cu/cv combination with the highest average correlation
    
    cu_opt = cu_cv_combination(comb_cu,index_RHO);
    cv_opt = cu_cv_combination(comb_cv,index_RHO);
    
    % save the optimized parameters
    save([permutation_folder '/opt_param.mat'],...
                'keep_in_data_x',...
                'keep_in_data_y',...
                'hold_out_data_x',...
                'hold_out_data_y',...
                'cu_opt',...
                'cv_opt',...
                'RHO_avg_collection',...
                'permutation_folder'); 
    
    %% Statistical Evaluation
    
    
    %train the model with the optimal cu/cv combination on the previous
    %train/test data set to get u and v
    [u_opt, v_opt, success_opt]=spls(keep_in_data_x, keep_in_data_y,cu_opt, cv_opt);

    %project the hold_out data set on the computed u and v vectors to
    %get the absolute correlations between these projections
    RHO_opt = abs(corr(hold_out_data_x*u_opt,hold_out_data_y*v_opt));

    % permute the keep_in data in one dimension so that the
    % relationship between these two views gets destroyed
    
%     Now comes the permutation step, where the order of the samples in one
%     view (in this case the Y matrix) gets permuted and the algorithm is
%     trained on the permuted data set with a destroyed relationship
%     between X and Y, this process is repeated B times

    for i=1:B
        temp = randperm(size(keep_in_data_y,1))';
        dp_txt_write(permutation_folder, ['perm_' num2str(i)], temp, '%d\n');
    end
    
    % write a bash script for hyperparameter optimization, which can
    % later be called
    permutation_bash = dp_bash_para_new(spls_standalone_path, permutation_folder, utilities_folder, 'permutation', mem_total, B, max_sim_jobs, queue_name);
    
    % submit the bash script for parallel computing of RHO_b values to the
    % queue
    cd(permutation_folder);
    system(['qsub ' permutation_bash]);
    
    % collect all files with RHO_b values, open them and store the RHO
    % values in a vector
    filecount = size(dir([permutation_folder '/RHO_b_*']),1);
    RHO_b_collection = nan(B,1);
    
    % see hyperparameter optimization
    while filecount < B
        filecount = size(dir([permutation_folder '/RHO_b_*.txt']),1);
    end
    
    for i=1:size(RHO_b_collection,1)
        if exist([permutation_folder '/RHO_b_',num2str(i),'.txt'],'file')
            RHO_b_collection(i,1) = dp_txtscan([permutation_folder '/RHO_b_' num2str(i) '.txt'], '%f');            
            delete([permutation_folder '/RHO_b_',num2str(i),'.txt']);
        else
            RHO_b_collection(i,1) = NaN;
        end
    end
    
    RHO_b_collection_completed = [];
    RHO_b_collection_completed = dp_perm_nan(RHO_b_collection, spls_standalone_path, permutation_folder, utilities_folder, mem_total, max_sim_jobs, queue_name);
    
    % test the following null hypothesis H_s: "There is no relationship
    % between the two views, therefore the correlation obtained with
    % the original data is not different from the correlation
    % obtained with the permuted data"

    % calculate how many times RHO_b was bigger than or equal to RHO_opt,
    % ie how many times the permutated data led to a higher correlation
    % using spls than the original data
    RHO_count_b = sum(RHO_b_collection_completed >= RHO_opt);

    % calculate the p value to test whether the null hypothesis is true
    p = ((1+RHO_count_b)/(B+1));
  
    % store all parameters of the w-th wide loop into a matrix
    opt_parameters(w,:) = {w cu_opt cv_opt u_opt v_opt success_opt RHO_opt p};
    
    save([analysis_folder '/opt_parameters_' num2str(ff) '_' num2str(w) '.mat'],'opt_parameters', 'RHO_b_collection_completed'); 
    
    end
        
    %% test the omnibus hypothesis using the p values

    % Omnibus hypothesis H_omni: "All the null hypotheses H_s are true" if
    % any of the p values are lower than the FWE-corrected threshold,
    % then the omnibus hypothesis can be rejected after that. 
    
    % Search for the statistically significant w-th iteration with the
    % lowest p value, if more than one iteration has the same p value,
    % select the one that has the highest RHO value => this is the final
    % w-th iteration with the corresponding cu/cv features and u/v scores
    
    % remove parameter combinations with nans as RHO_opt
    to_rmv=isnan([opt_parameters{:,opt_RHO}]');
    opt_parameters(to_rmv,:)=[];

    % check if p values are significant with FDR-correction
    [pvalues_FDR] = dp_FDR([opt_parameters{:,opt_p}]', FDRvalue);
    
    if isempty(pvalues_FDR)
        p_LV = 0;
    else
        
    logical_p_min = [opt_parameters{:,opt_p}]' == min([opt_parameters{:,opt_p}]);
    temp1 = opt_parameters(logical_p_min,:); % placeholder to make syntax easier
    
    if size(temp1,1) > 1 % if there is more than one p-value that is equal to the lowest p-value, ie there are two absolute minimal p-values, then continue in the loop
        RHO_max = [temp1{:,opt_RHO}]' == max([temp1{:,opt_RHO}]); % logical indexing the highest RHO value
        final_parameters(ff,:) = temp1(RHO_max,:); % take the row with the lowest p-value and the highest RHO value
    else
        final_parameters(ff,:) = temp1; % if there is just one absolute minimal p value, then take this row
    end

    %% Matrix Deflation => Projection Deflation

    % This method removes the covariance explained by u and v from X
    % and Y by projecting the data matrices onto the space spanned by the
    % corresponding weight vector, and subtracting this from the data:
    u = final_parameters{ff,opt_u}; % weight vector for X
    v = final_parameters{ff,opt_v}; % weight vector for Y
    [X,Y] = proj_def(X, Y, u, v);
    
    % the p value of this LV is updated from the previously set value of
    % zero to the actual p value of the LV, if this p value is lower than
    % the FWE_rate then the while loop continues after matrix deflation and
    % it keeps looking for the next possible significant LV. If the p value
    % is higher than the FWE_rate, then the while loop is discontinued and
    % the algorithm stops. Therefore, this function generates all possible
    % significant LVs and also the first non-significant LV, afterwards the
    % algorithm stops.
    p_LV = final_parameters{ff,opt_p}; 
    
    
    %% Projection of the original data onto the SPLS latent space
    
    % epsilon and omega are the projections of the original data onto the
    % SPLS latent space, spanned by the weight vector pair u and v of this
    % LV. Therefore each LV spans a specific latent space onto wihich the
    % samples can be projected. These latent spaces can then be used to
    % stratify patients.
    epsilon(ff,:) = X_original * final_parameters{ff,opt_u};
    omega(ff,:) = Y_original * final_parameters{ff,opt_v};

    % ff counts through the most outside LV while loop = counts the amount
    % of significant LVs
    ff = ff+1; 
    

end

%% clean up txt files and other files which are not needed anymore
files_to_clean = {'cu_', 'cv_'};
folder_to_clean = hyperopt_folder;
for i=1:size(files_to_clean,2)
    dp_cleanup(folder_to_clean, files_to_clean{i})
end

% after the LVs are computed, clean up the final parameters, epsilon and
% omega matrix by removing empty rows
final_parameters(ff:end,:) = [];
epsilon(ff:end,:) = [];
omega(ff:end,:) = [];

save([analysis_folder '/result.mat'], 'final_parameters', 'epsilon', 'omega');

end