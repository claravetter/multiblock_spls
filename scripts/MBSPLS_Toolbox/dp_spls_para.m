function [final_parameters, epsilon, omega] = dp_spls_para(folders_for_analysis, variables_for_analysis, parameters_for_analysis)
%
%   Implementation of SPLS algorithm and framework as described in Monteiro
%   et al. 2016, for details: doi:10.1016/j.jneumeth.2016.06.011
%   
% Inputs:
%
%     X, Y    - data matrices for MRI images (X) and behavioral data (Y) in
%               the form: samples x features. Use MRI data directly from
%               Neurominer output, the function will standardize both the
%               MRI data and the behavioral data, so that mean = 0 and std =
%               1.
% 
%     hp      - Ratio (in %) of hold out data, which is projected on both
%               the original data set with optimal cu/cv and on the
%               permuted dataset with optimal cu/cv. Default: 10 
% 
%     tp      - Ratio (in %) of the test data within the K-fold loops.
%               Default: 20 
%
%     B       - Number of Permutations for the statistical evaluation of
%               the discovered cu/cv combination. Default: 10000
% 
%     W       - Number of wide loops that define how many times the entire
%               SPLS procedure is repeated and how many final p values are
%               calculated. Default: 10
%
%     K       - Number of K splits for the training and testing of the held
%               in data for each cu/cv combination. Default: 100
%   
% 
% Outputs: 
% 
%     final_parameters  - Results matrix with rows = number of significant
%                       LVs and columns containing the following variables
%                       w, cu_opt, cv_opt, u_opt, v_opt, success_opt,
%                       RHO_opt, p. The last row contains the first LV that
%                       was not significant, it is for informative purposes
%                       only.
%
%     epsilon           - Matrix with rows = Number of LV and columns
%                       containing the projection of the weight vector u of
%                       the LV on X.
%
%     omega             - Matrix with rows = Number of LV and columns
%                       containing the projection of the weight vector v of
%                       the LV on Y. 
%
%
%   Version: 2018-04-12
%__________________________________________________________________________
%
% Written by David Popovic 
% Email: david.popovic@med.uni-muenchen.de
%

%% 1. Preparation of SPLS and assignment of matrices
% Please make sure that you have added the path for the folder containing
% the following function: spls.m, dp_spls.m, dp_partition_holdout.m,
% dp_perm.m
% standardize the data matrices X and Y for usage in SPLS. X_original and
% Y_original will later be used for projection onto the final SPLS latent
% space spanned by u and v

load(variables_for_analysis);

for i=1:size(folders_for_analysis,2)
    switch folders_for_analysis{1,i}
        case 'analysis_folder'
            analysis_folder = folders_for_analysis{2,i};
        case 'prep_folder'
            prep_folder = folders_for_analysis{2,i};
        case 'spls_toolbox_path'
            spls_toolbox_path = folders_for_analysis{2,i};
    end
end

for i=1:size(parameters_for_analysis,2)
    switch parameters_for_analysis{1,i}
        case 'hp'
            hp = parameters_for_analysis{2,i};
        case 'tp'
            tp = parameters_for_analysis{2,i};            
        case 'B'
            B = parameters_for_analysis{2,i};
        case 'W'
            W = parameters_for_analysis{2,i};
        case 'K'
            K = parameters_for_analysis{2,i};
        case 'RP'
            RP = parameters_for_analysis{2,i};
        case 'vmem_req'
            vmem_req = parameters_for_analysis{2,i};
        case 'max_sim_jobs'
            max_sim_jobs = parameters_for_analysis{2,i};
        case 'queue_name'
            queue_name = parameters_for_analysis{2,i};
        case 'vmem_req'
            vmem_req = parameters_for_analysis{2,i};
    end
end


X_original = zscore(X);
Y_original = zscore(Y);

% define the data matrices X and Y for input into the SPLS algorithm
X = X_original;
Y = Y_original;

% define the default values for variables tp, hp, B, W, K if no value for
% them is entered

% hp ratio
if ~exist('hp', 'var')
    hp = 10;
end

% tp ratio
if ~exist('tp', 'var')
    tp = 20;
end

% B
if ~exist('B', 'var')
    B = 10000;
end

% W
if ~exist('W', 'var')
    W = 10;
end

% K
if ~exist('K', 'var')
    K = 100;
end

% RP
if ~exist('RP', 'var')
    RP = 40;
end

%% 2. LV loop

% this loop is repeated for the next vector pair u and v as long as the
% previous vector pair was significant. If the previous vector was not
% significant, then the loop will stop

% pLV is a starting p value set to 0, so that the while loop starts. Later
% on it will be updated at the end of each LV iteration. Then it serves the
% purpose to stop the while loop as soon as an LV was not significant.
p_LV = 0; 

% FWE_rate corrects for multiple testing for the W amount of p values
FWE_rate = 0.05 / W; 

% ff counts the iteration through the LV loops, so that the results for each
% LV can be stored in separate rows
ff = 1;  

% cu_range and cv_range define the range for the grid search. The grid
% search is performed along 40 points from 1 to sqrt(number of variables)
% as proposed in Monteiro et al. 2016.
cu_range_points = linspace(1,sqrt(size(X,2)),RP);
cv_range_points = linspace(1,sqrt(size(Y,2)),RP);

% compile a matrix with separate row for all possible cu and cv
% combinations by taking cu and repeating every single element 40 times and
% then takin cv and repeating the entire vector 40 times
cu_range = repelem(cu_range_points,numel(cu_range_points));
cv_range = repmat(cv_range_points,1,numel(cv_range_points));
cu_cv_combination = [cu_range;cv_range]; % this matrix contains all the possible combinations between cu and cv

% define column names for matrices so that you can access them later by
% indexing
final_parameters_names = {'w', 'cu', 'cv', 'u', 'v', 'success', 'RHO', 'p'};
opt_parameters_names = {'w', 'cu', 'cv', 'u', 'v', 'success', 'RHO', 'p'}; % names of parameters for optimal cu/cv combination of the w-th loop
cu_cv_combination_names = {'cu';'cv'};

% preallocate placeholder matrices for higher speed
final_parameters = num2cell(nan(size(Y,2), numel(final_parameters_names))); % cell array to store final parameters for each LV iteration
epsilon = nan(size(final_parameters,1),size(X_original,1)); %projection of weight vector u on corresponding matrix X
omega = nan(size(final_parameters,1),size(Y_original,1)); % projection of weight vector v on corresponding matrix Y

% indices for later use
index_u = strcmp(opt_parameters_names,'u');
index_v = strcmp(opt_parameters_names,'v');
index_cu = strcmp(cu_cv_combination_names,'cu');
index_cv = strcmp(cu_cv_combination_names,'cv');
index_RHO = strcmp(opt_parameters_names, 'RHO');
index_p = strcmp(opt_parameters_names, 'p');  

% set up the parallel environment
if ~exist([analysis_folder '/hyperparameter'],'dir')
    mkdir([analysis_folder '/hyperparameter'])
    hyperparameter_folder = [analysis_folder '/hyperparameter'];
else
    hyperparameter_folder = [analysis_folder '/hyperparameter'];
end

if ~exist([analysis_folder '/permutation'],'dir')
    mkdir([analysis_folder '/permutation'])
    permutation_folder = [analysis_folder '/permutation'];
else
    permutation_folder = [analysis_folder '/permutation'];
end

for i=1:size(cu_range,2)
    FID = fopen([hyperparameter_folder '/cu_' num2str(i) '.txt'],'w');
    fprintf(FID,'%.4f',cu_range(i));
    fclose(FID);
end

for i=1:size(cv_range,2)
    FID = fopen([hyperparameter_folder '/cv_' num2str(i) '.txt'],'w');
    fprintf(FID,'%.4f',cv_range(i));
    fclose(FID);
end

% here starts the outer loop for each single LV, the next LV is only
% computed if the previous LV was significant (with FWE correction)
while p_LV <= FWE_rate
 
%% 3. Wide loop
% repeats the entire process W times to obtain W different final p values,
% which go into the omnibus hypothesis

% matrix to store the optimal parameters of the w-th loop
opt_parameters = num2cell(nan(W,numel(opt_parameters_names)));

    for w=1:W

        %% hyper-parameter optimisation
        % remove hp% of the data randomly and keep it in a hold-out
        % dataset, the rest is called the keep_in data for training/testing
        % the algorithm
        [hold_out_data_x, keep_in_data_x, hold_out_data_y, keep_in_data_y] = dp_partition_holdout(hp, X, Y);

        % save the hyperparameter information for optimization step
        save([hyperparameter_folder '/cucv.mat'],...
            'keep_in_data_x',...
            'keep_in_data_y',...
            'tp',...
            'K',...,
            'hyperparameter_folder');
                
        % write a bash script for hyperparameter optimization, which can
        % later be called
        no_jobs = size(cu_cv_combination,2);
        hyperopt_bash = dp_bash_para(spls_toolbox_path, hyperparameter_folder, 'hyperopt', vmem_req, no_jobs, max_sim_jobs, queue_name);
        
        % access the temp_RHO folder and initialize the bash script for
        % parallel computing of the hyperparameter optimization
        cd(hyperparameter_folder);
        system(['qsub ' hyperopt_bash]);

        % set a variable which counts the files in the RHO folder
        filecount = size(dir([hyperparameter_folder '/RHO_avg_*']),1);
        RHO_avg_collection = nan(size(cu_cv_combination,2),1);

        % use a while loop which recounts all the files in the RHO folder
        % until all computations from the parallelized hyperparameter
        % optimization step are completed and all files are saved in that
        % folder, the while loop will be exited as soon as all files are
        % saved in the folder
        while filecount < size(cu_cv_combination,2)
            filecount = size(dir([hyperparameter_folder '/RHO_avg_*']),1);
        end
        
        % load all RHO values which were saved in mat files into a
        % RHO_collection vector
        for i=1:size(cu_cv_combination,2)
            if exist([hyperparameter_folder '/RHO_avg_' num2str(i) '.mat'])
                load([hyperparameter_folder '/RHO_avg_',num2str(i),'.mat']);
                RHO_avg_collection(i,1) = RHO_avg;
                delete([hyperparameter_folder '/RHO_avg_',num2str(i),'.mat']); 
            else
                RHO_avg_collection(i,1) = NaN;
            end
        end


%     select the hyperparameter combination with the highest average
%     correlation, ie the one with the highest RHO_avg
        [~, index_RHO] = max(RHO_avg_collection); 
        
%     find the corresponding cu and cv values of the max RHO element, this
%     gives you the cu/cv combination with the highest average correlation
    
    cu_opt = cu_cv_combination(index_cu,index_RHO);
    cv_opt = cu_cv_combination(index_cv,index_RHO);

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

    % access the temp_B folder to store all RHO_b values which will be
    % computed in parallel in the next step
    
    for i=1:B
        temp = randperm(size(keep_in_data_y,1))';
        dlmwrite([permutation_folder '/perm_' num2str(i) '.txt'],temp);
    end
    
    % write a bash script for hyperparameter optimization, which can
    % later be called
    permutation_bash = dp_bash_para(spls_toolbox_path, permutation_folder, 'permutation', vmem_req, B, max_sim_jobs, queue_name);
        
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
        filecount = size(dir([permutation_folder '/RHO_b_*']),1);
    end
    
    for i=1:size(RHO_b_collection,1)
        if exist([permutation_folder '/RHO_b_' num2str(i) '.mat'],'file')==2
            load([permutation_folder '/RHO_b_',num2str(i),'.mat']);
            RHO_b_collection(i,1) = RHO_b;
            delete([permutation_folder '/RHO_b_',num2str(i),'.mat']);
        else
            RHO_b_collection(i,1) = NaN;
        end
    end
    
    % test the following null hypothesis H_s: "There is no relationship
    % between the two views, therefore the correlation obtained with
    % the original data is not different from the correlation
    % obtained with the permuted data"

    % calculate how many times RHO_b was bigger than or equal to RHO_opt,
    % ie how many times the permutated data led to a higher correlation
    % using spls than the original data
    RHO_count_b = sum(RHO_b_collection >= RHO_opt);

    % calculate the p value to test whether the null hypothesis is true
    p = ((1+RHO_count_b)/(B+1));
  
    % store all parameters of the w-th wide loop into a matrix
    opt_parameters(w,:) = {w cu_opt cv_opt u_opt v_opt success_opt RHO_opt p};
    
    
    end
    
    save([analysis_folder '/opt_parameters_' num2str(ff) '.mat'],'opt_parameters'); 
    
    %% test the omnibus hypothesis using the p values

    % Omnibus hypothesis H_omni: "All the null hypotheses H_s are true" if
    % any of the p values are lower than the FWE-corrected threshold,
    % then the omnibus hypothesis can be rejected after that. 
    
    % Search for the statistically significant w-th iteration with the
    % lowest p value, if more than one iteration has the same p value,
    % select the one that has the highest RHO value => this is the final
    % w-th iteration with the corresponding cu/cv features and u/v scores
    
    % remove parameter combinations with nans as RHO_opt

    if exist('opt_parameters_clean','var')
        clear opt_parameters_clean
    end
    
    nn=1;
    for i=1:size(opt_parameters,1)
        if ~isnan(cell2mat(opt_parameters(i,index_RHO)))
            opt_parameters_clean(nn,:) = opt_parameters(i,:);
            nn=nn+1;
        end
    end

    %create a logical indexing the lowest p-value(s)
    p_min = min(cell2mat(opt_parameters_clean(:,index_p))); 
    logical_p_min = cell2mat(opt_parameters_clean(:,index_p)) == p_min;
    temp1 = opt_parameters_clean(logical_p_min,:); % placeholder to make syntax easier
    
    if size(temp1,1) > 1 % if there is more than one p-value that is equal to the lowest p-value, ie there are two absolute minimal p-values, then continue in the loop
        RHO_max = cell2mat(temp1(:,index_RHO)) == max(cell2mat(temp1(:,index_RHO))); % logical indexing the highest RHO value
        final_parameters(ff,:) = temp1(RHO_max,:); % take the row with the lowest p-value and the highest RHO value
    else
        final_parameters(ff,:) = temp1; % if there is just one absolute minimal p value, then take this row
    end

    %% Matrix Deflation => Projection Deflation

    % This method removes the covariance explained by u and v from X
    % and Y by projecting the data matrices onto the space spanned by the
    % corresponding weight vector, and subtracting this from the data:

    u = final_parameters{ff,index_u}; % weight vector for X
    v = final_parameters{ff,index_v}; % weight vector for Y
    X = X - (X*u)*u'; % the effect of u on X is removed by projection deflation
    Y = Y - (Y*v)*v'; % the effect of v on Y is removed by projection deflation
    
    % the p value of this LV is updated from the previously set value of
    % zero to the actual p value of the LV, if this p value is lower than
    % the FWE_rate then the while loop continues after matrix deflation and
    % it keeps looking for the next possible significant LV. If the p value
    % is higher than the FWE_rate, then the while loop is discontinued and
    % the algorithm stops. Therefore, this function generates all possible
    % significant LVs and also the first non-significant LV, afterwards the
    % algorithm stops.
    p_LV = final_parameters{ff,index_p}; 
    
    
    %% Projection of the original data onto the SPLS latent space
    
    % epsilon and omega are the projections of the original data onto the
    % SPLS latent space, spanned by the weight vector pair u and v of this
    % LV. Therefore each LV spans a specific latent space onto wihich the
    % samples can be projected. These latent spaces can then be used to
    % stratify patients.
    epsilon(ff,:) = X_original * final_parameters{ff,index_u};
    omega(ff,:) = Y_original * final_parameters{ff,index_v};

    % ff counts through the most outside LV while loop = counts the amount
    % of significant LVs
    ff = ff+1; 
    

end

%% clean up txt files and other files which are not needed anymore
files_to_clean = {'cu_', 'cv_'};
folder_to_clean = hyperparameter_folder;
for i=1:size(files_to_clean,2)
    dp_cleanup(folder_to_clean, files_to_clean{i})
end

% after the LVs are computed, clean up the final parameters, epsilon and
% omega matrix by removing empty rows
final_parameters(ff:size(final_parameters,1),:) = [];
epsilon(ff:size(final_parameters,1),:) = [];
omega(ff:size(final_parameters,1),:) = [];

save([analysis_folder '/result.mat'], 'final_parameters', 'epsilon', 'omega');

end
