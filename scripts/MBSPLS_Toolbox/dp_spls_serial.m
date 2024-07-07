function [final_parameters, epsilon, omega] = dp_spls_serial(X, Y, grid_density, analysis_folder)
%
%   Implementation of SPLS algorithm and framework as described in Monteiro
%   et al. 2016, for details: doi:10.1016/j.jneumeth.2016.06.011
%   
% Inputs:
%
%     X, Y    - data matrices for MRI images (X) and behavioral data (Y) in
%               the form: samples x features. Use MRI data directly from
%               Neurominer output, the function will standardize both the
%               MRI data and the behavioraldata, so that mean = 0 and std =
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
%   Version: 2018-03-05
%__________________________________________________________________________

% Written by David Popovic 
% Email: david.popovic@med.uni-muenchen.de

%% 1. Preparation of SPLS and assignment of matrices
% Please make sure that you have added the path for the folder containing
% the following function: spls.m, dp_spls.m, dp_partition_holdout.m,
% dp_perm.m

% standardize the data matrices X and Y for usage in SPLS. X_original and
% Y_original will later be used for projection onto the final SPLS latent
% space spanned by u and v

X_original = zscore(X);
Y_original = zscore(Y);

% define the data matrices X and Y for input into the SPLS algorithm
X = X_original;
Y = Y_original;

% define the default values for variables tp, hp, B, W, K if no value for
% them is entered

% hp ratio
hp = 10;
tp = 20;
B = 5000;
W = 10;
K = 100;

%% 2. LV loop

% this loop is repeated for the next vector pair u and v as long as the
% previous vector pair was significant. If the previous vector was not
% significant, then the loop will stop

% pLV is a starting p value set to 0, so that the while loop starts. Later
% on it will be updated at the end of each LV iteration. Then it serves the
% purpose to stop the while loop as soon as an LV was not significant.
p_LV = 0; 

% FDR_rate corrects for multiple testing for the W amount of p values
FDRvalue = 0.05; % significance threshold for FDR testing, Default: 0.05
pvalue_FDR = 1;

% ff counts the iteration through the LV loops, so that the results for each
% LV can be stored in separate rows
ff = 1;  

% cu_range and cv_range define the range for the grid search. The grid
% search is performed along 40 points from 1 to sqrt(number of variables)
% as proposed in Monteiro et al. 2016.
cu_range_points = linspace(1,sqrt(size(X,2)),grid_density);
cu_range_points(1)=[]; cu_range_points(end)=[];
cv_range_points = linspace(1,sqrt(size(Y,2)),grid_density);
cv_range_points(1)=[]; cv_range_points(end)=[];

% compile a matrix with separate row for all possible cu and cv
% combinations by taking cu and repeating every single element 40 times and
% then takin cv and repeating the entire vector 40 times
cu_range = repelem(cu_range_points,numel(cu_range_points));
cv_range = repmat(cv_range_points,1,numel(cv_range_points));

spls_output_names = {'u', 'v', 'success'}; % output of spls function
spls_output_k_names = {'u','v','success', 'RHO'}; % output of each k iteration
final_parameters_names = {'w', 'cu', 'cv', 'u', 'v', 'success', 'RHO', 'p'};
final_parameters = num2cell(zeros(size(Y,2), numel(final_parameters_names))); % cell array to store final parameters for each LV iteration

% names of parameters for optimal cu/cv combination of the w-th loop
opt_parameters_names = {'w', 'cu', 'cv', 'u', 'v', 'success', 'RHO', 'p'};

% indices for later use
opt_u = strcmp(opt_parameters_names,'u');
opt_v = strcmp(opt_parameters_names,'v');
opt_RHO = strcmp(opt_parameters_names,'RHO');
opt_p = strcmp(opt_parameters_names, 'p');  

while p_LV <= pvalue_FDR

    
%% 3. Wide loop
% repeats the entire process W times to obtain W different final p values,
% which go into the omnibus hypothesis

% matrix to store the optimal parameters of the w-th loop
opt_parameters = num2cell(zeros(W,numel(opt_parameters_names)));

    for w=1:W

        %% hyper-parameter optimisation
        % remove hp% of the data randomly and keep it in a hold-out
        % dataset, the rest is called the keep_in data for training/testing
        % the algorithm
        
        [hold_out_data_x, keep_in_data_x, hold_out_data_y, keep_in_data_y] = dp_partition_holdout(hp, X, Y);

        % create an inner loop with k splits and a tp test perecentage set
        % within each split

        % RHO_avg is an l-by-ll matrix where all average RHO values, where
        %the 1600 cu/cv combinations are stored
        RHO_avg = zeros(numel(cu_range)); 
        
        for l=1:numel(cu_range)
            cu = cu_range(l);
            for ll=1:numel(cv_range)
                cv = cv_range(ll);
                spls_output_k_results = num2cell(zeros(K, numel(spls_output_names)));
                for k=1:K
                    % separate the keep_in data into training and test data
                    % according to the chosen test percentage tp
                    [test_data_x, training_data_x, test_data_y, training_data_y] = dp_partition_holdout(tp, keep_in_data_x, keep_in_data_y);
                    
                    %perform SPLS on the training data using the current cu/cv combination
                    [u,v,success] = spls_suppressed_display(training_data_x,training_data_y,cu,cv); 
                    
                    %compute the correlation between the projections of the
                    %training and test matrices onto the SPLS latent space
                    %spanned by the weight vectors
                    RHO = abs(corr(test_data_x*u,test_data_y*v)); 
                    
                    %store the results into a temporary cell array
                    temp = {u v success RHO}; 
                    
                    %store the relevant results into a output folder for
                    %the k-th iteration
                    for i=1:numel(temp)
                        spls_output_k_results{k,i} = temp{i};
                    end
                end
                %index for RHO value
                index_RHO = strcmp(spls_output_k_names,'RHO'); 
                
                %computing the mean value of all RHO calculated for the k
                %iterations for one specific cu/cv combination
                RHO_avg(l,ll) = mean(cell2mat(spls_output_k_results(:,index_RHO))); 
            end
        end


        % select the hyperparameter combination with the highest average
        % correlation, ie the one with the highest RHO_avg
        
        %convert RHO matrix into a column vector and then take the max
        %element from that vector, take the M value and the I index of this
        %element
        [~, I_RHO] = max(RHO_avg(:)); 
        
        %get the row and the column of this element in the RHO matrix
        [I_RHO_row, I_RHO_column] = ind2sub(size(RHO_avg),I_RHO); 

        % find the corresponding cu and cv values for the rows and vectors
        % of the max RHO element, this gives you the cu/cv combination with
        % the highest average correlation

        cu_opt = cu_range(I_RHO_row);
        cv_opt = cv_range(I_RHO_column);


        %% Statistical Evaluation

        %train the model with the optimal cu/cv combination on the previous
        %train/test data set to get u and v
        [u_opt, v_opt, success_opt]=spls(keep_in_data_x, keep_in_data_y,cu_opt, cv_opt);

        %project the hold_out data set on the computed u and v vectors to
        %get the absolute correlations between these projections
        RHO_opt = abs(corr(hold_out_data_x*u_opt,hold_out_data_y*v_opt));

        % permute the keep_in data in one dimension so that the
        % relationship between these two views gets destroyed

        % collection of RHO values for the B times permutated dataset 
        RHO_b_collection = zeros(B,1); 
        
        for b=1:B

            % use the permuted perm_data_y and the regular keep_in_data_x
            % and compute u_b and v_b for the permutated data using spls
            [u_b, v_b, ~]=spls_suppressed_display(keep_in_data_x, dp_perm(keep_in_data_y), cu_opt, cv_opt);

            % compute the absolute correlation between the hold_out data
            % and the permuted u_b and v_b
            RHO_b = abs(corr(hold_out_data_x*u_b,hold_out_data_y*v_b));

            % store the calculated RHO_b value in a vector
            RHO_b_collection(b,1) = RHO_b;
        end

        % test the following null hypothesis H_s: "There is no relationship
        % between the two views, therefore the correlation obtained with
        % the original data is not different from the correlation
        % obtained with the permuted data"

        % calculate how many times RHO_b was bigger than RHO_opt, ie how
        % many times the permutated data led to a higher correlation using
        % spls than the original data
        RHO_count_b = sum(RHO_b_collection >= RHO_opt);

        % calculate the p value to test whether the null hypothesis is true
        p = ((1+RHO_count_b)/(B+1));

        %store all parameters of the wide w-th wide loop into a matrix
        opt_parameters(w,:) = {w cu_opt cv_opt u_opt v_opt success_opt RHO_opt p};
        save([analysis_folder '/opt_parameters_' num2str(ff) '_' num2str(w) '.mat'],'opt_parameters', 'RHO_b_collection'); 

    end

    %% test the omnibus hypothesis using the p values

    % Omnibus hypothesis H_omni: "All the null hypothesis H_s are true" if
    % any of the W p values are lower than the FWE-corrected threshold,
    % then the omnibus hypothesis can be rejected after that, search for
    % the statistically significant W iteration with the lowest p value, if
    % more than one iteration have the same p value, select the one that
    % has the highest RHO value => this is the final W iteration with the
    % corresponding cu/cv features and u/v scores

    % remove parameter combinations with nans as RHO_opt
    to_rmv=isnan([opt_parameters{:,opt_RHO}]');
    opt_parameters(to_rmv,:)=[];
    
    % calculate FDR-corrected p value
    pvalue_FDR = dp_FDR([opt_parameters{:,opt_p}], FDRvalue);

    %create a logical indexing the lowest p-value(s)
    [~, I_pmin] = min([opt_parameters{:,opt_p}]);
    temp1 = opt_parameters(I_pmin,:); % placeholder to make syntax easier
    
    if size(temp1,1) > 1 % if there is more than one p-value that is equal to the lowest p-value, ie there are two absolute minimal p-values, then continue in the loop
        [~, I_Rhomax] = max([temp1{:,opt_RHO}]); % logical indexing the highest RHO value
        final_parameters(ff,:) = temp1(I_Rhomax,:); % take the row with the lowest p-value and the highest RHO value
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
    % the FDR_rate then the while loop continues after matrix deflation and
    % it keeps looking for the next possible significant LV. If the p value
    % is higher than the FDR_rate, then the while loop is discontinued and
    % the algorithm stops. Therefore, this function generates all possible
    % significant LVs and also the first non-significant LV, afterwards the
    % algorithm stops.
    p_LV = final_parameters{ff,opt_p}; 

end

% clean up the final parameters table and remove empty rows
final_parameters(ff:size(final_parameters,1),:) = [];


%% Projection of the weight vectors onto the SPLS latent space

epsilon = zeros(size(final_parameters,1),size(X_original,1)); %projection of weight vector u on corresponding matrix X
omega = zeros(size(final_parameters,1),size(Y_original,1)); % projection of weight vector v on corresponding matrix Y

for i = 1:size(final_parameters,1)
    epsilon(i,:) = X_original * final_parameters{i,index_u};
    omega(i,:) = Y_original * final_parameters{i,index_v};
end

save([analysis_folder '/result.mat'], 'final_parameters', 'epsilon', 'omega');

end

