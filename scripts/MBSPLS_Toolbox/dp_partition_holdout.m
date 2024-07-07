%% function for hold-out splits

function [hold_out_X, keep_in_X, hold_out_Y, keep_in_Y, hold_out_Z, keep_in_Z] = dp_partition_holdout(percentage_holdout, data_X, data_Y, data_Z)

if size(data_X,1) == size(data_Y,1)
    data = [data_X, data_Y];
else
    fprintf('The number of samples/rows in the matrices dont match!');
end

if exist('data_Z','var')
    if size(data_Z,1)==size(data,1)
        data=[data, data_Z];
    else
        fprintf('The number of samples/rows in the matrices dont match!');
    end
end
    
a = size(data,1); % amount of samples
aa = 1:size(data,1); %vector counting from 1 to a number of samples
hold_out = sort(randperm(a, round(percentage_holdout/100*a))); % create percentage_holdout/100*a random integers between 1 and a for the houldout dataset
keep_in = aa(~ismember(aa,hold_out)); % keep_in data which is used for training/test = the rest of the data

% index the data matrix and extract the hold_out and the keep_in data
hold_out_data = data(hold_out,:);
keep_in_data = data(keep_in,:);

%create hold_out and keep_in data sets for X and Y
hold_out_X = hold_out_data(:, 1:size(data_X,2));
hold_out_data(:, 1:size(data_X,2)) = [];
keep_in_X = keep_in_data(:, 1:size(data_X,2));
keep_in_data(:, 1:size(data_X,2)) = [];

hold_out_Y = hold_out_data(:, 1:size(data_Y,2));
hold_out_data(:, 1:size(data_Y,2)) = [];
keep_in_Y = keep_in_data(:, 1:size(data_Y,2));
keep_in_data(:, 1:size(data_Y,2)) = [];

if exist('data_Z','var')
    hold_out_Z = hold_out_data(:, 1:size(data_Z,2));
    hold_out_data(:, 1:size(data_Z,2)) = [];
    keep_in_Z = keep_in_data(:, 1:size(data_Z,2));
    keep_in_data(:, 1:size(data_Z,2)) = [];
end


end

