function data_perm = dp_perm(data)
perm = randperm(size(data,1)); %create a random permutation of the numbers 1 to a

% use the permuted numbers to reshuffle the rows of the data set
data_perm = data(perm,:);

end