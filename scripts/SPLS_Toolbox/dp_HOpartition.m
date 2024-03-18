%% DP function for hold-out partitions

function partition = dp_HOpartition(n, labels)

for w=1:n
    temp_cv_test=[]; temp_cv_train=[];
    for i=1:size(unique(labels),1)
        temp = find(labels == i);
        temp_test = randperm(size(temp,1),round((size(temp,1)/n)));
        temp_cv_test = [temp_cv_test; sortrows(temp(temp_test))]; 
        temp(temp_test) = [];
        temp_cv_train = [temp_cv_train; temp];
    end
    partition.TestInd{1,w} = temp_cv_test;
    partition.TrainInd{1,w} = temp_cv_train;
end

end