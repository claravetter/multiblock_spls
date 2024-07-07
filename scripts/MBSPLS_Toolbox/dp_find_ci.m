%% function to find CI of questionnaire items

function CI = dp_find_ci(data_vector)

temp=opt_parameters(:,5);
for i=1:size(temp,1)
    for ii=1:size(opt_parameters{1,5},1)
        temp_ex(i,ii) = opt_parameters{i,5}(ii);
    end
end

% mean_temp_ex = btstrp;
for i=1:size(temp_ex,2)
    [temp_ci, temp_bootstat] = bootci(1000, @(x) mean(x), temp_ex(:,i));
    ci_collection(1:2,i) = temp_ci;
    bootstat_collection(1,i) = mean(temp_bootstat);
end

errors=[ci_collection(1,:)-bootstat_collection;ci_collection(2,:)-bootstat_collection];


end