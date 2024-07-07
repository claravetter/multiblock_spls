%% DP function for bootstrapping from opt_parameters
function [ci_collection_temp, bootstat_collection_temp] = dp_bootstrap_opt(opt_matrix, number_bootstrap, target)

switch target
    case 'u'
        n=4;
    case 'v'
        n=5;
end

temp=opt_matrix(:,n); temp_ex=[]; ci_collection_temp=[]; bootstat_collection_temp=[];
for i=1:size(temp,1)
    for ii=1:size(opt_matrix{1,n},1)
        temp_ex(i,ii) = opt_matrix{i,n}(ii);
    end
end

nn=1;
for i=1:size(temp_ex,2)
    [temp_ci, temp_bootstat] = bootci(number_bootstrap, @(x) mean(x), temp_ex(:,i));
    ci_collection_temp(i,1:2) = temp_ci;
    bootstat_collection_temp(i,1) = mean(temp_bootstat);
    nn=nn+1;
    if mod(nn,1000)==0
        disp([num2str(nn) ' iterations reached']);
    end
end

end