%% extended significance testing of optimal parameters
y=[];
for i=1:size(opt_parameters,1)
    y=[y,opt_parameters{i,5}(49,1)];
end

bt_stats = bootstrp(1000, @(x)[mean(x), std(x)],y);
% stats = bootstrp(1000, @(x)[x],y);
% test=reshape(stats,1,(size(stats,1)*size(stats,2)));
bt_mean = mean(stats(:,1));
bt_std = mean(stats(:,2));

KI = [bt_mean-1.96*bt_std/(sqrt(size(stats,1))), bt_mean+1.96*bt_std/(sqrt(size(stats,1)))];

