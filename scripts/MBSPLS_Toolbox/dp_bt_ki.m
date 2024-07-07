%% extended significance testing of optimal parameters
function KI = dp_bt_ki(y,n)
n=1000;
y=1:size(opt_parameters,1);
bt_stats = bootstrp(n, @(x)[dp_empty(x)],y);
matrix_names = parameters_names;
opt_u = strcmp(matrix_names, 'u');
opt_v = strcmp(matrix_names, 'v');
best_collection = num2cell(nan(n,size(opt_parameters,2)));

[v_collect, inv_log_v] = dp_inverse_array(opt_parameters(:,opt_v));
opt_parameters(:,opt_v) = v_collect;

for i=1:size(bt_stats,1)
    temp_matrix = opt_parameters(bt_stats(i,:),:);
    [best_collection(i,:), ~] = dp_find(temp_matrix, matrix_names);
end

% then use vector u and invert the opt parameters iterations like in v
temp_inv = []; temp_inv_new=[];
temp_inv = opt_parameters(inv_log_v, opt_u);
for i=1:size(temp_inv,1)
    temp_inv_new {i,1} = temp_inv{i}.*(-1);
end

opt_parameters(inv_log_v, opt_u) = temp_inv_new;

% compute CI for u and v
temp_v=[];
for i=1:size(best_collection,1)
    temp_v(i,:) = best_collection{i,opt_v};
end

mean_v = mean(temp_v,1);
std_v = std(temp_v,1);
percentiles_v = prctile(temp_v, [2.5, 97.5],1);

v_zero = output.final_parameters{1,5} == 0;
percentiles_v(:,v_zero) = 0;

for i=1:size(temp_v,2)
    f=figure();
    histogram(temp_v(:,i));
    saveas(f,['figure_' num2str(i), '.png']);
end

KI_v=[];
for i=1:size(temp,1)
    temp1 = dp_bt_ki(temp(i,:),1000);
    if sum(temp1>0)==2 || sum(temp1<0)==2
        temp2=1;
    else
        temp2=0;
    end
    KI_v(i,:)=[temp1, temp2];
end


bt_mean = mean(bt_stats(:,1));
bt_std = mean(bt_stats(:,2));

KI = [bt_mean-1.96*bt_std/(sqrt(size(bt_stats,1))), bt_mean+1.96*bt_std/(sqrt(size(bt_stats,1)))];

end