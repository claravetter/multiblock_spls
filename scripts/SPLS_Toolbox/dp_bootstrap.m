%% test
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

barwitherr(errors',bootstat_collection');
hold on


test1=randi(size(test,2),1000,1);
for i=1:size(test1,1)
    btstrp(i) = test(test1(i));
end

mean_bts = mean(btstrp);
std_bts = std(btstrp);
K_int = 1.96*std_bts/(sqrt(size(test1,1)))

%% test every non-zero item of u and v for significance using bootstrapping and 95%KI

% first use vector v to find inverted variables
[v_collect, inv_log_v] = dp_inverse_array(opt_parameters(:,opt_v));
opt_parameters(:,opt_v) = v_collect;

% then use vector u and invert the opt parameters iterations like in v
temp_inv = []; temp_inv_new=[];
temp_inv = opt_parameters(inv_log_v, opt_u);
for i=1:size(temp_inv,1)
    temp_inv_new {i,1} = temp_inv{i}.*(-1);
end

opt_parameters(inv_log_v, opt_u) = temp_inv_new;


% compute CI for u and v
temp=[];
for i=1:size(opt_parameters,1)
    temp(:,i) = opt_parameters{i,opt_v};
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

temp=[];
for i=1:size(opt_parameters,1)
    temp(:,i) = opt_parameters{i,opt_u};
end

KI_u=[];
for i=1:size(temp,1)
    temp1 = dp_bt_ki(temp(i,:),1000);
    if sum(temp1>0)==2 || sum(temp1<0)==2
        temp2=1;
    else
        temp2=0;
    end
    KI_u(i,:)=[temp1, temp2];
end

KI_collection(ff,:) = {KI_u, KI_v};
KI_collection_names = {'KI_u', 'KI_v'};