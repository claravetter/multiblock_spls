%% prelim data

log_nan = isnan(cell2mat(output.final_parameters(:,1)));

output.final_parameters = output.final_parameters(~log_nan,:);

save('prelim.mat', 'input', 'output', 'setup');

IN.results_path = '/volume/HCStress/Analysis/07-Nov-2018/temporary/prelim.mat';
IN.overall_analysis = 'Stress';
IN.specific = {'atlas'}; % 'behavior', 'atlas', 'images', 'detailed'
dp_visualize_data(IN);

%% testing for deflation

X=randi(100,10,20);
Y=randi(100,30,20);
Z=[X;Y];
u=randi(5,20,1);

X_new = X - (X*u)*u';
Y_new = Y - (Y*u)*u';
Z_new = Z - (Z*u)*u';

sum(sum(Z_new-[X_new;Y_new]))

