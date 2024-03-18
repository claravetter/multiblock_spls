%% DP function for one k split

function [RHO, u, v, epsilon, omega] = dp_ICV_spls(training_data_x,training_data_y,test_data_x, test_data_y, cu, cv, correlation_method)

%perform SPLS on the training data using the current cu/cv combination
[u_temp, v_temp, ~] = spls_suppressed_display(training_data_x,training_data_y,cu,cv);

%compute the correlation between the projections of the training and
%test matrices onto the SPLS latent space spanned by the weight vectors
epsilon_temp = test_data_x*u_temp;
omega_temp = test_data_y*v_temp;
RHO_temp = corr(epsilon_temp, omega_temp, 'Type', correlation_method);

f_invert = @(x)(-1*x);

if RHO_temp<0
    RHO = f_invert(RHO_temp);
    u = u_temp;
    v = f_invert(v_temp);
    epsilon = epsilon_temp;
    omega = test_data_y*v;
else
    RHO = RHO_temp;
    u = u_temp;
    v = v_temp;
    epsilon = epsilon_temp;
    omega = omega_temp;
end

end