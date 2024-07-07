%% DP function for projection of u and v onto test data
function [RHO, p, epsilon, omega, u, v] = dp_projection_ext(data_x, data_y, u, v, correlation_method)

epsilon = data_x*u;
omega = data_y*v;
[RHO, p] = corr(epsilon, omega, 'Type', correlation_method);

f_invert = @(x)(-1*x);

if RHO<0
    v = f_invert(v);
    omega = data_y*v;
    [RHO, p] = corr(epsilon, omega, 'Type', correlation_method);
end

if sum(v)<0
    u = f_invert(u);
    v = f_invert(v);
    epsilon = data_x*u;
    omega = data_y*v;
    [RHO, p] = corr(epsilon, omega, 'Type', correlation_method);
end

end