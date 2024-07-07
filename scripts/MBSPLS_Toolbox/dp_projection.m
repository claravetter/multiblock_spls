%% DP function for projection of u and v onto test data
function [RHO, epsilon, omega, u, v] = dp_projection(data_x, data_y, u, v, correlation_method)

epsilon = data_x*u;
omega = data_y*v;
RHO = corr(epsilon, omega, 'Type', correlation_method);

f_invert = @(x)(-1*x);

if RHO<0 % if correlation negative, invert one weight vector --> correlation positive
    v = f_invert(v);
    omega = data_y*v;
    RHO = corr(epsilon, omega, 'Type', correlation_method);
end

if sum(v)<0 % easier interpretation; v = phenotypic vector
    u = f_invert(u);
    v = f_invert(v);
    epsilon = data_x*u;
    omega = data_y*v;
    RHO = corr(epsilon, omega, 'Type', correlation_method);
end

end