%% DP function for scaling and deflating

function [IN_x, IN_y] = dp_deflatescale(IN_x, IN_y, u_collection, v_collection, scaling_method)
IN_s.method = scaling_method;
for dd=1:size(u_collection,2)
    [IN_x.train, IN_x.test] = dp_standardize_comb(IN_x.train, IN_x.test, IN_s);
    [IN_y.train, IN_y.test_s] = dp_standardize_comb(IN_y.train, IN_y.test, IN_s);
    
    [IN_x.train,IN_y.train] = proj_def(IN_x.train, IN_y.train, u_collection(:,dd), v_collection(:,dd));
    [IN_x.test,IN_y.test] = proj_def(IN_x.test, IN_y.test, u_collection(:,dd), v_collection(:,dd));
end

end