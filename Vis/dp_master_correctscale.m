%% DP master function for correcting and scaling matrices

function [OUT_x, OUT_y] = dp_master_correctscale(IN_x, IN_y, COV, cs_method, correction_target)

switch correction_target
    case 1
        OUT_x = dp_correctscale(IN_x,COV,cs_method);
        try COV.test = nan(size(COV.test,1),1);
            COV.train = nan(size(COV.train,1),1);
        catch
            COV = nan(size(COV,1),1);
        end
        OUT_y = dp_correctscale(IN_y,COV,cs_method);
    case 2
        OUT_y = dp_correctscale(IN_y,COV,cs_method);
        try COV.test = nan(size(COV.test,1),1);
            COV.train = nan(size(COV.train,1),1);
        catch
            COV = nan(size(COV,1),1);
        end
        OUT_x = dp_correctscale(IN_x,COV,cs_method);
    case 3
        OUT_x = dp_correctscale(IN_x,COV,cs_method);
        OUT_y = dp_correctscale(IN_y,COV,cs_method);
end

end