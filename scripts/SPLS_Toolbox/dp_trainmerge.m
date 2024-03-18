%% DP function for training/retraining and merging training results

function [RHO_opt, u_opt, v_opt, V_opt] = dp_trainmerge(RHO_collection, u_collection, weights_u_collection, v_collection, weights_v_collection, V_collection, weights_V_collection, type_merge)

switch type_merge
    case 'mean'
        RHO_opt = mean(RHO_collection);
        u_opt = mean(u_collection,2);
        v_opt = mean(v_collection,2);
        V_opt = mean(V_collection,3);
        if V_opt==0
            V_opt = false;
        end
    case 'median'
        RHO_opt = median(RHO_collection);
        u_opt = median(u_collection,2);
        v_opt = median(v_collection,2);
        V_opt = median(V_collection,3);
    case 'weighted_mean'
        u_opt = wmean(u_collection, weights_u_collection, 2)';
        v_opt = wmean(v_collection, weights_v_collection, 2)';
        try V_opt = wmean(V_collection, weights_V_collection, 3);
        catch
            V_opt = false;
        end
        RHO_opt = wmean(RHO_collection, RHO_collection);
    case 'best'
        [RHO_opt,I_max] = max(RHO_collection);
        u_opt = u_collection(:,I_max);
        v_opt = v_collection(:,I_max);
        try V_opt = V_collection(:,:,I_max);
        catch
            V_opt = false;
        end
end

end