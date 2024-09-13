
function [OUT] = cv_master_correctscale(IN, COV, cs_method, correction_target)

if ~isfield(IN, 'knnimpute_k')
    k = 5;
else
    k = IN.knnimpute_k;
end

if isfield(IN, 'train')
    for i = 1:length(IN.train) % number matrices
        % check if block needs imputation
        % if any(isnan(IN.train{i}))
        %     % do imputation 
        %     IN.train{i} = knnimpute(IN.train{i}, k);
        % end
        % if any(isnan(IN.test{i}))
        %     % do imputation 
        %     IN.test{i} = knnimpute(IN.test{i}, k);
        % end

        % if correction_target(i) == 1 % check if current block needs correction
        IN_temp.train = IN.train{i};
        IN_temp.test = IN.test{i};

        if size(COV,2)>1
            COV_temp = COV{i}; % select covariates for this block
        else
            COV_temp = COV;
        end

        % if any(isnan(COV_temp.train))
        %     % do imputation
        %     COV_temp.train = knnimpute(COV_temp.train, k);
        % end
        % if any(isnan(COV_temp.test))
        %     % do imputation
        %     COV_temp.test = knnimpute(COV_temp.test, k);
        % end

        if size(cs_method,2) > 1 % can be deleted later; cs_method needs to be created for each matrix elsewhere
            cs_method_temp = cs_method{i};
        else
            cs_method_temp = cs_method;
        end
        if correction_target(i) == 1
            OUT_temp = dp_correctscale(IN_temp,COV_temp,cs_method_temp);
            OUT.train{i} = OUT_temp.train;
            OUT.test{i} = OUT_temp.test;
        else
            OUT.train{i} = IN_temp.train;
            OUT.test{i} = IN_temp.test;
        end


    end
else
    for i = 1:length(IN)
        IN_temp = IN{i};

        % if any(isnan(IN_temp))
        %     % do imputation
        %     IN_temp = knnimpute(IN_temp, k);
        % end

        if size(COV,2)>1
            COV_temp = COV{i}; % select covariates for this block
        else
            COV_temp = COV;
        end

        % if any(isnan(COV_temp.train))
        %     % do imputation
        %     COV_temp.train = knnimpute(COV_temp.train, k);
        % end
        % if any(isnan(COV_temp.test))
        %     % do imputation
        %     COV_temp.test = knnimpute(COV_temp.test, k);
        % end

        if size(cs_method,2) > 1 % can be deleted later; cs_method needs to be created for each matrix elsewhere
            cs_method_temp = cs_method{i};
        else
            cs_method_temp = cs_method{1};
        end
        if correction_target(i) == 1
            OUT_temp = dp_correctscale(IN_temp,COV_temp,cs_method_temp);
            OUT{i} = OUT_temp;
        else
            OUT{i} = IN_temp;
        end

    end
end



end