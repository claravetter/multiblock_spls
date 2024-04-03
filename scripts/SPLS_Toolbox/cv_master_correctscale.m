
function [OUT] = cv_master_correctscale(IN, COV, cs_method, correction_target)



if isfield(IN, 'train')
    for i = 1:length(IN.train)
        IN_temp.train = IN.train{i};
        IN_temp.test = IN.test{i};

        if correction_target(i) == 1
            OUT_temp = dp_correctscale(IN_temp,COV,cs_method);
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
        if correction_target(i) == 1
            OUT_temp = dp_correctscale(IN_temp,COV,cs_method);
            OUT{i} = OUT_temp;    
        else
            OUT{i} = IN_temp;
        end

    end
end



end