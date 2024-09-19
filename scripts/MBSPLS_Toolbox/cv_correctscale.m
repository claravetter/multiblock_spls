function OUT = cv_correctscale(DAT, COV, MET)

IN_s = struct;
IN_s.method = MET.method;
if ~isempty(MET.correction_subgroup) 
    if ~all(MET.subgroup_train == 0)
        IN_s.subgroup_train = MET.subgroup_train;
    else
        IN_s.subgroup_train = [];
    end
    if ~all(MET.subgroup_test == 0)
        IN_s.subgroup_test = MET.subgroup_test;
    else
        IN_s.subgroup_test = [];
    end
else
    IN_s.subgroup_train = [];
    IN_s.subgroup_test = [];
end

if isnumeric(DAT) % only operate on one dataset
    if ~any(isnan(COV), 'all') && ~isempty(COV)
        if any(sum(COV,1)==0) || any(range(COV)==0)
            log_remove = sum(COV,1)==0;
            COV(:,log_remove) = [];
        end
        if any(sum(DAT,1)==0) || any(range(DAT)==0)
            log_remove = double(sum(DAT,1)==0) + double(range(DAT)==0)>0;
            DAT_save = DAT;
            DAT(:,log_remove) = [];
        end
        [DAT_s, ~] = dp_standardize(DAT, IN_s);
        [COV_s, ~] = dp_standardize(COV, IN_s);
        
        IN_c.TrCovars = COV_s;
        [DAT_sc, ~] = nk_PartialCorrelationsObj(DAT_s, IN_c);
        [OUT, ~] = dp_standardize(DAT_sc, IN_s);
        if exist('DAT_save', 'var')
            DAT_save(:,~log_remove) = OUT;
            OUT = DAT_save;
        end
    else
        if any(sum(DAT,1)==0) || any(range(DAT)==0)
            log_remove = double(sum(DAT,1)==0) + double(range(DAT)==0)>0;
            DAT_save = DAT;
            DAT(:,log_remove) = [];
        end
        [OUT, ~] = dp_standardize(DAT, IN_s);
        if exist('DAT_save', 'var')
            DAT_save(:,~log_remove) = OUT;
            OUT = DAT_save;
        end
    end
    
elseif isstruct(DAT) % perform on train and apply to test dataset
    if ~any(isnan(COV.test), 'all') && ~isempty(COV.test)
        % standardization of data and covariates
        [COV.train_s, COV.test_s] = dp_standardize_comb(COV.train, COV.test, IN_s);
        [DAT.train_s, DAT.test_s] = dp_standardize_comb(DAT.train, DAT.test, IN_s);
        
        % correction of data
        [DAT.train_sc, DAT.test_sc] = dp_corrections(DAT.train_s, DAT.test_s, COV.train_s, COV.test_s, IN_s.subgroup_train, IN_s.subgroup_test);
        
        % standardization of data
        [OUT.train, OUT.test] = dp_standardize_comb(DAT.train_sc, DAT.test_sc, IN_s);
    else
        [OUT.train, OUT.test] = dp_standardize_comb(DAT.train, DAT.test, IN_s);
    end
end

end