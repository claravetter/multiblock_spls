%% DP function for combined correction and scaling training and testing data

function OUT = dp_correctscale_multi(DAT, COV, MET)

IN_s = struct;
IN_s.method = MET;

if isnumeric(DAT) % only operate on one dataset
    if ~any(isnan(COV))
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
        
        switch corr_option
            case 1 % use partial correlations to correct data
                IN_c.TrCovars = COV_s;
                [DAT_sc, ~] = nk_PartialCorrelationsObj(DAT_s, IN_c);
                
            case 2
                IN_c.S = DAT_s;
                IN_c.G = COV_s;
                [DAT_sc, ~, ~] = nk_PerfAdjForCovarsUsingPCAObj(DAT_s, IN_c);
        end
        
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
    if ~any(isnan(COV.test))
        % standardization of data and covariates
        [COV.train_s, COV.test_s] = dp_standardize_comb(COV.train, COV.test, IN_s);
        [DAT.train_s, DAT.test_s] = dp_standardize_comb(DAT.train, DAT.test, IN_s);
        
        % correction of data
        [DAT.train_sc, DAT.test_sc] = dp_corrections(DAT.train_s, DAT.test_s, COV.train_s, COV.test_s);
        
        % standardization of data
        [OUT.train, OUT.test] = dp_standardize_comb(DAT.train_sc, DAT.test_sc, IN_s);
    else
        [OUT.train, OUT.test] = dp_standardize_comb(DAT.train, DAT.test, IN_s);
    end
end

end