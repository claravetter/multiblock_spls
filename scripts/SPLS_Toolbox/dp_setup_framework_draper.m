%% DP function to create ML framework

function CV = dp_setup_framework_draper(IN)

switch IN.type % 1 = nested cross-validation, 2 = random hold-out splits, 3 = LSOV, 4 = random split-half
    case 1
        
        for ob=1:IN.OB
            outer_indices = crossvalind('Kfold',IN.labels,IN.OF);
            for oo=1:IN.OF
                CV.cv_outer_indices.TrainInd{ob,oo} = find(outer_indices~=oo);
                CV.cv_outer_indices.TestInd{ob,oo} = find(outer_indices==oo);
            end
        end
        for ob=1:IN.OB
            for w=1:IN.OF
                CV.cv_outer_test_labels{ob,w} = IN.labels(CV.cv_outer_indices.TestInd{ob,w},:);
                CV.cv_outer_train_labels{ob,w} = IN.labels(CV.cv_outer_indices.TrainInd{ob,w},:);
                for ib=1:IN.IB
                    inner_indices = crossvalind('Kfold',CV.cv_outer_train_labels{ob,w},IN.IF);
                    for ii=1:IN.IF
                        CV.cv_inner_indices{ob,w}.TrainInd{ib,ii} = find(inner_indices~=ii);
                        CV.cv_inner_indices{ob,w}.TestInd{ib,ii} = find(inner_indices==ii);
                    end
                end
            end
        end
    case 2
        try
            CV.cv_outer_indices = struct;
            CV.cv_outer_indices = dp_HOpartition(IN.OF, IN.labels);
        catch
            disp('Something went wrong with the hold-out partitions in the outer fold. Please check your label data.');
        end
        try
            for w=1:IN.OF
                CV.cv_outer_test_labels{1,w} = IN.labels(CV.cv_outer_indices.TestInd{1,w},:);
                CV.cv_outer_train_labels{1,w} = IN.labels(CV.cv_outer_indices.TrainInd{1,w},:);
                CV.cv_inner_indices{1,w} = dp_HOpartition(IN.IF, CV.cv_outer_train_labels{1,w});
            end
        catch
            disp('Something went wrong with the hold-out partitions in the inner fold. Please check your label data.');
        end
    case 3
        CV.cv_outer_indices = struct;
        CV.cv_outer_indices = dp_LSOVpartition(IN.labels);
        %         if ~isfield(IN, 'additional_NCV')
        %             IN.additional_NCV = 1;
        %         end
        if isempty(IN.sublabels)
            % LSOV on inner loop
            for w=1:size(IN.labels,2)
                CV.cv_outer_test_labels{1,w} = IN.labels(CV.cv_outer_indices.TestInd{1,w},:);
                CV.cv_outer_train_labels{1,w} = IN.labels(CV.cv_outer_indices.TrainInd{1,w},:);
                CV.cv_inner_indices{1,w} = dp_LSOVpartition(CV.cv_outer_train_labels{1,w});
            end
        else% CV on inner loop
            for w=1:size(IN.labels,2)
                CV.cv_outer_test_labels{1,w} = IN.sublabels(CV.cv_outer_indices.TestInd{1,w},:);
                CV.cv_outer_train_labels{1,w} = IN.sublabels(CV.cv_outer_indices.TrainInd{1,w},:);
                CV.cv_inner_indices{1,w} = nk_CVpartition2(IN.IB, IN.IF, CV.cv_outer_train_labels{1,w});
            end
        end
    case 4
        % use randperm to split the sample in half, then do NCV in the half
        CV.cv_outer_indices = struct;
        CV.cv_outer_indices = dp_RSHpartition(1, IN.labels);
        
        for w=1:IN.OF
            CV.cv_outer_train_labels{1,w} = IN.labels(CV.cv_outer_indices.TrainInd{1,w},:);
            CV.cv_inner_indices{1,w} = nk_CVpartition2(IN.IB, IN.IF, CV.cv_outer_train_labels{1,w});
        end
end

end