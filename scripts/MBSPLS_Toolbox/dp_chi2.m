%% playground

function [collection, collection_names, h, p, stats]=dp_chi2(obs, type)

switch type
    case 'absolute'
        sum_temp=sum(obs,2);
        size_temp=size(obs,1);
        
        bins = 1:size(sum_temp,1);
        [h.all,p.all,stats.all]=chi2gof(bins, 'Frequency', sum_temp, 'Expected', ones(size_temp,1)*(1/size_temp)*sum(sum_temp));
        
        h.detailed=nan(size_temp,size_temp);
        p.detailed=nan(size_temp,size_temp);
        
        for i=1:size_temp
            for ii=(i+1):size_temp
                temp_vector = [sum_temp(i); sum_temp(ii)];
                bins = 1:size(temp_vector,1);
                [~,p.detailed(i,ii),~]=chi2gof(bins, 'Frequency',temp_vector, 'Expected', ones(size(temp_vector,1),1)*(1/size(temp_vector,1))*sum(temp_vector));
            end
        end
        
        p.p_value_FDR = dp_FDR(p.detailed(~isnan(p.detailed)), 0.05);
        
        h.detailed=nan(size_temp,size_temp);
        p.detailed=nan(size_temp,size_temp);
        
        for i=1:size_temp
            for ii=(i+1):size_temp
                temp_vector = [sum_temp(i); sum_temp(ii)];
                bins = 1:size(temp_vector,1);
                [h.detailed(i,ii),p.detailed(i,ii),stats.detailed(i,ii)]=chi2gof(bins, 'Frequency',temp_vector, 'Expected', ones(size(temp_vector,1),1)*(1/size(temp_vector,1))*sum(temp_vector), 'Alpha', p.p_value_FDR);
            end
        end
        
        nn=1;
        for i=1:size(p.detailed,1)
            for ii=1:size(p.detailed,2)
                if ~isnan(p.detailed(i,ii))
                    collection(nn,:) = [i, ii, p.detailed(i,ii), stats.detailed(i,ii).chi2stat];
                    nn=nn+1;
                end
            end
        end
        
        collection_names = {'group 1', 'group 2', 'p', 'chi2'};
    case 'ratio'
        
        sum_temp=sum(obs,2)/sum(sum(obs,2));
        size_temp=size(obs,2);
        
        bins = 1:size(sum_temp,1);
        [h.all,p.all,stats.all]=chi2gof(bins, 'Frequency', sum_temp, 'Expected', ones(size_temp,1)*(1/size_temp)*sum(sum_temp));
        
        h.detailed=nan(size_temp,size_temp);
        p.detailed=nan(size_temp,size_temp);
        
        for i=1:size_temp
            for ii=(i+1):size_temp
                temp_vector = [sum_temp(i); sum_temp(ii)];
                bins = 1:size(temp_vector,1);
                [~,p.detailed(i,ii),~]=chi2gof(bins, 'Frequency',temp_vector, 'Expected', ones(size(temp_vector,1),1)*(1/size(temp_vector,1))*sum(temp_vector));
            end
        end
        
        p.p_value_FDR = dp_FDR(p.detailed(~isnan(p.detailed)), 0.05);
        
        h.detailed=nan(size_temp,size_temp);
        p.detailed=nan(size_temp,size_temp);
        
        for i=1:size_temp
            for ii=(i+1):size_temp
                temp_vector = [sum_temp(i); sum_temp(ii)];
                bins = 1:size(temp_vector,1);
                [h.detailed(i,ii),p.detailed(i,ii),stats.detailed(i,ii)]=chi2gof(bins, 'Frequency',temp_vector, 'Expected', ones(size(temp_vector,1),1)*(1/size(temp_vector,1))*sum(temp_vector), 'Alpha', p.p_value_FDR);
            end
        end
        
        nn=1;
        for i=1:size(p.detailed,1)
            for ii=1:size(p.detailed,2)
                if ~isnan(p.detailed(i,ii))
                    collection(nn,:) = [i, ii, p.detailed(i,ii), stats.detailed(i,ii).chi2stat];
                    nn=nn+1;
                end
            end
        end
        
end

end
