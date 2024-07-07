%% function for parallel CI computing

function dp_CI_100sets(i, analysis_folder)

dp_txt_write(analysis_folder, ['init_' i],'initialized','%s \n');

FID=fopen([analysis_folder '/bts_', i '.txt'], 'r');
formatSpec = '%f';
sizeA = [10 Inf];
input_matrix = fscanf(FID, formatSpec, sizeA);
fclose(FID);

for ii=1:size(input_matrix,2)
    temp = input_matrix(:,ii);
    [temp_upper, ~] = bootci(100, @(x) mean(x)+(1.96*std(x)/sqrt(numel(x))), temp);
    [temp_lower, ~] = bootci(100, @(x) mean(x)-(1.96*std(x)/sqrt(numel(x))), temp);
    log_pos =  [temp_upper;temp_lower]>0;
    log_zero = [temp_upper;temp_lower]==0;
    log_neg =  [temp_upper;temp_lower]<0;
    if sum(log_pos)>0 & sum(log_neg)>0
        log_unreliable(1,ii) = true;
    elseif sum(log_zero)>0
        log_unreliable(1,ii) = true;
    else
        log_unreliable(1,ii) = false;
    end
    
    if mod(ii,100) == 0
        FID = fopen([analysis_folder, '/init_' i '.txt'], 'a');
        fprintf(FID, '%d \n', ii);
        fclose(FID);
    end
    
end

dp_txt_write(analysis_folder, ['ci_', i], log_unreliable', '%d\n');
delete([analysis_folder '/bts_', i '.txt']);

end