%% DP function for bootci

function log_unreliable = dp_bootci(input_matrix)

for ii=1:size(input_matrix,2)
    temp = input_matrix(:,ii);
    [temp_upper, ~] = bootci(1000, @(x) mean(x)+1.96*std(x)/sqrt(numel(x)), temp);
    [temp_lower, ~] = bootci(1000, @(x) mean(x)-1.96*std(x)/sqrt(numel(x)), temp);
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
    
end

end