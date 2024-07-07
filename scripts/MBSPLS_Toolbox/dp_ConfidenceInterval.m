%% compute confidence interval from dataset

function CI = dp_ConfidenceInterval(IN)

data = IN.data;
if isfield(IN, 'sided')
    sided = IN.sided;
else
    sided = 2;
end

SEM = nanstd(data)/sqrt(length(data));            % Standard Error
switch sided
    case 1
        ts = tinv([0  0.95],length(data)-1);         % T-Score
    case 2
        ts = tinv([0.025  0.975],length(data)-1);         % T-Score
end
CI = nanmean(data) + ts*SEM;                      % Confidence Intervals

end