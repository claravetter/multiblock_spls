%% DP decision function for finding optimal threshold for voxels

decision = 'weights';
switch decision
    case 'voxels'
        n=2;
    case 'weights'
        n=4;
end

decision1 = 'negative';
switch decision1
    case 'positive'
        nn=1;
    case 'negative'
        nn=2;
end
        
f=figure();
bar(abs(cell2mat(output.regions.count.(decision){3,nn}(:,n))));
set(gca, 'XTickLabel', output.regions.count.(decision){3,nn}(:,1), 'XTick', 1:numel(output.regions.count.(decision){3,nn}(:,n)));
title([decision, ' ', decision1]);

f=figure();
bar(abs(cell2mat(output.regions.count.weights{3,1}(:,4))));
set(gca, 'XTickLabel', output.regions.count.weights{3,1}(:,1), 'XTick', 1:numel(output.regions.count.weights{3,1}(:,4)));

temp = cell2mat(output.regions.count.(decision){3,nn}(:,n));
for i=1:(size(temp,1)-1)
    deltas(i)=(abs(temp(i))-abs(temp(i+1)))/abs(temp(i));
end

[val, ind] = max(deltas);


