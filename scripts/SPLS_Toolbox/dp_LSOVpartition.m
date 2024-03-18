%% DP function for LSOV partitions
function partition = dp_LSOVpartition(sites)
% sites needs to be coded binary with each site being represented by one
% columns

for r=1:size(sites,2)
   if sum(sites(:,r))==0
       log_full(1,r)=false;
   else
       log_full(1,r)=true;
   end
end

sites = sites(:,log_full);

for i=1:size(sites,2)
    partition.TestInd{1,i} = find(sites(:,i)==1);
    partition.TrainInd{1,i} = find(sites(:,i)==0);
end

end