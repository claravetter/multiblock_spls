%% preprocessing for analysis

%  [sY, IN] = nk_PartialCorrelationsObj(Y, IN);
%  '/opt/NM/NeuroMiner_Release';

load('/volume/HCStress/Analysis/06-Jul-2018/DP_CISS_RSA_allgroups_T0_80P_raw.mat');

%  tIN.G                   = IN.G(trainind,:);  % covariates
%  tIN.DR.PercMode         = 1; %DR.PercMode   1 (dims specified in DR.dims),  2,(dims * number of subjects /0) 3 (specified in options, forget it)
%  tIN.DR.dims             = 0; % will take all components normally
%  tIN.DR.RedMode          = 'PCA'; % the other option is 'SparsePCA'but need options too
%  tIN.DR.DRsoft           = 1; % the other option is 0, uses different toolbox % Dimensionality reduction toolbox = 0
%  %     tIN.corrthresh  = 0.3; % cutoff for identification (def: 0.3)
%  [sYtr,tINtrans]         = RemoveCovWithPCA(sYtr,tIN);
%  [Y,IN] = nk_PerfAdjForCovarsUsingPCAObj(Y,IN,Y);
 
% rescale all columns from 0 to 1 to have equal effects
temp = [Age_Sex_for_analysis, sites_for_analysis];
for i=1:size(temp,2)
    norm_temp(:,i) = (temp(:,i)-min(temp(:,i)))/(max(temp(:,i))-min(temp(:,i)));
end

%% correct using PCA

 tIN.G = norm_temp;
 tIN.DR.PercMode         = 1; %DR.PercMode   1 (dims specified in DR.dims),  2,(dims * number of subjects /0) 3 (specified in options, forget it)
 tIN.DR.dims             = 0; % will take all components normally
 tIN.DR.RedMode          = 'PCA'; % the other option is 'SparsePCA'but need options too
 tIN.DR.DRsoft           = 1; % the other option is 0, uses different toolbox % Dimensionality reduction toolbox = 0
 
 X_t = MRI_for_analysis;
 Y_t = [data_for_analysis];
 [XNew,INnew] = nk_PerfAdjForCovarsUsingPCAObj(X_t,tIN,X_t);
 [YNew,INnew] = nk_PerfAdjForCovarsUsingPCAObj(Y_t,tIN,Y_t);
 
 X=XNew; Y=YNew;
 Y_names = [data_selected_names];
 save('/volume/HCStress/Analysis/06-Jul-2018/DP_CISS_RSA_allgroups_T0_80P_correctedforagesexsites.mat', 'X', 'Y', 'Y_names');

%% correct using normal correction

% rescale all columns from 0 to 1 to have equal effects
temp = [Age_Sex_HC_BOGsorted(:,1), Sites_HC_BOGsorted];
for i=1:size(temp,2)
    norm_temp(:,i) = (temp(:,i)-min(temp(:,i)))/(max(temp(:,i))-min(temp(:,i)));
end

IN.G = norm_temp;

[sX, IN] = nk_PartialCorrelationsObj(X, IN);
[sY, IN] = nk_PartialCorrelationsObj([Y, Age_Sex_HC_BOGsorted(:,2)], IN);
 
 
 
%  Test = zeros(215,1);
%  
%  for i = 1:7
%     Test = Test + i*(Y(:,26+i));
%  end
%   
%  [~,idx] = sort(Test);
%  
%  figure
%  subplot(1,2,1);imagesc(X(idx,:));
%  subplot(1,2,2);imagesc(YNew(idx,:));
 