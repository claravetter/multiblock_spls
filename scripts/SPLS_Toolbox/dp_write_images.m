%% script for writing images
MaskPath = '/volume/HCStress/Data/MRI/HC_215_average_optthr.nii';
Mask = spm_vol(MaskPath);
MaskData = spm_read_vols(Mask);

GMaskPath = '/opt/GeneralTools/BrainTemplates/Masks/abs_g_pronia_calib_study_MRI_sMRI_vbm8_mwrp1_2016-10-04_03-23-29.nii';
GMask = spm_vol(GMaskPath);
Mask = spm_vol(GMask);
mask_normal = spm_read_vols(Mask);
% new_mask_normal = mask_normal.*0;
Mask75 = MaskData(MaskData >= prctile(MaskData(MaskData ~= 0),75))';
MaskData_gmask = spm_read_vols(Mask75);

for i=1:size(final_parameters,1)
NewImg = MaskData_gmask.*0;
NewImg(MaskData_gmask>0) = final_parameters{i,4};
Mask.fname = ['/volume/HCStress/Analysis/09-Jun-2018/DP_CISS_HC_BO_5k_sets_FDR_75_gmask_2/LV',num2str(i),'.nii'];
Mask.dt = [4,0];
spm_write_vol(Mask,NewImg);
end


for i=1:size(final_parameters,1)
NewImg = MaskData.*0;
NewImg(MaskData>=1) = final_parameters{i,4};
Mask.fname = ['/volume/HCStress/Analysis/09-Jun-2018/DP_CISS_HC_BO_5k_sets_FDR/LV',num2str(i),'.nii'];
Mask.dt = [4,0];
spm_write_vol(Mask,NewImg);
end

info_brains = nan(size(final_parameters,1),8);
for i=1:size(final_parameters,1)
    info_brains(i,1) = min(final_parameters{i,4}(final_parameters{i,4}<0));
    info_brains(i,2) = max(final_parameters{i,4}(final_parameters{i,4}<0));
    info_brains(i,3) = sum(final_parameters{i,4}<0);
    info_brains(i,4) = mean(final_parameters{i,4}(final_parameters{i,4}<0));
    info_brains(i,5) = min(final_parameters{i,4}(final_parameters{i,4}>=0));
    info_brains(i,6) = max(final_parameters{i,4}(final_parameters{i,4}>=0));
    info_brains(i,7) = sum(final_parameters{i,4}>=0);
    info_brains(i,8) = mean(final_parameters{i,4}(final_parameters{i,4}>=0));    
end