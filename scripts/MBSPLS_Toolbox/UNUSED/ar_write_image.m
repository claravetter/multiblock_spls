%% test Anne

mask_path = NM.brainmask{1}; %as nii

ImgData = spm_vol(mask_path);
ImgDataData = spm_read_vols(ImgData);

ImgDataData(ImgDataData >= 0) = final_parameters{1,4}(:,1)';

ImgData.fname = '/volume/DP_FEF/Analysis/17-Apr-2018/LV1BrainWeights.nii';

spm_write_vol(ImgData,ImgDataData);