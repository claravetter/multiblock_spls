%% script to resample images

path_brains = pwd; %'/volume/HCStress/Analysis/Resilience/DP_CISS_RSA_HC_ROD_CHR_ROP_Diag_634_GM_80PI_20GD_correct10_10_XY_SC_3/images';

brains_selected = 2;
voxsiz=[1 1 1];
% name = '/volume/HCStress/Data/MRI/DP_CISS_636_average_optthr.nii';
for b=1:brains_selected
    V = spm_vol([path_brains, '/brain_LV_final_', num2str(b), '.nii']);
%     V = spm_vol(name);
    for i=1:numel(V)
        bb        = spm_get_bbox(V(i));
        VV(1:2)   = V(i);
        VV(1).mat = spm_matrix([bb(1,:) 0 0 0 voxsiz])*spm_matrix([-1 -1 -1]);
        VV(1).dim = ceil(VV(1).mat \ [bb(2,:) 1]' - 0.1)';
        VV(1).dim = VV(1).dim(1:3);
        spm_reslice(VV,struct('mean',false,'which',1,'interp',1)); % 1 for linear
    end
end