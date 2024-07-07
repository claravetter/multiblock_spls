% prep
cd('/volume/projects/CV_gs_PLS/ScrFun/multiblock_spls/')

addpath(genpath('/opt/NM/NeuroMiner_Current/'))
addpath(genpath('/volume/projects/CV_gs_PLS/ScrFun/multiblock_spls/scripts/MBSPLS_Toolbox'))

[matlab_version,d] = version; 

mbspls_version = 'Dev_2024';

mainflag = 1; 
hyperoptflag = 0;
permutationflag = 0; 
bootstrapflag = 0; 

slurmflag = 1; 
sgeflag = 1; 


if hyperoptflag
    mcc('/volume/projects/CV_gs_PLS/ScrFun/multiblock_spls/scripts/SPLS_Toolbox/cv_ICV_hyperopt_csv.m');
end

if permutationflag
    mcc('/volume/projects/CV_gs_PLS/ScrFun/multiblock_spls/scripts/SPLS_Toolbox/cv_ICV_permutation_csv.m');
end

if bootstrapflag
    mcc('/volume/projects/CV_gs_PLS/ScrFun/multiblock_spls/scripts/SPLS_Toolbox/cv_ICV_bootstrap_csv.m';)
end

if mainflag
    if slurmflag
        mcc('/volume/projects/CV_gs_PLS/ScrFun/multiblock_spls/scripts/SPLS_Toolbox/cv_gspls_SLURM_standalone_Dev_2024.m');
    end
    if sgeflag
        mcc('/volume/projects/CV_gs_PLS/ScrFun/multiblock_spls/scripts/SPLS_Toolbox/cv_gspls_standalone_Dev_2024.m');
    end
end

if hyperoptflag && slurmflag
    sys('cp -r ')
end
