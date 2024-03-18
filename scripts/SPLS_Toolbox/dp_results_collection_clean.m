%% DP function to analyze results in varying grid densities
% addpath(genpath('/opt/NM/NeuroMiner_Elessar_1.10_beorn/'));
% addpath(genpath('/opt/NM/NeuroMiner_Elessar_1.11_beorn/'));
% rmpath(genpath('/opt/NM/NeuroMiner_Elessar_1.11_beorn/'));
addpath(genpath('/opt/NM/NeuroMiner_1.1/NeuroMiner_1.1/'));
% rmpath(genpath('/opt/NM/NeuroMiner_1.1/NeuroMiner_1.1/'));
% addpath(genpath('/opt/NM/NeuroMiner_1.2'));
% rmpath(genpath('/opt/NM/NeuroMiner_1.2'));
% addpath(genpath('/opt/NM/test_versions/NeuroMiner_May16_2022/'));
addpath('/opt/SPM/spm12_v6685_cat12_r1207/');
addpath(genpath('/volume/DP_FEF/ScrFun/ScriptsRepository'));
% addpath(genpath('/opt/PRONIASoftware/Developpment/DataAllCron/PRONIA_JuliaCode'));

results_paths = {'/volume/projects/DP_Mimics/Plexus/Analysis/11-Jan-2024/CDP_Plexus_noDDBPD_phys_cograw244_HCSSD_NCV55_100100_100BS_5000P/final_results/result.mat'};

IN.atlases_chosen = {'yeo', 'buckner', 'brainnetome', 'cerebellum', 'YB17', 'YB7', 'YB8'};
IN.atlas_directory = '/volume/HCStress/Data/MRI/Atlases';
IN.plot_all = false;
IN.specific = {'atlas', 'behavior', 'images'}; % 'atlas', 'behavior', 'images', 'sociodemographic', 'medication_wais', 'outcome', 'medication_wais'
IN.SD_selection = {'BDI', 'PANSS', 'GAF', 'NEO', 'WHO', 'neurocognition'}; %'BDI', 'GAF', 'NEO', 'QOL'
IN.type = 1; % 1:new, 2:old
IN.analysis_origin = 2; % 1: PRONIA, 2: other

for i=1:numel(results_paths)
    
    IN.results_path = results_paths{i};
    load(IN.results_path)
%     input.X = input.MRI;
    max_iter = 1;
    
    if isfield(output, 'validation_results')
        max_iter=2;
    end
    
    for a=1:max_iter
        
        switch a
            
            case 1
                IN.results_path = results_paths{i};
                IN.results_path = strrep(IN.results_path, '.mat', '_final_vis.mat');
                save(IN.results_path, 'input', 'output', 'setup');
            case 2
                IN.results_path = results_paths{i};
                IN.results_path = strrep(IN.results_path, '.mat', '_validation_vis.mat');
                save(IN.results_path, 'input', 'output', 'setup');
        end
        
        if contains(IN.results_path, 'CTQ') & ~contains(IN.results_path, 'immune', 'IgnoreCase', true)
            IN.overall_analysis = 'Stress';
        elseif contains(IN.results_path, 'CISS')
            IN.overall_analysis = 'Resilience';
        elseif contains(IN.results_path, 'WSS')
            IN.overall_analysis = 'Schizotypy';
        elseif contains(IN.results_path, 'Mimics')
            IN.overall_analysis = 'Mimics';
        elseif contains(IN.results_path, 'immune', 'IgnoreCase', true)
            IN.overall_analysis = 'immune';
        elseif contains(IN.results_path, 'MDD', 'IgnoreCase', true)
            IN.overall_analysis = 'MDD_Trauma';
        elseif contains(IN.results_path, 'FF', 'IgnoreCase', true)
            IN.overall_analysis = 'FF';
        elseif contains(IN.results_path, 'COPE', 'IgnoreCase', true)
            IN.overall_analysis = 'COPE';
        else
            disp('Something''s wrong!');
        end
        
        dp_visualize_data_Devmod(IN);
        
    end
    
end

% compute adjusted p values (FDR)
% output = dp_fdr_posthoc_adjust(results_paths{1});


