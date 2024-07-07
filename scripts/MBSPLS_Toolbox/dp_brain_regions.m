%% function for brain regions

function dp_brain_regions(results_path, overall_analysis)

switch overall_analysis
    case 'Stress'
        overall_folder = '/volume/HCStress/Analysis/Stress';
    case 'Resilience'
        overall_folder = '/volume/HCStress/Analysis/Resilience';
end

%% visualize results
load(results_path);
collection_folder = [overall_folder, '/', input.name];
mkdir(collection_folder);

switch size(input.behavior,1)
    case 626
        hammers_path_full = '/volume/HCStress/Data/MRI/Hammers_Atlas/Hammers_mith-n30r95_spm12_full__626_NM_X.mat';
    case 627
        hammers_path_full = '/volume/HCStress/Data/MRI/Hammers_Atlas/Hammers_mith-n30r95_spm12_full__627_NM_X.mat';
    case 631
        hammers_path_full = '/volume/HCStress/Data/MRI/Hammers_Atlas/Hammers_mith-n30r95_spm12_full__631_NM_X.mat';
    case 634
        hammers_path_full = '/volume/HCStress/Data/MRI/Hammers_Atlas/Hammers_mith-n30r95_spm12_full__634_NM_X.mat';
end

load(hammers_path_full);
load('/volume/HCStress/Data/MRI/Hammers_Atlas/Hammers_mith-n30-ancillary-data.tar/Hammers_mith-n30-ancillary-data/Hammers_mith-n30-ancillary-data/labels_hammerssmith_n30_ancillary_data.mat');
output.regions.count_names = {'region_number', 'n_voxels', 'region_name', 'median_weights', 'mean_weights'};
output.regions.count = struct('raw', [], 'voxels', [], 'weights', []);
fields = fieldnames(output.regions.count);
for i=1:size(fields,1)
    FID = fopen([collection_folder, '/brain_regions_hammers_', fields{i}, '.txt'], 'w');
    fprintf(FID, [strrep(input.name,'_',' '), '\n', fields{i} ' sorted']);
    fclose(FID);
end

for i=1:size(output.final_parameters,1)
    output.regions.log{i,1} = output.final_parameters{i,4}>0;
    output.regions.log{i,2} = output.final_parameters{i,4}<0;
    output.regions.log{i,3} = output.final_parameters{i,4}~=0;
    for ii=1:size(output.regions.log,2)
        output.regions.sum{i,ii} = atlas_hammers_full_for_analysis((atlas_hammers_full_for_analysis~=0)' & output.regions.log{i,ii});
        [C, ~, ic] = unique(output.regions.sum{i,ii});
        a_counts = accumarray(ic, 1);
        output.regions.count.raw{i,ii} = [num2cell(C'), num2cell(a_counts), labels_regions_hammers((C'),2)];
        for iii=1:size(output.regions.count.raw{i,ii},1)
            output.regions.count.raw{i,ii}{iii,4} = median(abs(output.final_parameters{i,4}(output.regions.log{i,ii} & (atlas_hammers_full_for_analysis==output.regions.count.raw{i,ii}{iii,1})')));
            output.regions.count.raw{i,ii}{iii,5} = mean(abs(output.final_parameters{i,4}(output.regions.log{i,ii} & (atlas_hammers_full_for_analysis==output.regions.count.raw{i,ii}{iii,1})')));
        end
        output.regions.count.voxels{i,ii} = flipud(sortrows((output.regions.count.raw{i,ii}), find(strcmp(output.regions.count_names, 'n_voxels'))));
        output.regions.count.weights{i,ii} = flipud(sortrows((output.regions.count.raw{i,ii}), find(strcmp(output.regions.count_names, 'mean_weights'))));
    end
    
    temp_add = [];
    temp_print = output.regions.count_names([1,2,4,5,3]);
    for f=1:size(temp_print,2)
        %         try output.regions.count_names{f} = num2str(temp_print{f});
        %         catch ME
        %         end
        temp_add = [temp_add, sprintf('\t'), temp_print{f}];
    end
    
    
    fields = fieldnames(output.regions.count);
    for ff=1:size(fields,1)
        FID = fopen([collection_folder, '/brain_regions_hammers_', fields{ff}, '.txt'], 'a');
        fprintf(FID, ['\n \n \n Latent Variable ' num2str(i), '\n \n', temp_add, '\n']);
        fclose(FID);
    end
    
    fields = fieldnames(output.regions.count);
    for ff=1:size(fields,1)
        for r=1:size(output.regions.count.(fields{ff}){i,3},1)
            temp_print = output.regions.count.(fields{ff}){i,3}(r,[1,2,4,5,3]);
            temp_add = [];
            for ii=1:size(temp_print,2)
                try temp_print{ii} = num2str(temp_print{ii});
                catch
                end
                temp_add = [temp_add, sprintf('\t'), temp_print{ii}];
            end
            FID = fopen([collection_folder, '/brain_regions_hammers_', fields{ff}, '.txt'], 'a');
            fprintf(FID, ['\n' temp_add]);
            fclose(FID);
        end
    end
end

save(results_path, 'input', 'output', 'setup');



end