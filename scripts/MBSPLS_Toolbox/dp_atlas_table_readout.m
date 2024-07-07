%% DP script to write out voxels as table for paper
function atlas_table_readouts = dp_atlas_table_readout(IN)

% load output file
atlases = IN.atlases(contains(IN.atlases, {'brainnetome', 'cerebellum'}, 'IgnoreCase', true));

for a=1:size(atlases,2)
    if contains(atlases{a}, 'brainnetome', 'IgnoreCase', true)
        indices_temp = load('/volume/HCStress/Data/MRI/Atlases/Brainnetome_Atlas/brainnetome_indices.mat');
        indices_temp_names = fieldnames(indices_temp);
        labels_regions = indices_temp.(indices_temp_names{1});
        
        numbers_region_matrix = [[1:2:size(labels_regions,1)]',[2:2:size(labels_regions,1)]'];
        
        % create algorithm to sort all LV voxels in cerebrum_matrix template
        fields = fieldnames(IN.atlas_readouts);
        temp_names = IN.atlas_readouts.vector_1.brainnetome.positive.Properties.VariableNames;
        
        for i=1:(size(fields,1)-1)
            % create empty double template
            atlas_table_readouts.table_collection.(['LV_', num2str(i)]).cerebrum_matrix = nan(size(numbers_region_matrix,1), 2*size(numbers_region_matrix,2));
            fields1 = fieldnames(IN.atlas_readouts.(fields{i}).brainnetome);
            for ii=1:size(fields1,1)
                temp_vector = IN.atlas_readouts.(fields{i}).brainnetome.(fields1{ii}){:, (strcmp(temp_names, 'region_number')+strcmp(temp_names, 'voxel_percentage_weighted')>0)};
                for iii=1:size(temp_vector,1)
                    [row, col] = find(ismember(numbers_region_matrix, temp_vector{iii,1}));
                    if col == 1
                        col = ii*col;
                    elseif col == 2
                        col = ii+col;
                    end
                    atlas_table_readouts.table_collection.(['LV_', num2str(i)]).cerebrum_matrix(row, col) = temp_vector{iii,2};
                end
            end
        end
        
        % get absolute voxel values
        temp_atlas = load(atlases{a});
        temp_atlas_names = fieldnames(temp_atlas);
        atlas_for_analysis = round(temp_atlas.(temp_atlas_names{1}));
        
        [C_atlas, ~, ic_atlas] = unique(atlas_for_analysis);
        counts_atlas = accumarray(ic_atlas, 1);
        
        test = ismember(numbers_region_matrix, C_atlas);
        counts_atlas(1)=[];
        nn=1;
        for i=1:size(numbers_region_matrix,1)
            for ii=1:2
                if test(i,ii)
                    atlas_table_readouts.complete_voxels_regions(i,ii)=counts_atlas(nn);
                    nn=nn+1;
                else
                    atlas_table_readouts.complete_voxels_regions(i,ii)=NaN;
                end
            end
        end
        
        atlas_table_readouts.complete_voxels_regions_mean = nanmean(atlas_table_readouts.complete_voxels_regions,2);
        
        % get extra cerebellum data
        % create template to sort voxels in
        
    elseif contains(atlases{a}, 'cerebellum', 'IgnoreCase', true)
        load('/volume/HCStress/Data/MRI/Atlases/Cerebellum-MNIflirt-MRICroN/cerebellum_vermis_indices.mat');
        
        temp_names = IN.atlas_readouts.vector_1.cerebellum.positive.Properties.VariableNames;
        
        % first compute the main cerebellum hemispheres using
        % cerebellumhemispheres_indices
        
        for i=1:(size(fields,1)-1)
            % create empty double template
            atlas_table_readouts.table_collection.(['LV_', num2str(i)]).cerebellum_matrix = nan(size(cerebellumhemispheres_indices,1), 2*size(cerebellumhemispheres_indices,2));
            fields1 = fieldnames(IN.atlas_readouts.(fields{i}).cerebellum);
            for ii=1:size(fields1,1)
                try temp_vector = IN.atlas_readouts.(fields{i}).cerebellum.(fields1{ii}){:, (strcmp(temp_names, 'region_number')+strcmp(temp_names, 'voxel_percentage_weighted')>0)};
                    for iii=1:size(temp_vector,1)
                        [row, col] = find(ismember(cerebellumhemispheres_indices, temp_vector{iii,1}));
                        if col == 1
                            col = ii*col;
                        elseif col == 2
                            col = ii+col;
                        end
                        atlas_table_readouts.table_collection.(['LV_', num2str(i)]).cerebellum_matrix(row, col) = temp_vector{iii,2};
                    end
                end
            end
            
            atlas_table_readouts.table_collection.(['LV_', num2str(i)]).vermis_matrix = nan(size(vermis_indices,1), 2*size(vermis_indices,2));
            fields1 = fieldnames(IN.atlas_readouts.(fields{i}).cerebellum);
            for ii=1:size(fields1,1)
                try temp_vector = IN.atlas_readouts.(fields{i}).cerebellum.(fields1{ii}){:, (strcmp(temp_names, 'region_number')+strcmp(temp_names, 'voxel_percentage_weighted')>0)};
                    for iii=1:size(temp_vector,1)
                        [row, col] = find(ismember(vermis_indices, temp_vector{iii,1}));
                        atlas_table_readouts.table_collection.(['LV_', num2str(i)]).vermis_matrix(row, ii) = temp_vector{iii,2};
                    end
                end
            end
            
            
        end
        
        % get absolute voxel values for cerebellum
        temp_atlas = load(atlases{a});
        temp_atlas_names = fieldnames(temp_atlas);
        atlas_for_analysis = round(temp_atlas.(temp_atlas_names{1}));
        
        [C_atlas, ~, ic_atlas] = unique(atlas_for_analysis);
        counts_atlas = accumarray(ic_atlas, 1);
        
        test = ismember(cerebellumhemispheres_indices, C_atlas);
        counts_atlas(1)=[];
        nn=1;
        for i=1:size(cerebellumhemispheres_indices,1)
            for ii=1:2
                if test(i,ii)
                    atlas_table_readouts.complete_voxels_cerebellum(i,ii)=counts_atlas(nn);
                    nn=nn+1;
                else
                    atlas_table_readouts.complete_voxels_cerebellum(i,ii)=NaN;
                end
            end
        end
        
        atlas_table_readouts.complete_voxels_cerebellum_mean = nanmean(atlas_table_readouts.complete_voxels_cerebellum,2);
        
        % get absolute voxel values for vermis
        test = ismember(vermis_indices, C_atlas);
        try counts_atlas(1)=[];
%         catch
%             counts_atlas(1)=[];
        end
        nn=1;
        for i=1:size(vermis_indices,1)
            %     for ii=1:2
            if test(i)
                atlas_table_readouts.complete_voxels_vermis(i,1)=counts_atlas(nn);
                nn=nn+1;
            else
                atlas_table_readouts.complete_voxels_vermis(i,1)=NaN;
            end
            %     end
        end
    end
end

LV_temp = fieldnames(atlas_table_readouts.table_collection);
temp_coll_bl = [round(atlas_table_readouts.complete_voxels_regions_mean); round(atlas_table_readouts.complete_voxels_cerebellum_mean)]; 
temp_coll_ul = atlas_table_readouts.complete_voxels_vermis;
temp_coll_bl_names = {'Total voxels'}; 
temp_coll_ul_names = {'Total voxels'};
bl_names = {'left', 'right'};
pn_names = {'pos', 'neg'};
bl_pn_names = {};
for i=1:size(bl_names,2)
    for ii=1:size(pn_names,2)
        bl_pn_names = [bl_pn_names, [bl_names{i}, '_', pn_names{ii}]];
    end
end

for i=1:size(LV_temp,1)
    temp_coll_bl = [temp_coll_bl, [atlas_table_readouts.table_collection.(LV_temp{i}).cerebrum_matrix; atlas_table_readouts.table_collection.(LV_temp{i}).cerebellum_matrix]];
    temp_coll_ul = [temp_coll_ul, atlas_table_readouts.table_collection.(LV_temp{i}).vermis_matrix];
    temp_coll_bl_names = [temp_coll_bl_names, cellfun(@(x) append(x, ['_', LV_temp{i}]), bl_pn_names, 'UniformOutput', false)];
    temp_coll_ul_names = [temp_coll_ul_names, cellfun(@(x) append(x, ['_', LV_temp{i}]), pn_names, 'UniformOutput', false)];
end

load('/volume/HCStress/Data/MRI/Atlases/Cerebellum-MNIflirt-MRICroN/cerebellum_names.mat');

reduced_labels_regions = strrep(labels_regions(1:2:end, 2), ' L', '');

atlas_table_readouts.complete_table_cerebrum_cerebellum = array2table(temp_coll_bl, 'VariableNames', temp_coll_bl_names, 'RowNames', [reduced_labels_regions; cerebellum_hemispheres_names]);
atlas_table_readouts.complete_table_vermis = array2table(temp_coll_ul, 'VariableNames', temp_coll_ul_names, 'RowNames', cerebellum_vermis_names);

writetable(atlas_table_readouts.complete_table_cerebrum_cerebellum, [IN.a_folder, '/table_readouts.xlsx'],  'WriteRowNames', true, 'Sheet', 'Cerebrum_Cerebellum');
writetable(atlas_table_readouts.complete_table_vermis, [IN.a_folder, '/table_readouts.xlsx'], 'WriteRowNames', true, 'Sheet', 'Vermis');

end

