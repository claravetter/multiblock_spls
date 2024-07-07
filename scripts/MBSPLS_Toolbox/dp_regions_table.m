%% DP script to write out voxels as table for paper

% load output file
results_path = '/volume/HCStress/Analysis/02-Jul-2020/CISS_636_IQRadd_NCV55_single_folds_bestmerge_noval_min10_2020_5000AUC_Dev/final_results/result_final_vis.mat';
load(results_path);

readouts = {'brainnetome', 'cerebellum'};

for r=1:size(readouts,2)
    
    switch readouts{r}
        case 'brainnetome'
            fields = load('/volume/HCStress/Data/MRI/Atlases/Brainnetome_Atlas/brainnetome_indices.mat');
            atlas_path_full = '/volume/HCStress/Data/MRI/Atlases/Brainnetome_Atlas/brainnetome_3mm_636_CISS_NM_X.mat';
        case 'cerebellum'
            load('/volume/HCStress/Data/MRI/Atlases/Cerebellum-MNIflirt-MRICroN/cerebellum_vermis_indices.mat');
            atlas_path_full = '/volume/HCStress/Data/MRI/Cerebellum-MNIflirt-MRICroN/Cerebellum-MNIflirt_649_CTQ_X.mat';
            
    end
    
    % create template to sort voxels in
    fields_temp = fieldnames(fields);
    labels_regions = fields.(fields_temp{1});
    
    numbers_region_matrix = [[1:2:size(labels_regions,1)]',[2:2:size(labels_regions,1)]'];
    
    % create algorithm to sort all LV voxels in final_matrix template
    fields = fieldnames(output.atlas_readouts);
    temp_names = output.atlas_readouts.vector_1.(readouts{r}).positive.Properties.VariableNames;
    
    for i=1:(size(fields,1)-1)
        % create empty double template
        output.regions.table_collection.(['LV_', num2str(i)]).final_matrix = nan(size(numbers_region_matrix,1), 2*size(numbers_region_matrix,2));
        fields1 = fieldnames(output.atlas_readouts.(fields{i}).(readouts{r}));
        for ii=1:size(fields1,1)
            temp_vector = output.atlas_readouts.(fields{i}).(readouts{r}).(fields1{ii}){:, (strcmp(temp_names, 'region_number')+strcmp(temp_names, 'voxel_percentage_weighted')>0)};
            for iii=1:size(temp_vector,1)
                [row, col] = find(ismember(numbers_region_matrix, temp_vector{iii,1}));
                switch col
                    case 1
                        col_final = ii*col;
                    case 2
                        col_final = ii+col;
                end
                output.regions.table_collection.(['LV_', num2str(i)]).final_matrix(row, col_final) = temp_vector{iii,2};
            end
        end
    end
    
    % get absolute voxel values
    load(atlas_path_full);
    a_number = 1;
    atlas_for_analysis = round(MRI_for_analysis(a_number,:));
    [C_atlas, ~, ic_atlas] = unique(atlas_for_analysis);
    counts_atlas = accumarray(ic_atlas, 1);
    
    test = ismember(numbers_region_matrix, C_atlas);
    counts_atlas(1)=[];
    nn=1;
    for i=1:size(numbers_region_matrix,1)
        for ii=1:2
            if test(i,ii)
                output.regions.complete_voxels_regions(i,ii)=counts_atlas(nn);
                nn=nn+1;
            else
                output.regions.complete_voxels_regions(i,ii)=NaN;
            end
        end
    end
    
    output.regions.complete_voxels_regions_mean = nanmean(output.regions.complete_voxels_regions,2);
    
    % get extra cerebellum data
    % create template to sort voxels in
    fields = output.regions.cerebellum.count.voxels;
    temp_names = output.regions.cerebellum.count_names;
    
    % first compute the main cerebellum hemispheres using
    % cerebellumhemispheres_indices
    
    for i=1:(size(fields,1)-1)
        % create empty double template
        output.regions.table_collection.(['LV_', num2str(i)]).cerebellum_matrix = nan(size(cerebellumhemispheres_indices,1), 2*size(cerebellumhemispheres_indices,2));
        
        for ii=1:(size(fields,2)-1)
            try temp_vector = fields{i,ii}(:, (strcmp(temp_names, 'region_number')+strcmp(temp_names, 'voxel_percentage')>0));
                for iii=1:size(temp_vector,1)
                    [row, col] = find(ismember(cerebellumhemispheres_indices, temp_vector{iii,1}));
                    switch col
                        case 1
                            col_final = ii*col;
                        case 2
                            col_final = ii+col;
                    end
                    output.regions.table_collection.(['LV_', num2str(i)]).cerebellum_matrix(row, col_final) = temp_vector{iii,2};
                end
            end
        end
        
        output.regions.table_collection.(['LV_', num2str(i)]).vermis_matrix = nan(size(vermis_indices,1), 2*size(vermis_indices,2));
        
        for ii=1:(size(fields,2)-1)
            try temp_vector = fields{i,ii}(:, (strcmp(temp_names, 'region_number')+strcmp(temp_names, 'voxel_percentage')>0));
                for iii=1:size(temp_vector,1)
                    [row, ~] = find(ismember(vermis_indices, temp_vector{iii,1}));
                    output.regions.table_collection.(['LV_', num2str(i)]).vermis_matrix(row, ii) = temp_vector{iii,2};
                end
            end
        end
        
        
    end
    
    
    % get absolute voxel values for cerebellum
    load(atlas_path_full);
    a_number = 1;
    atlas_for_analysis = round(MRI_for_analysis(a_number,:));
    [C_atlas, ~, ic_atlas] = unique(atlas_for_analysis);
    counts_atlas = accumarray(ic_atlas, 1);
    
    test = ismember(cerebellumhemispheres_indices, C_atlas);
    counts_atlas(1)=[];
    nn=1;
    for i=1:size(cerebellumhemispheres_indices,1)
        for ii=1:2
            if test(i,ii)
                output.regions.complete_voxels_cerebellum(i,ii)=counts_atlas(nn);
                nn=nn+1;
            else
                output.regions.complete_voxels_cerebellum(i,ii)=NaN;
            end
        end
    end
    
    output.regions.complete_voxels_cerebellum_mean = nanmean(output.regions.complete_voxels_cerebellum,2);
    
    % get absolute voxel values for vermis
    atlas_path_full = '/volume/HCStress/Data/MRI/Cerebellum-MNIflirt-MRICroN/Cerebellum-MNIflirt_649_CTQ_X.mat';
    load(atlas_path_full);
    a_number = 1;
    atlas_for_analysis = round(MRI_for_analysis(a_number,:));
    [C_atlas, ~, ic_atlas] = unique(atlas_for_analysis);
    counts_atlas = accumarray(ic_atlas, 1);
    
    test = ismember(vermis_indices, C_atlas);
    counts_atlas(1)=[];
    nn=1;
    for i=1:size(vermis_indices,1)
        %     for ii=1:2
        if test(i)
            output.regions.complete_voxels_vermis(i,1)=counts_atlas(nn);
            nn=nn+1;
        else
            output.regions.complete_voxels_vermis(i,1)=NaN;
        end
        %     end
    end
    
    save(results_path, 'input', 'setup', 'output');
    
    
end

