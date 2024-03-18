%% DP function to find atlases

function atlases = dp_atlas_find(IN)

atlases_dir = dir(IN.atlas_directory);
dirs_names = {atlases_dir.name};
atlases={};
for i=1:size(IN.atlases_chosen,2)
    log_find = contains(dirs_names, IN.atlases_chosen{i}, 'IgnoreCase', true);
    atlas_dir = [IN.atlas_directory, '/', dirs_names{log_find}];
    spec_atlas_dir = dir([atlas_dir, '/*', num2str(IN.sample_size), '*.mat']);
    spec_names = {spec_atlas_dir.name};
    spec_atlas_name = spec_names{sum([contains(spec_names, IN.var_names, 'IgnoreCase', true); contains(spec_names, 'X')],1)==2};
    atlases{1,i} = [atlas_dir, '/', spec_atlas_name];
end

end