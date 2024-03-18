%% DP function for cleaning files from a specified folder

function dp_cleanup_files(folder, file_ID)

temp_dir = dir([folder, '/', file_ID, '*']);

for i=1:size(temp_dir,1)
    delete(temp_dir(i).name);
end


end