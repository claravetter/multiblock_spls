%% DP range file

function dp_range_file(output_folder, R)
dlmwrite([output_folder 'range_file.txt'],R)
end