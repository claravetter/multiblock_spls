%% script for scanning txt files

function [output_file] = dp_txtscan(filepath, filetype)

fileID_temp = fopen(filepath);
output_file = fscanf(fileID_temp, filetype);
fclose(fileID_temp);

end