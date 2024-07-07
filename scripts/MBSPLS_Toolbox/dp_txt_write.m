%% DP script for writing data to txt files
function dp_txt_write(folder, name, data, format)

FID = fopen([folder '/' name '.txt'],'w');
fprintf(FID, format, data);
fclose(FID);

end