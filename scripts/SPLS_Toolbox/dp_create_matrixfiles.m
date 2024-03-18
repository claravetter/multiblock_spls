%% Create X and Y files with dimensions

function dp_create_matrixfiles(output_folder,X,Y)

output_files = {X Y size(X,2) size(Y,2)};
output_files_names = {'X' 'Y' 'columns_X' 'columns_Y'};

for i =1:size(output_files,2)
    M = output_files{i};
    dlmwrite([output_folder output_files_names{i} '_file.txt'],M,'delimiter','\t');
end

end