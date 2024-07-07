%% Clean up for SPLS parallel script

function dp_cleanup(folder, filecore)

filecount = size(dir([folder '/' filecore '*']),1);

for i=1:filecount
    delete([folder '/' filecore '*']);
end

end