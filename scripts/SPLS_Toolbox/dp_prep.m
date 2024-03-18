%% prep file for spls parallel environment

function [variable_path] = dp_prep(folder, variable, variable_name)

if ischar(variable)
    FID = fopen([folder '/' variable_name '.txt'],'w');
    fprintf(FID,'%s',variable);
    fclose(FID);
elseif isnumeric(variable)
    FID = fopen([folder '/' variable_name '.txt'],'w');
    fprintf(FID,'%f',variable);
    fclose(FID); 
end

variable_path = [folder '/' variable_name '.txt'];

end