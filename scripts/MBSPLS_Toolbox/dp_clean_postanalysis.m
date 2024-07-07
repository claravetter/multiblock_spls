%% loop to clean up keep_in_partition, opt_param

folders = dir;

for i=1:numel(folders)
    if folders(i).isdir
        if numel(folders(i).name) > 5
            subfolders=dir(folders(i).name); 
            for ii=1:numel(subfolders)
                if(subfolders(ii).isdir)
                    if numel(subfolders(ii).name) > 5
                        delete([folders(i).name, '/', subfolders(ii).name, '/hyperopt/keep_in_partition.mat']);
                        delete([folders(i).name, '/', subfolders(ii).name, '/permutation/opt_param.mat']);
                        if ~exist([folders(i).name, '/', subfolders(ii).name, '/final_results/result.mat'])
                            rmdir([folders(i).name, '/', subfolders(ii).name],'s');
                            disp([folders(i).name, '/', subfolders(ii).name, ' was deleted']);
                        end
                    end
                end
            end
        end
    end
end