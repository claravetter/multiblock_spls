%% DP function for bootstrapping from opt_parameters
function ci_collection = dp_bootstrap_opt_para(spls_standalone_path, analysis_folder, type_analysis, mem_total, max_sim_jobs, queue_name, matrix, size_sets_bootstrap)

% set up splits for bootstrapping
rest_bootstrap = mod(size(matrix,2),size_sets_bootstrap);
if rest_bootstrap>0
    bootstrap_sets = ((size(matrix,2) - rest_bootstrap)/size_sets_bootstrap)+1;
else
    bootstrap_sets = ((size(matrix,2) - rest_bootstrap)/size_sets_bootstrap);
end

interval_range= 1:size_sets_bootstrap:size(matrix,2);
bootstrap_interval = 1:bootstrap_sets;

for i=1:bootstrap_sets
    try
        temp1 = matrix(:,interval_range(bootstrap_interval(i)):(interval_range(bootstrap_interval(i+1))-1));
        dp_txt_write(analysis_folder, ['bts_', num2str(i)], temp1, '%.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f\n');
    catch
        temp1 = matrix(:,interval_range(bootstrap_interval(i)):end);
        dp_txt_write(analysis_folder, ['bts_', num2str(i)], temp1, '%.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f\n');
    end
end

CI_bash = dp_bash_parallel(spls_standalone_path, analysis_folder, type_analysis, mem_total, max_sim_jobs, queue_name, bootstrap_sets, size_sets_bootstrap);
cd(analysis_folder);
system(['qsub ' CI_bash]);

mydir = size(dir([analysis_folder '/ci_*.txt']),1);
ci_collection = [];

if bootstrap_sets>1
    while mydir<bootstrap_sets;
        if exist([analysis_folder '/ci_' num2str(bootstrap_sets) '.txt'])
            mydir_initiated = size(dir([analysis_folder '/init_*.txt']),1);
            while mydir_initiated < bootstrap_sets
                for i=1:bootstrap_sets
                    if ~exist([analysis_folder '/init_' num2str(i) '.txt'])
                        dp_CI_1000sets(num2str(i), analysis_folder)
                        mydir_initiated = size(dir([analysis_folder '/init_*.txt']),1);
                    end
                end
            end
        end
        mydir = size(dir([analysis_folder '/ci_*.txt']),1);
    end
else
    while mydir<bootstrap_sets;
        mydir = size(dir([analysis_folder '/ci_*.txt']),1);
    end
end

for i=1:bootstrap_sets
    FID=fopen([analysis_folder '/ci_', num2str(i) '.txt'], 'r');
    formatSpec = '%d';
    sizeA = [1 Inf];
    ci_temp = fscanf(FID, formatSpec, sizeA)';
    ci_collection = [ci_collection; ci_temp];
    fclose(FID);
    delete([analysis_folder '/ci_', num2str(i) '.txt']);
    delete([analysis_folder '/init_' num2str(i) '.txt']);
end

end