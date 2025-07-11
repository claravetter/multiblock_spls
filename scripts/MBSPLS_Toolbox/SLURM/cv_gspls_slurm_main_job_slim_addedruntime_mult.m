%% Script to create bash files

function script_name = cv_gspls_slurm_main_job_slim_addedruntime_mult(datafile)

load(datafile);

FID = fopen([setup.analysis_folder, '/' input.name '.sh'],'w');

fprintf(FID,['# /bin/bash \n',...
    '# \n',...
    '#SBATCH --chdir=' analysis_folder '/ #output directory \n',...
    '#SBATCH --output=' setup.analysis_folder '/$SLURM_ARRAY_TASK_ID-output-main.txt #output directory \n',...
    '#SBATCH --job-name=' input.name ' # Name of the job \n',...
    '#SBATCH --mem=' num2str(setup.mem_request), 'GB \n',...
    '#SBATCH --partition=jobs-matlab \n',...
    '#SBATCH --account=core-psy \n' ...
    '#SBATCH --nodes=1 \n' ...
    'datafile=' datafile '\n',...
    'export MCR_CACHE_ROOT=', setup.cache_path, ' \n',...
    setup.spls_standalone_path '/cv_gspls_standalone_Dev_2024/', setup.compilation_subpath, '/run_cv_gspls_standalone_Dev_2024.sh /opt/matlab/', setup.matlab_version, ' $datafile']);

fclose(FID);

script_name = [setup.analysis_folder, '/' input.name '.sh'];

end