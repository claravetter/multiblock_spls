%% DP function to create folders for analysis

function [permutation_folder, hyperopt_folder, bootstrap_folder, detailed_results, final_results] = dp_create_folders(temp_folder, final_folder, analysis_name)

permutation_folder = [temp_folder, '/' analysis_name, '/permutation']; % folder for permutation testing
mkdir(permutation_folder);
hyperopt_folder = [temp_folder, '/' analysis_name,  '/hyperopt']; % folder for permutation testing
mkdir(hyperopt_folder);
bootstrap_folder = [temp_folder, '/' analysis_name,  '/bootstrap']; % folder for bootstrapping
mkdir(bootstrap_folder);

detailed_results = [final_folder '/detailed_results'];
mkdir(detailed_results);
final_results = [final_folder '/final_results'];
mkdir(final_results);

end