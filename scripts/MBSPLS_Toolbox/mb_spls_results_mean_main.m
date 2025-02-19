function mb_spls_results_mean_main(filepath, varargin)
% MB_SPLS_RESULTS_MAIN Generate analysis reports, tables, and figures for mbSPLS results.
%
%   MB_SPLS_RESULTS_MAIN(FILEPATH) processes the mbSPLS results stored in the
%   specified FILEPATH. It generates summary tables, figures (bar plots and
%   heatmaps), and a report in PDF format. 
%
%   MB_SPLS_RESULTS_MAIN(FILEPATH, 'PARAMETER', VALUE, ...) allows customization
%   by specifying parameter-value pairs.
%
%   INPUT:
%       FILEPATH: Path to the mbSPLS results file (MAT-file).
%
%   PARAMETERS (Name-Value pairs):
%       'maxFeatures' (numeric)    : Maximum number of features to include in bar plots (default: []).
%       'report_flag' (logical)    : Generate a PDF report (default: true).
%       'flip_flag' (logical)      : Flip feature weights for interpretation (default: false).
%
%   OUTPUT:
%       The function generates the following:
%         - Summary tables (Excel files) saved in a "Tables" directory.
%         - Figures (bar plots and heatmaps) saved in a "Figures" directory.
%         - A report (PDF) summarizing the analysis, saved in a "Reports" directory.
%
%   EXAMPLE:
%       % Run the analysis on a results file:
%       mb_spls_results_main('path/to/result.mat', 'maxFeatures', 10, 'flip_flag', true);
%
%   See also: MBSPLS_RESULTS_LATENT_SCORES, MBSPLS_RESULTS_FIGURES, MBSPLS_RESULTS_BOOTSTRAP_PRUNING
%
%   AUTHORS: 
%       Clara Vetter <clara.vetter@med.uni-muenchen.de>
%       Clara Weyer <clara.weyer@med.uni-muenchen.de>
%
%   DATE:
%       2024-12-10
%
%   COPYRIGHT:
%       © 2024 Section for Precision Psychiatry. All rights reserved.

p = inputParser;

addRequired(p, 'filepath', @ischar)
addParameter(p, 'maxFeatures', [], @isnumeric);
addParameter(p, 'report_flag', true, @islogical);
addParameter(p, 'flip_flag', false, @islogical);
% addParameter(p, 'figureFormat', 'png', @ischar);

parse(p, filepath, varargin{:});
maxFeatures = p.Results.maxFeatures;
report_flag = p.Results.report_flag;
flip_flag = p.Results.flip_flag;

% figureFormat = p.Results.figureFormat; 

% LOAD FILEPATH
if iscellstr(filepath)
    filepath = char(filepath);
    data = load(filepath);
elseif ~iscell(filepath) && (ischar(filepath) || istring(filepath))
    data = load(filepath);
end

% RETRIEVE PATH TO ANALYSIS FOLDER
[folderpath, ~, ~] = fileparts(filepath);

% CREATE FOLDER FOR TABLES [IN ANALYSIS FOLDER]
if flip_flag
    Path2Tables = fullfile(folderpath, 'Tables_Flipped_Weights');
else
    Path2Tables = fullfile(folderpath, 'Tables');
end
if ~exist(Path2Tables)
    mkdir(Path2Tables)
end

% CREATE FOLDER FOR FIGURES [IN ANALYSIS FOLDER]
if flip_flag
    Path2Figures = fullfile(folderpath, 'Figures_Flipped_Weights');
else
    Path2Figures = fullfile(folderpath, 'Figures');
end
if ~exist(Path2Figures)
    mkdir(Path2Figures)
end

% CREATE FOLDER FOR REPORT
% path2test = '/opt/PrecisionCodeRep/SPLS_Toolbox/mbSPLS/3_Results/';
% Path2Report = fullfile(path2test, 'Reports');
if flip_flag
    Path2Report = fullfile(folderpath, 'Reports_Flipped_Weights');
else
    Path2Report = fullfile(folderpath, 'Reports');
end
if ~exist(Path2Report)
    mkdir(Path2Report)
end

% CHECK WHETHER NUMBER OF TOP FEATURES FOR FIGURES WAS DEFINED (BARPLOT)
% if isempty(varargin)
%     maxFeatures = [];
% else
%     maxFeatures = varargin{1};
% end

%% BEFORE BOOTSTRAPPING
% CORRECTION
switch data.input.type_correction
    case {'correct', 'corrected'}
        correct_log = true;
    case {'uncorrected','uncorrect'}
        correct_log = false;
end

% LOOP THROUGH MATRICES
for matrix_idx = 1:numel(data.input.Xs)
    % CREATE EMPTY TABLE
    T = table();
    % LOOP THROUGH LATENT VARIABLES
    for lv_idx = 1:height(data.output.final_parameters)
        % SAVE FEATURE NAMES IN TABLE
        T.VariableName = data.input.Xs_feature_names{1, matrix_idx}.';
        temp_weights = data.output.final_parameters{lv_idx, 3}{1, matrix_idx}; 
        if flip_flag
            f_invert = @(x)(-1*x);
            temp_weights= f_invert(temp_weights);
        end
        T.(['LV', num2str(lv_idx)]) = temp_weights;
    end
    % SAVE TABLE AS EXCEL FILE (SHEET)
    
    writetable(T, fullfile(Path2Tables, 'LV_results.xlsx'), 'Sheet', data.input.Xs_names{matrix_idx})
    clear T
end
clear matrix_idx lv_idx


if report_flag
    import mlreportgen.report.*
    import mlreportgen.dom.*

    % [folderpath, ~, ~] = fileparts(filepath);

    % folderpath = '/opt/PrecisionCodeRep/SPLS_Toolbox/mbSPLS/3_Results';

    % Define the path where you want to save the report
    outputFileName = [date, '_Report.pdf'];

    % Combine the path and filename
    fullFilePath = fullfile(Path2Report, outputFileName);

    % // REPORT
    rpt = Report(fullFilePath, 'pdf');

    % // REPORT: TITLE PAGE
    titlepg = TitlePage;
    titlepg.Title = 'REPORT';
    % titlepg.Subtitle = 'Analysis of Experimental Data';
    titlepg.Author = '  ';
    titlepg.PubDate = '   ';
    % Get the folder of the current function
    current_folder = fileparts(mfilename('fullpath'));

    % Dynamically search for the image file in the folder and its subfolders
    image_file = 'mb_spls_logo_small.png';
    file_info = dir(fullfile(current_folder, '**', image_file)); % Recursively search subfolders

    if ~isempty(file_info)
        titlepg.Image = fullfile(file_info(1).folder, file_info(1).name); % Use the first match
    end

    % titlepg.Image = '/opt/PrecisionCodeRep/SPLS_Toolbox/mbSPLS/3_Results/mb_spls_logo_small.png'; % Optional
    add(rpt, titlepg);
    desc = Paragraph(['ANALYSIS:', newline, data.input.name]);
    desc.Style = {HAlign('center')};
    add(rpt, desc);
    for i = 1:12
        add(rpt, Paragraph(' '));
    end

    desc = Paragraph(['© ',num2str(year(datetime('now'))), 'Section for Precision Psychiatry. All rights reserved.']);
    desc.Style = {HAlign('center')};
    add(rpt, desc);

    % // REPORT: TABLE OF CONTENTS
    toc = TableOfContents;
    add(rpt, toc);

    % // REPORT: CHAPTER 1 (SETUP)
    ch1 = Chapter('Title', 'Setup');
    vars = {'date', 'analysis_folder', 'data_folder', 'standalone_version', ...
        'mbspls_standalone_path', 'matlab_version'};
    mb_spls_results_pdf_add_fields(ch1, data.setup, 'Settings', vars);
    add(rpt, ch1);

    % // REPORT: CHAPTER 1 (INPUT)
    % Input Data
    ch2 = Chapter('Title', 'Input');
    vars = {'Xs', 'covariates'};
    mb_spls_results_pdf_add_fields(ch2, data.input, 'Input Data', vars); clear vars
    add(ch2, Paragraph(' '))

    % Machine Learning Framework
    vars = {'framework', 'outer_folds', 'inner_folds', ...
        'permutation_testing', 'bootstrap_testing',...
        'optimization_strategy', 'density', 'correlation_method', ...
        'mult_test', 'statistical_testing'};

    mb_spls_results_pdf_add_fields(ch2, data.input, 'Machine Learning Framework', vars); clear vars
    add(rpt, ch2);

    % // REPORT: CHAPTER 3 (RESULTS)
    ch3 = Chapter('Title', 'Results');
    add(ch3, Paragraph(' '))

    T = table([], [], [], [], 'VariableNames', {'LV', 'P value', 'Frobenius Norm', 'R2'});

    for lv_idx = 1:height(data.output.final_parameters)
        T.LV(lv_idx) = lv_idx;
        T.("P value")(lv_idx) = data.output.final_parameters{lv_idx, matches(data.output.parameters_names, 'p')};
        T.("Frobenius Norm")(lv_idx) = data.output.final_parameters{lv_idx, matches(data.output.parameters_names, 'RHO')};
        T.R2(lv_idx) = 1;
    end

    mb_spls_results_pdf_table(ch3, T)
    add(rpt, ch3)
end

% LATENT SCORES
[LS] = mb_spls_results_latent_scores(data.input, data.output, correct_log, [], Path2Tables, flip_flag);
% [AUTOMATICALLY SAVES TABLE]

% FIGURES: HEATMAP [LATENT SCORES]
mb_spls_results_figures(LS, [], 'heatmap', Path2Figures, flip_flag)
clear LS

% FIGURES: BARPLOTS [LATENT VARIABLES]
mb_spls_results_mean_figures(data, [], 'barplot', Path2Figures, maxFeatures,flip_flag)
clear input output setup clear data

%% WITH BOOTSTRAPPING
% boot_options = {'BS', 'CI'};
% for ii=1:numel(boot_options)
%     % GET RESULTS AFTER BOOTSTRAPPING
%     % [boot_results_file, input, output] = cv_cw_mbspls_bootstrap_pruning(filepath, boot_options{ii});
%     [BS_input, BS_output] = mb_spls_results_bootstrap_pruning(filepath, boot_options{ii});
% 
%     % LOOP THROUGH MATRICES
%     for matrix_idx = 1:numel(BS_input.Xs)
%         % CREATE EMPTY TABLE
%         T = table();
%         % LOOP THROUGH LATENT VARIABLES
%         for lv_idx = 1:height(BS_output.final_parameters)
%             T.VariableName = BS_input.Xs_feature_names{1, matrix_idx}.';
%             temp_weights = BS_output.final_parameters{lv_idx, matches(BS_output.parameters_names, 'weights')}{1, matrix_idx}; 
%             if flip_flag
%                 f_invert = @(x)(-1*x);
%                 temp_weights= f_invert(temp_weights);
%             end
%             T.(['LV', num2str(lv_idx)]) = temp_weights;
%         end
% 
%         % SAVE AS EXCEL FILE
%         writetable(T, fullfile(Path2Tables, ['LV_results_', boot_options{ii}, '.xlsx']),  'Sheet', BS_input.Xs_names{matrix_idx})
%         clear T
%     end
% 
%     % LATENT SCORES
%     [LS] = mb_spls_results_latent_scores(BS_input, BS_output, correct_log, boot_options{ii}, Path2Tables, flip_flag);
%     % [LS] = cv_cw_spls_get_latent_scores(boot_results_file, correct_log, boot_options{ii});
% 
%     % FIGURES: HEATMAP [LATENT SCORES]
%     mb_spls_results_figures(LS, boot_options{ii}, 'heatmap', Path2Figures)
% 
%     % FIGURES: BARPLOTS [LATENT VARIABLES]
%     BS_data.input = BS_input; BS_data.output = BS_output;
%     mb_spls_results_figures(BS_data, boot_options{ii}, 'barplot', Path2Figures, maxFeatures, flip_flag)
%     clear boot_results_file LS BS_input BS_output BS_data
% end

