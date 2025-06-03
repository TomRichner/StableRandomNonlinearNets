function save_data_figs_mfiles(folder_path, folder_name, note_string, save_data, varargin)
%SAVE_DATA_FIGS_MFILES Saves current M-files, open figures, and optionally .mat files to a specified location.
%
%   save_data_figs_mfiles(FOLDER_PATH, FOLDER_NAME, NOTE_STRING, SAVE_DATA)
%   Creates a directory structure within FOLDER_PATH. The main directory
%   will be named FOLDER_NAME. If FOLDER_NAME is empty or not provided,
%   a name based on the current timestamp and NOTE_STRING will be generated.
%   Inside this main directory, 'code' and 'figs' subdirectories are created.
%   All .m files from the current working directory are copied to the 'code'
%   subdirectory. All currently open figures are saved as .fig and .png
%   files into the 'figs' subdirectory.
%   If SAVE_DATA is true, all .mat files in the current working directory
%   are also copied to the 'code' subdirectory.
%
%   Inputs:
%       FOLDER_PATH   - String. The base path where the data should be saved.
%       FOLDER_NAME   - String. The name for the main save folder. If empty,
%                       a name is generated using timestamp and note_string.
%       NOTE_STRING   - String. A descriptive note to include in the folder
%                       name if FOLDER_NAME is not provided.
%       SAVE_DATA     - Logical. If true, also copy .mat files to 'code' folder.
%       ADDITIONAL_VARS_CELL - Optional. Cell array of strings, where each
%                              string is the name of a variable in the caller's
%                              workspace to save into 'additional_workspace_data.mat'.

% Input Handling
if nargin < 1 || isempty(folder_path)
    error('FOLDER_PATH is a required input.');
end
if nargin < 2
    folder_name = ''; % Default to empty if not provided
end
if nargin < 3
    note_string = ''; % Default to empty if not provided
end
if nargin < 4
    save_data = false; % Default to false if not provided
end
additional_vars_cell = {}; % Default to empty cell
if nargin > 4 % Check if the fifth argument is provided
    % The fifth argument from the call was passed directly.
    % In the previous response, I called it additional_vars_to_save when calling.
    % Let's assume the actual fifth argument passed is the cell array.
    if iscellstr(varargin{1}) %#ok<ISCLSTR>
        additional_vars_cell = varargin{1};
    elseif ~isempty(varargin{1})
        warning('Fifth argument to save_data_figs_mfiles should be a cell array of variable name strings. Ignoring.');
    end
end

% Generate folder name if not provided
if isempty(folder_name)
    time_string = datestr(now, 'yy-mm-dd_HH_MM_SS');
    if ~isempty(note_string)
        folder_name = [time_string '_' note_string];
    else
        folder_name = time_string;
    end
end

% Construct full save path
save_folder = fullfile(folder_path, folder_name);

% Create directories
code_folder = fullfile(save_folder, 'code');
figs_folder = fullfile(save_folder, 'figs');

if ~exist(save_folder, 'dir')
    mkdir(save_folder);
    disp(['Created folder: ' save_folder]);
else
    disp(['Folder already exists: ' save_folder]);
end

if ~exist(code_folder, 'dir')
    mkdir(code_folder);
    disp(['Created subfolder: ' code_folder]);
end

if ~exist(figs_folder, 'dir')
    mkdir(figs_folder);
    disp(['Created subfolder: ' figs_folder]);
end

% Save m-files
disp('Copying .m files...');
current_dir = pwd; % Get current working directory
list_m_files = dir(fullfile(current_dir, '*.m'));
copied_count = 0;
for i_m_file = 1:length(list_m_files)
    source_file = fullfile(current_dir, list_m_files(i_m_file).name);
    [status, msg, msgID] = copyfile(source_file, code_folder, 'f');
    if status
        copied_count = copied_count + 1;
    else
        warning('Failed to copy %s: %s (%s)', list_m_files(i_m_file).name, msg, msgID);
    end
end
disp(['Copied ' num2str(copied_count) ' .m files to ' code_folder]);

% Optionally save .mat files
if save_data
    disp('Copying .mat files...');
    list_mat_files = dir(fullfile(current_dir, '*.mat'));
    mat_copied_count = 0;
    for i_mat_file = 1:length(list_mat_files)
        source_file = fullfile(current_dir, list_mat_files(i_mat_file).name);
        [status, msg, msgID] = copyfile(source_file, code_folder, 'f');
        if status
            mat_copied_count = mat_copied_count + 1;
        else
            warning('Failed to copy %s: %s (%s)', list_mat_files(i_mat_file).name, msg, msgID);
        end
    end
    disp(['Copied ' num2str(mat_copied_count) ' .mat files to ' code_folder]);
end

% Save figures
disp('Saving open figures...');
fig_handles = findall(groot, 'Type', 'figure'); % Use groot to get all figure handles
saved_count = 0;

if isempty(fig_handles)
    disp('No open figures found to save.');
else
    for i = 1:length(fig_handles)
        fig = fig_handles(i);
        fig_number_str = num2str(fig.Number);
        base_filename = fullfile(figs_folder, ['figure_' fig_number_str]);

        try
            % % Save as .fig
            savefig(fig, [base_filename '.fig']);
            % Save as .png
            print(fig, [base_filename '.png'], '-dpng', '-r150');

            set(fig, 'Renderer', 'painters');
            saveas(fig, base_filename, 'svg');

            saved_count = saved_count + 1;
            
        catch ME
            warning('Failed to save figure %s: %s', fig_number_str, ME.message);
        end
    end
    disp(['Saved ' num2str(saved_count) ' figures to ' figs_folder]);
end

% Save additional variables from caller's workspace if specified
if ~isempty(additional_vars_cell)
    additional_data_filename = fullfile(save_folder, 'additional_workspace_data.mat');
    try
        vars_to_save_struct = struct();
        valid_vars_found = false;
        for i_var = 1:length(additional_vars_cell)
            var_name = additional_vars_cell{i_var};
            if ischar(var_name) && evalin('caller', ['exist(''' var_name ''', ''var'')'])
                vars_to_save_struct.(var_name) = evalin('caller', var_name);
                valid_vars_found = true;
            else
                warning('Variable "%s" not found in caller workspace or not a string. Skipping.', var_name);
            end
        end
        
        if valid_vars_found
            save(additional_data_filename, '-struct', 'vars_to_save_struct');
            disp(['Saved additional workspace variables to ' additional_data_filename]);
        else
            disp('No valid additional workspace variables found to save.');
        end
    catch ME_vars
        warning(ME_vars.identifier, 'Failed to save additional workspace variables: %s', ME_vars.message);
    end
end

disp('Operation complete.');

end % End of function
