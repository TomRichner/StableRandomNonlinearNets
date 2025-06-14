function save_data_figs_mfiles(input_folder_path, output_folder_path, folder_name, note_string, save_mat_files, save_m_files, save_open_figs, varargin)
%SAVE_DATA_FIGS_MFILES Saves M-files, open figures, and optionally .mat files to a specified location.
%
%   save_data_figs_mfiles(input_folder_path, output_folder_path, FOLDER_NAME, NOTE_STRING, save_mat_files, SAVE_M_FILES, SAVE_OPEN_FIGS, ADDITIONAL_VARS_CELL)
%   Creates a directory structure within output_folder_path. The main directory
%   will be named FOLDER_NAME. If FOLDER_NAME is empty or not provided,
%   a name based on the current timestamp and NOTE_STRING will be generated.
%   Inside this main directory, 'code' and 'figs' subdirectories are created.
%
%   If SAVE_M_FILES is true, all .m files from input_folder_path
%   are copied to the 'code' subdirectory.
%
%   If SAVE_OPEN_FIGS is true, all currently open figures are saved as .fig,
%   .png, and .svg files into the 'figs' subdirectory.
%
%   If save_mat_files is true, all .mat files in input_folder_path
%   are also copied to the 'code' subdirectory.
%
%   Inputs:
%       input_folder_path    - String. The path to copy .m and .mat files from.
%       output_folder_path   - String. The base path where the data should be saved.
%       FOLDER_NAME          - String. The name for the main save folder. If empty,
%                              a name is generated using timestamp and note_string.
%       NOTE_STRING          - String. A descriptive note to include in the folder
%                              name if FOLDER_NAME is not provided.
%       save_mat_files       - Logical. If true, copy .mat files to 'code' folder.
%       save_m_files         - Logical. If true, copy .m files. Defaults to true.
%       SAVE_OPEN_FIGS       - Logical. If true, save open figures. Defaults to true.
%       ADDITIONAL_VARS_CELL - Optional. Cell array of strings, where each
%                              string is the name of a variable in the caller's
%                              workspace to save into 'additional_workspace_data.mat'.

% Input Handling
if nargin < 1 || isempty(input_folder_path)
    error('input_folder_path is a required input. You can use pwd for the current directory.');
end
if nargin < 2 || isempty(output_folder_path)
    error('output_folder_path is a required input.');
end
if nargin < 3
    folder_name = ''; % Default to empty if not provided
end
if nargin < 4
    note_string = ''; % Default to empty if not provided
end
if nargin < 5
    save_mat_files = false; % Default to false if not provided
end
if nargin < 6
    save_m_files = true; % Default to true
end
if nargin < 7
    save_open_figs = true; % Default to true
end
additional_vars_cell = {}; % Default to empty cell
if nargin > 7 % Check if the eighth argument is provided
    if iscellstr(varargin{1}) %#ok<ISCLSTR>
        additional_vars_cell = varargin{1};
    elseif ~isempty(varargin{1})
        warning('Eighth argument to save_data_figs_mfiles should be a cell array of variable name strings. Ignoring.');
    end
end

if ~exist(input_folder_path, 'dir')
    error('Provided input_folder_path does not exist or is not a directory: %s', input_folder_path);
end

% Generate folder name
time_string = datestr(now, 'yy-mm-dd_HH_MM_SS');
name_parts = {};
if ~isempty(folder_name)
    name_parts{end+1} = folder_name;
end
if ~isempty(note_string)
    name_parts{end+1} = note_string;
end
name_parts{end+1} = time_string;
folder_name = strjoin(name_parts, '_');

% Construct full save path
save_folder = fullfile(output_folder_path, folder_name);

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
if save_m_files
    disp(['Copying .m files from ' input_folder_path '...']);
    list_m_files = dir(fullfile(input_folder_path, '*.m'));
    copied_count = 0;
    for i_m_file = 1:length(list_m_files)
        source_file = fullfile(input_folder_path, list_m_files(i_m_file).name);
        [status, msg, msgID] = copyfile(source_file, code_folder, 'f');
        if status
            copied_count = copied_count + 1;
        else
            warning('Failed to copy %s: %s (%s)', list_m_files(i_m_file).name, msg, msgID);
        end
    end
    disp(['Copied ' num2str(copied_count) ' .m files to ' code_folder]);
end

% Optionally save .mat files
if save_mat_files
    disp(['Copying .mat files from ' input_folder_path '...']);
    list_mat_files = dir(fullfile(input_folder_path, '*.mat'));
    mat_copied_count = 0;
    for i_mat_file = 1:length(list_mat_files)
        source_file = fullfile(input_folder_path, list_mat_files(i_mat_file).name);
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
if save_open_figs
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
