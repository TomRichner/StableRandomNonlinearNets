function copyCurrentMFile(destination)
    % COPYCURRENTMFILE Copies the currently running M-file to a specified destination.
    %
    %   copyCurrentMFile(destination) identifies the currently executing M-file
    %   (or its caller) and copies it to the provided destination path.
    %
    %   Inputs:
    %       destination - String specifying the target folder or full file path.
    %
    %   Example:
    %       copyCurrentMFile('C:\Backups\myfile_copy.m');
    
    % Get the call stack
    stack = dbstack('-completenames');
    
    if length(stack) < 1
        disp('No current M-file detected (e.g., called from command line). Nothing to copy.');
        return;
    end
    
    % If there's a caller, use the caller's file; otherwise, use the current one
    if length(stack) > 1
        sourceFile = stack(2).file;  % Caller's file
    else
        sourceFile = stack(1).file;  % This function itself (if called directly)
    end
    
    % Perform the copy
    try
        copyfile(sourceFile, destination);
        disp(['Successfully copied ' sourceFile ' to ' destination]);
    catch err
        disp('Error during copy:');
        disp(err.message);
    end
end