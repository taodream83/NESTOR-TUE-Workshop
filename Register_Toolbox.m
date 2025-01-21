% This script will recursively generate path variables for all folders inside the functions folder.
% This way you always have access to all toolbox functions.
% Matlab uses a similar method for registering its toolboxes.
% The generated paths are put in the user startup.m file.
% This allows to have different toolbox versions and locations per user

% Start from a blank workspace
clear all; %#ok
fclose('all'); % Make sure any stuck file descriptors are closed
fprintf('Registering toolbox\n');

% Generate full functions folder path by getting full current script path,
% removing script name and adding folder name
scriptPath = mfilename('fullpath');
toolboxPath = scriptPath(1:end-length(mfilename));
functionsPath = strcat(toolboxPath,'functions');

% Determine startup file location
startupFilePathOld = fullfile(userpath,'startup.m');
startupFilePathNew = fullfile(userpath,'startupNew.m');

% Open new file for writing
fNew = fopen(startupFilePathNew,'w');
% Open file for reading if it exists (2 means the file is found)
if exist(startupFilePathOld,'file') == 2
    % Open old file for reading
    fOld = fopen(startupFilePathOld,'r');
    while true
        % Read lines until -1 denotes end of file
        line = fgetl(fOld); % Use fgetl instead of fgets to ignore end of line characters
        if line == -1
            break;
            % Don't copy previous entries of toolbox path generator
        elseif ~(contains(line,'addpath(') && (contains(line,'opticalcommtoolbox') || contains(line,'userpath')))
            % Copy existing line to new file, use explicit newline for consistency
            fprintf(fNew,'%s\n',line);
        end
    end
    fclose(fOld);
end

% Print new location for toolbox to file
fprintf(fNew,'addpath(genpath(''%s''));\n',functionsPath);
% Make sure userpath is on top of the path variable for compiled mex files
% to be recognized instead of the toolbox .m files.
fprintf(fNew,'addpath(userpath);\n');
fclose(fNew);

% Execute startup to apply settings for current session
startupNew;

% Copy the new file over the old file, if anything fails before this, the old file is still in place
movefile(startupFilePathNew,startupFilePathOld);

% Check for availability of necessary toolboxes
V=ver;
if (~any(strcmp({V.Name},'Signal Processing Toolbox')))
    simplewarning('Missing MATLAB Toolbox for full functionality: Signal Processing Toolbox');
end
if (~any(strcmp({V.Name},'Communications Toolbox')))
    simplewarning('Missing MATLAB Toolbox for full functionality: Communications Toolbox');
end
if (~any(strcmp({V.Name},'DSP System Toolbox')))
    simplewarning('Missing MATLAB Toolbox for full functionality: DSP System Toolbox');
end
if (~any(strcmp({V.Name},'MATLAB Coder')))
    simplewarning('Missing MATLAB Toolbox for full functionality: MATLAB Coder');
end
if (~any(strcmp({V.Name},'Parallel Computing Toolbox')))
    simplewarning('Missing MATLAB Toolbox for full functionality: Parallel Computing Toolbox');
end

% Inform the user that registering is finished
fprintf('Done\n');
