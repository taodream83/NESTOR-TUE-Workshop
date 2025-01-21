% This script tries to compile all C/C++ files as defined in this script.
% Before compiling it will test for available compilers and provides an
% error with instructions when absent.
% On Windows operating systems, this will automatically open the add-on
% explorer for installing the MinGW-w64 compiler as supported by Matlab.
% On Mac and Linux it will prompt the user to install the corresponding
% compiler themselves.
% Mex files will be placed in the "userpath" folder. Which is often the
% Documents/MATLAB folder on Windows and Linux and is unique per user.
% This allows compilation for different users/OS's from the same toolbox
% folder without interfering.

%% Clear all to reset compilation cache
clear all; %#ok
fprintf('Compiling toolbox\n')

%% Some logic to give feedback to the user if problems with compiler arise
try
    % Test for available compiler
    mex('-setup','C++');
catch ME
    if ispc
        % Prompt the user to install MinGW-w64 (Windows only)
        if isMATLABReleaseOlderThan("R2024a")
            matlab.internal.language.introspective.showAddon('ML_MINGW');
        else
            matlab.internal.addons.showAddon('ML_MINGW')
        end
        error('Please install the MinGW-w64 C/C++ compiler from the add-on manager (opens automatically in separate window). For more options, visit https://www.mathworks.com/support/compilers');
    else
        % Use the default error message for Linux and Mac
        error(ME.message);
    end
end

%% Compile all ESS functions
mex(fullfile('functions','Probabilistic_Shaping','ESS_C','ESS_generateTrellis.cpp'),fullfile('functions','Probabilistic_Shaping','ESS_C','ESS_functions.cpp'),'-output',fullfile(userpath,'ESS_generateTrellis'));
mex(fullfile('functions','Probabilistic_Shaping','ESS_C','ESS_encode.cpp'),fullfile('functions','Probabilistic_Shaping','ESS_C','ESS_functions.cpp'),'-output',fullfile(userpath,'ESS_encode'));
mex(fullfile('functions','Probabilistic_Shaping','ESS_C','ESS_decode.cpp'),fullfile('functions','Probabilistic_Shaping','ESS_C','ESS_functions.cpp'),'-output',fullfile(userpath,'ESS_decode'));

%% Compile all CCDM functions
V=ver; hascodertoolbox = any(strcmp({V.Name},'MATLAB Coder'));
if hascodertoolbox
    cfg = coder.config('mex');
    cfg.SaturateOnIntegerOverflow = false;
    cfg.ExtrinsicCalls = false;
    cfg.IntegrityChecks = false;
    cfg.ResponsivenessChecks = false;
    codegen('CCDM_encode.m','-config',cfg,'-d',fullfile('functions','Probabilistic_Shaping','CCDM_C'),'-args',{coder.typeof(0,[Inf,1]),coder.typeof(0),coder.typeof(0),coder.typeof(0,[1,Inf])},'-d',fullfile(userpath,'codegen'),'-o',fullfile(userpath,'CCDM_encode'));
    codegen('CCDM_decode.m','-config',cfg,'-d',fullfile('functions','Probabilistic_Shaping','CCDM_C'),'-args',{coder.typeof(0,[Inf,1]),coder.typeof(0),coder.typeof(0),coder.typeof(0,[1,Inf])},'-d',fullfile(userpath,'codegen'),'-o',fullfile(userpath,'CCDM_decode'));
    rmdir(fullfile(userpath,'codegen'),'s');
else % This assumes a suitable compiler has already been installed in earlier steps
        if isMATLABReleaseOlderThan("R2024a")
            matlab.internal.language.introspective.showAddon('ME');
        else
            matlab.internal.addons.showAddon('ME')
        end
    error('Please install Matlab Coder toolbox (opens automatically in separate window).')
end

%% Compilation is finished
fprintf('Done\n');
