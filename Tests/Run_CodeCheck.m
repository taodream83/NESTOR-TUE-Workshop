%% Run_CodeCheck
% Automatically present feedback for Matlab warnings and errors
% Uses the 'codecheck' function on all '.m' files in the toolbox
% Outputs feedback to 'codecheck_summary.xml' if run from GitLab

%% Clean up workspace
clear;
close all;

%% Get toolbox base location and list all '.m' files
toolboxPath = fileparts(fileparts(mfilename('fullpath')));
filepaths = dir(fullfile('..', fullfile('**', '*.m')));

%% Create XML document structure
docNode = com.mathworks.xml.XMLUtils.createDocument('testsuites');
testsuites = docNode.getDocumentElement;
testsuite = docNode.createElement('testsuite');
testsuites.appendChild(testsuite);

%% Loop over all found files and count those containing warnings
numWarningFiles = 0;
for n = 1:length(filepaths)
    %% Append each test case to the XML file
    testcase = docNode.createElement('testcase');
    testcase.setAttribute('name', filepaths(n).name);
    testcase.setAttribute('classname', 'CodeCheck');
    testsuite.appendChild(testcase);

    %% Capture output from 'checkcode' in a string
    file = fullfile(filepaths(n).folder, filepaths(n).name);
    outText = evalc('checkcode(file)');

    %% Skip if no warnings occurred
    if ~isempty(outText)
        numWarningFiles = numWarningFiles + 1;

        %% Print feedback to Matlab command window
        beginSentence = ['Warnings occurred in ', fullfile(filepaths(n).folder(length(toolboxPath)+2:end), filepaths(n).name)];
        fprintf('%s\n%s', beginSentence, outText);

        %% Append error element to the XML test case
        error = docNode.createElement('error');
        testcase.appendChild(error);

        %% Put error contents in the XML (remove hyperlinks with regexprep)
        endSentence = 'Put %#ok behind a line to suppress a warning.';
        message = [repmat('=', 1, 80), newline, beginSentence, newline, newline, regexprep(outText, '<.*?>', ''), endSentence, newline, repmat('=', 1, 80)];
        error.appendChild(docNode.createTextNode(message));
    end
end

%% Only actually write the XML file if called from GitLab
if batchStartupOptionUsed
    testsuite.setAttribute('name', 'CodeCheck');
    testsuite.setAttribute('failures', num2str(numWarningFiles));
    testsuite.setAttribute('tests', num2str(length(filepaths)));

    xmlwrite('codecheck_summary.xml', docNode);
end