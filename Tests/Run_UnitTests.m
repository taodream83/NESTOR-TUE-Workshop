%% Master file from where all unit test scripts are called
% Full documentation available on:
% - https://mathworks.com/help/matlab/script-based-unit-tests.html
% - https://mathworks.com/help/matlab/class-based-unit-tests.html

%% Basic script based unit test execution
% - Test scripts are run using: runner.run(suite), with suite = testsuite('ScriptName/ScriptFolder')
% - If a test in a section fails, the section is aborted and continues testing the next sections
% - The function returns a 'TestResult' type variable, can be converted with table() for a nicer overview

%% Basic script based unit test rules:
% - Each script must start or end with 'test' (case insensitive)
% - Each unit test should be in a separate script section (%%)
% - Text following a new %% will define the name of that test, if empty Matlab will assign a name
% - Variables defined before the first section are shared between all tests and reset after each section, other variables are only shared within a section
% - Any assert error before the first %% will prevent all tests from running in a script
% - When the script is run as a test, variables defined in a section cannot be accessed in other sections

%% JUnit XML report
% - Matlab support outputting test summaries to JUnit style XML files for use by GitLab
% - This is only executed for automated tests called from GitLab

%% Clean up workspace
clear;
close all;

%% Import necessary plugins
import matlab.unittest.TestRunner
import matlab.unittest.TestSuite
import matlab.unittest.plugins.XMLPlugin

%% Define which tests to run in what order and define test runner
suites = [
    testsuite('UnitTests/TestCompilation.m') ...
    testsuite('UnitTests/TestExamples.m') ...
    ];
runner = TestRunner.withTextOutput;

%% If run from GitLab, force local functions and output XML summary
if batchStartupOptionUsed
    restoredefaultpath;
    toolboxPath = fileparts(fileparts(mfilename('fullpath')));
    functionsPath = fullfile(toolboxPath, 'functions');
    addpath(genpath(functionsPath));

    % Set current folder as userpath, since this doesn't exist when running
    % as windows SYSTEM account. Needed for mex compilation.
    userpath(pwd);
    addpath(userpath);

    xmlFile = 'unittests_summary.xml';
    plugin = XMLPlugin.producingJUnitFormat(xmlFile);
    runner.addPlugin(plugin);
end

%% Run unit tests and collect results in a table
result = table(runner.run(suites));

%% Print result summary
fprintf('\nTotals:\n\t%d Passed, %d Failed, %d Incomplete.\n\t%.4f seconds testing time.\n', sum(result.Passed), sum(result.Failed), sum(result.Incomplete), sum(result.Duration));
