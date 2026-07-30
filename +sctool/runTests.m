function result = runTests(varargin)
%RUNTESTS  Run the SC Toolbox regression test suite.
%
%   sctool.runTests()              run all tests
%   sctool.runTests('disk')        run only testDisk
%   sctool.runTests('disk', 'ann') run a subset
%
%   Valid names: 'disk', 'ext', 'hpl', 'strip', 'rect', 'crdisk', 'ann',
%                'vderiv'

CLASS_MAP = struct( ...
    'disk',   'testDisk', ...
    'ext',    'testExterior', ...
    'hpl',    'testHalfplane', ...
    'strip',  'testStrip', ...
    'rect',   'testRectangle', ...
    'crdisk', 'testCRDisk', ...
    'ann',    'testAnnulus', ...
    'vderiv', 'testVertexDerivative');

if isempty(varargin)
    names = fieldnames(CLASS_MAP);
else
    names = varargin;
    for i = 1:numel(names)
        if ~isfield(CLASS_MAP, names{i})
            error('Unknown test suite ''%s''. Valid names: %s', ...
                  names{i}, strjoin(fieldnames(CLASS_MAP)', ', '));
        end
    end
end

% Ensure tests/ directory is on the path
here    = fileparts(which('lapsolve'));
testdir = fullfile(here, 'tests');
added   = false;
if ~any(strcmp(strsplit(path, pathsep), testdir))
    addpath(testdir);
    added = true;
end
cleanup = onCleanup(@() localRmpath(added, testdir)); %#ok<NASGU>

% Build suite
suite = [];
for i = 1:numel(names)
    cls   = CLASS_MAP.(names{i});
    suite = [suite, matlab.unittest.TestSuite.fromClass(meta.class.fromName(cls))]; %#ok<AGROW>
end

% Run
runner = matlab.unittest.TestRunner.withTextOutput('OutputDetail', ...
    matlab.unittest.Verbosity.Detailed);
result = runner.run(suite);

% Summary
nPass = sum([result.Passed]);
nFail = sum([result.Failed]);
nSkip = sum([result.Incomplete]);
fprintf('\n%d passed, %d failed, %d incomplete\n', nPass, nFail, nSkip);

if nargout == 0
    clear result
end
end

function localRmpath(added, testdir)
if added; rmpath(testdir); end
end
