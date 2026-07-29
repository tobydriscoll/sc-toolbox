function generateGoldens(which)
%GENERATEGOLDENS  Generate golden values for C++ translation testing.
%
%   generateGoldens()           regenerate all golden files
%   generateGoldens('gaussj')   regenerate one group
%
%   Valid group names:
%     'gaussj', 'scqdata', 'scangle_scfix',
%     'diskmap_private', 'hplmap_private', 'extermap_private',
%     'stripmap_private', 'rectmap_private', 'crdiskmap_private',
%     'annulus_private', 'nesolve'
%
%   Output: tests/cpp/goldens/<group>.mat, each containing a struct array
%   'cases' with fields: desc, inputs, outputs, tol.
%
%   Prerequisites: SC Toolbox on the MATLAB path.

ALL = {'gaussj', 'scqdata', 'scangle_scfix', ...
       'diskmap_private', 'hplmap_private', 'extermap_private', ...
       'stripmap_private', 'rectmap_private', 'crdiskmap_private', ...
       'annulus_private', 'nesolve', 'isinpoly', 'polygon'};

if nargin < 1
    which = ALL;
elseif ischar(which)
    which = {which};
end

for i = 1:numel(which)
    if ~ismember(which{i}, ALL)
        error('Unknown group ''%s''. Valid groups: %s', which{i}, strjoin(ALL, ', '));
    end
end

outDir = fullfile(fileparts(mfilename('fullpath')), 'goldens');
if ~exist(outDir, 'dir'), mkdir(outDir); end

genDir = fullfile(fileparts(mfilename('fullpath')), 'generators');
addpath(genDir);
cleanup = onCleanup(@() rmpath(genDir));

for i = 1:numel(which)
    fprintf('--- %s ---\n', which{i});
    % Some private functions (e.g. sctool.findz0, used by deinvmap/
    % crinvmap/etc.) draw from the global RNG (rand(1)) to pick fallback
    % search directions. Reseed before every group so golden values are
    % reproducible regardless of which other groups ran earlier in the
    % same MATLAB session.
    rng(0);
    feval(['gen_' which{i}], outDir);
end
fprintf('\nAll goldens written to %s\n', outDir);
end
