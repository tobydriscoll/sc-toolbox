classdef testCppGoldens < matlab.unittest.TestCase
%TESTCPPGOLDENS  Verify that MATLAB golden values are self-consistent.
%
%   This test class recomputes every golden case and confirms that the
%   current MATLAB implementation matches the stored values to within tol.
%   It is also the acceptance specification for the C++ port: a C++
%   implementation must produce results matching these goldens to within
%   the same tolerances.
%
%   Prerequisites: run generateGoldens() once to populate tests/cpp/goldens/.
%
%   Run from MATLAB (SC Toolbox on path):
%       result = matlab.unittest.TestRunner.withTextOutput().run( ...
%                  matlab.unittest.TestSuite.fromClass(?testCppGoldens));

    properties
        goldenDir  = ''
        hasGoldens = false
    end

    methods (TestClassSetup)
        function locateGoldens(tc)
            here = fileparts(mfilename('fullpath'));
            tc.goldenDir = fullfile(here, 'goldens');
            tc.hasGoldens = exist(tc.goldenDir, 'dir') == 7;
        end
    end

    % ------------------------------------------------------------------
    % gaussj
    % ------------------------------------------------------------------
    methods (Test)
        function testGaussj(tc)
            tc.assumeTrue(tc.hasGoldens, 'goldens/ dir missing — run generateGoldens() first');
            f = fullfile(tc.goldenDir, 'gaussj.mat');
            tc.assumeTrue(exist(f,'file')==2, 'gaussj.mat missing');
            s = load(f);
            for k = 1:numel(s.cases)
                c = s.cases(k);
                [z, w] = sctool.gaussj(c.inputs.n, c.inputs.alf, c.inputs.bet);
                tc.verifyEqual(z, c.outputs.z, 'AbsTol', c.tol, c.desc);
                tc.verifyEqual(w, c.outputs.w, 'AbsTol', c.tol, c.desc);
            end
        end

        % ------------------------------------------------------------------
        % scqdata
        % ------------------------------------------------------------------
        function testScqdata(tc)
            tc.assumeTrue(tc.hasGoldens, 'goldens/ dir missing');
            f = fullfile(tc.goldenDir, 'scqdata.mat');
            tc.assumeTrue(exist(f,'file')==2, 'scqdata.mat missing');
            s = load(f);
            for k = 1:numel(s.cases)
                c = s.cases(k);
                qdat = sctool.scqdata(c.inputs.beta, c.inputs.nqpts);
                tc.verifyEqual(qdat, c.outputs.qdat, 'AbsTol', c.tol, c.desc);
            end
        end

        % ------------------------------------------------------------------
        % scangle and scfix
        % ------------------------------------------------------------------
        function testScangle(tc)
            tc.assumeTrue(tc.hasGoldens, 'goldens/ dir missing');
            f = fullfile(tc.goldenDir, 'scangle_scfix.mat');
            tc.assumeTrue(exist(f,'file')==2, 'scangle_scfix.mat missing');
            s = load(f);
            for k = 1:numel(s.cases_angle)
                c = s.cases_angle(k);
                beta = sctool.scangle(c.inputs.w);
                tc.verifyEqual(beta, c.outputs.beta, 'AbsTol', c.tol, c.desc);
            end
        end

        function testScfix(tc)
            tc.assumeTrue(tc.hasGoldens, 'goldens/ dir missing');
            f = fullfile(tc.goldenDir, 'scangle_scfix.mat');
            tc.assumeTrue(exist(f,'file')==2, 'scangle_scfix.mat missing');
            s = load(f);
            for k = 1:numel(s.cases_scfix)
                c = s.cases_scfix(k);
                [wf, bf] = scfix(c.inputs.type, c.inputs.w, c.inputs.beta);
                tc.verifyEqual(wf, c.outputs.w,    'AbsTol', c.tol, c.desc);
                tc.verifyEqual(bf, c.outputs.beta,  'AbsTol', c.tol, c.desc);
            end
        end

        % ------------------------------------------------------------------
        % diskmap private
        % ------------------------------------------------------------------
        function testDiskPrivate(tc)
            tc.assumeTrue(tc.hasGoldens, 'goldens/ dir missing');
            f = fullfile(tc.goldenDir, 'diskmap_private.mat');
            tc.assumeTrue(exist(f,'file')==2, 'diskmap_private.mat missing');
            s = load(f);

            opt = sctool.scmapopt('trace', 0, 'tol', 1e-12);

            for k = 1:numel(s.cases_dparam)
                c = s.cases_dparam(k);
                [z, cv, qdat] = diskmap.private_('dparam', c.inputs.w, c.inputs.beta, [], opt);
                tc.verifyEqual(z,    c.outputs.z,    'AbsTol', c.tol, c.desc);
                tc.verifyEqual(cv,   c.outputs.c,    'AbsTol', c.tol, c.desc);
                tc.verifyEqual(qdat, c.outputs.qdat, 'AbsTol', c.tol, c.desc);
            end
            for k = 1:numel(s.cases_dquad)
                c = s.cases_dquad(k);
                I = diskmap.private_('dquad', c.inputs.z1, c.inputs.z2, c.inputs.sing1, ...
                          c.inputs.z, c.inputs.beta, c.inputs.qdat);
                tc.verifyEqual(I, c.outputs.I, 'AbsTol', c.tol, c.desc);
            end
            for k = 1:numel(s.cases_dderiv)
                c = s.cases_dderiv(k);
                fp = diskmap.private_('dderiv', c.inputs.zp, c.inputs.z, c.inputs.beta, c.inputs.c);
                tc.verifyEqual(fp, c.outputs.fp, 'AbsTol', c.tol, c.desc);
            end
            for k = 1:numel(s.cases_dmap)
                c = s.cases_dmap(k);
                wp = diskmap.private_('dmap', c.inputs.zp, c.inputs.z, c.inputs.c, [], c.inputs.beta, c.inputs.qdat);
                tc.verifyEqual(wp, c.outputs.wp, 'AbsTol', c.tol, c.desc);
            end
            for k = 1:numel(s.cases_dinvmap)
                c = s.cases_dinvmap(k);
                zp = diskmap.private_('dinvmap', c.inputs.wp, c.inputs.w, c.inputs.beta, ...
                             c.inputs.z, c.inputs.c, c.inputs.qdat);
                tc.verifyEqual(zp, c.outputs.zp, 'AbsTol', c.tol, c.desc);
            end
        end

        % ------------------------------------------------------------------
        % rectmap private (elliptic functions)
        % ------------------------------------------------------------------
        function testRectPrivate(tc)
            tc.assumeTrue(tc.hasGoldens, 'goldens/ dir missing');
            f = fullfile(tc.goldenDir, 'rectmap_private.mat');
            tc.assumeTrue(exist(f,'file')==2, 'rectmap_private.mat missing');
            s = load(f);

            for k = 1:numel(s.cases_ellipkkp)
                c = s.cases_ellipkkp(k);
                [K, Kp] = rectmap.private_('ellipkkp', c.inputs.L);
                tc.verifyEqual(K,  c.outputs.K,  'AbsTol', c.tol, c.desc);
                tc.verifyEqual(Kp, c.outputs.Kp, 'AbsTol', c.tol, c.desc);
            end
            for k = 1:numel(s.cases_ellipjc)
                c = s.cases_ellipjc(k);
                [sn, cn, dn] = rectmap.private_('ellipjc', c.inputs.u, c.inputs.L);
                tc.verifyEqual(sn, c.outputs.sn, 'AbsTol', c.tol, c.desc);
                tc.verifyEqual(cn, c.outputs.cn, 'AbsTol', c.tol, c.desc);
                tc.verifyEqual(dn, c.outputs.dn, 'AbsTol', c.tol, c.desc);
            end
            for k = 1:numel(s.cases_rmap)
                c = s.cases_rmap(k);
                wp = rectmap.private_('rmap', c.inputs.zp, c.inputs.w, c.inputs.beta, ...
                          c.inputs.z, c.inputs.c, c.inputs.L, c.inputs.qdat);
                tc.verifyEqual(wp, c.outputs.wp, 'AbsTol', c.tol, c.desc);
            end
            for k = 1:numel(s.cases_rinvmap)
                c = s.cases_rinvmap(k);
                zp = rectmap.private_('rinvmap', c.inputs.wp, c.inputs.w, c.inputs.beta, ...
                             c.inputs.z, c.inputs.c, c.inputs.L, c.inputs.qdat);
                tc.verifyEqual(zp, c.outputs.zp, 'AbsTol', c.tol, c.desc);
            end
        end

        % ------------------------------------------------------------------
        % nesolve components
        % ------------------------------------------------------------------
        function testNesolve(tc)
            tc.assumeTrue(tc.hasGoldens, 'goldens/ dir missing');
            f = fullfile(tc.goldenDir, 'nesolve.mat');
            tc.assumeTrue(exist(f,'file')==2, 'nesolve.mat missing');
            s = load(f);

            for k = 1:numel(s.cases_neqrdcmp)
                c = s.cases_neqrdcmp(k);
                [M, M1, M2, sing] = sctool.neqrdcmp(c.inputs.A);
                tc.verifyEqual(M,    c.outputs.M,    'AbsTol', c.tol, c.desc);
                tc.verifyEqual(M1,   c.outputs.M1,   'AbsTol', c.tol, c.desc);
                tc.verifyEqual(M2,   c.outputs.M2,   'AbsTol', c.tol, c.desc);
                tc.verifyEqual(sing, c.outputs.sing, 'AbsTol', c.tol, c.desc);
            end
            for k = 1:numel(s.cases_nechdcmp)
                c = s.cases_nechdcmp(k);
                [L, mu] = sctool.nechdcmp(c.inputs.H, c.inputs.maxoffl);
                tc.verifyEqual(L,  c.outputs.L,  'AbsTol', c.tol, c.desc);
                tc.verifyEqual(mu, c.outputs.mu, 'AbsTol', c.tol, c.desc);
            end
            f2d = @(x) [x(1)^2 - 1; x(2)^2 - 4];
            f3d = @(x) [x(1) + x(2) - 1; x(2) + x(3) - 2; x(1)*x(3) - 0.5];
            for k = 1:numel(s.cases_nesolve)
                c = s.cases_nesolve(k);
                tc.verifyEqual(c.outputs.termcode, 1, c.desc);
                if strcmp(c.inputs.system, '2d')
                    f = f2d;
                else
                    f = f3d;
                end
                [xf, termcode] = sctool.nesolve(f, c.inputs.x0, c.inputs.details);
                % Verify residual is small, not bit-exact (solver path may vary)
                tc.verifyLessThan(norm(f(xf), inf), 1e-8, c.desc);
                tc.verifyEqual(termcode, 1, c.desc);
            end
        end

        % ------------------------------------------------------------------
        % nefdjac
        % ------------------------------------------------------------------
        function testNefdjac(tc)
            tc.assumeTrue(tc.hasGoldens, 'goldens/ dir missing');
            f = fullfile(tc.goldenDir, 'nesolve.mat');
            tc.assumeTrue(exist(f,'file')==2, 'nesolve.mat missing');
            s = load(f);

            fvec = @(x) [x(1)^2 + x(2) - 2; x(1) - x(2)^2 + 1];
            Sx = ones(2,1);
            details = zeros(16,1);
            for k = 1:numel(s.cases_nefdjac)
                c = s.cases_nefdjac(k);
                [J, ~] = sctool.nefdjac(fvec, c.inputs.F, c.inputs.x, Sx, details, 0);
                tc.verifyEqual(J, c.outputs.J, 'AbsTol', c.tol, c.desc);
            end
        end

    end % methods (Test)

end
