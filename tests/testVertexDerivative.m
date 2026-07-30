classdef testVertexDerivative < matlab.unittest.TestCase
%TESTVERTEXDERIVATIVE  Finite-difference baseline for df/dw_i (shape derivative).
%
%   Computes the derivative of the SC rectangle map f with respect to each
%   polygon vertex w_i at fixed source-domain points, using central finite
%   differences.  The perturbed map is re-solved (warm-started) after each
%   perturbation.
%
%   Three evaluation-point classes are tested:
%
%     Interior   — three points strictly inside the rectangle
%     Edge       — one point near each of the four rectangle edges
%                  (slightly inside so points remain in all perturbed rects)
%     Corners    — one point near each of the four rectangle corners
%                  (slightly inside for the same reason)
%
%   Each test verifies FD convergence by checking that step h and step h/2
%   agree to a relative+absolute tolerance.  These tests are the numerical
%   baseline for a future analytic implementation via the implicit function
%   theorem.
%
%   Output convention of vertexDerivFD:
%     dfdwre(k,i)  =  df/d(Re w_i)  at zp(k)   [complex-valued]
%     dfdwim(k,i)  =  df/d(Im w_i)  at zp(k)   [complex-valued]

    properties
        m        % base rectmap
        p        % polygon as stored in the map (post-scfix ordering)
        K        % rectangle half-width   (corners at Re = ±K)
        Kp       % rectangle height       (corners at Im = 0 and Im = Kp)
        zRect    % 4 rectangle corner prevertices [K; K+iKp; -K+iKp; -K]

        zInterior       % 3 strictly interior evaluation points
        zEdge           % 4 points, one near each edge (slightly inside)
        zCorners        % 4 points, one near each corner (slightly inside)
        zPrevertex      % exact non-corner prevertex positions
        zNearPrevertex  % non-corner prevertices offset slightly into interior
    end

    methods (TestClassSetup)
        function setup(testCase)
            opt        = sctool.scmapopt('trace', 0, 'tol', 1e-10);
            p0         = polygon([4, 2i, -2+4i, -3, -3-1i, 2-2i]);
            testCase.m = rectmap(p0, 1:4, opt);

            % Use the polygon stored in the map (post-scfix) so that vertex
            % indices are consistent with the continuation call.
            testCase.p = polygon(testCase.m);

            zr           = rectangle(testCase.m);   % [K; K+iKp; -K+iKp; -K]
            testCase.K   = max(real(zr));
            testCase.Kp  = max(imag(zr));
            testCase.zRect = zr;

            K  = testCase.K;
            Kp = testCase.Kp;

            testCase.zInterior = [ 0.30*K + 0.40*Kp*1i, ...
                                  -0.50*K + 0.60*Kp*1i, ...
                                   0.10*K + 0.20*Kp*1i ];

            % Near-edge points: offset 2–5% inside each side so that they
            % remain inside the rectangle when the polygon is perturbed
            % (K and Kp shift by O(h) = O(1e-4) per perturbation).
            testCase.zEdge = [  0.50*K  + 0.02*Kp*1i, ...   % near bottom (Im=0)
                                0.95*K  + 0.50*Kp*1i, ...   % near right  (Re=K)
                                0       + 0.98*Kp*1i, ...   % near top    (Im=Kp)
                               -0.95*K  + 0.50*Kp*1i ];     % near left   (Re=-K)

            % Near-corner points: 5% inside along both dimensions from each
            % rectangle corner.  The exact corners are singular points of the
            % SC map in z; their position shifts with every vertex perturbation,
            % so evaluating at a fixed exact-corner point gives unreliable FD.
            d = 0.05;
            testCase.zCorners = [  (1-d)*K + d*Kp*1i, ...          % near K
                                   (1-d)*K + (1-d)*Kp*1i, ...      % near K+iKp
                                  -(1-d)*K + (1-d)*Kp*1i, ...      % near -K+iKp
                                  -(1-d)*K + d*Kp*1i ];            % near -K

            % Non-corner prevertex positions: boundary points where the SC map
            % has a z-singularity (like rectangle corners, but on the edges).
            zAll  = testCase.m.prevertex;
            cidx  = corners(testCase.m);
            ncidx = setdiff((1:numel(zAll))', cidx(:));
            testCase.zPrevertex = zAll(ncidx);

            % Near-prevertex: same positions offset 2% of the perpendicular
            % rectangle dimension into the interior.
            eps_K  = 0.02 * K;
            eps_Kp = 0.02 * Kp;
            edge_tol = 1e-8 * max(K, Kp);
            zPV = testCase.zPrevertex;
            zNPV = zPV;
            for k = 1:numel(zPV)
                z = zPV(k);
                if     abs(imag(z))      < edge_tol,  zNPV(k) = z + 1i*eps_Kp;  % bottom
                elseif abs(imag(z) - Kp) < edge_tol,  zNPV(k) = z - 1i*eps_Kp;  % top
                elseif abs(real(z) - K)  < edge_tol,  zNPV(k) = z - eps_K;       % right
                else,                                  zNPV(k) = z + eps_K;       % left
                end
            end
            testCase.zNearPrevertex = zNPV;
        end
    end

    % -----------------------------------------------------------------
    methods (Access = private)

        function [dfdwre, dfdwim] = vertexDerivFD(testCase, zp, h)
        %VERTEXDERIVFD  Central FD approximation of df/dw_i for all vertices.
        %
        %   [dfdwre, dfdwim] = vertexDerivFD(testCase, zp, h)
        %
        %   For each polygon vertex index i (1 … n):
        %     dfdwre(k,i)  ≈  ∂f/∂(Re w_i)  at  zp(k)
        %     dfdwim(k,i)  ≈  ∂f/∂(Im w_i)  at  zp(k)
        %
        %   polygon(wp) recomputes the turning angles from the perturbed
        %   vertices, which is the correct perturbation for shape sensitivity.
        %
        %   evalRaw bypasses the isinpoly domain check so that near-edge and
        %   near-corner points are evaluated even when the perturbed map's
        %   rectangle boundary has shifted slightly.

            p    = testCase.p;
            m    = testCase.m;
            n    = length(p);
            zp   = zp(:);
            nz   = numel(zp);
            w    = vertex(p);

            dfdwre = zeros(nz, n);
            dfdwim = zeros(nz, n);

            for i = 1:n
                % ∂f / ∂(Re w_i)
                wp = w;  wp(i) = wp(i) + h;
                wm = w;  wm(i) = wm(i) - h;
                fp = evalRaw(rectmap(m, polygon(wp)), zp);
                fm = evalRaw(rectmap(m, polygon(wm)), zp);
                dfdwre(:, i) = (fp - fm) / (2*h);

                % ∂f / ∂(Im w_i)
                wp = w;  wp(i) = wp(i) + 1i*h;
                wm = w;  wm(i) = wm(i) - 1i*h;
                fp = evalRaw(rectmap(m, polygon(wp)), zp);
                fm = evalRaw(rectmap(m, polygon(wm)), zp);
                dfdwim(:, i) = (fp - fm) / (2*h);
            end
        end

    end

    % -----------------------------------------------------------------
    methods (Test)

        % --- Analytic vs FD baseline ---

        function testAnalyticInterior(testCase)
        %TESTANALYTICINTERIOR  vertexDeriv matches FD at interior points.
            h = 1e-4;
            [fdre, fdim] = testCase.vertexDerivFD(testCase.zInterior, h);
            [dre, dim]   = vertexDeriv(testCase.m, testCase.zInterior);
            testCase.verifyEqual(dre, fdre, 'AbsTol', 1e-5, 'RelTol', 1e-2, ...
                'Interior: vertexDeriv dfdwre mismatch with FD');
            testCase.verifyEqual(dim, fdim, 'AbsTol', 1e-5, 'RelTol', 1e-2, ...
                'Interior: vertexDeriv dfdwim mismatch with FD');
        end

        function testAnalyticEdge(testCase)
        %TESTANALYTICEDGE  vertexDeriv matches FD at near-edge points.
            h = 1e-4;
            [fdre, fdim] = testCase.vertexDerivFD(testCase.zEdge, h);
            [dre, dim]   = vertexDeriv(testCase.m, testCase.zEdge);
            testCase.verifyEqual(dre, fdre, 'AbsTol', 1e-5, 'RelTol', 1e-2, ...
                'Edge: vertexDeriv dfdwre mismatch with FD');
            testCase.verifyEqual(dim, fdim, 'AbsTol', 1e-5, 'RelTol', 1e-2, ...
                'Edge: vertexDeriv dfdwim mismatch with FD');
        end

        function testAnalyticNearPrevertex(testCase)
        %TESTANALYTICNEARPREVERTEX  vertexDeriv matches FD at near-prevertex points.
            h = 1e-4;
            [fdre, fdim] = testCase.vertexDerivFD(testCase.zNearPrevertex, h);
            [dre, dim]   = vertexDeriv(testCase.m, testCase.zNearPrevertex);
            testCase.verifyEqual(dre, fdre, 'AbsTol', 1e-5, 'RelTol', 1e-2, ...
                'NearPrevertex: vertexDeriv dfdwre mismatch with FD');
            testCase.verifyEqual(dim, fdim, 'AbsTol', 1e-5, 'RelTol', 1e-2, ...
                'NearPrevertex: vertexDeriv dfdwim mismatch with FD');
        end

        % --- FD convergence baseline (unchanged) ---

        function testInteriorFDConverges(testCase)
        %TESTINTERIORFDCONVERGES  FD shape derivatives converge at interior pts.
            h = 1e-4;
            [d1re, d1im] = testCase.vertexDerivFD(testCase.zInterior, h);
            [d2re, d2im] = testCase.vertexDerivFD(testCase.zInterior, h/2);
            testCase.verifyEqual(d1re, d2re, 'AbsTol', 1e-6, 'RelTol', 1e-3, ...
                'Interior dfdwre: step h and h/2 disagree');
            testCase.verifyEqual(d1im, d2im, 'AbsTol', 1e-6, 'RelTol', 1e-3, ...
                'Interior dfdwim: step h and h/2 disagree');
        end

        function testEdgeFDConverges(testCase)
        %TESTEDGEFDCONVERGES  FD shape derivatives converge at near-edge pts.
            h = 1e-4;
            [d1re, d1im] = testCase.vertexDerivFD(testCase.zEdge, h);
            [d2re, d2im] = testCase.vertexDerivFD(testCase.zEdge, h/2);
            testCase.verifyEqual(d1re, d2re, 'AbsTol', 1e-6, 'RelTol', 1e-3, ...
                'Edge dfdwre: step h and h/2 disagree');
            testCase.verifyEqual(d1im, d2im, 'AbsTol', 1e-6, 'RelTol', 1e-3, ...
                'Edge dfdwim: step h and h/2 disagree');
        end

        function testCornerFDConverges(testCase)
        %TESTCORNERFDCONVERGES  FD shape derivatives converge at near-corner pts.
            h = 1e-4;
            [d1re, d1im] = testCase.vertexDerivFD(testCase.zCorners, h);
            [d2re, d2im] = testCase.vertexDerivFD(testCase.zCorners, h/2);
            testCase.verifyEqual(d1re, d2re, 'AbsTol', 1e-6, 'RelTol', 1e-3, ...
                'Corner dfdwre: step h and h/2 disagree');
            testCase.verifyEqual(d1im, d2im, 'AbsTol', 1e-6, 'RelTol', 1e-3, ...
                'Corner dfdwim: step h and h/2 disagree');
        end

        function testNearPrevertexFDConverges(testCase)
        %TESTNEARPREVERTEXFDCONVERGES  FD at near-prevertex positions.
        %   Each point is offset 2% of the perpendicular rectangle dimension
        %   into the interior from the corresponding non-corner prevertex.
        %
        %   Note: evaluating at exact non-corner prevertex positions produces the
        %   same power-law divergence as exact rectangle corners.  When a polygon
        %   vertex is perturbed, the prevertex position shifts by O(h), placing
        %   the fixed evaluation point near (but not at) the singularity of the
        %   perturbed map.  The FD then scales as h^(alpha-1) where alpha is the
        %   polygon turning angle, which diverges for alpha < 1.  The ratio
        %   d(h/2)/d(h) -> 2^alpha (~1.35 for the test polygon) rather than 1.
        %   Offsetting 2% into the interior eliminates this effect.
            h = 1e-4;
            [d1re, d1im] = testCase.vertexDerivFD(testCase.zNearPrevertex, h);
            [d2re, d2im] = testCase.vertexDerivFD(testCase.zNearPrevertex, h/2);
            testCase.verifyEqual(d1re, d2re, 'AbsTol', 1e-6, 'RelTol', 1e-3, ...
                'Near-prevertex dfdwre: step h and h/2 disagree');
            testCase.verifyEqual(d1im, d2im, 'AbsTol', 1e-6, 'RelTol', 1e-3, ...
                'Near-prevertex dfdwim: step h and h/2 disagree');
        end

    end

end
