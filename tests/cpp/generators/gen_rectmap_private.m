function gen_rectmap_private(outDir)
%GEN_RECTMAP_PRIVATE  Golden values for @rectmap/private functions:
%   ellipkkp, ellipjc, r2strip, rparam, rderiv, rmap, rinvmap.
%
%   Actual signatures (see @rectmap/private/*.m):
%     rparam(w,beta,cnr,z0,options)        -> [z,c,L,qdat]
%     rderiv(zp,z,beta,c,L,zs)             -> fprime   (zs optional)
%     rmap(zp,w,beta,z,c,L,qdat)           -> wp        (note: w precedes z)
%     rinvmap(wp,w,beta,z,c,L,qdat)        -> zp        (L precedes qdat)

opt = sctool.scmapopt('trace', 0, 'tol', 1e-12);

% --- ellipkkp ---
cases_ellipkkp = struct('desc',{},'inputs',{},'outputs',{},'tol',{});
Lvals = [0.1, 0.3, 0.5, 1/sqrt(2), 0.7, 0.9, 0.99];
for k = 1:numel(Lvals)
    L = Lvals(k);
    [K, Kp] = rectmap.private_('ellipkkp', L);
    cases_ellipkkp(end+1) = struct( ...
        'desc',    sprintf('ellipkkp L=%.4g', L), ...
        'inputs',  struct('L', L), ...
        'outputs', struct('K', K, 'Kp', Kp), ...
        'tol',     1e-13); %#ok<AGROW>
end

% --- ellipjc ---
cases_ellipjc = struct('desc',{},'inputs',{},'outputs',{},'tol',{});
L = 0.5;
[K, ~] = rectmap.private_('ellipkkp', L);
uSpecs = { ...
    0.3,         'real argument'; ...
    0.3 + 0.2i,  'complex argument'; ...
    0 + 0.5i,    'imaginary argument'; ...
    K/2,         'real = K/2'; ...
    K/2 + 0.3i,  'complex near quarter-period'; ...
};
for k = 1:size(uSpecs, 1)
    u = uSpecs{k,1};
    [sn, cn, dn] = rectmap.private_('ellipjc', u, L);
    cases_ellipjc(end+1) = struct( ...
        'desc',    ['ellipjc L=0.5 ' uSpecs{k,2}], ...
        'inputs',  struct('u', u, 'L', L), ...
        'outputs', struct('sn', sn, 'cn', cn, 'dn', dn), ...
        'tol',     1e-13); %#ok<AGROW>
end

% --- rparam + rderiv + rmap + rinvmap ---
% Polygons/corner choices copied from tests/fixturePolygons.m's 'hex6' and
% 'L' rect fixtures (known to converge cleanly in the regression suite) —
% an ad hoc L-shape with corners [1,2,5,6] failed to converge ("Nonlinear
% equations solver did not terminate normally").
rectSpecs = { ...
    polygon([4, 2i, -2+4i, -3, -3-1i, 2-2i]), 1:4,      'hex6 / corners 1234'; ...
    polygon([1i, -1+1i, -1-1i, 1-1i, 1, 0]),  [2,4,5,1], 'L-shape / corners 2451'; ...
};

cases_rparam  = struct('desc',{},'inputs',{},'outputs',{},'tol',{});
cases_rderiv  = struct('desc',{},'inputs',{},'outputs',{},'tol',{});
cases_rmap    = struct('desc',{},'inputs',{},'outputs',{},'tol',{});
cases_rinvmap = struct('desc',{},'inputs',{},'outputs',{},'tol',{});

for k = 1:size(rectSpecs,1)
    p       = rectSpecs{k,1};
    corners = rectSpecs{k,2};
    try
        m = rectmap(p, corners, opt);
    catch
        fprintf('  rectmap skipped for %s\n', rectSpecs{k,3});
        continue
    end
    w    = vertex(m.polygon);
    beta = angle(m.polygon) - 1;
    z    = m.prevertex;
    c    = m.constant;
    L    = m.stripL;
    qdat = m.qdata;

    % rectmap.m's constructor stores `m.polygon` AFTER applying scfix to
    % (w,beta,corner) -- so `vertex(m.polygon)`/`angle(m.polygon)-1` are the
    % POST-scfix w/beta, but the ORIGINAL (pre-scfix) `corners` local
    % variable is never updated. Pairing them (as a naive extraction would)
    % mismatches whenever scfix actually renumbers (confirmed directly
    % against MATLAB: corners=[2,4,5,1] for the L-shape case becomes
    % corner2=[1,3,4,6] post-scfix, and calling rparam with the original
    % [2,4,5,1] against this w/beta does NOT reproduce the recorded z/c/L).
    % Recompute the post-scfix corners directly so this golden case is
    % self-consistent, mirroring the deparam/dparam generator fixes above.
    [~, ~, corners_post] = scfix('r', vertex(p), angle(p) - 1, corners);

    cases_rparam(end+1) = struct( ...
        'desc',    ['rparam / ' rectSpecs{k,3}], ...
        'inputs',  struct('w', w, 'beta', beta, 'corners', corners_post), ...
        'outputs', struct('z', z, 'c', c, 'L', L, 'qdat', qdat), ...
        'tol',     1e-9); %#ok<AGROW>

    % rderiv and rmap: points in the rectangle [0,K] x [0,Kp] spanned by z
    K  = max(real(z));
    Kp = max(imag(z));
    zp = [K/4 + Kp/4*1i, K*3/4 + Kp/2*1i].';
    fp = rectmap.private_('rderiv', zp, z, beta, c, L);
    wp = rectmap.private_('rmap', zp, w, beta, z, c, L, qdat);
    cases_rderiv(end+1) = struct( ...
        'desc',    ['rderiv / ' rectSpecs{k,3}], ...
        'inputs',  struct('zp', zp, 'z', z, 'beta', beta, 'c', c, 'L', L), ...
        'outputs', struct('fp', fp), ...
        'tol',     1e-11); %#ok<AGROW>
    cases_rmap(end+1) = struct( ...
        'desc',    ['rmap / ' rectSpecs{k,3}], ...
        'inputs',  struct('zp', zp, 'w', w, 'beta', beta, 'z', z, 'c', c, 'L', L, 'qdat', qdat), ...
        'outputs', struct('wp', wp), ...
        'tol',     1e-11); %#ok<AGROW>

    % rinvmap
    zp_inv = rectmap.private_('rinvmap', wp, w, beta, z, c, L, qdat);
    cases_rinvmap(end+1) = struct( ...
        'desc',    ['rinvmap / ' rectSpecs{k,3}], ...
        'inputs',  struct('wp', wp, 'w', w, 'beta', beta, 'z', z, 'c', c, 'L', L, 'qdat', qdat), ...
        'outputs', struct('zp', zp_inv), ...
        'tol',     1e-10); %#ok<AGROW>
end

save(fullfile(outDir, 'rectmap_private.mat'), ...
    'cases_ellipkkp', 'cases_ellipjc', ...
    'cases_rparam', 'cases_rderiv', 'cases_rmap', 'cases_rinvmap');
fprintf('  rectmap_private: ellipkkp=%d ellipjc=%d rparam=%d rderiv=%d rmap=%d rinvmap=%d\n', ...
    numel(cases_ellipkkp), numel(cases_ellipjc), numel(cases_rparam), ...
    numel(cases_rderiv), numel(cases_rmap), numel(cases_rinvmap));
end
