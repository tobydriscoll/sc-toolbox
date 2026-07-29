function gen_extermap_private(outDir)
%GEN_EXTERMAP_PRIVATE  Golden values for @extermap/private functions:
%   deparam, dequad, dederiv, demap, deinvmap.

opt = sctool.scmapopt('trace', 0, 'tol', 1e-12);

% Exterior maps require outward-oriented polygons. The unit square was
% replaced with an asymmetric quadrilateral: deparam's parameter problem
% for the square is fully symmetric (all beta equal), making its Jacobian
% singular/degenerate at the solution -- confirmed against MATLAB directly
% that *consecutive* deparam(w,beta,...) calls with bit-identical inputs
% return different (rotated/mirrored) but equally-valid prevertex orderings
% from run to run, i.e. genuinely non-reproducible MATLAB output, not a
% generator or C++-port bug. The asymmetric quadrilateral below has a
% unique, well-conditioned solution and reproduces identically across runs.
polys = { ...
    polygon([0, 2i, 2+3i, 3]),     'asymmetric quad (CCW exterior)'; ...
    polygon([0, 1, 2, 2+1i, 1i]),  'pentagon exterior'; ...
};

cases_deparam  = struct('desc',{},'inputs',{},'outputs',{},'tol',{});
cases_dequad   = struct('desc',{},'inputs',{},'outputs',{},'tol',{});
cases_dederiv  = struct('desc',{},'inputs',{},'outputs',{},'tol',{});
cases_demap    = struct('desc',{},'inputs',{},'outputs',{},'tol',{});
cases_deinvmap = struct('desc',{},'inputs',{},'outputs',{},'tol',{});

for k = 1:size(polys,1)
    p = polys{k,1};

    % Build deparam's literal input the same way extermap.m's constructor
    % does: `w_de = flipud(vertex(poly)); beta_de = 1-flipud(angle(poly));`
    % (note: 1-alpha, the *negation* of the usual alpha-1 beta, needed so
    % that sccheck's orientation-sum check passes for the clockwise-
    % traversed exterior polygon), then scfix.
    %
    % IMPORTANT: this generator calls deparam directly (via
    % extermap.private_) rather than constructing `extermap(p,opt)` and
    % extracting prevertex/constant/qdata from the resulting object.
    % Verified directly against MATLAB that those two routes are NOT
    % equivalent here: calling deparam with bit-identical (w,beta,tol,opt)
    % reproducibly returns a *different* (but equally valid -- satisfies
    % the same residual equations to the same tolerance) prevertex
    % configuration depending on whether the call goes through
    % extermap(p,opt)'s constructor or is made directly/via feval, in BOTH
    % the n=4 and n=5 test cases here. This was tracked down to MATLAB's
    % own deparam/nesolve internals (not this generator or scfix): the
    % parameter problem's Jacobian is sufficiently close to degenerate that
    % which of several valid roots nesolve converges to is sensitive to
    % low-level floating-point execution-path differences between
    % interpreted/feval and compiled call contexts, even with identical
    % inputs (confirmed bit-for-bit via hex dumps of w/beta/tol). The
    % direct/feval call path is fully deterministic and reproducible on its
    % own (10+ repeated calls, identical result every time), so it is used
    % here for a stable, well-defined golden value -- exactly mirroring how
    % dparam/hpparam's goldens call scfix+dparam/hpparam directly rather
    % than extracting from a constructed diskmap/hplmap object.
    try
        w_de = flipud(vertex(p));
        beta_de = 1 - flipud(angle(p));
        [w_de, beta_de] = scfix('de', w_de, beta_de);
        [z, c, qdat] = extermap.private_('deparam', w_de, beta_de, [], opt);
    catch
        fprintf('  extermap skipped for %s\n', polys{k,2});
        continue
    end

    % Un-negate/un-flip back to the normal (alpha-1) convention that
    % dequad/dederiv/demap/deinvmap expect, matching
    % `poly = polygon(flipud(w),1-flipud(beta))` in extermap.m.
    w = flipud(w_de);
    beta = 1 - flipud(beta_de);

    cases_deparam(end+1) = struct( ...
        'desc',    ['deparam / ' polys{k,2}], ...
        'inputs',  struct('w', w_de, 'beta', beta_de), ...
        'outputs', struct('z', z, 'c', c, 'qdat', qdat), ...
        'tol',     1e-10); %#ok<AGROW>

    % dequad: arc on unit circle from one prevertex toward interior
    finIdx = find(~isinf(z));
    if numel(finIdx) >= 2
        za = z(finIdx(1));
        zb = mean([z(finIdx(1)), z(finIdx(2))]);
        zb = zb / abs(zb);  % snap to circle
        I  = extermap.private_('dequad', za, zb, finIdx(1), z, beta, qdat);
        cases_dequad(end+1) = struct( ...
            'desc',    ['dequad / ' polys{k,2}], ...
            'inputs',  struct('z1',za,'z2',zb,'sing1',finIdx(1),'z',z,'beta',beta,'qdat',qdat), ...
            'outputs', struct('I', I), ...
            'tol',     1e-11); %#ok<AGROW>
    end

    % dederiv: points outside the unit disk, scaled well clear of it so
    % their forward images are unambiguously outside the polygon (close-in
    % points can map close to an edge and make deinvmap's automatic
    % z0-search/ODE-continuation converge poorly).
    zp = 8 * [1.2+0.3i, -1.5+0.8i, 2-0.5i].';
    fp = extermap.private_('dederiv', zp, z, beta, c);
    cases_dederiv(end+1) = struct( ...
        'desc',    ['dederiv / ' polys{k,2}], ...
        'inputs',  struct('zp', zp, 'z', z, 'beta', beta, 'c', c), ...
        'outputs', struct('fp', fp), ...
        'tol',     1e-11); %#ok<AGROW>

    % demap
    wp = extermap.private_('demap', zp, w, beta, z, c, qdat);
    cases_demap(end+1) = struct( ...
        'desc',    ['demap / ' polys{k,2}], ...
        'inputs',  struct('zp', zp, 'w', w, 'beta', beta, 'z', z, 'c', c, 'qdat', qdat), ...
        'outputs', struct('wp', wp), ...
        'tol',     1e-11); %#ok<AGROW>

    % deinvmap: points close to the unit circle (the natural domain for the
    % ODE-continuation z0 search). Far-field points like the dederiv/demap
    % ones above land deinvmap's Newton iteration in a non-convergent regime
    % (verified against MATLAB itself: even with tol=1e-12/maxiter=80 the
    % residual stalls around 0.02-0.3 for those far points), so this case
    % uses its own well-conditioned test points and a tightened tolerance/
    % iteration budget to get a tightly-converged, reproducible golden value.
    zp_di = [1.05+0.1i, 1.2-0.1i, -1.1+0.2i].';
    wp_di = extermap.private_('demap', zp_di, w, beta, z, c, qdat);
    invopt = [0, 1e-12, 80];
    zp_inv = extermap.private_('deinvmap', wp_di, w, beta, z, c, qdat, [], invopt);
    cases_deinvmap(end+1) = struct( ...
        'desc',    ['deinvmap / ' polys{k,2}], ...
        'inputs',  struct('wp', wp_di, 'w', w, 'beta', beta, 'z', z, 'c', c, 'qdat', qdat), ...
        'outputs', struct('zp', zp_inv), ...
        'tol',     1e-10); %#ok<AGROW>
end

save(fullfile(outDir, 'extermap_private.mat'), ...
    'cases_deparam', 'cases_dequad', 'cases_dederiv', 'cases_demap', 'cases_deinvmap');
fprintf('  extermap_private: deparam=%d dequad=%d dederiv=%d demap=%d deinvmap=%d\n', ...
    numel(cases_deparam), numel(cases_dequad), numel(cases_dederiv), ...
    numel(cases_demap), numel(cases_deinvmap));
end
