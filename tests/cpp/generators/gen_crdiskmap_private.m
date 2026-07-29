function gen_crdiskmap_private(outDir)
%GEN_CRDISKMAP_PRIVATE  Golden values for @crdiskmap/private functions:
%   crparam, crquad, crderiv, crmap, crinvmap.
%
%   Actual signatures (see @crdiskmap/private/*.m) differ from a naive
%   dparam-style guess: crdiskmap works with crossratios (cr), an affine
%   table (aff), a fixed conformal center (wcfix), and a quadrilateral
%   graph (Q), rather than prevertices z and constant c directly.
%     crparam(w,beta,cr0,options)             -> [w,beta,cr,aff,Q,orig,qdat]
%     crquad(z1,sing1,z,beta,qdat)             -> I
%     crderiv(zp,beta,cr,aff,wcfix,Q)          -> fp
%     crmap(zp,w,beta,cr,aff,wcfix,Q,qdat)     -> wp
%     crinvmap(wp,w,beta,cr,aff,wcfix,Q,qdat)  -> zp

opt = sctool.scmapopt('trace', 0, 'tol', 1e-12);

% crdiskmap works best with quadrilaterals
polys = { ...
    polygon([0, 1, 1+1i, 1i]),              'unit square'; ...
    polygon([0, 1.5, 1.5+1i, 0.5+1i, 0+1i]), 'irregular quad'; ...
};

cases_crparam  = struct('desc',{},'inputs',{},'outputs',{},'tol',{});
cases_crquad   = struct('desc',{},'inputs',{},'outputs',{},'tol',{});
cases_crderiv  = struct('desc',{},'inputs',{},'outputs',{},'tol',{});
cases_crmap    = struct('desc',{},'inputs',{},'outputs',{},'tol',{});
cases_crinvmap = struct('desc',{},'inputs',{},'outputs',{},'tol',{});

for k = 1:size(polys,1)
    p = polys{k,1};
    try
        m = crdiskmap(p, opt);
    catch
        fprintf('  crdiskmap skipped for %s\n', polys{k,2});
        continue
    end

    % Extract internals from the solved map object
    w     = vertex(m.polygon);
    beta  = angle(m.polygon) - 1;
    cr    = m.crossratio;
    aff   = m.affine;
    Q     = m.qlgraph;
    qdat  = m.qdata;
    wcfix = m.center{2};

    cases_crparam(end+1) = struct( ...
        'desc',    ['crparam / ' polys{k,2}], ...
        'inputs',  struct('w', w, 'beta', beta), ...
        'outputs', struct('cr', cr, 'qdat', qdat), ...
        'tol',     1e-10); %#ok<AGROW>

    % crquad: integrate from a singular prevertex toward 0
    z = crdiskmap.private_('crembed', cr, Q, 1);
    I = crdiskmap.private_('crquad', z(1), 1, z, beta, qdat);
    cases_crquad(end+1) = struct( ...
        'desc',    ['crquad / ' polys{k,2}], ...
        'inputs',  struct('z1', z(1), 'sing1', 1, 'z', z, 'beta', beta, 'qdat', qdat), ...
        'outputs', struct('I', I), ...
        'tol',     1e-11); %#ok<AGROW>

    % Q (the quadrilateral graph) is a struct of purely numeric matrices
    % (qlvert, qledge, adjacent); flatten its fields into the case inputs
    % under a "Q" prefix so the plain-text golden exporter (which only
    % handles numeric/char fields, not nested structs) can serialize it.
    % This unblocks crderiv/crmap/crinvmap without needing to port
    % crqgraph/crtriang/crcdt (the triangulation that builds Q from
    % scratch) -- exactly mirroring how rderiv/rmap/rinvmap took z/c/L
    % directly instead of requiring rparam.
    Qfields = struct('Qqlvert', double(Q.qlvert), 'Qqledge', double(Q.qledge), ...
                      'Qadjacent', double(Q.adjacent));

    % crderiv / crmap / crinvmap.
    % NOTE: for a polygon with only one quadrilateral (n3==1, e.g. a plain
    % unit square), MATLAB's `[~,idx] = min(abs(zl))` in crderiv.m/crmap.m
    % silently misbehaves for multi-point ZP: zl is 1-by-m there, and
    % MATLAB's min() treats a 1-by-m matrix as a *vector*, collapsing to a
    % single scalar index instead of an m-vector of per-point embedding
    % choices. Downstream `mask = (idx==q)` then only touches one of the m
    % points, leaving the others at their zero-initialized default -- a
    % genuine latent bug in the original .m code for this degenerate case,
    % not a generator or C++-port issue. Side-step it here (rather than
    % "fix" the .m source) by using a single query point for n3==1
    % polygons and reserving multi-point batches for n3>1, where the bug
    % cannot trigger.
    if numel(cr) == 1
        zp = 0.1+0.2i;
    else
        zp = [0.1+0.2i, -0.15+0.32i].';
    end
    fp = crdiskmap.private_('crderiv', zp, beta, cr, aff, wcfix, Q);
    cases_crderiv(end+1) = struct( ...
        'desc',    ['crderiv / ' polys{k,2}], ...
        'inputs',  mergeStructs(struct('zp', zp, 'beta', beta, 'cr', cr, 'aff', aff, 'wcfix', wcfix), Qfields), ...
        'outputs', struct('fp', fp), ...
        'tol',     1e-11); %#ok<AGROW>

    % crmap
    wp = crdiskmap.private_('crmap', zp, w, beta, cr, aff, wcfix, Q, qdat);
    cases_crmap(end+1) = struct( ...
        'desc',    ['crmap / ' polys{k,2}], ...
        'inputs',  mergeStructs(struct('zp', zp, 'w', w, 'beta', beta, 'cr', cr, 'aff', aff, 'wcfix', wcfix, 'qdat', qdat), Qfields), ...
        'outputs', struct('wp', wp), ...
        'tol',     1e-11); %#ok<AGROW>

    % crinvmap
    zp_inv = crdiskmap.private_('crinvmap', wp, w, beta, cr, aff, wcfix, Q, qdat);
    cases_crinvmap(end+1) = struct( ...
        'desc',    ['crinvmap / ' polys{k,2}], ...
        'inputs',  mergeStructs(struct('wp', wp, 'w', w, 'beta', beta, 'cr', cr, 'aff', aff, 'wcfix', wcfix, 'qdat', qdat), Qfields), ...
        'outputs', struct('zp', zp_inv), ...
        'tol',     1e-10); %#ok<AGROW>
end

save(fullfile(outDir, 'crdiskmap_private.mat'), ...
    'cases_crparam', 'cases_crquad', 'cases_crderiv', 'cases_crmap', 'cases_crinvmap');
fprintf('  crdiskmap_private: crparam=%d crquad=%d crderiv=%d crmap=%d crinvmap=%d\n', ...
    numel(cases_crparam), numel(cases_crquad), numel(cases_crderiv), ...
    numel(cases_crmap), numel(cases_crinvmap));
end

function s = mergeStructs(a, b)
s = a;
fn = fieldnames(b);
for i = 1:numel(fn)
    s.(fn{i}) = b.(fn{i});
end
end
