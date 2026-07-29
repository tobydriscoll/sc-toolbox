function gen_crdiskmap_class(outDir)
%GEN_CRDISKMAP_CLASS  Golden values for the @crdiskmap class-level methods:
%   eval, evalinv, evaldiff, accuracy.

opt = sctool.scmapopt('trace', 0, 'tol', 1e-12);

% Same polygons as gen_crdiskmap_private.m (the unit square's n3==1
% quadrilateral graph is degenerate for some private functions, but the
% class-level methods here only ever use a single point/embedding per
% call, so it isn't an issue; kept for parity with the existing fixtures).
polys = { ...
    polygon([0, 1, 1+1i, 1i]),              'unit square'; ...
    polygon([0, 1.5, 1.5+1i, 0.5+1i, 0+1i]), 'irregular quad'; ...
};

cases_eval     = struct('desc',{},'inputs',{},'outputs',{},'tol',{});
cases_evalinv  = struct('desc',{},'inputs',{},'outputs',{},'tol',{});
cases_evaldiff = struct('desc',{},'inputs',{},'outputs',{},'tol',{});
cases_accuracy = struct('desc',{},'inputs',{},'outputs',{},'tol',{});

for k = 1:size(polys,1)
    p = polys{k,1};
    m = crdiskmap(p, opt);

    w     = vertex(m.polygon);
    beta  = angle(m.polygon) - 1;
    cr    = m.crossratio;
    aff   = m.affine;
    Q     = m.qlgraph;
    qdat  = m.qdata;
    wcfix = m.center{2};

    % Q (the quadrilateral graph) is a struct of plain numeric matrices;
    % flatten it the same way gen_crdiskmap_private.m does so the
    % plain-text exporter (which only handles numeric/char top-level
    % fields) can serialize it.
    Qfields = struct('Qqlvert', double(Q.qlvert), 'Qqledge', double(Q.qledge), ...
                      'Qadjacent', double(Q.adjacent));

    % NOTE: as documented in CPP_PLAN.md 9.3, MATLAB's crderiv.m/crmap.m
    % have a latent bug for n3==1 polygons (e.g. this unit square) where
    % batching multiple query points silently zeros out all but the first.
    % Use a single point there and reserve multi-point batches for n3>1.
    if numel(cr) == 1
        zp = 0.1+0.2i;
    else
        zp = [0.1+0.2i, -0.15+0.32i].';
    end
    wp = eval(m, zp);
    cases_eval(end+1) = struct( ...
        'desc',    ['eval / ' polys{k,2}], ...
        'inputs',  mergeStructs(struct('w', w, 'beta', beta, 'cr', cr, 'aff', aff, 'wcfix', wcfix, 'qdat', qdat, 'zp', zp), Qfields), ...
        'outputs', struct('wp', wp), ...
        'tol',     1e-10); %#ok<AGROW>

    zp_inv = evalinv(m, wp);
    cases_evalinv(end+1) = struct( ...
        'desc',    ['evalinv / ' polys{k,2}], ...
        'inputs',  mergeStructs(struct('w', w, 'beta', beta, 'cr', cr, 'aff', aff, 'wcfix', wcfix, 'qdat', qdat, 'wp', wp), Qfields), ...
        'outputs', struct('zp', zp_inv), ...
        'tol',     1e-9); %#ok<AGROW>

    fp = evaldiff(m, zp);
    cases_evaldiff(end+1) = struct( ...
        'desc',    ['evaldiff / ' polys{k,2}], ...
        'inputs',  mergeStructs(struct('beta', beta, 'cr', cr, 'aff', aff, 'wcfix', wcfix, 'zp', zp), Qfields), ...
        'outputs', struct('fp', fp), ...
        'tol',     1e-10); %#ok<AGROW>

    acc = accuracy(m);
    cases_accuracy(end+1) = struct( ...
        'desc',    ['accuracy / ' polys{k,2}], ...
        'inputs',  mergeStructs(struct('w', w, 'beta', beta, 'cr', cr, 'qdat', qdat), Qfields), ...
        'outputs', struct('acc', acc), ...
        'tol',     1e-9); %#ok<AGROW>
end

save(fullfile(outDir, 'crdiskmap_class.mat'), ...
    'cases_eval', 'cases_evalinv', 'cases_evaldiff', 'cases_accuracy');
fprintf('  crdiskmap_class: eval=%d evalinv=%d evaldiff=%d accuracy=%d\n', ...
    numel(cases_eval), numel(cases_evalinv), numel(cases_evaldiff), numel(cases_accuracy));
end

function s = mergeStructs(a, b)
s = a;
fn = fieldnames(b);
for i = 1:numel(fn)
    s.(fn{i}) = b.(fn{i});
end
end
