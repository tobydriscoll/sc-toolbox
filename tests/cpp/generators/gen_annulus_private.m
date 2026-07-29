function gen_annulus_private(outDir)
%GEN_ANNULUS_PRIVATE  Golden values for @annulusmap/private functions:
%   qinit, and the full annulusmap eval/evalinv pipeline (constructor +
%   zdsc/wdsc), used to validate the C++ AnnulusMap port end to end.

% Doubly-connected polygon pairs: [outer polygon, inner polygon]
pairSpecs = { ...
    polygon([0, 2, 2+2i, 2i]),   polygon([0.5+0.5i, 1.5+0.5i, 1.5+1.5i, 0.5+1.5i]), ...
      'square annulus (outer 2x2, inner 1x1)'; ...
};

cases_anneval    = struct('desc',{},'inputs',{},'outputs',{},'tol',{});
cases_annevalinv = struct('desc',{},'inputs',{},'outputs',{},'tol',{});
cases_qinit      = struct('desc',{},'inputs',{},'outputs',{},'tol',{});

for k = 1:size(pairSpecs,1)
    outerPoly = pairSpecs{k,1};
    innerPoly = pairSpecs{k,2};
    desc      = pairSpecs{k,3};
    try
        m = annulusmap(outerPoly, innerPoly);
    catch ME
        fprintf('  annulusmap skipped for %s: %s\n', desc, ME.message);
        continue
    end

    % qinit
    nptq   = 8;
    qwork  = annulusmap.private_('qinit', m, nptq);
    f      = annulusmap.private_('golden_fields', m);
    cases_qinit(end+1) = struct( ...
        'desc',    ['qinit / ' desc], ...
        'inputs',  struct('M', f.M, 'N', f.N, 'ALFA0', f.ALFA0, 'ALFA1', f.ALFA1, 'nptq', nptq), ...
        'outputs', struct('qwork_size', size(qwork)), ...
        'tol',     1e-14); %#ok<AGROW>

    % eval and evalinv on points in the annular region (u is the inner
    % radius of the canonical annulus directly -- NOT a log-radius -- since
    % |w1(k)| == u by construction in xwtran.m).
    u      = f.u;
    rInner = u + 0.05 * (1 - u);
    rOuter = u + 0.7  * (1 - u);
    thetas = (0:3) * pi/2;
    zp = [rInner * exp(1i*thetas), rOuter * exp(1i*thetas)].';

    wOuter = vertex(outerPoly); aOuter = angle(outerPoly);
    wInner = vertex(innerPoly); aInner = angle(innerPoly);

    wp = eval(m, zp);
    cases_anneval(end+1) = struct( ...
        'desc',    ['annulusmap eval / ' desc], ...
        'inputs',  struct('wOuter', wOuter, 'aOuter', aOuter, 'wInner', wInner, 'aInner', aInner, 'zp', zp), ...
        'outputs', struct('wp', wp), ...
        'tol',     1e-9); %#ok<AGROW>

    zp_inv = evalinv(m, wp);
    cases_annevalinv(end+1) = struct( ...
        'desc',    ['annulusmap evalinv / ' desc], ...
        'inputs',  struct('wOuter', wOuter, 'aOuter', aOuter, 'wInner', wInner, 'aInner', aInner, 'wp', wp), ...
        'outputs', struct('zp', zp_inv), ...
        'tol',     1e-7); %#ok<AGROW>
end

save(fullfile(outDir, 'annulus_private.mat'), ...
    'cases_qinit', 'cases_anneval', 'cases_annevalinv');
fprintf('  annulus_private: qinit=%d eval=%d evalinv=%d\n', ...
    numel(cases_qinit), numel(cases_anneval), numel(cases_annevalinv));
end
