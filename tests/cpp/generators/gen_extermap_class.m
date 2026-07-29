function gen_extermap_class(outDir)
%GEN_EXTERMAP_CLASS  Golden values for the @extermap class-level methods:
%   eval, evalinv, evaldiff, accuracy.
%
%   Uses the same asymmetric-quad/pentagon test polygons as
%   gen_extermap_private.m (the unit square was found there to make
%   deparam's parameter problem non-deterministic across calling contexts;
%   see CPP_PLAN.md 9.4).

opt = sctool.scmapopt('trace', 0, 'tol', 1e-12);

polys = { ...
    polygon([0, 2i, 2+3i, 3]),     'asymmetric quad (CCW exterior)'; ...
    polygon([0, 1, 2, 2+1i, 1i]),  'pentagon exterior'; ...
};

cases_eval     = struct('desc',{},'inputs',{},'outputs',{},'tol',{});
cases_evalinv  = struct('desc',{},'inputs',{},'outputs',{},'tol',{});
cases_evaldiff = struct('desc',{},'inputs',{},'outputs',{},'tol',{});
cases_accuracy = struct('desc',{},'inputs',{},'outputs',{},'tol',{});

for k = 1:size(polys,1)
    p = polys{k,1};

    % Build deparam's literal input directly (not via the extermap(p,opt)
    % constructor), exactly as gen_extermap_private.m now does, to avoid
    % the direct-vs-constructor non-determinism documented in CPP_PLAN.md
    % 9.4.
    w_de = flipud(vertex(p));
    beta_de = 1 - flipud(angle(p));
    [w_de, beta_de] = scfix('de', w_de, beta_de);
    [z, c, qdat] = extermap.private_('deparam', w_de, beta_de, [], opt);

    % Normal (alpha-1) convention, matching m.polygon's vertex()/angle()-1
    % (and what eval/evalinv/evaldiff/accuracy each re-derive internally
    % via flipud(vertex)/flipud(1-angle)). NOTE: angle(poly) = 1-flipud(beta_de)
    % (the *alpha* convention), so beta = angle(poly)-1 = -flipud(beta_de),
    % not 1-flipud(beta_de) -- that's alpha, not beta.
    w = flipud(w_de);
    beta = -flipud(beta_de);

    % eval.m's documented domain is the unit disk's interior (abs(zp)<=1+eps,
    % per its own validity filter) -- the exterior map's domain is still the
    % disk, just with a different SC integrand than diskmap's, so unlike the
    % dederiv/demap-private test points below (which probe the *target*
    % side and intentionally use far-field points), eval's own test points
    % must stay inside the unit disk.
    zp_eval = [0.3+0.2i, -0.5+0.1i, 0.1-0.6i].';
    wp = extermap.private_('demap', zp_eval, w_de, beta_de, z, c, qdat);
    cases_eval(end+1) = struct( ...
        'desc',    ['eval / ' polys{k,2}], ...
        'inputs',  struct('w', w, 'beta', beta, 'z', z, 'c', c, 'qdat', qdat, 'zp', zp_eval), ...
        'outputs', struct('wp', wp), ...
        'tol',     1e-10); %#ok<AGROW>

    zp = 8 * [1.2+0.3i, -1.5+0.8i, 2-0.5i].';

    % deinvmap: well-conditioned points near the unit circle (see 9.1/9.4's
    % notes on deinvmap's poor conditioning for far-field points).
    zp_di = [1.05+0.1i, 1.2-0.1i, -1.1+0.2i].';
    wp_di = extermap.private_('demap', zp_di, w_de, beta_de, z, c, qdat);
    invopt = [0, 1e-12, 80];
    zp_inv = extermap.private_('deinvmap', wp_di, w_de, beta_de, z, c, qdat, [], invopt);
    cases_evalinv(end+1) = struct( ...
        'desc',    ['evalinv / ' polys{k,2}], ...
        'inputs',  struct('w', w, 'beta', beta, 'z', z, 'c', c, 'qdat', qdat, 'wp', wp_di), ...
        'outputs', struct('zp', zp_inv), ...
        'tol',     1e-10); %#ok<AGROW>

    fp = extermap.private_('dederiv', zp, z, beta_de, c);
    cases_evaldiff(end+1) = struct( ...
        'desc',    ['evaldiff / ' polys{k,2}], ...
        'inputs',  struct('z', z, 'beta', beta, 'c', c, 'zp', zp), ...
        'outputs', struct('fp', fp), ...
        'tol',     1e-10); %#ok<AGROW>

    mid = z(1)*exp(1i*angle(z(2)/z(1))/2);
    n = length(z);
    idx = (1:n)';
    idx = [idx idx([2:end 1])];
    dtheta = mod(angle(z(idx(:,2))./z(idx(:,1))),2*pi);
    midv = z(idx(:,1)).*exp(1i*dtheta/2);
    I = extermap.private_('dequad', z(idx(:,1)), midv, idx(:,1), z, beta_de, qdat) - ...
        extermap.private_('dequad', z(idx(:,2)), midv, idx(:,2), z, beta_de, qdat);
    acc = max(abs( c*I - diff(w_de([1:end 1])) ));
    cases_accuracy(end+1) = struct( ...
        'desc',    ['accuracy / ' polys{k,2}], ...
        'inputs',  struct('w', w, 'beta', beta, 'z', z, 'c', c, 'qdat', qdat), ...
        'outputs', struct('acc', acc), ...
        'tol',     1e-9); %#ok<AGROW>
end

save(fullfile(outDir, 'extermap_class.mat'), ...
    'cases_eval', 'cases_evalinv', 'cases_evaldiff', 'cases_accuracy');
fprintf('  extermap_class: eval=%d evalinv=%d evaldiff=%d accuracy=%d\n', ...
    numel(cases_eval), numel(cases_evalinv), numel(cases_evaldiff), numel(cases_accuracy));
end
