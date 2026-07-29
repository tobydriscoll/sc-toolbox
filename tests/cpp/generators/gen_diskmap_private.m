function gen_diskmap_private(outDir)
%GEN_DISKMAP_PRIVATE  Golden values for @diskmap/private functions:
%   dparam, dquad, dderiv, dmap, dimapfun, dinvmap.
%
%   Private functions are accessed via the diskmap.private_(name, ...)
%   static method (see @diskmap/diskmap.m), since MATLAB does not allow
%   adding class-folder private/ directories to the path directly.

opt = sctool.scmapopt('trace', 0, 'tol', 1e-12);

% Polygon library: square and L-shape
polys = { ...
    polygon([0, 1, 1+1i, 1i]),             'square'; ...
    polygon([0, 2, 2+1i, 1+1i, 1+2i, 2i]), 'L-shape'; ...
    polygon([0, 1, 1+0.5i, 0.5+1i, 0+1i]), 'irregular quad'; ...
};

% --- dparam ---
cases_dparam = struct('desc', {}, 'inputs', {}, 'outputs', {}, 'tol', {});
for k = 1:size(polys,1)
    p    = polys{k,1};
    w    = vertex(p);
    beta = angle(p) - 1;
    [w2, beta2] = scfix('d', w, beta);
    [z, c, qdat] = diskmap.private_('dparam', w2, beta2, [], opt);
    cases_dparam(end+1) = struct( ...
        'desc',    ['dparam / ' polys{k,2}], ...
        'inputs',  struct('w', w2, 'beta', beta2), ...
        'outputs', struct('z', z, 'c', c, 'qdat', qdat), ...
        'tol',     1e-10); %#ok<AGROW>
end

% --- dquad ---
cases_dquad = struct('desc', {}, 'inputs', {}, 'outputs', {}, 'tol', {});
for k = 1:size(polys,1)
    p    = polys{k,1};
    w    = vertex(p); beta = angle(p) - 1;
    [w2, beta2] = scfix('d', w, beta);
    [z, c, qdat] = diskmap.private_('dparam', w2, beta2, [], opt);
    n = length(z);

    % Singular-to-interior
    midArc = exp(1i * (angle(z(1)) + angle(z(2))) / 2) * 0.9;
    I = diskmap.private_('dquad', z(1), midArc, 1, z, beta2, qdat);
    cases_dquad(end+1) = struct( ...
        'desc',    ['dquad sing->interior / ' polys{k,2}], ...
        'inputs',  struct('z1', z(1), 'z2', midArc, 'sing1', 1, 'z', z, 'beta', beta2, 'qdat', qdat), ...
        'outputs', struct('I', I), ...
        'tol',     1e-11); %#ok<AGROW>

    % Regular interior-to-interior
    p1 = 0.3+0.2i; p2 = -0.1+0.5i;
    I2 = diskmap.private_('dquad', p1, p2, 0, z, beta2, qdat);
    cases_dquad(end+1) = struct( ...
        'desc',    ['dquad interior->interior / ' polys{k,2}], ...
        'inputs',  struct('z1', p1, 'z2', p2, 'sing1', 0, 'z', z, 'beta', beta2, 'qdat', qdat), ...
        'outputs', struct('I', I2), ...
        'tol',     1e-11); %#ok<AGROW>
end

% --- dabsquad (+sctool/dabsquad.m; used internally by dparam's residual fn) ---
cases_dabsquad = struct('desc', {}, 'inputs', {}, 'outputs', {}, 'tol', {});
for k = 1:size(polys,1)
    p    = polys{k,1};
    w    = vertex(p); beta = angle(p) - 1;
    [w2, beta2] = scfix('d', w, beta);
    [z, c, qdat] = diskmap.private_('dparam', w2, beta2, [], opt);

    % Singular-to-interior-point, along the unit circle (the "mid" point
    % dpfun.m uses): same scenario dabsquad is actually called with.
    midArc = exp(1i * (angle(z(1)) + angle(z(2))) / 2) * 0.9;
    I = sctool.dabsquad(z(1), midArc, 1, z, beta2, qdat);
    cases_dabsquad(end+1) = struct( ...
        'desc',    ['dabsquad sing->interior / ' polys{k,2}], ...
        'inputs',  struct('z1', z(1), 'z2', midArc, 'sing1', 1, 'z', z, 'beta', beta2, 'qdat', qdat), ...
        'outputs', struct('I', I), ...
        'tol',     1e-11); %#ok<AGROW>

    % Two adjacent singularities (both endpoints on the unit circle, no
    % singularity flag at either end -- mirrors dpfun's two-sided sum).
    midArc2 = exp(1i * (angle(z(2)) + angle(z(3))) / 2);
    I2a = sctool.dabsquad(z(2), midArc2, 1, z, beta2, qdat);
    I2b = sctool.dabsquad(z(3), midArc2, 2, z, beta2, qdat);
    cases_dabsquad(end+1) = struct( ...
        'desc',    ['dabsquad two-sided / ' polys{k,2}], ...
        'inputs',  struct('z1', [z(2);z(3)], 'z2', [midArc2;midArc2], 'sing1', [1;2], ...
                           'z', z, 'beta', beta2, 'qdat', qdat), ...
        'outputs', struct('I', [I2a;I2b]), ...
        'tol',     1e-11); %#ok<AGROW>
end

% --- dderiv ---
cases_dderiv = struct('desc', {}, 'inputs', {}, 'outputs', {}, 'tol', {});
for k = 1:size(polys,1)
    p    = polys{k,1};
    w    = vertex(p); beta = angle(p) - 1;
    [w2, beta2] = scfix('d', w, beta);
    [z, c, qdat] = diskmap.private_('dparam', w2, beta2, [], opt);

    zp = [0.1+0.2i, -0.3+0.15i, 0.05-0.4i].';
    fp = diskmap.private_('dderiv', zp, z, beta2, c);
    cases_dderiv(end+1) = struct( ...
        'desc',    ['dderiv / ' polys{k,2}], ...
        'inputs',  struct('zp', zp, 'z', z, 'beta', beta2, 'c', c), ...
        'outputs', struct('fp', fp), ...
        'tol',     1e-11); %#ok<AGROW>
end

% --- dmap ---
cases_dmap = struct('desc', {}, 'inputs', {}, 'outputs', {}, 'tol', {});
for k = 1:size(polys,1)
    p    = polys{k,1};
    w    = vertex(p); beta = angle(p) - 1;
    [w2, beta2] = scfix('d', w, beta);
    [z, c, qdat] = diskmap.private_('dparam', w2, beta2, [], opt);

    zp = [0.1+0.2i, -0.3+0.15i, 0.05-0.4i, 0+0i].';
    wp = diskmap.private_('dmap', zp, w2, beta2, z, c, qdat);
    cases_dmap(end+1) = struct( ...
        'desc',    ['dmap / ' polys{k,2}], ...
        'inputs',  struct('zp', zp, 'w', w2, 'beta', beta2, 'z', z, 'c', c, 'qdat', qdat), ...
        'outputs', struct('wp', wp), ...
        'tol',     1e-11); %#ok<AGROW>
end

% --- dinvmap ---
cases_dinvmap = struct('desc', {}, 'inputs', {}, 'outputs', {}, 'tol', {});
for k = 1:size(polys,1)
    p    = polys{k,1};
    w    = vertex(p); beta = angle(p) - 1;
    [w2, beta2] = scfix('d', w, beta);
    [z, c, qdat] = diskmap.private_('dparam', w2, beta2, [], opt);

    zp_ref = [0.1+0.2i, -0.3+0.15i, 0.05-0.4i].';
    wp_ref = diskmap.private_('dmap', zp_ref, w2, beta2, z, c, qdat);
    zp_inv = diskmap.private_('dinvmap', wp_ref, w2, beta2, z, c, qdat);
    cases_dinvmap(end+1) = struct( ...
        'desc',    ['dinvmap / ' polys{k,2}], ...
        'inputs',  struct('wp', wp_ref, 'w', w2, 'beta', beta2, 'z', z, 'c', c, 'qdat', qdat), ...
        'outputs', struct('zp', zp_inv), ...
        'tol',     1e-10); %#ok<AGROW>
end

save(fullfile(outDir, 'diskmap_private.mat'), ...
    'cases_dparam', 'cases_dquad', 'cases_dabsquad', 'cases_dderiv', 'cases_dmap', 'cases_dinvmap');
fprintf('  diskmap_private: dparam=%d dquad=%d dabsquad=%d dderiv=%d dmap=%d dinvmap=%d\n', ...
    numel(cases_dparam), numel(cases_dquad), numel(cases_dabsquad), numel(cases_dderiv), ...
    numel(cases_dmap), numel(cases_dinvmap));
end
