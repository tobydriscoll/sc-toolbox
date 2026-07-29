function gen_composite(outDir)
%GEN_COMPOSITE  Golden values for the @composite class (eval and inv).
%
%   Each case pins a composition built from a @diskmap on a fixed polygon
%   and one or two @moebius maps, so the C++ side can rebuild every member
%   from the stored polygon vertices/angles and Moebius coefficients.
%
%   Groups:
%     cases_eval     forward evaluation of the composition
%     cases_inv      evaluation of inv(composite), i.e. the reversed chain
%                    of inverted members
%     cases_flatten  member count after nesting one composite inside another
%
%   The 'order' field records which member comes first:
%     'map_then_mob'      diskmap, then moebius
%     'mob_then_map'      moebius (disk -> disk), then diskmap
%     'mob_map_mob'       three members

opt = sctool.scmapopt('trace', 0, 'tol', 1e-12);

% A disk automorphism, so it can legally precede a diskmap: (z-a)/(1-conj(a)z).
a = 0.3 - 0.2i;
blaschke = [-a, 1, 1, -conj(a)];
% A generic (non-automorphism) Moebius map, used only after a diskmap.
generic = [1+1i, 2, 3, -1i];

polys = { ...
    polygon([0, 1, 1+1i, 1i]),             'square'; ...
    polygon([0, 2, 2+1i, 1+1i, 1+2i, 2i]), 'L-shape'; ...
};

cases_eval    = struct('desc',{},'inputs',{},'outputs',{},'tol',{});
cases_inv     = struct('desc',{},'inputs',{},'outputs',{},'tol',{});
cases_flatten = struct('desc',{},'inputs',{},'outputs',{},'tol',{});

for k = 1:size(polys,1)
    p = polys{k,1};
    m = diskmap(p, opt);
    w    = vertex(m.polygon);
    beta = angle(m.polygon) - 1;

    Mb = moebius(blaschke);
    Mg = moebius(generic);

    % Interior disk points, kept away from the boundary so that the inverse
    % chain (which runs evalinv on the diskmap) stays well conditioned.
    zp = [0.1+0.2i, -0.3+0.15i, 0.05-0.4i, 0+0i].';

    orders = { ...
        'map_then_mob', composite(m, Mg); ...
        'mob_then_map', composite(Mb, m); ...
        'mob_map_mob',  composite(Mb, m, Mg); ...
    };

    for j = 1:size(orders,1)
        f  = orders{j,2};
        wp = eval(f, zp);

        cases_eval(end+1) = struct( ...
            'desc',    ['eval / ' polys{k,2} ' / ' orders{j,1}], ...
            'inputs',  struct('w', w, 'beta', beta, 'coeff_pre', blaschke(:), ...
                              'coeff_post', generic(:), 'order', orders{j,1}, 'zp', zp), ...
            'outputs', struct('wp', wp), ...
            'tol',     1e-10); %#ok<AGROW>

        fi = inv(f);
        zp_inv = eval(fi, wp);
        cases_inv(end+1) = struct( ...
            'desc',    ['inv / ' polys{k,2} ' / ' orders{j,1}], ...
            'inputs',  struct('w', w, 'beta', beta, 'coeff_pre', blaschke(:), ...
                              'coeff_post', generic(:), 'order', orders{j,1}, 'wp', wp), ...
            'outputs', struct('zp', zp_inv), ...
            'tol',     1e-8); %#ok<AGROW>
    end

    % A composite passed to composite() is flattened into its members.
    f2 = composite(Mb, m);
    fn = composite(f2, Mg);
    cases_flatten(end+1) = struct( ...
        'desc',    ['flatten / ' polys{k,2} ' / composite(composite(mob,map), mob)'], ...
        'inputs',  struct('w', w, 'beta', beta, 'coeff_pre', blaschke(:), ...
                          'coeff_post', generic(:), 'zp', zp), ...
        'outputs', struct('nmaps', numel(members(fn)), 'wp', eval(fn, zp)), ...
        'tol',     1e-10); %#ok<AGROW>
end

save(fullfile(outDir, 'composite.mat'), 'cases_eval', 'cases_inv', 'cases_flatten');
fprintf('  composite: eval=%d inv=%d flatten=%d\n', ...
    numel(cases_eval), numel(cases_inv), numel(cases_flatten));
end
