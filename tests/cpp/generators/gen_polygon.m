function gen_polygon(outDir)
%GEN_POLYGON  Golden values for the @polygon constructor (vertex/angle/isinf).

cases = struct('desc', {}, 'inputs', {}, 'outputs', {}, 'tol', {});

% Each row: {w (column), alpha or [] to auto-compute, desc}.
specs = { ...
    [0,1,1+1i,1i], [], 'unit square, already CCW'; ...
    [0,1i,1+1i,1], [], 'unit square, CW input gets reversed'; ...
    [1i,-1i,Inf], [1.5,0.5,-1], 'unbounded triangle, explicit angles'; ...
    [0,2,2+1i,1+1i,1+2i,2i], [], 'L-shape, already CCW'; ...
};

for k = 1:size(specs, 1)
    w = specs{k,1}(:);
    alphaSpec = specs{k,2};
    if isempty(alphaSpec)
        p = polygon(w);
    else
        p = polygon(w, alphaSpec(:));
    end
    inp = struct('w', w);
    if ~isempty(alphaSpec)
        inp.alpha = alphaSpec(:);
    else
        inp.alpha = zeros(0,1);
    end
    out = struct('vertex', p.vertex, 'angle', p.angle, 'isinf', double(isinf(p)));
    cases(end+1) = struct( ...
        'desc',    specs{k,3}, ...
        'inputs',  inp, ...
        'outputs', out, ...
        'tol',     1e-12); %#ok<AGROW>
end

cases_polygon = cases;
save(fullfile(outDir, 'polygon.mat'), 'cases_polygon');
fprintf('  polygon: %d cases\n', numel(cases_polygon));
end
