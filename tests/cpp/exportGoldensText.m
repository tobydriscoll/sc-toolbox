function exportGoldensText(groupFiles, outDir)
%EXPORTGOLDENSTEXT  Convert tests/cpp/goldens/<name>.mat files into a
%   simple whitespace-delimited text format that the C++ test suite can
%   read without any MATLAB-file-format library (no libmatio dependency).
%
%   exportGoldensText({'gaussj','scqdata','scangle_scfix'})
%   exportGoldensText(groupFiles, outDir)
%
%   Format (one .gold file per input .mat file):
%     GROUP <variable-name>
%     CASE
%     DESC <free text, rest of line>
%     TOL <scalar>
%     INPUT <fieldname> <rows> <cols> <R|C>
%     <rows lines of rows*cols (or 2*rows*cols for C) space-separated doubles>
%     OUTPUT <fieldname> <rows> <cols> <R|C>
%     <...>
%     ENDCASE
%     ENDGROUP
%
%   R fields list one value per element (row-major). C fields list two
%   values (real, imag) per element (row-major). Inf/-Inf/NaN are written
%   as literal tokens, which C++'s strtod/std::stod parse natively.

if nargin < 2
    outDir = fullfile(fileparts(mfilename('fullpath')), 'goldens_text');
end
if ~exist(outDir, 'dir'), mkdir(outDir); end
goldenDir = fullfile(fileparts(mfilename('fullpath')), 'goldens');

if ischar(groupFiles), groupFiles = {groupFiles}; end

for gi = 1:numel(groupFiles)
    name = groupFiles{gi};
    matFile = fullfile(goldenDir, [name '.mat']);
    if ~exist(matFile, 'file')
        fprintf('  skip %s (missing %s)\n', name, matFile);
        continue
    end
    s = load(matFile);
    fid = fopen(fullfile(outDir, [name '.gold']), 'w');
    varNames = fieldnames(s);
    nCases = 0;
    for vi = 1:numel(varNames)
        vname = varNames{vi};
        cases = s.(vname);
        if ~isstruct(cases), continue, end
        fprintf(fid, 'GROUP %s\n', vname);
        for k = 1:numel(cases)
            c = cases(k);
            fprintf(fid, 'CASE\n');
            fprintf(fid, 'DESC %s\n', c.desc);
            fprintf(fid, 'TOL %.17g\n', c.tol);
            writeFieldStruct(fid, 'INPUT', c.inputs);
            writeFieldStruct(fid, 'OUTPUT', c.outputs);
            fprintf(fid, 'ENDCASE\n');
            nCases = nCases + 1;
        end
        fprintf(fid, 'ENDGROUP\n');
    end
    fclose(fid);
    fprintf('  %s: %d cases written to %s.gold\n', name, nCases, name);
end
end

function writeFieldStruct(fid, tag, st)
fields = fieldnames(st);
for fi = 1:numel(fields)
    fname = fields{fi};
    v = st.(fname);
    if ischar(v)
        fprintf(fid, '%s %s STR\n%s\n', tag, fname, v);
    elseif isnumeric(v) || islogical(v)
        v = double(v);
        [r, c] = size(v);
        isComplex = ~isreal(v);
        fprintf(fid, '%s %s %d %d %s\n', tag, fname, r, c, ternary(isComplex, 'C', 'R'));
        for row = 1:r
            rowvals = v(row, :);
            if isComplex
                parts = cell(1, 2*c);
                parts(1:2:end) = arrayfun(@(x) sprintf('%.17g', real(x)), rowvals, 'UniformOutput', false);
                parts(2:2:end) = arrayfun(@(x) sprintf('%.17g', imag(x)), rowvals, 'UniformOutput', false);
            else
                parts = arrayfun(@(x) sprintf('%.17g', x), rowvals, 'UniformOutput', false);
            end
            fprintf(fid, '%s\n', strjoin(parts, ' '));
        end
    else
        fprintf(fid, '%s %s SKIP\n', tag, fname);
        warning('exportGoldensText: skipping non-numeric/non-char field %s (class %s)', fname, class(v));
    end
end
end

function out = ternary(cond, a, b)
if cond, out = a; else, out = b; end
end
