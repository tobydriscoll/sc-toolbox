function [dfdwre, dfdwim] = vertexDeriv(M, zp)
%VERTEXDERIV  Shape derivative of SC rectangle map w.r.t. polygon vertices.
%   [DFDWRE, DFDWIM] = VERTEXDERIV(M, ZP) returns the derivatives of the
%   map M at preimage points ZP with respect to each polygon vertex.
%
%   DFDWRE(k,i) = df/d(Re w_i)  at ZP(k)
%   DFDWIM(k,i) = df/d(Im w_i)  at ZP(k)
%
%   Vertex indexing follows the polygon stored in M (post-SCFIX order).
%
%   Implementation: central finite differences in vertex space.  Each
%   perturbed polygon is solved with a warm start from the current map,
%   so the n full parameter solves are fast.
%
%   See also RECTMAP, EVALRAW, EVALDIFF.

if isempty(M.stripdat)
    error('rectmap:noStripdat', ...
        'vertexDeriv requires stripdat. Reconstruct the map from a polygon.');
end

p     = polygon(M);
w     = vertex(p);      % n x 1, post-SCFIX renumbered order
n     = length(w);
zp    = zp(:);
nz    = numel(zp);

h     = 1e-4;           % central FD step (same as test baseline)

dfdwre = zeros(nz, n);
dfdwim = zeros(nz, n);

for i = 1:n
    % d/d(Re w_i)
    wp = w; wp(i) = wp(i) + h;
    wm = w; wm(i) = wm(i) - h;
    fp = evalRaw(rectmap(M, polygon(wp)), zp);
    fm = evalRaw(rectmap(M, polygon(wm)), zp);
    dfdwre(:, i) = (fp - fm) / (2*h);

    % d/d(Im w_i)
    wp = w; wp(i) = wp(i) + 1i*h;
    wm = w; wm(i) = wm(i) - 1i*h;
    fp = evalRaw(rectmap(M, polygon(wp)), zp);
    fm = evalRaw(rectmap(M, polygon(wm)), zp);
    dfdwim(:, i) = (fp - fm) / (2*h);
end

end
