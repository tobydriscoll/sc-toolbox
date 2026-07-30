function wp = evalRaw(M, zp)
%EVALRAW Evaluate SC rectangle map at points, bypassing domain check.
%   EVALRAW(M,ZP) is like EVAL but does not restrict to points inside the
%   source rectangle.  Needed for boundary and corner evaluation.
%
%   See also EVAL, RECTMAP.

p    = polygon(M);
beta = angle(p) - 1;
wp   = rmap(zp, vertex(p), beta, M.prevertex, M.constant, M.stripL, M.qdata);
