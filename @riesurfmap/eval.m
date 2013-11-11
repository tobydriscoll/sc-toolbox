function wp = eval(M,zp,tol)
%EVAL Evaluate Schwarz-Christoffel disk map at points.
%   EVAL(M,ZP) evaluates the Schwarz-Christoffel map M at the points ZP
%   in the unit disk. The default tolerance of M is used.
%   
%   EVAL(M,ZP,TOL) attempts to give an answer accurate to TOL. If TOL is
%   less than the accuracy of M, this is unlikely to be met.
%   
%   See also DISKMAP, EVALINV.

%   Copyright 1998 by Toby Driscoll.
%   $Id: eval.m 198 2002-09-13 18:44:12Z driscoll $

if nargin < 3
  qdata = M.qdata;
else 
  qdata = tol;
end

p = polygon(M);
wp = NaN*zp;
idx = abs(zp) <= 1+eps;
zb = M.prebranch;
wp(idx) = rsmap(zp(idx),vertex(p),angle(p)-1,M.prevertex,zb,M.constant,qdata);
